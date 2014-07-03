module histFldsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: histFldsMod
!
! !DESCRIPTION:
! Module containing initialization of clm history fields and files
! This is the module that the user must modify in order to add new
! history fields or modify defaults associated with existing history
! fields.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varcon, only: dzsoi_decomp
  use clm_varctl  , only : iulog
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public hist_initFlds ! Build master field list of all possible history
                       ! file fields
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 03/2003
! heald (11/28/06)
! F. Li and S. Levis (11/06/12)
!EOP
!------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hist_initFlds
!
! !INTERFACE:
  subroutine hist_initFlds()
!
! !DESCRIPTION:
! Build master field list of all possible fields in a history file.
! Each field has associated with it a ``long\_name'' netcdf attribute that
! describes what the field is, and a ``units'' attribute. A subroutine is
! called to add each field to the masterlist.
!
! !USES:
    use clmtype
    use clm_varcon , only : spval
    use clm_varpar , only : maxpatch_glcmec
    use clm_atmlnd , only : clm_a2l
    use clm_varctl , only : create_glacier_mec_landunit, &
                            use_c13, use_c14
#if (defined LCH4)
    use clm_atmlnd , only : clm_l2a
    use histFileMod, only : hist_wrtch4diag
    use ch4varcon  , only : allowlakeprod
#endif
    use clm_glclnd , only : clm_s2x
    use histFileMod, only : hist_add_subscript, hist_addfld1d, hist_addfld2d, &
                            hist_printflds
    use surfrdMod  , only : crop_prog
    use shr_megan_mod  , only : shr_megan_linkedlist, shr_megan_megcomp_t, shr_megan_megcomps_n

    use clm_varpar , only :  ndecomp_cascade_transitions, ndecomp_pools, nlevdecomp, nlevdecomp_full
!
! !ARGUMENTS:
    implicit none

    type(shr_megan_megcomp_t), pointer :: meg_cmp
    integer :: imeg
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 03/2003
! Mariana Vertenstein: Updated interface to create history fields 10/2003
!
!EOP
!-----------------------------------------------------------------------
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    integer :: k,l, ii, jj
    character(24) :: fieldname
    character(100) :: longname
    character(8) :: vr_suffix
    character(10) :: active

    ! Determine what subscripts to add
    ! (uncomment the following call and modify it appropriately)

    ! call hist_add_subscript(subname='subscript_name', subdim=subscript_dim)

    ! NOTE: make a field not appear on the primary history tape by default -
    ! add the keyword to default='inactive' to the call to addfld_1d or addfld_2d

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp .gt. 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif


    call hist_addfld1d (fname='H2OSFC',  units='mm',  &
         avgflag='A', long_name='surface water depth', &
         ptr_col=clm3%g%l%c%cws%h2osfc)

    call hist_addfld1d (fname='FH2OSFC',  units='unitless',  &
         avgflag='A', long_name='fraction of ground covered by surface water', &
         ptr_col=clm3%g%l%c%cps%frac_h2osfc)

    call hist_addfld1d (fname='QH2OSFC',  units='mm/s',  &
         avgflag='A', long_name='surface water runoff', &
         ptr_col=clm3%g%l%c%cwf%qflx_h2osfc_surf)

    call hist_addfld1d (fname='TH2OSFC',  units='K',  &
         avgflag='A', long_name='surface water temperature', &
         ptr_col=clm3%g%l%c%ces%t_h2osfc)

    call hist_addfld1d (fname='QFLOOD',  units='mm/s',  &
         avgflag='A', long_name='runoff from river flooding', &
         ptr_lnd=clm_a2l%forc_flood)

    call hist_addfld1d (fname='TWS',  units='mm',  &
         avgflag='A', long_name='total water storage', &
         ptr_lnd=clm3%g%tws)

    call hist_addfld1d (fname='VOLR',  units='m3',  &
         avgflag='A', long_name='river channel water storage', &
         ptr_lnd=clm_a2l%volr)


    ! Snow properties
    ! These will be vertically averaged over the snow profile

    call hist_addfld1d (fname='SNOW_DEPTH',  units='m',  &
         avgflag='A', long_name='snow height of snow covered area', &
         ptr_col=clm3%g%l%c%cps%snow_depth, c2l_scale_type='urbanf')!, default='inactive')

    call hist_addfld1d (fname='SNOWDP',  units='m',  &
         avgflag='A', long_name='gridcell mean snow height', &
         ptr_col=clm3%g%l%c%cps%snowdp, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSNO',  units='unitless',  &
         avgflag='A', long_name='fraction of ground covered by snow', &
         ptr_col=clm3%g%l%c%cps%frac_sno, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSNO_EFF',  units='unitless',  &
         avgflag='A', long_name='effective fraction of ground covered by snow', &
         ptr_col=clm3%g%l%c%cps%frac_sno_eff, c2l_scale_type='urbanf')!, default='inactive')

    ! Temperatures

    call hist_addfld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m)

    call hist_addfld1d (fname='TSA_U', units='K',  &
         avgflag='A', long_name='Urban 2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_u, set_nourb=spval)

    call hist_addfld1d (fname='TSA_R', units='K',  &
         avgflag='A', long_name='Rural 2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_r, set_spec=spval)

    call hist_addfld1d(fname='TBUILD', units='K',  &
         avgflag='A', long_name='internal urban building temperature', &
         ptr_lunit=clm3%g%l%lps%t_building, set_nourb=spval, l2g_scale_type='unity')

    call hist_addfld1d (fname='TREFMNAV', units='K',  &
         avgflag='A', long_name='daily minimum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_min)

    call hist_addfld1d (fname='TREFMXAV', units='K',  &
         avgflag='A', long_name='daily maximum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_max)

    call hist_addfld1d (fname='TREFMNAV_U', units='K',  &
         avgflag='A', long_name='Urban daily minimum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_min_u, set_nourb=spval)

    call hist_addfld1d (fname='TREFMXAV_U', units='K',  &
         avgflag='A', long_name='Urban daily maximum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_max_u, set_nourb=spval)

    call hist_addfld1d (fname='TREFMNAV_R', units='K',  &
         avgflag='A', long_name='Rural daily minimum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_min_r, set_spec=spval)

    call hist_addfld1d (fname='TREFMXAV_R', units='K',  &
         avgflag='A', long_name='Rural daily maximum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_max_r, set_spec=spval)

    call hist_addfld1d (fname='TV', units='K',  &
         avgflag='A', long_name='vegetation temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_veg)

    call hist_addfld1d (fname='TV24', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 24hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%t_veg24, default='inactive')

    call hist_addfld1d (fname='TV240', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 240hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%t_veg240, default='inactive')

    call hist_addfld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=clm3%g%l%c%ces%t_grnd, c2l_scale_type='urbans')

    call hist_addfld1d (fname='TG_U', units='K',  &
         avgflag='A', long_name='Urban ground temperature', &
         ptr_col=clm3%g%l%c%ces%t_grnd_u, set_nourb=spval, c2l_scale_type='urbans')

    call hist_addfld1d (fname='TG_R', units='K',  &
         avgflag='A', long_name='Rural ground temperature', &
         ptr_col=clm3%g%l%c%ces%t_grnd_r, set_spec=spval)

    call hist_addfld1d (fname='HCSOI',  units='MJ/m2',  &
         avgflag='A', long_name='soil heat content', &
         ptr_col=clm3%g%l%c%ces%hc_soi, set_lake=spval, set_urb=spval, l2g_scale_type='veg')

    call hist_addfld1d (fname='HC',  units='MJ/m2',  &
         avgflag='A', long_name='heat content of soil/snow/lake', &
         ptr_col=clm3%g%l%c%ces%hc_soisno, set_urb=spval)

    call hist_addfld2d (fname='TSOI',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%ces%t_soisno, l2g_scale_type='veg')

    call hist_addfld2d (fname='TSOI_ICE',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (ice landunits only)', &
         ptr_col=clm3%g%l%c%ces%t_soisno, l2g_scale_type='ice')

    call hist_addfld1d (fname='TSOI_10CM',  units='K', &
         avgflag='A', long_name='soil temperature in top 10cm of soil', &
         ptr_col=clm3%g%l%c%ces%t_soi_10cm, set_urb=spval)

    call hist_addfld2d (fname='TLAKE',  units='K', type2d='levlak', &
         avgflag='A', long_name='lake temperature', &
         ptr_col=clm3%g%l%c%ces%t_lake)

    ! New lake fields
    call hist_addfld2d (fname='LAKEICEFRAC',  units='unitless', type2d='levlak', &
         avgflag='A', long_name='lake layer ice mass fraction', &
         ptr_col=clm3%g%l%c%cws%lake_icefrac)

         ! This will be more useful than LAKEICEFRAC for many users.
    call hist_addfld1d (fname='LAKEICETHICK', units='m', &
         avgflag='A', long_name='thickness of lake ice (including physical expansion on freezing)', &
         ptr_col=clm3%g%l%c%cws%lake_icethick, set_nolake=spval)

    call hist_addfld1d (fname='TKE1',  units='W/(mK)', &
         avgflag='A', long_name='top lake level eddy thermal conductivity', &
         ptr_col=clm3%g%l%c%cps%savedtke1)

    call hist_addfld1d (fname='EFLX_GRND_LAKE', units='W/m^2', &
         avgflag='A', long_name='net heat flux into lake/snow surface, excluding light transmission', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_grnd_lake, set_nolake=spval)

    call hist_addfld1d (fname='RAM_LAKE', units='s/m', &
         avgflag='A', long_name='aerodynamic resistance for momentum (lakes only)', &
         ptr_pft=clm3%g%l%c%p%pps%ram1_lake, set_nolake=spval, default='inactive')

    call hist_addfld1d (fname='RH_LEAF', units='fraction', &
         avgflag='A', long_name='fractional humidity at leaf surface', &
         ptr_pft=clm3%g%l%c%p%pps%rh_leaf, set_spec=spval, default='inactive')

    call hist_addfld1d (fname='RHAF', units='fraction', &
         avgflag='A', long_name='fractional humidity of canopy air', &
         ptr_pft=clm3%g%l%c%p%pps%rhaf, set_spec=spval, default='inactive')

    call hist_addfld1d (fname='UST_LAKE', units='m/s', &
         avgflag='A', long_name='friction velocity (lakes only)', &
         ptr_col=clm3%g%l%c%cps%ust_lake, set_nolake=spval, default='inactive')

    call hist_addfld1d (fname='Z0MG', units='m', &
         avgflag='A', long_name='roughness length over ground, momentum', &
         ptr_col=clm3%g%l%c%cps%z0mg, default='inactive')

    call hist_addfld1d (fname='Z0HG', units='m', &
         avgflag='A', long_name='roughness length over ground, sensible heat', &
         ptr_col=clm3%g%l%c%cps%z0hg, default='inactive')

    call hist_addfld1d (fname='Z0QG', units='m', &
         avgflag='A', long_name='roughness length over ground, latent heat', &
         ptr_col=clm3%g%l%c%cps%z0qg, default='inactive')

    ! End new lake fields

    ! Allow active layer fields to be optionally output even if not running CN
#ifndef CN
    call hist_addfld1d (fname='ALT', units='m', &
         avgflag='A', long_name='current active layer thickness', &
         ptr_col=clm3%g%l%c%cps%alt, default='inactive')

    call hist_addfld1d (fname='ALTMAX', units='m', &
         avgflag='A', long_name='maximum annual active layer thickness', &
         ptr_col=clm3%g%l%c%cps%altmax, default='inactive')

    call hist_addfld1d (fname='ALTMAX_LASTYEAR', units='m', &
         avgflag='A', long_name='maximum prior year active layer thickness', &
         ptr_col=clm3%g%l%c%cps%altmax_lastyear, default='inactive')
#endif

    ! Specific humidity

    call hist_addfld1d (fname='Q2M', units='kg/kg',  &
         avgflag='A', long_name='2m specific humidity', &
         ptr_pft=clm3%g%l%c%p%pes%q_ref2m)

    ! Relative humidity

    call hist_addfld1d (fname='RH2M', units='%',  &
         avgflag='A', long_name='2m relative humidity', &
         ptr_pft=clm3%g%l%c%p%pes%rh_ref2m)

    call hist_addfld1d (fname='RH2M_U', units='%',  &
         avgflag='A', long_name='Urban 2m relative humidity', &
         ptr_pft=clm3%g%l%c%p%pes%rh_ref2m_u, set_nourb=spval)

    call hist_addfld1d (fname='RH2M_R', units='%',  &
         avgflag='A', long_name='Rural 2m specific humidity', &
         ptr_pft=clm3%g%l%c%p%pes%rh_ref2m_r, set_spec=spval)

    ! Wind

    call hist_addfld1d (fname='U10', units='m/s', &
         avgflag='A', long_name='10-m wind', &
         ptr_pft=clm3%g%l%c%p%pps%u10_clm)
    call hist_addfld1d (fname='VA', units='m/s', &
         avgflag='A', long_name='atmospheric wind speed plus convective velocity', &
         ptr_pft=clm3%g%l%c%p%pps%va, default='inactive')

    ! Surface radiation

    call hist_addfld1d (fname='SABV', units='W/m^2',  &
         avgflag='A', long_name='solar rad absorbed by veg', &
         ptr_pft=clm3%g%l%c%p%pef%sabv, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='SABG', units='W/m^2',  &
         avgflag='A', long_name='solar rad absorbed by ground', &
         ptr_pft=clm3%g%l%c%p%pef%sabg, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='SABG_PEN', units='watt/m^2',  &
         avgflag='A', long_name='Rural solar rad penetrating top soil or snow layer', &
         ptr_pft=clm3%g%l%c%p%pef%sabg_pen, set_spec=spval)

    call hist_addfld1d (fname='FSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_d)

    call hist_addfld1d (fname='FSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_d)

    call hist_addfld1d (fname='FSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_i)

    call hist_addfld1d (fname='FSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_i)

    call hist_addfld1d (fname='FSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_d, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_d, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_i, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_i, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSDSVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_d_ln)

    call hist_addfld1d (fname='FSDSVILN', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_i_ln)

    call hist_addfld1d (fname='PARVEGLN', units='W/m^2',  &
         avgflag='A', long_name='absorbed par by vegetation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%parveg_ln)

    call hist_addfld1d (fname='FSDSNDLN', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_d_ln)

    call hist_addfld1d (fname='FSRVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_d_ln, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRNDLN', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_d_ln, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSA', units='W/m^2',  &
         avgflag='A', long_name='absorbed solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsa, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban absorbed solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsa_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural absorbed solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsa_r, set_spec=spval)

    call hist_addfld1d (fname='FSR', units='W/m^2',  &
         avgflag='A', long_name='reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr, c2l_scale_type='urbanf')

    ! Rename of FSR for Urban intercomparision project
    call hist_addfld1d (fname='SWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling shortwave radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FIRA', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_net, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FIRA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban net infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_net_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FIRA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural net infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_net_r, set_spec=spval)

    call hist_addfld1d (fname='FIRE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_out, c2l_scale_type='urbanf')

    ! Rename of FIRE for Urban intercomparision project
    call hist_addfld1d (fname='LWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling longwave radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_out, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='BUILDHEAT', units='W/m^2',  &
         avgflag='A', long_name='heat flux from urban building interior to walls and roof', &
         ptr_col=clm3%g%l%c%cef%eflx_building_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
         avgflag='A', long_name='urban air conditioning flux', &
         ptr_col=clm3%g%l%c%cef%eflx_urban_ac, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
         avgflag='A', long_name='urban heating flux', &
         ptr_col=clm3%g%l%c%cef%eflx_urban_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='TRAFFICFLUX', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from urban traffic', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_traffic_pft, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    call hist_addfld1d (fname='WASTEHEAT', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from heating/cooling sources of urban waste heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_wasteheat_pft, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='HEAT_FROM_AC', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux put into canyon due to heat removed from air conditioning', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_heat_from_ac_pft, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='Qanth', units='W/m^2',  &
         avgflag='A', long_name='anthropogenic heat flux', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_anthro, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    call hist_addfld1d (fname='Rnet', units='W/m^2',  &
         avgflag='A', long_name='net radiation', &
         ptr_pft=clm3%g%l%c%p%pef%netrad, c2l_scale_type='urbanf', &
         default='inactive')

    ! Solar zenith angle and solar declination angle

    call hist_addfld1d (fname='COSZEN', units='none', &
         avgflag='A', long_name='cosine of solar zenith angle', &
         ptr_col=clm3%g%l%c%cps%coszen, default='inactive')

    call hist_addfld1d (fname='DECL', units='radians', &
         avgflag='A', long_name='solar declination angle', &
         ptr_col=clm3%g%l%c%cps%decl, default='inactive')

    ! Surface energy fluxes

    call hist_addfld1d (fname='FCTR', units='W/m^2',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_vegt, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FCEV', units='W/m^2',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_vege, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGEV', units='W/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_grnd, c2l_scale_type='urbanf') 

    call hist_addfld1d (fname='FSH_NODYNLNDUSE', units='W/m^2',  &
         avgflag='A', long_name='sensible heat not including correction for land use change', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSH', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_lnd=clm3%g%gef%eflx_sh_totg)

    call hist_addfld1d (fname='FSH_U', units='W/m^2',  &
         avgflag='A', long_name='Urban sensible heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSH_R', units='W/m^2',  &
         avgflag='A', long_name='Rural sensible heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot_r, set_spec=spval)

    call hist_addfld1d (fname='GC_HEAT1',  units='J/m^2',  &
         avgflag='A', long_name='initial gridcell total heat content', &
         ptr_lnd=clm3%g%ges%gc_heat1)

    call hist_addfld1d (fname='GC_HEAT2',  units='J/m^2',  &
         avgflag='A', long_name='post land cover change total heat content', &
         ptr_lnd=clm3%g%ges%gc_heat2, default='inactive')

    call hist_addfld1d (fname='EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_lnd=clm3%g%gef%eflx_dynbal)

    call hist_addfld1d (fname='Qh', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='Qle', units='W/m^2',  &
         avgflag='A', long_name='total evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='EFLX_LH_TOT', units='W/m^2', &
         avgflag='A', long_name='total latent heat flux [+ to atm]', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_tot, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='EFLX_LH_TOT_U', units='W/m^2',  &
         avgflag='A', long_name='Urban total evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='EFLX_LH_TOT_R', units='W/m^2',  &
         avgflag='A', long_name='Rural total evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_tot_r, set_spec=spval)

    call hist_addfld1d (fname='Qstor', units='W/m^2',  &
         avgflag='A', long_name='storage heat flux (includes snowmelt)', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='FSH_V', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from veg', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSH_G', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGR', units='W/m^2',  &
!         avgflag='A', long_name='heat flux into soil/snow including snow melt', &
         avgflag='A', long_name='heat flux into soil/snow including snow melt and lake / snow light transmission', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGR_U', units='W/m^2',  &
         avgflag='A', long_name='Urban heat flux into soil/snow including snow melt', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FGR_R', units='W/m^2',  &
!         avgflag='A', long_name='Rural heat flux into soil/snow including snow melt', &
         avgflag='A', long_name='Rural heat flux into soil/snow including snow melt and snow light transmission', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd_r, set_spec=spval)

    call hist_addfld1d (fname='FSM',  units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=clm3%g%l%c%cef%eflx_snomelt, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSM_U',  units='W/m^2',  &
         avgflag='A', long_name='Urban snow melt heat flux', &
         ptr_col=clm3%g%l%c%cef%eflx_snomelt_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSM_R',  units='W/m^2',  &
         avgflag='A', long_name='Rural snow melt heat flux', &
         ptr_col=clm3%g%l%c%cef%eflx_snomelt_r, set_spec=spval)

    call hist_addfld1d (fname='FGR12',  units='W/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=clm3%g%l%c%cef%eflx_fgr12, set_lake=spval)

    call hist_addfld2d (fname='FGR_SOIL_R', units='watt/m^2', type2d='levgrnd', &
         avgflag='A', long_name='Rural downward heat flux at interface below each soil layer', &
         ptr_col=clm3%g%l%c%cef%eflx_fgr, set_spec=spval, default='inactive')

    call hist_addfld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_pft=clm3%g%l%c%p%pmf%taux)

    ! Rename of TAUX for Urban intercomparision project (when U=V)
    call hist_addfld1d (fname='Qtau', units='kg/m/s^2',  &
         avgflag='A', long_name='momentum flux', &
         ptr_pft=clm3%g%l%c%p%pmf%taux, default='inactive')

    call hist_addfld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_pft=clm3%g%l%c%p%pmf%tauy)

    ! Vegetation phenology

    call hist_addfld1d (fname='ELAI', units='m^2/m^2', &
          avgflag='A', long_name='exposed one-sided leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%elai)

    call hist_addfld1d (fname='ESAI', units='m^2/m^2', &
          avgflag='A', long_name='exposed one-sided stem area index', &
         ptr_pft=clm3%g%l%c%p%pps%esai)

    call hist_addfld1d (fname='LAISUN', units='none', &
         avgflag='A', long_name='sunlit projected leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%laisun, set_urb=0._r8)

    call hist_addfld1d (fname='LAISHA', units='none', &
         avgflag='A', long_name='shaded projected leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%laisha, set_urb=0._r8)

    call hist_addfld1d (fname='TLAI', units='none', &
         avgflag='A', long_name='total projected leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%tlai)

    call hist_addfld1d (fname='TSAI', units='none', &
         avgflag='A', long_name='total projected stem area index', &
         ptr_pft=clm3%g%l%c%p%pps%tsai)

    ! Canopy physiology

    call hist_addfld1d (fname='RSSUN', units='s/m',  &
         avgflag='M', long_name='sunlit leaf stomatal resistance', &
         ptr_pft=clm3%g%l%c%p%pps%rssun, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='RSSHA', units='s/m',  &
         avgflag='M', long_name='shaded leaf stomatal resistance', &
         ptr_pft=clm3%g%l%c%p%pps%rssha, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='BTRAN', units='unitless',  &
         avgflag='A', long_name='transpiration beta factor', &
         ptr_pft=clm3%g%l%c%p%pps%btran, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='FPSN', units='umol/m2s',  &
         avgflag='A', long_name='photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%fpsn, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='FPSN_WC', units='umol/m2s',  &
         avgflag='A', long_name='Rubisco-limited photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%fpsn_wc, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='FPSN_WJ', units='umol/m2s',  &
         avgflag='A', long_name='RuBP-limited photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%fpsn_wj, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='FPSN_WP', units='umol/m2s',  &
         avgflag='A', long_name='Product-limited photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%fpsn_wp, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='DSTFLXT', units='kg/m2/s',  &
         avgflag='A', long_name='total surface dust emission', &
         ptr_pft=clm3%g%l%c%p%pdf%flx_mss_vrt_dst_tot, set_lake=0._r8, set_urb=0._r8)
    call hist_addfld1d (fname='DPVLTRB1', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 1', &
         ptr_pft=clm3%g%l%c%p%pdf%vlc_trb_1, default='inactive')
    call hist_addfld1d (fname='DPVLTRB2', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 2', &
         ptr_pft=clm3%g%l%c%p%pdf%vlc_trb_2, default='inactive')
    call hist_addfld1d (fname='DPVLTRB3', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 3', &
         ptr_pft=clm3%g%l%c%p%pdf%vlc_trb_3, default='inactive')
    call hist_addfld1d (fname='DPVLTRB4', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 4', &
         ptr_pft=clm3%g%l%c%p%pdf%vlc_trb_4, default='inactive')

    ! for MEGAN emissions diagnositics
    if (shr_megan_megcomps_n>0) then
       
       ! loop over megan compounds
       meg_cmp => shr_megan_linkedlist
       do while(associated(meg_cmp))
          imeg = meg_cmp%index

          call hist_addfld1d ( fname='MEG_'//trim(meg_cmp%name), units='kg/m2/sec',  &
               avgflag='A', long_name='MEGAN flux', &
               ptr_pft=clm3%g%l%c%p%pvf%meg(imeg)%flux_out, set_lake=0._r8, set_urb=0._r8 )

          meg_cmp => meg_cmp%next_megcomp
       enddo
       
       call hist_addfld1d (fname='VOCFLXT', units='moles/m2/sec',  &
            avgflag='A', long_name='total VOC flux into atmosphere', &
            ptr_pft=clm3%g%l%c%p%pvf%vocflx_tot, set_lake=0._r8, set_urb=0._r8)

       call hist_addfld1d (fname='GAMMA', units='non',  &
            avgflag='A', long_name='total gamma for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%gamma_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAL', units='non',  &
            avgflag='A', long_name='gamma L for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%gammaL_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAT', units='non',  &
            avgflag='A', long_name='gamma T for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%gammaT_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAP', units='non',  &
            avgflag='A', long_name='gamma P for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%gammaP_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAA', units='non',  &
            avgflag='A', long_name='gamma A for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%gammaA_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAS', units='non',  &
            avgflag='A', long_name='gamma S for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%gammaS_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAC', units='non',  &
            avgflag='A', long_name='gamma C for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%gammaC_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='EOPT', units='non',  &
            avgflag='A', long_name='Eopt coefficient for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%Eopt_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='TOPT', units='non',  &
            avgflag='A', long_name='topt coefficient for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%topt_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='ALPHA', units='non',  &
            avgflag='A', long_name='alpha coefficient for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%alpha_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='CP', units='non',  &
            avgflag='A', long_name='cp coefficient for VOC calc', &
            ptr_pft=clm3%g%l%c%p%pvf%cp_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR', &
            ptr_pft=clm3%g%l%c%p%pvf%paru_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR24_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR (24 hrs)', &
            ptr_pft=clm3%g%l%c%p%pvf%par24u_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR240_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR (240 hrs)', &
            ptr_pft=clm3%g%l%c%p%pvf%par240u_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR', &
            ptr_pft=clm3%g%l%c%p%pvf%para_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR24_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR (24 hrs)', &
            ptr_pft=clm3%g%l%c%p%pvf%par24a_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR240_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR (240 hrs)', &
            ptr_pft=clm3%g%l%c%p%pvf%par240a_out, set_lake=0._r8, default='inactive')

    endif

    call hist_addfld1d (fname='FSUN24', units='K',  &
         avgflag='A', long_name='fraction sunlit (last 24hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsun24, default='inactive')

    call hist_addfld1d (fname='FSUN240', units='K',  &
         avgflag='A', long_name='fraction sunlit (last 240hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsun240, default='inactive')

    call hist_addfld1d (fname='FSI24', units='K',  &
         avgflag='A', long_name='indirect radiation (last 24hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsi24, default='inactive')

    call hist_addfld1d (fname='FSI240', units='K',  &
         avgflag='A', long_name='indirect radiation (last 240hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsi240, default='inactive')

    call hist_addfld1d (fname='FSD24', units='K',  &
         avgflag='A', long_name='direct radiation (last 24hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsd24, default='inactive')

    call hist_addfld1d (fname='FSD240', units='K',  &
         avgflag='A', long_name='direct radiation (last 240hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsd240, default='inactive')

    ! Hydrology

    call hist_addfld1d (fname='SoilAlpha',  units='unitless',  &
         avgflag='A', long_name='factor limiting ground evap', &
         ptr_col=clm3%g%l%c%cws%soilalpha, set_urb=spval)

    call hist_addfld1d (fname='SoilAlpha_U',  units='unitless',  &
         avgflag='A', long_name='urban factor limiting ground evap', &
         ptr_col=clm3%g%l%c%cws%soilalpha_u, set_nourb=spval)

    call hist_addfld1d (fname='FCOV',  units='unitless',  &
         avgflag='A', long_name='fractional impermeable area', &
         ptr_col=clm3%g%l%c%cws%fcov, l2g_scale_type='veg')
    call hist_addfld1d (fname='FSAT',  units='unitless',  &
         avgflag='A', long_name='fractional area with water table at surface', &
         ptr_col=clm3%g%l%c%cws%fsat, l2g_scale_type='veg')
    call hist_addfld1d (fname='ZWT',  units='m',  &
         avgflag='A', long_name='water table depth (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%zwt, l2g_scale_type='veg')

    call hist_addfld1d (fname='INT_SNOW',  units='mm',  &
         avgflag='A', long_name='accumulated swe (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%int_snow, l2g_scale_type='veg')
    call hist_addfld1d (fname='FROST_TABLE',  units='m',  &
         avgflag='A', long_name='frost table depth (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%frost_table, l2g_scale_type='veg')
    call hist_addfld1d (fname='ZWT_PERCH',  units='m',  &
         avgflag='A', long_name='perched water table depth (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%zwt_perched, l2g_scale_type='veg')
    call hist_addfld1d (fname='QDRAI_PERCH',  units='mm/s',  &
         avgflag='A', long_name='perched wt drainage', &
         ptr_col=clm3%g%l%c%cwf%qflx_drain_perched, c2l_scale_type='urbanf')
    call hist_addfld1d (fname='QDRAI_XS',  units='mm/s',  &
         avgflag='A', long_name='saturation excess drainage', &
         ptr_col=clm3%g%l%c%cwf%qflx_rsub_sat, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='WA',  units='mm',  &
         avgflag='A', long_name='water in the unconfined aquifer (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%wa, l2g_scale_type='veg')

    call hist_addfld1d (fname='QCHARGE',  units='mm/s',  &
         avgflag='A', long_name='aquifer recharge rate (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%qcharge, l2g_scale_type='veg')

#ifdef LCH4
    if (hist_wrtch4diag) then
       active = "active"
    else
       active = "inactive"
    end if
#else
    active = "inactive"
#endif
    call hist_addfld2d (fname='SMP',  units='mm', type2d='levgrnd',  &
         avgflag='A', long_name='soil matric potential (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%smp_l, set_spec=spval, l2g_scale_type='veg', default=active)

    call hist_addfld2d (fname='HK',  units='mm/s', type2d='levgrnd',  &
         avgflag='A', long_name='hydraulic conductivity (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%hk_l, set_spec=spval, l2g_scale_type='veg', default='inactive')

    call hist_addfld1d (fname='H2OSNO',  units='mm',  &
         avgflag='A', long_name='snow depth (liquid water)', &
         ptr_col=clm3%g%l%c%cws%h2osno, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='ERRH2OSNO',  units='mm',  &
         avgflag='A', long_name='imbalance in snow depth (liquid water)', &
         ptr_col=clm3%g%l%c%cws%errh2osno, c2l_scale_type='urbanf')

    ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
    ! any given time step but only if there is at least one snow layer (for all landunits 
    ! except lakes).  h2osno also includes snow that is part of the soil column (an 
    ! initial snow layer is only created if h2osno > 10mm). Also note that monthly average
    ! files of snow_sources and snow sinks must be weighted by number of days in the month to 
    ! diagnose, for example, an annual value of the change in h2osno. 

    call hist_addfld1d (fname='SNOW_SOURCES',  units='mm/s',  &
         avgflag='A', long_name='snow sources (liquid water)', &
         ptr_col=clm3%g%l%c%cws%snow_sources, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='SNOW_SINKS',  units='mm/s',  &
         avgflag='A', long_name='snow sinks (liquid water)', &
         ptr_col=clm3%g%l%c%cws%snow_sinks, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='H2OCAN', units='mm',  &
         avgflag='A', long_name='intercepted water', &
         ptr_pft=clm3%g%l%c%p%pws%h2ocan, set_lake=0._r8)

    call hist_addfld2d (fname='H2OSOI',  units='mm3/mm3', type2d='levgrnd', &
         avgflag='A', long_name='volumetric soil water (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%h2osoi_vol, l2g_scale_type='veg')

    call hist_addfld2d (fname='SOILLIQ',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='soil liquid water (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%h2osoi_liq, l2g_scale_type='veg')

    call hist_addfld2d (fname='SOILICE',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='soil ice (vegetated landunits only)', &
         ptr_col=clm3%g%l%c%cws%h2osoi_ice, l2g_scale_type='veg')

    call hist_addfld1d (fname='SOILWATER_10CM',  units='kg/m2', &
         avgflag='A', long_name='soil liquid water + ice in top 10cm of soil (veg landunits only)', &
         ptr_col=clm3%g%l%c%cws%h2osoi_liqice_10cm, set_urb=spval, l2g_scale_type='veg', &
         set_lake=spval)

    call hist_addfld1d (fname='SNOWLIQ',  units='kg/m2',  &
         avgflag='A', long_name='snow liquid water', &
         ptr_col=clm3%g%l%c%cws%snowliq)

    call hist_addfld1d (fname='SNOWICE',  units='kg/m2',  &
         avgflag='A', long_name='snow ice', &
         ptr_col=clm3%g%l%c%cws%snowice)

    call hist_addfld1d (fname='QTOPSOIL',  units='mm/s',  &
         avgflag='A', long_name='water input to surface', &
         ptr_col=clm3%g%l%c%cwf%qflx_top_soil, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='QINFL',  units='mm/s',  &
         avgflag='A', long_name='infiltration', &
         ptr_col=clm3%g%l%c%cwf%qflx_infl, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QOVER',  units='mm/s',  &
         avgflag='A', long_name='surface runoff', &
         ptr_col=clm3%g%l%c%cwf%qflx_surf, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRGWL',  units='mm/s',  &
         avgflag='A', long_name='surface runoff at glaciers (liquid only), wetlands, lakes', &
         ptr_col=clm3%g%l%c%cwf%qflx_qrgwl, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNWCPLIQ', units='mm H2O/s', &
         avgflag='A', long_name='excess rainfall due to snow capping', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_snwcp_liq, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='QSNWCPICE_NODYNLNDUSE', units='mm H2O/s', &
         avgflag='A', &
         long_name='excess snowfall due to snow capping not including correction for land use change', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_snwcp_ice, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNWCPICE',  units='mm/s',  &
         avgflag='A', long_name='excess snowfall due to snow capping', &
         ptr_lnd=clm3%g%gwf%qflx_snwcp_iceg)

    call hist_addfld1d (fname='QDRAI',  units='mm/s',  &
         avgflag='A', long_name='sub-surface drainage', &
         ptr_col=clm3%g%l%c%cwf%qflx_drain, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF_NODYNLNDUSE',  units='mm/s',  &
         avgflag='A', &
         long_name='total liquid runoff (does not include QSNWCPICE) not including correction for land use change', &
         ptr_col=clm3%g%l%c%cwf%qflx_runoff, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF',  units='mm/s',  &
         avgflag='A', long_name='total liquid runoff (does not include QSNWCPICE)', &
         ptr_lnd=clm3%g%gwf%qflx_runoffg)

    call hist_addfld1d (fname='GC_LIQ1',  units='mm',  &
         avgflag='A', long_name='initial gridcell total liq content', &
         ptr_lnd=clm3%g%gws%gc_liq1)

    call hist_addfld1d (fname='GC_LIQ2',  units='mm',  &  
         avgflag='A', long_name='post landuse change gridcell total liq content', &              
         ptr_lnd=clm3%g%gws%gc_liq2, default='inactive')     

    call hist_addfld1d (fname='QFLX_LIQ_DYNBAL',  units='mm/s',  &  
         avgflag='A', long_name='liq dynamic land cover change conversion runoff flux', &              
         ptr_lnd=clm3%g%gwf%qflx_liq_dynbal)     

    call hist_addfld1d (fname='GC_ICE1',  units='mm',  &  
         avgflag='A', long_name='initial gridcell total ice content', &              
         ptr_lnd=clm3%g%gws%gc_ice1)     

    call hist_addfld1d (fname='GC_ICE2',  units='mm',  &  
         avgflag='A', long_name='post land cover change total ice content', &              
         ptr_lnd=clm3%g%gws%gc_ice2, default='inactive')

    call hist_addfld1d (fname='QFLX_ICE_DYNBAL',  units='mm/s',  &
         avgflag='A', long_name='ice dynamic land cover change conversion runoff flux', &                                   
         ptr_lnd=clm3%g%gwf%qflx_ice_dynbal)

    call hist_addfld1d (fname='QRUNOFF_U', units='mm/s',  &
         avgflag='A', long_name='Urban total runoff', &
         ptr_col=clm3%g%l%c%cwf%qflx_runoff_u, set_nourb=spval, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF_R', units='mm/s',  &
         avgflag='A', long_name='Rural total runoff', &
         ptr_col=clm3%g%l%c%cwf%qflx_runoff_r, set_spec=spval)

    call hist_addfld1d (fname='QINTR', units='mm/s',  &
         avgflag='A', long_name='interception', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_prec_intr, set_lake=0._r8)

    call hist_addfld1d (fname='QDRIP', units='mm/s',  &
         avgflag='A', long_name='throughfall', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_prec_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNOMELT',  units='mm/s',  &
         avgflag='A', long_name='snow melt', &
         ptr_col=clm3%g%l%c%cwf%qflx_snow_melt, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNOFRZ', units='kg/m2/s', &
         avgflag='A', long_name='column-integrated snow freezing rate', &
         ptr_col=clm3%g%l%c%cwf%qflx_snofrz_col, default='inactive', &
         set_lake=spval, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSOIL', units='mm/s',  &
         avgflag='A', long_name= &
         'Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_soi, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QVEGE', units='mm/s',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_can, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QVEGT', units='mm/s',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_tran_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QIRRIG', units='mm/s', &
         avgflag='A', long_name='water added through irrigation', &
         ptr_col=clm3%g%l%c%cwf%qflx_irrig, set_lake=0._r8)

    if (create_glacier_mec_landunit) then

       call hist_addfld1d (fname='QICE',  units='mm/s',  &
            avgflag='A', long_name='ice growth/melt', &
            ptr_col=clm3%g%l%c%cwf%qflx_glcice, set_noglcmec=spval)

       call hist_addfld1d (fname='QICE_FRZ',  units='mm/s',  &
            avgflag='A', long_name='ice growth', &
            ptr_col=clm3%g%l%c%cwf%qflx_glcice_frz, set_noglcmec=spval)

       call hist_addfld1d (fname='QICE_MELT',  units='mm/s',  &
            avgflag='A', long_name='ice melt', &
            ptr_col=clm3%g%l%c%cwf%qflx_glcice_melt, set_noglcmec=spval)

       call hist_addfld1d (fname='gris_mask',  units='unitless',  &
            avgflag='A', long_name='Greenland mask', &
            ptr_gcell=clm3%g%gris_mask)

       call hist_addfld1d (fname='gris_area',  units='km^2',  &
            avgflag='A', long_name='Greenland ice area', &
            ptr_gcell=clm3%g%gris_area)

       call hist_addfld1d (fname='aais_mask',  units='unitless',  &
            avgflag='A', long_name='Antarctic mask', &
            ptr_gcell=clm3%g%aais_mask)

       call hist_addfld1d (fname='aais_area',  units='km^2',  &
            avgflag='A', long_name='Antarctic ice area', &
            ptr_gcell=clm3%g%aais_area)

   endif

    ! Water and energy balance checks

    call hist_addfld1d (fname='ERRSOI',  units='W/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=clm3%g%l%c%cebal%errsoi)

    call hist_addfld1d (fname='ERRSEB',  units='W/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_pft=clm3%g%l%c%p%pebal%errseb)

    call hist_addfld1d (fname='ERRSOL',  units='W/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_pft=clm3%g%l%c%p%pebal%errsol, set_urb=spval)

    call hist_addfld1d (fname='ERRH2O', units='mm',  &
         avgflag='A', long_name='total water conservation error', &
         ptr_col=clm3%g%l%c%cwbal%errh2o)

    ! Atmospheric forcing

    call hist_addfld1d (fname='RAIN', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain', &
         ptr_lnd=clm_a2l%forc_rain)

    call hist_addfld1d (fname='SNOW', units='mm/s',  &
         avgflag='A', long_name='atmospheric snow', &
         ptr_lnd=clm_a2l%forc_snow)

    call hist_addfld1d (fname='TBOT', units='K',  &
         avgflag='A', long_name='atmospheric air temperature', &
         ptr_lnd=clm_a2l%forc_t)

    call hist_addfld1d (fname='THBOT', units='K',  &
         avgflag='A', long_name='atmospheric air potential temperature', &
         ptr_lnd=clm_a2l%forc_th)

    call hist_addfld1d (fname='WIND', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_lnd=clm_a2l%forc_wind)

    ! Rename of WIND for Urban intercomparision project
    call hist_addfld1d (fname='Wind', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_gcell=clm_a2l%forc_wind, default = 'inactive')

    call hist_addfld1d (fname='Tair', units='K',  &
         avgflag='A', long_name='atmospheric air temperature', &
         ptr_gcell=clm_a2l%forc_t, default='inactive')

    call hist_addfld1d (fname='PSurf', units='Pa',  &
         avgflag='A', long_name='surface pressure', &
         ptr_gcell=clm_a2l%forc_pbot, default='inactive')

    call hist_addfld1d (fname='Rainf', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain', &
         ptr_gcell=clm_a2l%forc_rain, default='inactive')

    call hist_addfld1d (fname='SWdown', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_gcell=clm_a2l%forc_solar, default='inactive')

    call hist_addfld1d (fname='LWdown', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation', &
         ptr_gcell=clm_a2l%forc_lwrad, default='inactive')

    call hist_addfld1d (fname='RH', units='%',  &
         avgflag='A', long_name='atmospheric relative humidity', &
         ptr_gcell=clm_a2l%forc_rh, default='inactive')

    call hist_addfld1d (fname='QBOT', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity', &
         ptr_lnd=clm_a2l%forc_q)

    ! Rename of QBOT for Urban intercomparision project
    call hist_addfld1d (fname='Qair', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity', &
         ptr_lnd=clm_a2l%forc_q, default='inactive')

    call hist_addfld1d (fname='ZBOT', units='m',  &
         avgflag='A', long_name='atmospheric reference height', &
         ptr_lnd=clm_a2l%forc_hgt)

    call hist_addfld1d (fname='FLDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation', &
         ptr_lnd=clm_a2l%forc_lwrad)

    call hist_addfld1d (fname='FSDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_lnd=clm_a2l%forc_solar)

    call hist_addfld1d (fname='PCO2', units='Pa',  &
         avgflag='A', long_name='atmospheric partial pressure of CO2', &
         ptr_lnd=clm_a2l%forc_pco2)

    call hist_addfld1d (fname='PBOT', units='Pa',  &
         avgflag='A', long_name='atmospheric pressure', &
         ptr_lnd=clm_a2l%forc_pbot)

#if (defined CNDV) || (defined CROP)
    active = "active"
#else
    active = "inactive"
#endif
    call hist_addfld1d (fname='T10', units='K',  &
         avgflag='A', long_name='10-day running mean of 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t10, default=active)


#if (defined CNDV)
    call hist_addfld1d (fname='TDA', units='K',  &
         avgflag='A', long_name='daily average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pdgvs%t_mo)

    call hist_addfld1d (fname='AGDD', units='K',  &
         avgflag='A', long_name='growing degree-days base 5C', &
         ptr_pft=clm3%g%l%c%p%pdgvs%agdd)
#endif

#if (defined LCH4)

    call hist_addfld1d (fname='PCH4', units='Pa',  &
         avgflag='A', long_name='atmospheric partial pressure of CH4', &
         ptr_lnd=clm_a2l%forc_pch4)

    call hist_addfld1d (fname='WTGQ', units='m/s',  &
         avgflag='A', long_name='surface tracer conductance', &
         ptr_col=clm3%g%l%c%cps%pps_a%grnd_ch4_cond)

    call hist_addfld1d (fname='FINUNDATED', units='unitless', &
         avgflag='A', long_name='fractional inundated area of vegetated columns', &
         ptr_col=clm3%g%l%c%cws%finundated)

    call hist_addfld1d (fname='CH4_SURF_DIFF_SAT', units='mol/m2/s',  &
         avgflag='A', long_name='diffusive surface CH4 flux for inundated / lake area; (+ to atm)', &
         ptr_col=clm3%g%l%c%cch4%ch4_surf_diff_sat)

    call hist_addfld1d (fname='CH4_SURF_DIFF_UNSAT', units='mol/m2/s',  &
         avgflag='A', long_name='diffusive surface CH4 flux for non-inundated area; (+ to atm)', &
         ptr_col=clm3%g%l%c%cch4%ch4_surf_diff_unsat)

    call hist_addfld1d (fname='CH4_EBUL_TOTAL_SAT', units='mol/m2/s',  &
         avgflag='A', long_name='ebullition surface CH4 flux; (+ to atm)', &
         ptr_col=clm3%g%l%c%cch4%ch4_ebul_total_sat, default='inactive')

    call hist_addfld1d (fname='CH4_EBUL_TOTAL_UNSAT', units='mol/m2/s',  &
         avgflag='A', long_name='ebullition surface CH4 flux; (+ to atm)', &
         ptr_col=clm3%g%l%c%cch4%ch4_ebul_total_unsat, default='inactive')

    call hist_addfld1d (fname='CH4_SURF_EBUL_SAT', units='mol/m2/s',  &
         avgflag='A', long_name='ebullition surface CH4 flux for inundated / lake area; (+ to atm)', &
         ptr_col=clm3%g%l%c%cch4%ch4_surf_ebul_sat)

    call hist_addfld1d (fname='CH4_SURF_EBUL_UNSAT', units='mol/m2/s',  &
         avgflag='A', long_name='ebullition surface CH4 flux for non-inundated area; (+ to atm)', &
         ptr_col=clm3%g%l%c%cch4%ch4_surf_ebul_unsat)

    call hist_addfld1d (fname='CH4_SURF_AERE_SAT', units='mol/m2/s',  &
         avgflag='A', long_name='aerenchyma surface CH4 flux for inundated area; (+ to atm)', &
         ptr_col=clm3%g%l%c%cch4%ch4_surf_aere_sat)

    call hist_addfld1d (fname='CH4_SURF_AERE_UNSAT', units='mol/m2/s',  &
         avgflag='A', long_name='aerenchyma surface CH4 flux for non-inundated area; (+ to atm)', &
         ptr_col=clm3%g%l%c%cch4%ch4_surf_aere_unsat)

    call hist_addfld1d (fname='TOTCOLCH4', units='gC/m2',  &
         avgflag='A', long_name='total belowground CH4, (0 for non-lake special landunits)', &
         ptr_col=clm3%g%l%c%cch4%totcolch4)

    call hist_addfld2d (fname='CONC_CH4_SAT', units='mol/m3', type2d='levgrnd', &
         avgflag='A', long_name='CH4 soil Concentration for inundated / lake area', &
         ptr_col=clm3%g%l%c%cch4%conc_ch4_sat)
                                                                       
    call hist_addfld2d (fname='CONC_CH4_UNSAT', units='mol/m3', type2d='levgrnd', &
         avgflag='A', long_name='CH4 soil Concentration for non-inundated area', &
         ptr_col=clm3%g%l%c%cch4%conc_ch4_unsat)
                                                                       
! Begin if hist_wrtch4diag block
    if (hist_wrtch4diag) then
       call hist_addfld2d (fname='CH4_PROD_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil production for inundated / lake area', &
            ptr_col=clm3%g%l%c%cch4%ch4_prod_depth_sat)
   
       call hist_addfld2d (fname='CH4_PROD_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil production for non-inundated area', &
            ptr_col=clm3%g%l%c%cch4%ch4_prod_depth_unsat)
   
       call hist_addfld2d (fname='CH4_OXID_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil oxidation for inundated / lake area', &
            ptr_col=clm3%g%l%c%cch4%ch4_oxid_depth_sat)
   
       call hist_addfld2d (fname='CH4_OXID_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil oxidation for non-inundated area', &
            ptr_col=clm3%g%l%c%cch4%ch4_oxid_depth_unsat)
   
       call hist_addfld2d (fname='CH4_AERE_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil aerenchyma loss for inundated / lake area (including transpiration flux if activated)', &
            ptr_col=clm3%g%l%c%cch4%ch4_aere_depth_sat)
   
       call hist_addfld2d (fname='CH4_AERE_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil aerenchyma loss for non-inundated area (including transpiration flux if activated)', &
            ptr_col=clm3%g%l%c%cch4%ch4_aere_depth_unsat)
   
       call hist_addfld2d (fname='O2_AERE_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='O2 aerenchyma diffusion into soil for inundated / lake area', &
            ptr_col=clm3%g%l%c%cch4%o2_aere_depth_sat)
   
       call hist_addfld2d (fname='O2_AERE_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='O2 aerenchyma diffusion into soil for non-inundated area', &
            ptr_col=clm3%g%l%c%cch4%o2_aere_depth_unsat)
   
       call hist_addfld2d (fname='O2_DECOMP_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='O2 consumption from HR and AR for inundated / lake area', &
            ptr_col=clm3%g%l%c%cch4%o2_decomp_depth_sat)
   
       call hist_addfld2d (fname='O2_DECOMP_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='O2 consumption from HR and AR for non-inundated area', &
            ptr_col=clm3%g%l%c%cch4%o2_decomp_depth_unsat)

       call hist_addfld2d (fname='CH4_TRAN_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil loss from transpiration for inundated / lake area', &
            ptr_col=clm3%g%l%c%cch4%ch4_tran_depth_sat)
   
       call hist_addfld2d (fname='CH4_TRAN_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil loss from transpiration for non-inundated area', &
            ptr_col=clm3%g%l%c%cch4%ch4_tran_depth_unsat)
   
       call hist_addfld2d (fname='CH4_EBUL_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil ebullition for inundated / lake area', &
            ptr_col=clm3%g%l%c%cch4%ch4_ebul_depth_sat)
   
       call hist_addfld2d (fname='CH4_EBUL_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil ebullition for non-inundated area', &
            ptr_col=clm3%g%l%c%cch4%ch4_ebul_depth_unsat)
   
       call hist_addfld2d (fname='O2STRESS_UNSAT', units='unitless', type2d='levgrnd',  &
            avgflag='A', long_name='Ratio of oxygen available to demanded for inundated / lake area', &
            ptr_col=clm3%g%l%c%cch4%o2stress_unsat)
   
       call hist_addfld2d (fname='O2STRESS_SAT', units='unitless', type2d='levgrnd',  &
            avgflag='A', long_name='Ratio of oxygen available to demanded for non-inundated area', &
            ptr_col=clm3%g%l%c%cch4%o2stress_sat)
   
       call hist_addfld2d (fname='CH4STRESS_UNSAT', units='unitless', type2d='levgrnd',  &
            avgflag='A', long_name='Ratio of methane available to total potential sink for inundated / lake area', &
            ptr_col=clm3%g%l%c%cch4%ch4stress_unsat)
   
       call hist_addfld2d (fname='CH4STRESS_SAT', units='unitless', type2d='levgrnd',  &
            avgflag='A', long_name='Ratio of methane available to total potential sink for non-inundated area', &
            ptr_col=clm3%g%l%c%cch4%ch4stress_sat)
   
       ! Lake diagnostics
       if (allowlakeprod) then
          call hist_addfld2d (fname='CH4_PROD_DEPTH_LAKE', units='mol/m3/s', type2d='levgrnd', &
               avgflag='A', long_name='CH4 production in each soil layer, lake col. only', &
               ptr_col=clm3%g%l%c%cch4%ch4_prod_depth_sat)
      
          call hist_addfld2d (fname='CONC_CH4_LAKE', units='mol/m3', type2d='levgrnd', &
               avgflag='A', long_name='CH4 Concentration each soil layer, lake col. only', &
               ptr_col=clm3%g%l%c%cch4%conc_ch4_sat)
      
          call hist_addfld2d (fname='CONC_O2_LAKE', units='mol/m3', type2d='levgrnd', &
               avgflag='A', long_name='O2 Concentration each soil layer, lake col. only', &
               ptr_col=clm3%g%l%c%cch4%conc_o2_sat)
      
          call hist_addfld1d (fname='CH4_SURF_DIFF_LAKE', units='mol/m2/s',  &
               avgflag='A', long_name='diffusive surface CH4 flux, lake col. only (+ to atm)', &
               ptr_col=clm3%g%l%c%cch4%ch4_surf_diff_sat)
      
          call hist_addfld1d (fname='CH4_SURF_EBUL_LAKE', units='mol/m2/s',  &
               avgflag='A', long_name='ebullition surface CH4 flux, lake col. only (+ to atm)', &
               ptr_col=clm3%g%l%c%cch4%ch4_surf_ebul_sat)
      
          call hist_addfld2d (fname='CH4_OXID_DEPTH_LAKE', units='mol/m2/s', type2d='levgrnd',  &
               avgflag='A', long_name='CH4 oxidation in each soil layer, lake col. only', &
               ptr_col=clm3%g%l%c%cch4%ch4_oxid_depth_sat)
       end if

       call hist_addfld2d (fname='LAYER_SAT_LAG', units='unitless', type2d='levgrnd',  &
            avgflag='A', long_name='lagged saturation status of layer in unsat. zone', &
            ptr_col=clm3%g%l%c%cch4%layer_sat_lag)

#ifdef CN

       call hist_addfld_decomp (fname='FPHR'//trim(vr_suffix), units='unitless', type2d='levdcmp', &
            avgflag='A', long_name='fraction of potential HR due to N limitation', &
            ptr_col=clm3%g%l%c%cch4%fphr)

       call hist_addfld1d (fname='ANNAVG_FINRW', units='unitless',  &
            avgflag='A', long_name='annual average respiration-weighted FINUNDATED', &
            ptr_col=clm3%g%l%c%cch4%annavg_finrw)

       call hist_addfld1d (fname='SIF', units='unitless',  &
            avgflag='A', long_name='seasonal inundation factor calculated for sat. CH4 prod. (non-lake)', &
            ptr_col=clm3%g%l%c%cch4%sif)
#endif
! CN

    end if
! End hist_wrtch4diag block if.

    call hist_addfld2d (fname='CONC_O2_SAT', units='mol/m3', type2d='levgrnd', &
         avgflag='A', long_name='O2 soil Concentration for inundated / lake area', &
         ptr_col=clm3%g%l%c%cch4%conc_o2_sat)
                                                                       
    call hist_addfld2d (fname='CONC_O2_UNSAT', units='mol/m3', type2d='levgrnd', &
         avgflag='A', long_name='O2 soil Concentration for non-inundated area', &
         ptr_col=clm3%g%l%c%cch4%conc_o2_unsat)

    call hist_addfld1d (fname='FCH4', units='kgC/m2/s', &
         avgflag='A', long_name='Gridcell surface CH4 flux to atmosphere (+ to atm)', &
         ptr_lnd=clm_l2a%flux_ch4)

    call hist_addfld1d (fname='FCH4TOCO2', units='gC/m2/s', &
         avgflag='A', long_name='Gridcell oxidation of CH4 to CO2', &
         ptr_lnd=clm3%g%gch4%ch4co2f)

    call hist_addfld1d (fname='CH4PROD', units='gC/m2/s', &
         avgflag='A', long_name='Gridcell total production of CH4', &
         ptr_lnd=clm3%g%gch4%ch4prodg)

    call hist_addfld1d (fname='NEM', units='gC/m2/s', &
         avgflag='A', long_name='Gridcell net adjustment to NEE passed to atm. for methane production', &
         ptr_lnd=clm3%g%gch4%nem)

    call hist_addfld1d (fname='FCH4_DFSAT', units='kgC/m2/s',  &
         avgflag='A', long_name='CH4 additional flux due to changing fsat, vegetated landunits only', &
         ptr_col=clm3%g%l%c%cch4%ch4_dfsat_flux)

    call hist_addfld1d (fname='ZWT_CH4_UNSAT', units='m',  &
         avgflag='A', long_name='depth of water table for methane production used in non-inundated area', &
         ptr_col=clm3%g%l%c%cch4%zwt_ch4_unsat)

    call hist_addfld1d (fname='QOVER_LAG', units='mm/s',  &
         avgflag='A', long_name='time-lagged surface runoff for soil columns', &
         ptr_col=clm3%g%l%c%cch4%qflx_surf_lag)

    call hist_addfld1d (fname='FINUNDATED_LAG', units='unitless',  &
         avgflag='A', long_name='time-lagged inundated fraction of vegetated columns', &
         ptr_col=clm3%g%l%c%cch4%finundated_lag)

    if (allowlakeprod) then
       call hist_addfld2d (fname='LAKE_SOILC', units='gC/m3', type2d='levgrnd', &
            avgflag='A', long_name='Soil carbon under lakes', &
            ptr_col=clm3%g%l%c%cch4%lake_soilc)
    end if

#endif
! CH4

#if (defined CN)
    call hist_addfld2d (fname='SOILPSI', units='MPa', type2d='levgrnd', &
         avgflag='A', long_name='soil water potential in each soil layer', &
         ptr_col=clm3%g%l%c%cps%soilpsi)
#endif

#if (defined CN)
    ! add history fields for all CN variables, always set as default='inactive'

    if ( crop_prog )then

       call hist_addfld1d (fname='A5TMIN', units='K',  &
            avgflag='A', long_name='5-day running mean of min 2-m temperature', &
            ptr_pft=clm3%g%l%c%p%pes%a5tmin, default='inactive')

       call hist_addfld1d (fname='A10TMIN', units='K',  &
            avgflag='A', long_name='10-day running mean of min 2-m temperature', &
            ptr_pft=clm3%g%l%c%p%pes%a10tmin, default='inactive')

    end if
    
    !-------------------------------
    ! C state variables - native to PFT 
    !-------------------------------
    ! add history fields for all CLAMP CN variables

    if (crop_prog) then
       call hist_addfld1d (fname='GRAINC', units='gC/m^2', &
            avgflag='A', long_name='grain C', &
            ptr_pft=clm3%g%l%c%p%pcs%grainc)
    end if

    call hist_addfld1d (fname='WOODC', units='gC/m^2', &
             avgflag='A', long_name='wood C', &
             ptr_pft=clm3%g%l%c%p%pcs%woodc)
    
    call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
         avgflag='A', long_name='leaf C', &
         ptr_pft=clm3%g%l%c%p%pcs%leafc)

    call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='leaf C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%leafc_storage, default='inactive')

    call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
         avgflag='A', long_name='leaf C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%leafc_xfer, default='inactive')

    call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
         avgflag='A', long_name='fine root C', &
         ptr_pft=clm3%g%l%c%p%pcs%frootc)

    call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='fine root C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%frootc_storage, default='inactive')

    call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
         avgflag='A', long_name='fine root C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%frootc_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
         avgflag='A', long_name='live stem C', &
         ptr_pft=clm3%g%l%c%p%pcs%livestemc)

    call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='live stem C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%livestemc_storage, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
         avgflag='A', long_name='live stem C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%livestemc_xfer, default='inactive')

    call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
         avgflag='A', long_name='dead stem C', &
         ptr_pft=clm3%g%l%c%p%pcs%deadstemc)

    call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='dead stem C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%deadstemc_storage, default='inactive')

    call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
         avgflag='A', long_name='dead stem C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%deadstemc_xfer, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
         avgflag='A', long_name='live coarse root C', &
         ptr_pft=clm3%g%l%c%p%pcs%livecrootc)

    call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='live coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%livecrootc_storage, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
         avgflag='A', long_name='live coarse root C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%livecrootc_xfer, default='inactive')

    call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
         avgflag='A', long_name='dead coarse root C', &
         ptr_pft=clm3%g%l%c%p%pcs%deadcrootc)

    call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='dead coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%deadcrootc_storage,  default='inactive')

    call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
         avgflag='A', long_name='dead coarse root C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%deadcrootc_xfer, default='inactive')

    call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='growth respiration storage', &
         ptr_pft=clm3%g%l%c%p%pcs%gresp_storage, default='inactive')

    call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
         avgflag='A', long_name='growth respiration transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%gresp_xfer, default='inactive')

    call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
         avgflag='A', long_name='temporary photosynthate C pool', &
         ptr_pft=clm3%g%l%c%p%pcs%cpool)

    call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
         avgflag='A', long_name='temporary photosynthate C pool', &
         ptr_pft=clm3%g%l%c%p%pcs%xsmrpool)

    call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
         avgflag='A', long_name='pft-level sink for C truncation', &
         ptr_pft=clm3%g%l%c%p%pcs%pft_ctrunc)

    call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
         avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
         ptr_pft=clm3%g%l%c%p%pcs%dispvegc)

    call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
         avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
         ptr_pft=clm3%g%l%c%p%pcs%storvegc)

    call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
         avgflag='A', long_name='total vegetation carbon, excluding cpool', &
         ptr_pft=clm3%g%l%c%p%pcs%totvegc)

    call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
         avgflag='A', long_name='total pft-level carbon, including cpool', &
         ptr_pft=clm3%g%l%c%p%pcs%totpftc)

    if ( use_c13 ) then
       !-------------------------------
       ! C13 state variables - native to PFT 
       !-------------------------------
       
       call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C', &
            ptr_pft=clm3%g%l%c%p%pc13s%leafc)
       
       call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C storage', &
            ptr_pft=clm3%g%l%c%p%pc13s%leafc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C transfer', &
            ptr_pft=clm3%g%l%c%p%pc13s%leafc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C', &
            ptr_pft=clm3%g%l%c%p%pc13s%frootc)
       
       call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C storage', &
            ptr_pft=clm3%g%l%c%p%pc13s%frootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C transfer', &
            ptr_pft=clm3%g%l%c%p%pc13s%frootc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C', &
            ptr_pft=clm3%g%l%c%p%pc13s%livestemc)
       
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C storage', &
            ptr_pft=clm3%g%l%c%p%pc13s%livestemc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C transfer', &
            ptr_pft=clm3%g%l%c%p%pc13s%livestemc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C', &
            ptr_pft=clm3%g%l%c%p%pc13s%deadstemc)
       
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C storage', &
            ptr_pft=clm3%g%l%c%p%pc13s%deadstemc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C transfer', &
            ptr_pft=clm3%g%l%c%p%pc13s%deadstemc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C', &
            ptr_pft=clm3%g%l%c%p%pc13s%livecrootc)
       
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C storage', &
            ptr_pft=clm3%g%l%c%p%pc13s%livecrootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C transfer', &
            ptr_pft=clm3%g%l%c%p%pc13s%livecrootc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C', &
            ptr_pft=clm3%g%l%c%p%pc13s%deadcrootc)
       
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C storage', &
            ptr_pft=clm3%g%l%c%p%pc13s%deadcrootc_storage,  default='inactive')
       
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C transfer', &
            ptr_pft=clm3%g%l%c%p%pc13s%deadcrootc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 growth respiration storage', &
            ptr_pft=clm3%g%l%c%p%pc13s%gresp_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 growth respiration transfer', &
            ptr_pft=clm3%g%l%c%p%pc13s%gresp_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
            avgflag='A', long_name='C13 temporary photosynthate C pool', &
            ptr_pft=clm3%g%l%c%p%pc13s%cpool)
       
       call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
            avgflag='A', long_name='C13 temporary photosynthate C pool', &
            ptr_pft=clm3%g%l%c%p%pc13s%xsmrpool)
       
       call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
            avgflag='A', long_name='C13 pft-level sink for C truncation', &
            ptr_pft=clm3%g%l%c%p%pc13s%pft_ctrunc)
       
       call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
            ptr_pft=clm3%g%l%c%p%pc13s%dispvegc)
       
       call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
            ptr_pft=clm3%g%l%c%p%pc13s%storvegc)
       
       call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
            ptr_pft=clm3%g%l%c%p%pc13s%totvegc)
       
       call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total pft-level carbon, including cpool', &
            ptr_pft=clm3%g%l%c%p%pc13s%totpftc)
    endif

    if ( use_c14 ) then
       !-------------------------------
       ! C14 state variables - native to PFT 
       !-------------------------------
       
       call hist_addfld1d (fname='C14_LEAFC', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C', &
            ptr_pft=clm3%g%l%c%p%pc14s%leafc)
       
       call hist_addfld1d (fname='C14_LEAFC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C storage', &
            ptr_pft=clm3%g%l%c%p%pc14s%leafc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_LEAFC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C transfer', &
            ptr_pft=clm3%g%l%c%p%pc14s%leafc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_FROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C', &
            ptr_pft=clm3%g%l%c%p%pc14s%frootc)
       
       call hist_addfld1d (fname='C14_FROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C storage', &
            ptr_pft=clm3%g%l%c%p%pc14s%frootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_FROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C transfer', &
            ptr_pft=clm3%g%l%c%p%pc14s%frootc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVESTEMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C', &
            ptr_pft=clm3%g%l%c%p%pc14s%livestemc)
       
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C storage', &
            ptr_pft=clm3%g%l%c%p%pc14s%livestemc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C transfer', &
            ptr_pft=clm3%g%l%c%p%pc14s%livestemc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_DEADSTEMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C', &
            ptr_pft=clm3%g%l%c%p%pc14s%deadstemc)
       
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C storage', &
            ptr_pft=clm3%g%l%c%p%pc14s%deadstemc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C transfer', &
            ptr_pft=clm3%g%l%c%p%pc14s%deadstemc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVECROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C', &
            ptr_pft=clm3%g%l%c%p%pc14s%livecrootc)
       
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C storage', &
            ptr_pft=clm3%g%l%c%p%pc14s%livecrootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C transfer', &
            ptr_pft=clm3%g%l%c%p%pc14s%livecrootc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_DEADCROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C', &
            ptr_pft=clm3%g%l%c%p%pc14s%deadcrootc)
       
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C storage', &
            ptr_pft=clm3%g%l%c%p%pc14s%deadcrootc_storage,  default='inactive')
       
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C transfer', &
            ptr_pft=clm3%g%l%c%p%pc14s%deadcrootc_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_GRESP_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 growth respiration storage', &
            ptr_pft=clm3%g%l%c%p%pc14s%gresp_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_GRESP_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 growth respiration transfer', &
            ptr_pft=clm3%g%l%c%p%pc14s%gresp_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL', units='gC14/m^2', &
            avgflag='A', long_name='C14 temporary photosynthate C pool', &
            ptr_pft=clm3%g%l%c%p%pc14s%cpool)
       
       call hist_addfld1d (fname='C14_XSMRPOOL', units='gC14/m^2', &
            avgflag='A', long_name='C14 temporary photosynthate C pool', &
            ptr_pft=clm3%g%l%c%p%pc14s%xsmrpool)
       
       call hist_addfld1d (fname='C14_PFT_CTRUNC', units='gC14/m^2', &
            avgflag='A', long_name='C14 pft-level sink for C truncation', &
            ptr_pft=clm3%g%l%c%p%pc14s%pft_ctrunc)
       
       call hist_addfld1d (fname='C14_DISPVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 displayed veg carbon, excluding storage and cpool', &
            ptr_pft=clm3%g%l%c%p%pc14s%dispvegc)
       
       call hist_addfld1d (fname='C14_STORVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 stored vegetation carbon, excluding cpool', &
            ptr_pft=clm3%g%l%c%p%pc14s%storvegc)
       
       call hist_addfld1d (fname='C14_TOTVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total vegetation carbon, excluding cpool', &
            ptr_pft=clm3%g%l%c%p%pc14s%totvegc)
       
       call hist_addfld1d (fname='C14_TOTPFTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total pft-level carbon, including cpool', &
            ptr_pft=clm3%g%l%c%p%pc14s%totpftc)
    endif

    !-------------------------------
    ! C state variables - native to column
    !-------------------------------
     ! add history fields for all CLAMP CN variables
     call hist_addfld1d (fname='SOILC', units='gC/m^2', &
          avgflag='A', long_name='soil C', &
          ptr_col=clm3%g%l%c%ccs%totsomc)

     call hist_addfld1d (fname='LITTERC', units='gC/m^2', &
          avgflag='A', long_name='litter C', &
          ptr_col=clm3%g%l%c%ccs%totlitc)

    
    do l  = 1, ndecomp_pools
       if ( nlevdecomp_full .gt. 1 ) then
          data2dptr => clm3%g%l%c%ccs%decomp_cpools_vr(:,:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='gC/m^3',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif


       data1dptr => clm3%g%l%c%ccs%decomp_cpools(:,l)
       fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
       longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
       
       call hist_addfld1d (fname=fieldname, units='gC/m^2', &
            avgflag='A', long_name=longname, &
            ptr_col=data1dptr)

       if ( nlevdecomp_full .gt. 1 ) then
          data1dptr => clm3%g%l%c%ccs%decomp_cpools_1m(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_1m'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
          call hist_addfld1d (fname=fieldname, units='gC/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default = 'inactive')
       endif

    end do

    if ( nlevdecomp_full .gt. 1 ) then
       call hist_addfld1d (fname='TOTLITC_1m', units='gC/m^2', &
            avgflag='A', long_name='total litter carbon to 1 meter depth', &
            ptr_col=clm3%g%l%c%ccs%totlitc_1m)
       
       call hist_addfld1d (fname='TOTSOMC_1m', units='gC/m^2', &
            avgflag='A', long_name='total soil organic matter carbon to 1 meter depth', &
            ptr_col=clm3%g%l%c%ccs%totsomc_1m)
    end if

    call hist_addfld1d (fname='COL_CTRUNC', units='gC/m^2',  &
         avgflag='A', long_name='column-level sink for C truncation', &
         ptr_col=clm3%g%l%c%ccs%col_ctrunc)

    call hist_addfld_decomp (fname='NFIXATION_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for biological N fixation', &
         ptr_col=clm3%g%l%c%cps%nfixation_prof, default='inactive')

    call hist_addfld_decomp (fname='NDEP_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for atmospheric N  deposition', &
         ptr_col=clm3%g%l%c%cps%ndep_prof, default='inactive')

    call hist_addfld1d (fname='ALT', units='m', &
         avgflag='A', long_name='current active layer thickness', &
         ptr_col=clm3%g%l%c%cps%alt)

    call hist_addfld1d (fname='ALTMAX', units='m', &
         avgflag='A', long_name='maximum annual active layer thickness', &
         ptr_col=clm3%g%l%c%cps%altmax)

    call hist_addfld1d (fname='ALTMAX_LASTYEAR', units='m', &
         avgflag='A', long_name='maximum prior year active layer thickness', &
         ptr_col=clm3%g%l%c%cps%altmax_lastyear)

    call hist_addfld_decomp (fname='SOM_ADV_COEF', units='m/s',  type2d='levdcmp', &
         avgflag='A', long_name='advection term for vertical SOM translocation', &
         ptr_col=clm3%g%l%c%cps%som_adv_coef, default='inactive')

    call hist_addfld_decomp (fname='SOM_DIFFUS_COEF', units='m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='diffusion coefficient for vertical SOM translocation', &
         ptr_col=clm3%g%l%c%cps%som_diffus_coef, default='inactive')

    call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs', &
         ptr_col=clm3%g%l%c%ccs%seedc)
    
    call hist_addfld1d (fname='TOTLITC', units='gC/m^2', &
         avgflag='A', long_name='total litter carbon', &
         ptr_col=clm3%g%l%c%ccs%totlitc)
    
    call hist_addfld1d (fname='TOTSOMC', units='gC/m^2', &
         avgflag='A', long_name='total soil organic matter carbon', &
         ptr_col=clm3%g%l%c%ccs%totsomc)
    
    call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
         avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool', &
         ptr_col=clm3%g%l%c%ccs%totecosysc)
    
    call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
         avgflag='A', long_name='total column carbon, incl veg and cpool', &
         ptr_col=clm3%g%l%c%ccs%totcolc)

    call hist_addfld1d (fname='PROD10C', units='gC/m^2', &
         avgflag='A', long_name='10-yr wood product C', &
         ptr_col=clm3%g%l%c%ccs%prod10c)

    call hist_addfld1d (fname='PROD100C', units='gC/m^2', &
         avgflag='A', long_name='100-yr wood product C', &
         ptr_col=clm3%g%l%c%ccs%prod100c)

    call hist_addfld1d (fname='TOTPRODC', units='gC/m^2', &
         avgflag='A', long_name='total wood product C', &
         ptr_col=clm3%g%l%c%ccs%totprodc)
    
     call hist_addfld1d (fname='FUELC', units='gC/m^2', &
         avgflag='A', long_name='fuel load', &
         ptr_col=clm3%g%l%c%ccs%fuelc)
    
    
    if ( use_c13 ) then
       !-------------------------------
       ! C13 state variables - native to column
       !-------------------------------
       
       do l = 1, ndecomp_pools
          if ( nlevdecomp_full .gt. 1 ) then
             data2dptr => clm3%g%l%c%cc13s%decomp_cpools_vr(:,:,l)
             fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr)
          endif
          data1dptr => clm3%g%l%c%cc13s%decomp_cpools(:,l)
          fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC13/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       end do
       
       call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
            avgflag='A', long_name='C13 pool for seeding new PFTs', &
            ptr_col=clm3%g%l%c%ccs%seedc)
       
       call hist_addfld1d (fname='C13_COL_CTRUNC', units='gC13/m^2',  &
            avgflag='A', long_name='C13 column-level sink for C truncation', &
            ptr_col=clm3%g%l%c%cc13s%col_ctrunc)
       
       call hist_addfld1d (fname='C13_TOTLITC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total litter carbon', &
            ptr_col=clm3%g%l%c%cc13s%totlitc)
       
       call hist_addfld1d (fname='C13_TOTSOMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total soil organic matter carbon', &
            ptr_col=clm3%g%l%c%cc13s%totsomc)

       if ( nlevdecomp_full .gt. 1 ) then
          call hist_addfld1d (fname='C13_TOTLITC_1m', units='gC13/m^2', &
               avgflag='A', long_name='C13 total litter carbon to 1 meter', &
               ptr_col=clm3%g%l%c%cc13s%totlitc_1m)
          
          call hist_addfld1d (fname='C13_TOTSOMC_1m', units='gC13/m^2', &
               avgflag='A', long_name='C13 total soil organic matter carbon to 1 meter', &
               ptr_col=clm3%g%l%c%cc13s%totsomc_1m)
       endif

       call hist_addfld1d (fname='C13_TOTECOSYSC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=clm3%g%l%c%cc13s%totecosysc)
       
       call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total column carbon, incl veg and cpool', &
            ptr_col=clm3%g%l%c%cc13s%totcolc)
       
       call hist_addfld1d (fname='C13_PROD10C', units='gC13/m^2', &
            avgflag='A', long_name='C13 10-yr wood product C', &
            ptr_col=clm3%g%l%c%cc13s%prod10c)
       
       call hist_addfld1d (fname='C13_PROD100C', units='gC13/m^2', &
            avgflag='A', long_name='C13 100-yr wood product C', &
            ptr_col=clm3%g%l%c%cc13s%prod100c)
       
       call hist_addfld1d (fname='C13_TOTPRODC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total wood product C', &
            ptr_col=clm3%g%l%c%cc13s%totprodc)
    endif
    
    if ( use_c14 ) then
       !-------------------------------
       ! C14 state variables - native to column
       !-------------------------------
       
       
       do l = 1, ndecomp_pools
          if ( nlevdecomp_full .gt. 1 ) then
             data2dptr => clm3%g%l%c%cc14s%decomp_cpools_vr(:,:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr)
          endif
          
          data1dptr => clm3%g%l%c%cc14s%decomp_cpools(:,l)
          fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC14/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
          
          if ( nlevdecomp_full .gt. 1 ) then
             data1dptr => clm3%g%l%c%cc14s%decomp_cpools_1m(:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_1m'
             longname =  'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr, default='inactive')
          endif
          
       end do
       
       call hist_addfld1d (fname='C14_SEEDC', units='gC14/m^2', &
            avgflag='A', long_name='C14 pool for seeding new PFTs', &
            ptr_col=clm3%g%l%c%cc14s%seedc)
       
       call hist_addfld1d (fname='C14_COL_CTRUNC', units='gC14/m^2', &
            avgflag='A', long_name='C14 column-level sink for C truncation', &
            ptr_col=clm3%g%l%c%cc14s%col_ctrunc)
       
       call hist_addfld1d (fname='C14_TOTLITC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total litter carbon', &
            ptr_col=clm3%g%l%c%cc14s%totlitc)
       
       call hist_addfld1d (fname='C14_TOTSOMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total soil organic matter carbon', &
            ptr_col=clm3%g%l%c%cc14s%totsomc)

       if ( nlevdecomp_full .gt. 1 ) then       
          call hist_addfld1d (fname='C14_TOTLITC_1m', units='gC14/m^2', &
               avgflag='A', long_name='C14 total litter carbon to 1 meter', &
               ptr_col=clm3%g%l%c%cc14s%totlitc_1m)
          
          call hist_addfld1d (fname='C14_TOTSOMC_1m', units='gC14/m^2', &
               avgflag='A', long_name='C14 total soil organic matter carbon to 1 meter', &
               ptr_col=clm3%g%l%c%cc14s%totsomc_1m)
       endif

       call hist_addfld1d (fname='C14_TOTECOSYSC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=clm3%g%l%c%cc14s%totecosysc)
       
       call hist_addfld1d (fname='C14_TOTCOLC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total column carbon, incl veg and cpool', &
            ptr_col=clm3%g%l%c%cc14s%totcolc)
       
       call hist_addfld1d (fname='C14_PROD10C', units='gC14/m^2', &
            avgflag='A', long_name='C14 10-yr wood product C', &
            ptr_col=clm3%g%l%c%cc14s%prod10c)
       
       call hist_addfld1d (fname='C14_PROD100C', units='gC14/m^2', &
            avgflag='A', long_name='C14 100-yr wood product C', &
            ptr_col=clm3%g%l%c%cc14s%prod100c)
       
       call hist_addfld1d (fname='C14_TOTPRODC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total wood product C', &
            ptr_col=clm3%g%l%c%cc14s%totprodc)
    endif

    !-------------------------------
    ! N state variables - native to PFT
    !-------------------------------

    if (crop_prog) then
       call hist_addfld1d (fname='GRAINN', units='gN/m^2', &
            avgflag='A', long_name='grain N', &
            ptr_pft=clm3%g%l%c%p%pns%grainn)
    end if

    call hist_addfld1d (fname='LEAFN', units='gN/m^2', &
         avgflag='A', long_name='leaf N', &
         ptr_pft=clm3%g%l%c%p%pns%leafn)

    call hist_addfld1d (fname='LEAFN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='leaf N storage', &
         ptr_pft=clm3%g%l%c%p%pns%leafn_storage, default='inactive')

    call hist_addfld1d (fname='LEAFN_XFER', units='gN/m^2', &
         avgflag='A', long_name='leaf N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%leafn_xfer, default='inactive')

    call hist_addfld1d (fname='FROOTN', units='gN/m^2', &
         avgflag='A', long_name='fine root N', &
         ptr_pft=clm3%g%l%c%p%pns%frootn)

    call hist_addfld1d (fname='FROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='fine root N storage', &
         ptr_pft=clm3%g%l%c%p%pns%frootn_storage, default='inactive')

    call hist_addfld1d (fname='FROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='fine root N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%frootn_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN', units='gN/m^2', &
         avgflag='A', long_name='live stem N', &
         ptr_pft=clm3%g%l%c%p%pns%livestemn)

    call hist_addfld1d (fname='LIVESTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live stem N storage', &
         ptr_pft=clm3%g%l%c%p%pns%livestemn_storage, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live stem N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%livestemn_xfer, default='inactive')

    call hist_addfld1d (fname='DEADSTEMN', units='gN/m^2', &
         avgflag='A', long_name='dead stem N', &
         ptr_pft=clm3%g%l%c%p%pns%deadstemn)

    call hist_addfld1d (fname='DEADSTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead stem N storage', &
         ptr_pft=clm3%g%l%c%p%pns%deadstemn_storage, default='inactive')

    call hist_addfld1d (fname='DEADSTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead stem N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%deadstemn_xfer, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N', &
         ptr_pft=clm3%g%l%c%p%pns%livecrootn)

    call hist_addfld1d (fname='LIVECROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N storage', &
         ptr_pft=clm3%g%l%c%p%pns%livecrootn_storage, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%livecrootn_xfer, default='inactive')

    call hist_addfld1d (fname='DEADCROOTN', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N', &
         ptr_pft=clm3%g%l%c%p%pns%deadcrootn)

    call hist_addfld1d (fname='DEADCROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N storage', &
         ptr_pft=clm3%g%l%c%p%pns%deadcrootn_storage, default='inactive')

    call hist_addfld1d (fname='DEADCROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%deadcrootn_xfer, default='inactive')

    call hist_addfld1d (fname='RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='plant pool of retranslocated N', &
         ptr_pft=clm3%g%l%c%p%pns%retransn)

    call hist_addfld1d (fname='NPOOL', units='gN/m^2', &
         avgflag='A', long_name='temporary plant N pool', &
         ptr_pft=clm3%g%l%c%p%pns%npool, default='inactive')

    call hist_addfld1d (fname='PFT_NTRUNC', units='gN/m^2', &
         avgflag='A', long_name='pft-level sink for N truncation', &
         ptr_pft=clm3%g%l%c%p%pns%pft_ntrunc)

    call hist_addfld1d (fname='DISPVEGN', units='gN/m^2', &
         avgflag='A', long_name='displayed vegetation nitrogen', &
         ptr_pft=clm3%g%l%c%p%pns%dispvegn)

    call hist_addfld1d (fname='STORVEGN', units='gN/m^2', &
         avgflag='A', long_name='stored vegetation nitrogen', &
         ptr_pft=clm3%g%l%c%p%pns%storvegn)

    call hist_addfld1d (fname='TOTVEGN', units='gN/m^2', &
         avgflag='A', long_name='total vegetation nitrogen', &
         ptr_pft=clm3%g%l%c%p%pns%totvegn)

    call hist_addfld1d (fname='TOTPFTN', units='gN/m^2', &
         avgflag='A', long_name='total PFT-level nitrogen', &
         ptr_pft=clm3%g%l%c%p%pns%totpftn)

    !-------------------------------
    ! N state variables - native to column
    !-------------------------------

    do l  = 1, ndecomp_pools
       if ( nlevdecomp_full .gt. 1 ) then
          data2dptr => clm3%g%l%c%cns%decomp_npools_vr(:,:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif
       data1dptr => clm3%g%l%c%cns%decomp_npools(:,l)
       fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N'
       longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N'
       call hist_addfld1d (fname=fieldname, units='gN/m^2', &
            avgflag='A', long_name=longname, &
            ptr_col=data1dptr)

       if ( nlevdecomp_full .gt. 1 ) then
          data1dptr => clm3%g%l%c%cns%decomp_npools_1m(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N_1m'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N to 1 meter'
          call hist_addfld1d (fname=fieldname, units='gN/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default = 'inactive')
       endif
    end do


    if ( nlevdecomp_full .gt. 1 ) then
       call hist_addfld1d (fname='SMINN', units='gN/m^2', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=clm3%g%l%c%cns%sminn)
       
       call hist_addfld1d (fname='TOTLITN_1m', units='gN/m^2', &
            avgflag='A', long_name='total litter N to 1 meter', &
            ptr_col=clm3%g%l%c%cns%totlitn_1m)
       
       call hist_addfld1d (fname='TOTSOMN_1m', units='gN/m^2', &
            avgflag='A', long_name='total soil organic matter N to 1 meter', &
            ptr_col=clm3%g%l%c%cns%totsomn_1m)
    endif

    call hist_addfld1d (fname='COL_NTRUNC', units='gN/m^2',  &
         avgflag='A', long_name='column-level sink for N truncation', &
         ptr_col=clm3%g%l%c%cns%col_ntrunc)


#ifdef NITRIF_DENITRIF

    call hist_addfld_decomp (fname='SMIN_NO3'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil mineral NO3 (vert. res.)', &
         ptr_col=clm3%g%l%c%cns%smin_no3_vr)

    call hist_addfld_decomp (fname='SMIN_NH4'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil mineral NH4 (vert. res.)', &
         ptr_col=clm3%g%l%c%cns%smin_nh4_vr)

    if ( nlevdecomp_full .gt. 1 ) then
       call hist_addfld1d (fname='SMIN_NO3', units='gN/m^2', &
            avgflag='A', long_name='soil mineral NO3', &
            ptr_col=clm3%g%l%c%cns%smin_no3)
       
       call hist_addfld1d (fname='SMIN_NH4', units='gN/m^2', &
            avgflag='A', long_name='soil mineral NH4', &
            ptr_col=clm3%g%l%c%cns%smin_nh4)
    endif

    call hist_addfld_decomp (fname='SMINN'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil mineral N', &
         ptr_col=clm3%g%l%c%cns%sminn_vr, default = 'inactive')

#else
    call hist_addfld_decomp (fname='SMINN'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil mineral N', &
         ptr_col=clm3%g%l%c%cns%sminn_vr)

#endif
    call hist_addfld1d (fname='TOTLITN', units='gN/m^2', &
         avgflag='A', long_name='total litter N', &
         ptr_col=clm3%g%l%c%cns%totlitn)

    call hist_addfld1d (fname='TOTSOMN', units='gN/m^2', &
         avgflag='A', long_name='total soil organic matter N', &
         ptr_col=clm3%g%l%c%cns%totsomn)

    call hist_addfld1d (fname='TOTECOSYSN', units='gN/m^2', &
         avgflag='A', long_name='total ecosystem N', &
         ptr_col=clm3%g%l%c%cns%totecosysn)

    call hist_addfld1d (fname='TOTCOLN', units='gN/m^2', &
         avgflag='A', long_name='total column-level N', &
         ptr_col=clm3%g%l%c%cns%totcoln)

    call hist_addfld1d (fname='SEEDN', units='gN/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs ', &
         ptr_col=clm3%g%l%c%cns%seedn)

    call hist_addfld1d (fname='PROD10N', units='gN/m^2', &
         avgflag='A', long_name='10-yr wood product N', &
         ptr_col=clm3%g%l%c%cns%prod10n)

    call hist_addfld1d (fname='PROD100N', units='gN/m^2', &
         avgflag='A', long_name='100-yr wood product N', &
         ptr_col=clm3%g%l%c%cns%prod100n)

    call hist_addfld1d (fname='TOTPRODN', units='gN/m^2', &
         avgflag='A', long_name='total wood product N', &
         ptr_col=clm3%g%l%c%cns%totprodn)

    !-------------------------------
    ! C flux variables - native to PFT
    !-------------------------------

     ! add history fields for all CLAMP CN variables

    if (crop_prog) then
       call hist_addfld1d (fname='GRAINC_TO_FOOD', units='gC/m^2/s', &
            avgflag='A', long_name='grain C to food', &
            ptr_pft=clm3%g%l%c%p%pcf%grainc_to_food)
     end if

     call hist_addfld1d (fname='WOODC_ALLOC', units='gC/m^2/s', &
          avgflag='A', long_name='wood C allocation', &
          ptr_pft=clm3%g%l%c%p%pcf%woodc_alloc)

     call hist_addfld1d (fname='WOODC_LOSS', units='gC/m^2/s', &
          avgflag='A', long_name='wood C loss', &
          ptr_pft=clm3%g%l%c%p%pcf%woodc_loss)

     call hist_addfld1d (fname='LEAFC_LOSS', units='gC/m^2/s', &
          avgflag='A', long_name='leaf C loss', &
          ptr_pft=clm3%g%l%c%p%pcf%leafc_loss)

     call hist_addfld1d (fname='LEAFC_ALLOC', units='gC/m^2/s', &
          avgflag='A', long_name='leaf C allocation', &
          ptr_pft=clm3%g%l%c%p%pcf%leafc_alloc)

     call hist_addfld1d (fname='FROOTC_LOSS', units='gC/m^2/s', &
          avgflag='A', long_name='fine root C loss', &
          ptr_pft=clm3%g%l%c%p%pcf%frootc_loss)

     call hist_addfld1d (fname='FROOTC_ALLOC', units='gC/m^2/s', &
          avgflag='A', long_name='fine root C allocation', &
          ptr_pft=clm3%g%l%c%p%pcf%frootc_alloc)

    call hist_addfld1d (fname='PSNSUN', units='umolCO2/m^2/s', &
         avgflag='A', long_name='sunlit leaf photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%psnsun)

    call hist_addfld1d (fname='PSNSHA', units='umolCO2/m^2/s', &
         avgflag='A', long_name='shaded leaf photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%psnsha)

    call hist_addfld1d (fname='M_LEAFC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_GRESP_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter, default='inactive')

     call hist_addfld1d (fname='M_LEAFC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_to_fire, default='inactive')
 
   call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_xfer_to_fire, default='inactive')

   call hist_addfld1d (fname='M_LIVESTEMC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_to_fire, default='inactive')
 
    call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_storage_to_fire, default='inactive')

   call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_xfer_to_fire, default='inactive')

   call hist_addfld1d (fname='M_DEADSTEMC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_to_fire, default='inactive')
 
  call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_storage_to_fire, default='inactive')

   call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_fire, default='inactive')

   call hist_addfld1d (fname='M_FROOTC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_to_fire, default='inactive')

   call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_xfer_to_fire, default='inactive')

   call hist_addfld1d (fname='M_LIVEROOTC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_to_fire, default='inactive')

   call hist_addfld1d (fname='M_LIVEROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVEROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADROOTC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_to_fire, default='inactive')

   call hist_addfld1d (fname='M_DEADROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_GRESP_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_GRESP_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_xfer_to_fire, default='inactive')
  
    call hist_addfld1d (fname='M_LEAFC_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_to_litter_fire, default='inactive')

! add by F. Li and S. Levis
   call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_storage_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C transfer fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter_fire, default='inactive')

   call hist_addfld1d (fname='M_LIVESTEMC_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_to_litter_fire, default='inactive')
 
    call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter_fire, default='inactive')

   call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C transfer fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter_fire, default='inactive')
	 
 call hist_addfld1d (fname='M_LIVESTEMC_TO_DEADSTEMC_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C fire mortality to dead stem C', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_to_deadstemc_fire, default='inactive')

   call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_to_litter_fire, default='inactive')
 
  call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter_fire, default='inactive')

   call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C transfer fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter_fire, default='inactive')

   call hist_addfld1d (fname='M_FROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_to_litter_fire, default='inactive')

   call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_storage_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C transfer fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter_fire, default='inactive')

   call hist_addfld1d (fname='M_LIVEROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live root C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_to_litter_fire, default='inactive')

   call hist_addfld1d (fname='M_LIVEROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live root C storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter_fire, default='inactive')
	 
    call hist_addfld1d (fname='M_LIVEROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live root C transfer fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter_fire, default='inactive')
	 
   call hist_addfld1d (fname='M_LIVEROOTC_TO_DEADROOTC_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live root C fire mortality to dead root C', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_to_deadcrootc_fire, default='inactive')


    call hist_addfld1d (fname='M_DEADROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead root C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_to_litter_fire, default='inactive')

   call hist_addfld1d (fname='M_DEADROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead root C storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead root C transfer fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter_fire, default='inactive')

 call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter_fire,  default='inactive')

    call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration storage fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_storage_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_GRESP_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration transfer fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter_fire, default='inactive')   

    call hist_addfld1d (fname='LEAFC_XFER_TO_LEAFC', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%leafc_xfer_to_leafc, default='inactive')

    call hist_addfld1d (fname='FROOTC_XFER_TO_FROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%frootc_xfer_to_frootc, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC_XFER_TO_LIVESTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%livestemc_xfer_to_livestemc, default='inactive')

    call hist_addfld1d (fname='DEADSTEMC_XFER_TO_DEADSTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%deadstemc_xfer_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC_XFER_TO_LIVECROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%livecrootc_xfer_to_livecrootc, default='inactive')

    call hist_addfld1d (fname='DEADCROOTC_XFER_TO_DEADCROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%deadcrootc_xfer_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='LEAFC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C litterfall', &
         ptr_pft=clm3%g%l%c%p%pcf%leafc_to_litter, default='inactive')

    call hist_addfld1d (fname='FROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C litterfall', &
         ptr_pft=clm3%g%l%c%p%pcf%frootc_to_litter, default='inactive')

    call hist_addfld1d (fname='LEAF_MR', units='gC/m^2/s', &
         avgflag='A', long_name='leaf maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%leaf_mr)

    call hist_addfld1d (fname='FROOT_MR', units='gC/m^2/s', &
         avgflag='A', long_name='fine root maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%froot_mr, default='inactive')

    call hist_addfld1d (fname='LIVESTEM_MR', units='gC/m^2/s', &
         avgflag='A', long_name='live stem maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%livestem_mr, default='inactive')

    call hist_addfld1d (fname='LIVECROOT_MR', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%livecroot_mr, default='inactive')

    call hist_addfld1d (fname='PSNSUN_TO_CPOOL', units='gC/m^2/s', &
         avgflag='A', long_name='C fixation from sunlit canopy', &
         ptr_pft=clm3%g%l%c%p%pcf%psnsun_to_cpool)

    call hist_addfld1d (fname='PSNSHADE_TO_CPOOL', units='gC/m^2/s', &
         avgflag='A', long_name='C fixation from shaded canopy', &
         ptr_pft=clm3%g%l%c%p%pcf%psnshade_to_cpool)

    call hist_addfld1d (fname='CPOOL_TO_LEAFC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to leaf C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_leafc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LEAFC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to leaf C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_leafc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_FROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to fine root C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_frootc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_FROOTC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to fine root C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_frootc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to live stem C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_livestemc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to live stem C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_livestemc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to dead stem C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to dead stem C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_deadstemc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_livecrootc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_livecrootc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_deadcrootc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_GRESP_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to growth respiration storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_gresp_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_LEAF_GR', units='gC/m^2/s', &
         avgflag='A', long_name='leaf growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_leaf_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LEAF_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='leaf growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_leaf_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_LEAF_GR', units='gC/m^2/s', &
         avgflag='A', long_name='leaf growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_leaf_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_FROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='fine root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_froot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_FROOT_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='fine root  growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_froot_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_FROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='fine root  growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_froot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LIVESTEM_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live stem growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_livestem_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LIVESTEM_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live stem growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_livestem_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_LIVESTEM_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live stem growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_livestem_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_DEADSTEM_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_deadstem_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_DEADSTEM_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_deadstem_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_DEADSTEM_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_deadstem_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LIVECROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_livecroot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LIVECROOT_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_livecroot_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_LIVECROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_livecroot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_DEADCROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_deadcroot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_DEADCROOT_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_deadcroot_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_DEADCROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_deadcroot_gr, default='inactive')

    call hist_addfld1d (fname='LEAFC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%leafc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='FROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%frootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%livestemc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='DEADSTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%deadstemc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%livecrootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='DEADCROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%deadcrootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='GRESP_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%gresp_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC_TO_DEADSTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C turnover', &
         ptr_pft=clm3%g%l%c%p%pcf%livestemc_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC_TO_DEADCROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C turnover', &
         ptr_pft=clm3%g%l%c%p%pcf%livecrootc_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='GPP', units='gC/m^2/s', &
         avgflag='A', long_name='gross primary production', &
         ptr_pft=clm3%g%l%c%p%pcf%gpp)

    call hist_addfld1d (fname='MR', units='gC/m^2/s', &
         avgflag='A', long_name='maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%mr)

    call hist_addfld1d (fname='CURRENT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='growth resp for new growth displayed in this timestep', &
         ptr_pft=clm3%g%l%c%p%pcf%current_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_GR', units='gC/m^2/s', &
         avgflag='A', long_name='growth resp for transfer growth displayed in this timestep', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_gr, default='inactive')

    call hist_addfld1d (fname='STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='growth resp for growth sent to storage for later display', &
         ptr_pft=clm3%g%l%c%p%pcf%storage_gr, default='inactive')

    call hist_addfld1d (fname='GR', units='gC/m^2/s', &
         avgflag='A', long_name='total growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%gr)

    call hist_addfld1d (fname='AR', units='gC/m^2/s', &
         avgflag='A', long_name='autotrophic respiration (MR + GR)', &
         ptr_pft=clm3%g%l%c%p%pcf%ar)

    call hist_addfld1d (fname='RR', units='gC/m^2/s', &
         avgflag='A', long_name='root respiration (fine root MR + total root GR)', &
         ptr_pft=clm3%g%l%c%p%pcf%rr)

    call hist_addfld1d (fname='NPP', units='gC/m^2/s', &
         avgflag='A', long_name='net primary production', &
         ptr_pft=clm3%g%l%c%p%pcf%npp)

    call hist_addfld1d (fname='AGNPP', units='gC/m^2/s', &
         avgflag='A', long_name='aboveground NPP', &
         ptr_pft=clm3%g%l%c%p%pcf%agnpp)

    call hist_addfld1d (fname='BGNPP', units='gC/m^2/s', &
         avgflag='A', long_name='belowground NPP', &
         ptr_pft=clm3%g%l%c%p%pcf%bgnpp)

    call hist_addfld1d (fname='LITFALL', units='gC/m^2/s', &
         avgflag='A', long_name='litterfall (leaves and fine roots)', &
         ptr_pft=clm3%g%l%c%p%pcf%litfall)

    call hist_addfld1d (fname='VEGFIRE', units='gC/m^2/s', &
         avgflag='A', long_name='pft-level fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%vegfire, default='inactive')

    call hist_addfld1d (fname='WOOD_HARVESTC', units='gC/m^2/s', &
         avgflag='A', long_name='wood harvest carbon (to product pools)', &
         ptr_pft=clm3%g%l%c%p%pcf%wood_harvestc)

    call hist_addfld1d (fname='PFT_FIRE_CLOSS', units='gC/m^2/s', &
         avgflag='A', long_name='total pft-level fire C loss for non-peat fires outside land-type converted region', &
         ptr_pft=clm3%g%l%c%p%pcf%pft_fire_closs)

    if ( use_c13 ) then
       !-------------------------------
       ! C13 flux variables - native to PFT
       !-------------------------------
       
       call hist_addfld1d (fname='C13_PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C13 sunlit leaf photosynthesis', &
            ptr_pft=clm3%g%l%c%p%pc13f%psnsun)
       
       call hist_addfld1d (fname='C13_PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C13 shaded leaf photosynthesis', &
            ptr_pft=clm3%g%l%c%p%pc13f%psnsha)
       
       call hist_addfld1d (fname='C13_M_LEAFC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_FROOTC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_gresp_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_gresp_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LEAFC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_FROOTC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_storage_to_fire,  default='inactive')
       
       call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C fire mortality to litter', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_to_litter_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C fire mortality to litter', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_to_litter_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_gresp_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%m_gresp_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C13_LEAFC_XFER_TO_LEAFC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%leafc_xfer_to_leafc, default='inactive')
       
       call hist_addfld1d (fname='C13_FROOTC_XFER_TO_FROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%frootc_xfer_to_frootc, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER_TO_LIVESTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%livestemc_xfer_to_livestemc, default='inactive')
       
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER_TO_DEADSTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%deadstemc_xfer_to_deadstemc, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER_TO_LIVECROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%livecrootc_xfer_to_livecrootc, default='inactive')
       
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER_TO_DEADCROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%deadcrootc_xfer_to_deadcrootc, default='inactive')
       
       call hist_addfld1d (fname='C13_LEAFC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C litterfall', &
            ptr_pft=clm3%g%l%c%p%pc13f%leafc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_FROOTC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C litterfall', &
            ptr_pft=clm3%g%l%c%p%pc13f%frootc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C13_LEAF_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%leaf_mr, default='inactive')
       
       call hist_addfld1d (fname='C13_FROOT_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%froot_mr, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVESTEM_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%livestem_mr, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVECROOT_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%livecroot_mr, default='inactive')
       
       call hist_addfld1d (fname='C13_PSNSUN_TO_CPOOL', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 C fixation from sunlit canopy', &
            ptr_pft=clm3%g%l%c%p%pc13f%psnsun_to_cpool)
       
       call hist_addfld1d (fname='C13_PSNSHADE_TO_CPOOL', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 C fixation from shaded canopy', &
            ptr_pft=clm3%g%l%c%p%pc13f%psnshade_to_cpool)
       
       call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to leaf C', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_leafc, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to leaf C storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_leafc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to fine root C', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_frootc, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to fine root C storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_frootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to live stem C', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_livestemc, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to live stem C storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_livestemc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to dead stem C', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_deadstemc, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to dead stem C storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_deadstemc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to live coarse root C', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_livecrootc, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to live coarse root C storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_livecrootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to dead coarse root C', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_deadcrootc, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to dead coarse root C storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_deadcrootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_TO_GRESP_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to growth respiration storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_gresp_storage, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_LEAF_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_leaf_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_LEAF_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_leaf_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_TRANSFER_LEAF_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%transfer_leaf_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_FROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_froot_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_FROOT_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root  growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_froot_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_TRANSFER_FROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root  growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%transfer_froot_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_livestem_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_livestem_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_TRANSFER_LIVESTEM_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%transfer_livestem_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_deadstem_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_deadstem_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_TRANSFER_DEADSTEM_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%transfer_deadstem_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_livecroot_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_livecroot_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_TRANSFER_LIVECROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%transfer_livecroot_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_deadcroot_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%cpool_deadcroot_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_TRANSFER_DEADCROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc13f%transfer_deadcroot_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_LEAFC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc13f%leafc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_FROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc13f%frootc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc13f%livestemc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc13f%deadstemc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc13f%livecrootc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc13f%deadcrootc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_GRESP_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc13f%gresp_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVESTEMC_TO_DEADSTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C turnover', &
            ptr_pft=clm3%g%l%c%p%pc13f%livestemc_to_deadstemc, default='inactive')
       
       call hist_addfld1d (fname='C13_LIVECROOTC_TO_DEADCROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C turnover', &
            ptr_pft=clm3%g%l%c%p%pc13f%livecrootc_to_deadcrootc, default='inactive')
       
       call hist_addfld1d (fname='C13_GPP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 gross primary production', &
            ptr_pft=clm3%g%l%c%p%pc13f%gpp)
       
       call hist_addfld1d (fname='C13_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%mr)
       
       call hist_addfld1d (fname='C13_CURRENT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth resp for new growth displayed in this timestep', &
            ptr_pft=clm3%g%l%c%p%pc13f%current_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_TRANSFER_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth resp for transfer growth displayed in this timestep', &
            ptr_pft=clm3%g%l%c%p%pc13f%transfer_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth resp for growth sent to storage for later display', &
            ptr_pft=clm3%g%l%c%p%pc13f%storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C13_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc13f%gr)
       
       call hist_addfld1d (fname='C13_AR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 autotrophic respiration (MR + GR)', &
            ptr_pft=clm3%g%l%c%p%pc13f%ar)
       
       call hist_addfld1d (fname='C13_RR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 root respiration (fine root MR + total root GR)', &
            ptr_pft=clm3%g%l%c%p%pc13f%rr)
       
       call hist_addfld1d (fname='C13_NPP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 net primary production', &
            ptr_pft=clm3%g%l%c%p%pc13f%npp)
       
       call hist_addfld1d (fname='C13_AGNPP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 aboveground NPP', &
            ptr_pft=clm3%g%l%c%p%pc13f%agnpp)
       
       call hist_addfld1d (fname='C13_BGNPP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 belowground NPP', &
            ptr_pft=clm3%g%l%c%p%pc13f%bgnpp)
       
       call hist_addfld1d (fname='C13_LITFALL', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 litterfall (leaves and fine roots)', &
            ptr_pft=clm3%g%l%c%p%pc13f%litfall, default='inactive')
       
       call hist_addfld1d (fname='C13_VEGFIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 pft-level fire loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%vegfire, default='inactive')
       
       call hist_addfld1d (fname='C13_PFT_FIRE_CLOSS', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total pft-level fire C loss', &
            ptr_pft=clm3%g%l%c%p%pc13f%pft_fire_closs)
    endif
    
    if ( use_c14 ) then
       !-------------------------------
       ! C14 flux variables - native to PFT
       !-------------------------------
       
       call hist_addfld1d (fname='C14_PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C14 sunlit leaf photosynthesis', &
            ptr_pft=clm3%g%l%c%p%pc14f%psnsun)
       
       call hist_addfld1d (fname='C14_PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C14 shaded leaf photosynthesis', &
            ptr_pft=clm3%g%l%c%p%pc14f%psnsha)
       
       call hist_addfld1d (fname='C14_M_LEAFC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_leafc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_FROOTC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_frootc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LEAFC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_leafc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_FROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_frootc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVESTEMC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livestemc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADSTEMC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadstemc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVECROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livecrootc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADCROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadcrootc_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LEAFC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_leafc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_FROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_frootc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVESTEMC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livestemc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADSTEMC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadstemc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVECROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livecrootc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADCROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadcrootc_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVESTEMC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livestemc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadstemc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVECROOTC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livecrootc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadcrootc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_GRESP_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration storage mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_gresp_storage_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_GRESP_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration transfer mortality', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_gresp_xfer_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LEAFC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_leafc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_FROOTC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_frootc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LEAFC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_leafc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_FROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_frootc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVESTEMC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livestemc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADSTEMC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadstemc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVECROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livecrootc_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADCROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadcrootc_storage_to_fire,  default='inactive')
       
       call hist_addfld1d (fname='C14_M_LEAFC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_leafc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_FROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_frootc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVESTEMC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livestemc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADSTEMC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadstemc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVECROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livecrootc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADCROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadcrootc_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVESTEMC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livestemc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadstemc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_LITTER_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C fire mortality to litter', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadstemc_to_litter_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_LIVECROOTC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_livecrootc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadcrootc_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_LITTER_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C fire mortality to litter', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_deadcrootc_to_litter_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_GRESP_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration storage fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_gresp_storage_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_M_GRESP_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration transfer fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%m_gresp_xfer_to_fire, default='inactive')
       
       call hist_addfld1d (fname='C14_LEAFC_XFER_TO_LEAFC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%leafc_xfer_to_leafc, default='inactive')
       
       call hist_addfld1d (fname='C14_FROOTC_XFER_TO_FROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%frootc_xfer_to_frootc, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER_TO_LIVESTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%livestemc_xfer_to_livestemc, default='inactive')
       
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER_TO_DEADSTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%deadstemc_xfer_to_deadstemc, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER_TO_LIVECROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%livecrootc_xfer_to_livecrootc, default='inactive')
       
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER_TO_DEADCROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C growth from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%deadcrootc_xfer_to_deadcrootc, default='inactive')
       
       call hist_addfld1d (fname='C14_LEAFC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C litterfall', &
            ptr_pft=clm3%g%l%c%p%pc14f%leafc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_FROOTC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C litterfall', &
            ptr_pft=clm3%g%l%c%p%pc14f%frootc_to_litter, default='inactive')
       
       call hist_addfld1d (fname='C14_LEAF_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%leaf_mr, default='inactive')
       
       call hist_addfld1d (fname='C14_FROOT_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%froot_mr, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVESTEM_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%livestem_mr, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVECROOT_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%livecroot_mr, default='inactive')
       
       call hist_addfld1d (fname='C14_PSNSUN_TO_CPOOL', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 C fixation from sunlit canopy', &
            ptr_pft=clm3%g%l%c%p%pc14f%psnsun_to_cpool)
       
       call hist_addfld1d (fname='C14_PSNSHADE_TO_CPOOL', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 C fixation from shaded canopy', &
            ptr_pft=clm3%g%l%c%p%pc14f%psnshade_to_cpool)
       
       call hist_addfld1d (fname='C14_CPOOL_TO_LEAFC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to leaf C', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_leafc, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_LEAFC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to leaf C storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_leafc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_FROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to fine root C', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_frootc, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_FROOTC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to fine root C storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_frootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_LIVESTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to live stem C', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_livestemc, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_LIVESTEMC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to live stem C storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_livestemc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_DEADSTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to dead stem C', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_deadstemc, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_DEADSTEMC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to dead stem C storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_deadstemc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_LIVECROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to live coarse root C', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_livecrootc, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_LIVECROOTC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to live coarse root C storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_livecrootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_DEADCROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to dead coarse root C', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_deadcrootc, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_DEADCROOTC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to dead coarse root C storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_deadcrootc_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_TO_GRESP_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to growth respiration storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_to_gresp_storage, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_LEAF_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_leaf_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_LEAF_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_leaf_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_TRANSFER_LEAF_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%transfer_leaf_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_FROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_froot_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_FROOT_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root  growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_froot_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_TRANSFER_FROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root  growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%transfer_froot_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_LIVESTEM_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_livestem_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_LIVESTEM_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_livestem_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_TRANSFER_LIVESTEM_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%transfer_livestem_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_DEADSTEM_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_deadstem_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_DEADSTEM_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_deadstem_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_TRANSFER_DEADSTEM_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%transfer_deadstem_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_LIVECROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_livecroot_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_LIVECROOT_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_livecroot_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_TRANSFER_LIVECROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%transfer_livecroot_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_DEADCROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_deadcroot_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_CPOOL_DEADCROOT_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root growth respiration to storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%cpool_deadcroot_storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_TRANSFER_DEADCROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root growth respiration from storage', &
            ptr_pft=clm3%g%l%c%p%pc14f%transfer_deadcroot_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_LEAFC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc14f%leafc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_FROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc14f%frootc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc14f%livestemc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc14f%deadstemc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc14f%livecrootc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc14f%deadcrootc_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_GRESP_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration shift storage to transfer', &
            ptr_pft=clm3%g%l%c%p%pc14f%gresp_storage_to_xfer, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVESTEMC_TO_DEADSTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C turnover', &
            ptr_pft=clm3%g%l%c%p%pc14f%livestemc_to_deadstemc, default='inactive')
       
       call hist_addfld1d (fname='C14_LIVECROOTC_TO_DEADCROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C turnover', &
            ptr_pft=clm3%g%l%c%p%pc14f%livecrootc_to_deadcrootc, default='inactive')
       
       call hist_addfld1d (fname='C14_GPP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 gross primary production', &
            ptr_pft=clm3%g%l%c%p%pc14f%gpp)
       
       call hist_addfld1d (fname='C14_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 maintenance respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%mr)
       
       call hist_addfld1d (fname='C14_CURRENT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth resp for new growth displayed in this timestep', &
            ptr_pft=clm3%g%l%c%p%pc14f%current_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_TRANSFER_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth resp for transfer growth displayed in this timestep', &
            ptr_pft=clm3%g%l%c%p%pc14f%transfer_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth resp for growth sent to storage for later display', &
            ptr_pft=clm3%g%l%c%p%pc14f%storage_gr, default='inactive')
       
       call hist_addfld1d (fname='C14_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total growth respiration', &
            ptr_pft=clm3%g%l%c%p%pc14f%gr)
       
       call hist_addfld1d (fname='C14_AR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 autotrophic respiration (MR + GR)', &
            ptr_pft=clm3%g%l%c%p%pc14f%ar)
       
       call hist_addfld1d (fname='C14_RR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 root respiration (fine root MR + total root GR)', &
            ptr_pft=clm3%g%l%c%p%pc14f%rr)
       
       call hist_addfld1d (fname='C14_NPP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 net primary production', &
            ptr_pft=clm3%g%l%c%p%pc14f%npp)
       
       call hist_addfld1d (fname='C14_AGNPP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 aboveground NPP', &
            ptr_pft=clm3%g%l%c%p%pc14f%agnpp)
       
       call hist_addfld1d (fname='C14_BGNPP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 belowground NPP', &
            ptr_pft=clm3%g%l%c%p%pc14f%bgnpp)
       
       call hist_addfld1d (fname='C14_LITFALL', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 litterfall (leaves and fine roots)', &
            ptr_pft=clm3%g%l%c%p%pc14f%litfall, default='inactive')
       
       call hist_addfld1d (fname='C14_VEGFIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 pft-level fire loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%vegfire, default='inactive')
       
       call hist_addfld1d (fname='C14_PFT_FIRE_CLOSS', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total pft-level fire C loss', &
            ptr_pft=clm3%g%l%c%p%pc14f%pft_fire_closs)
    endif

    !-------------------------------
    ! C flux variables - native to column 
    !-------------------------------
    ! add history fields for all CLAMP CN variables

    call hist_addfld1d (fname='CWDC_HR', units='gC/m^2/s', &
         avgflag='A', long_name='coarse woody debris C heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%cwdc_hr)

    call hist_addfld1d (fname='CWDC_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='coarse woody debris C loss', &
         ptr_col=clm3%g%l%c%ccf%cwdc_loss)

    call hist_addfld1d (fname='LITTERC_HR', units='gC/m^2/s', &
         avgflag='A', long_name='litter C heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%lithr)

    call hist_addfld1d (fname='LITTERC_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='litter C loss', &
         ptr_col=clm3%g%l%c%ccf%litterc_loss)

    call hist_addfld1d (fname='SOILC_HR', units='gC/m^2/s', &
         avgflag='A', long_name='soil C heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%somhr)

    call hist_addfld1d (fname='SOILC_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='soil C loss', &
         ptr_col=clm3%g%l%c%ccf%somhr)
   
! F. Li and S. Levis
   call hist_addfld1d (fname='LF_CONV_CFLUX', units='gC/m^2/s', &
         avgflag='A', long_name='conversion carbon due to BET and BDT area decreasing', &
         ptr_col=clm3%g%l%c%ccf%lf_conv_cflux)   
   
   call hist_addfld1d (fname='SOMC_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='C loss due to peat burning', &
         ptr_col=clm3%g%l%c%ccf%somc_fire)

 

  
    do k = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
          data1dptr => clm3%g%l%c%ccf%m_decomp_cpools_to_fire(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
          call hist_addfld1d (fname=fieldname, units='gC/m^2/s',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          if ( nlevdecomp_full .gt. 1 ) then
             data2dptr => clm3%g%l%c%ccf%m_decomp_cpools_to_fire_vr(:,:,k)
             fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'//trim(vr_suffix)
             longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
             call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s', type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       endif

       ! decomposition k
       data2dptr => clm3%g%l%c%ccf%decomp_k(:,:,k)
       fieldname = 'K_'//trim(decomp_cascade_con%decomp_pool_name_history(k))
       longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' potential loss coefficient'
       call hist_addfld_decomp (fname=fieldname, units='1/s',  type2d='levdcmp', &
            avgflag='A', long_name=longname, &
            ptr_col=data2dptr, default='inactive')
    end do
    
    do l = 1, ndecomp_cascade_transitions
       ! output the vertically integrated fluxes only as  default
       !-- HR fluxes (none from CWD)
       if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
          data1dptr => clm3%g%l%c%ccf%decomp_cascade_hr(:,l)
          ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
          ii = 0
          do jj = 1, ndecomp_cascade_transitions
             if ( decomp_cascade_con%cascade_donor_pool(jj) .eq. decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
          end do
          if ( ii .eq. 1 ) then
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR'
          else
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_'//&
             trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))
          endif
          longname =  'Het. Resp. from '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
          call hist_addfld1d (fname=fieldname, units='gC/m^2/s',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       endif
       !-- transfer fluxes (none from terminal pool, if present)
       if ( decomp_cascade_con%cascade_receiver_pool(l) .ne. 0 ) then
          data1dptr => clm3%g%l%c%ccf%decomp_cascade_ctransfer(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'C_TO_'//&
               trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'C'
          longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
               ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
          call hist_addfld1d (fname=fieldname, units='gC/m^2/s', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       endif
       
       ! output the vertically resolved fluxes 
       if ( nlevdecomp_full .gt. 1 ) then  
          !-- HR fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             data2dptr => clm3%g%l%c%ccf%decomp_cascade_hr_vr(:,:,l)
             ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
             ii = 0
             do jj = 1, ndecomp_cascade_transitions
                if ( decomp_cascade_con%cascade_donor_pool(jj) .eq. decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
             end do
             if ( ii .eq. 1 ) then
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR'//trim(vr_suffix)
             else
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_'//&
                  trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))//trim(vr_suffix)
             endif
             longname =  'Het. Resp. from '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
             call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
          !-- transfer fluxes (none from terminal pool, if present)
          if ( decomp_cascade_con%cascade_receiver_pool(l) .ne. 0 ) then
             data2dptr => clm3%g%l%c%ccf%decomp_cascade_ctransfer_vr(:,:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'C_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'C'//trim(vr_suffix)
             longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
             call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       end if
    end do

    call hist_addfld_decomp (fname='T_SCALAR', units='unitless',  type2d='levdcmp', &
         avgflag='A', long_name='temperature inhibition of decomposition', &
         ptr_col=clm3%g%l%c%ccf%t_scalar)

    call hist_addfld_decomp (fname='W_SCALAR', units='unitless',  type2d='levdcmp', &
         avgflag='A', long_name='Moisture (dryness) inhibition of decomposition', &
         ptr_col=clm3%g%l%c%ccf%w_scalar)

    call hist_addfld_decomp (fname='O_SCALAR', units='unitless', type2d='levdcmp', &
         avgflag='A', long_name='fraction by which decomposition is reduced due to anoxia', &
         ptr_col=clm3%g%l%c%ccf%o_scalar)

    call hist_addfld1d (fname='SOM_C_LEACHED', units='gC/m^2/s', &
         avgflag='A', long_name='total flux of C from SOM pools due to leaching', &
         ptr_col=clm3%g%l%c%ccf%som_c_leached)!, default='inactive')
    
    do k = 1, ndecomp_pools
       if ( .not. decomp_cascade_con%is_cwd(k) ) then
          data1dptr => clm3%g%l%c%ccf%decomp_cpools_leached(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_LEACHING'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C leaching loss'
          call hist_addfld1d (fname=fieldname, units='gC/m^2/s', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)!, default='inactive')
          
          data2dptr => clm3%g%l%c%ccf%decomp_cpools_transport_tendency(:,:,k)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TNDNCY_VERT_TRANSPORT'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C tendency due to vertical transport'
          call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr, default='inactive')
       endif
    end do


    call hist_addfld1d (fname='LITHR', units='gC/m^2/s', &
         avgflag='A', long_name='litter heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%lithr)

    call hist_addfld1d (fname='SOMHR', units='gC/m^2/s', &
         avgflag='A', long_name='soil organic matter heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%somhr)

    if ( nlevdecomp_full .gt. 1 ) then
       call hist_addfld2d (fname='HR_vr', units='gC/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='total vertically resolved heterotrophic respiration', &
            ptr_col=clm3%g%l%c%ccf%hr_vr)
    endif

    call hist_addfld1d (fname='HR', units='gC/m^2/s', &
         avgflag='A', long_name='total heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%hr)

    call hist_addfld1d (fname='SR', units='gC/m^2/s', &
         avgflag='A', long_name='total soil respiration (HR + root resp)', &
         ptr_col=clm3%g%l%c%ccf%sr)

    call hist_addfld1d (fname='ER', units='gC/m^2/s', &
         avgflag='A', long_name='total ecosystem respiration, autotrophic + heterotrophic', &
         ptr_col=clm3%g%l%c%ccf%er)

    call hist_addfld1d (fname='LITFIRE', units='gC/m^2/s', &
         avgflag='A', long_name='litter fire losses', &
         ptr_col=clm3%g%l%c%ccf%litfire, default='inactive')

    call hist_addfld1d (fname='SOMFIRE', units='gC/m^2/s', &
         avgflag='A', long_name='soil organic matter fire losses', &
         ptr_col=clm3%g%l%c%ccf%somfire, default='inactive')

    call hist_addfld1d (fname='TOTFIRE', units='gC/m^2/s', &
         avgflag='A', long_name='total ecosystem fire losses', &
         ptr_col=clm3%g%l%c%ccf%totfire, default='inactive')

    call hist_addfld1d (fname='NEP', units='gC/m^2/s', &
         avgflag='A', long_name='net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink', &
         ptr_col=clm3%g%l%c%ccf%nep)

    call hist_addfld1d (fname='NBP', units='gC/m^2/s', &
         avgflag='A', long_name='net biome production, includes fire, landuse, and harvest flux, positive for sink', &
         ptr_col=clm3%g%l%c%ccf%nbp)

    call hist_addfld1d (fname='NEE', units='gC/m^2/s', &
         avgflag='A', long_name='net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source', &
         ptr_col=clm3%g%l%c%ccf%nee)

    call hist_addfld1d (fname='COL_FIRE_CLOSS', units='gC/m^2/s', &
         avgflag='A', long_name='total column-level fire C loss for non-peat fires outside land-type converted region', &
         ptr_col=clm3%g%l%c%ccf%col_fire_closs)

    call hist_addfld1d (fname='DWT_SEEDC_TO_LEAF', units='gC/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level leaf', &
         ptr_col=clm3%g%l%c%ccf%dwt_seedc_to_leaf)

    call hist_addfld1d (fname='DWT_SEEDC_TO_DEADSTEM', units='gC/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level deadstem', &
         ptr_col=clm3%g%l%c%ccf%dwt_seedc_to_deadstem)

    call hist_addfld1d (fname='DWT_CONV_CFLUX', units='gC/m^2/s', &
         avgflag='A', long_name='conversion C flux (immediate loss to atm)', &
         ptr_col=clm3%g%l%c%ccf%dwt_conv_cflux)

    call hist_addfld1d (fname='DWT_PROD10C_GAIN', units='gC/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%ccf%dwt_prod10c_gain)

    call hist_addfld1d (fname='PROD10C_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='loss from 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%ccf%prod10c_loss)

    call hist_addfld1d (fname='DWT_PROD100C_GAIN', units='gC/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%ccf%dwt_prod100c_gain)

    call hist_addfld1d (fname='PROD100C_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='loss from 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%ccf%prod100c_loss)

    call hist_addfld_decomp (fname='DWT_FROOTC_TO_LITR_MET_C', units='gC/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_frootc_to_litr_met_c, default='inactive')

    call hist_addfld_decomp (fname='DWT_FROOTC_TO_LITR_CEL_C', units='gC/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_frootc_to_litr_cel_c, default='inactive')

    call hist_addfld_decomp (fname='DWT_FROOTC_TO_LITR_LIG_C', units='gC/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_frootc_to_litr_lig_c, default='inactive')

    call hist_addfld_decomp (fname='DWT_LIVECROOTC_TO_CWDC', units='gC/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_livecrootc_to_cwdc, default='inactive')

    call hist_addfld_decomp (fname='DWT_DEADCROOTC_TO_CWDC', units='gC/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_deadcrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='DWT_CLOSS', units='gC/m^2/s', &
         avgflag='A', long_name='total carbon loss from land cover conversion', &
         ptr_col=clm3%g%l%c%ccf%dwt_closs)

    call hist_addfld1d (fname='PRODUCT_CLOSS', units='gC/m^2/s', &
         avgflag='A', long_name='total carbon loss from wood product pools', &
         ptr_col=clm3%g%l%c%ccf%product_closs)

    call hist_addfld1d (fname='LAND_USE_FLUX', units='gC/m^2/s', &
         avgflag='A', long_name='total C emitted from land cover conversion and wood product pools', &
         ptr_col=clm3%g%l%c%ccf%landuseflux)

    call hist_addfld1d (fname='LAND_UPTAKE', units='gC/m^2/s', &
         avgflag='A', long_name='NEE minus LAND_USE_FLUX, negative for update', &
         ptr_col=clm3%g%l%c%ccf%landuptake)

    if ( use_c13 ) then
       !-------------------------------
       ! C13 flux variables - native to column 
       !-------------------------------
       
       
       do k = 1, ndecomp_pools
          if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
             data1dptr => clm3%g%l%c%cc13f%m_decomp_cpools_to_fire(:,k)
             fieldname = 'C13_M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'
             longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
             call hist_addfld1d (fname=fieldname, units='gC13/m^2',  &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr, default='inactive')
             
             if ( nlevdecomp_full .gt. 1 ) then
                data2dptr => clm3%g%l%c%cc13f%m_decomp_cpools_to_fire_vr(:,:,k)
                fieldname = 'C13_M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'//trim(vr_suffix)
                longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
                call hist_addfld_decomp (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                     avgflag='A', long_name=longname, &
                     ptr_col=data2dptr, default='inactive')
             end if
          endif
          
       end do
       
       do l = 1, ndecomp_cascade_transitions
          !-- HR fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             data2dptr => clm3%g%l%c%cc13f%decomp_cascade_hr_vr(:,:,l)
             ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
             ii = 0
             do jj = 1, ndecomp_cascade_transitions
                if ( decomp_cascade_con%cascade_donor_pool(jj) .eq. decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
             end do
             if ( ii .eq. 1 ) then
                fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR'//trim(vr_suffix)
             else
                fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_'//&
                     trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))//trim(vr_suffix)
             endif
             longname =  'C13 Het. Resp. from '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
             call hist_addfld_decomp (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
          !-- transfer fluxes (none from terminal pool, if present)
          if ( decomp_cascade_con%cascade_receiver_pool(l) .ne. 0 ) then
             data2dptr => clm3%g%l%c%cc13f%decomp_cascade_ctransfer_vr(:,:,l)
             fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'C_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'C'//trim(vr_suffix)
             longname =  'C13 decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
             call hist_addfld_decomp (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       end do
       
       
       call hist_addfld1d (fname='C13_LITHR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C litterfall to litter 3 C', &
            ptr_col=clm3%g%l%c%cc13f%lithr)
       
       call hist_addfld1d (fname='C13_SOMHR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 soil organic matter heterotrophic respiration', &
            ptr_col=clm3%g%l%c%cc13f%somhr)
       
       call hist_addfld1d (fname='C13_HR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total heterotrophic respiration', &
            ptr_col=clm3%g%l%c%cc13f%hr)
       
       call hist_addfld1d (fname='C13_SR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total soil respiration (HR + root resp)', &
            ptr_col=clm3%g%l%c%cc13f%sr)
       
       call hist_addfld1d (fname='C13_ER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total ecosystem respiration, autotrophic + heterotrophic', &
            ptr_col=clm3%g%l%c%cc13f%er)
       
       call hist_addfld1d (fname='C13_LITFIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 litter fire losses', &
            ptr_col=clm3%g%l%c%cc13f%litfire, default='inactive')
       
       call hist_addfld1d (fname='C13_SOMFIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 soil organic matter fire losses', &
            ptr_col=clm3%g%l%c%cc13f%somfire, default='inactive')
       
       call hist_addfld1d (fname='C13_TOTFIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total ecosystem fire losses', &
            ptr_col=clm3%g%l%c%cc13f%totfire, default='inactive')
       
       call hist_addfld1d (fname='C13_NEP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 net ecosystem production, excludes fire flux, positive for sink', &
            ptr_col=clm3%g%l%c%cc13f%nep)
       
       call hist_addfld1d (fname='C13_NEE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 net ecosystem exchange of carbon, includes fire flux, positive for source', &
            ptr_col=clm3%g%l%c%cc13f%nee)
       
       call hist_addfld1d (fname='C13_COL_FIRE_CLOSS', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total column-level fire C loss', &
            ptr_col=clm3%g%l%c%cc13f%col_fire_closs)
       
       call hist_addfld1d (fname='C13_DWT_SEEDC_TO_LEAF', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 seed source to PFT-level leaf', &
            ptr_col=clm3%g%l%c%cc13f%dwt_seedc_to_leaf)
       
       call hist_addfld1d (fname='C13_DWT_SEEDC_TO_DEADSTEM', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 seed source to PFT-level deadstem', &
            ptr_col=clm3%g%l%c%cc13f%dwt_seedc_to_deadstem)
       
       call hist_addfld1d (fname='C13_DWT_CONV_CFLUX', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 conversion C flux (immediate loss to atm)', &
            ptr_col=clm3%g%l%c%cc13f%dwt_conv_cflux)
       
       call hist_addfld1d (fname='C13_DWT_PROD10C_GAIN', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 addition to 10-yr wood product pool', &
            ptr_col=clm3%g%l%c%cc13f%dwt_prod10c_gain)
       
       call hist_addfld1d (fname='C13_PROD10C_LOSS', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 loss from 10-yr wood product pool', &
            ptr_col=clm3%g%l%c%cc13f%prod10c_loss)
       
       call hist_addfld1d (fname='C13_DWT_PROD100C_GAIN', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 addition to 100-yr wood product pool', &
            ptr_col=clm3%g%l%c%cc13f%dwt_prod100c_gain)
       
       call hist_addfld1d (fname='C13_PROD100C_LOSS', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 loss from 100-yr wood product pool', &
            ptr_col=clm3%g%l%c%cc13f%prod100c_loss)
       
       call hist_addfld_decomp (fname='C13_DWT_FROOTC_TO_LITR_MET_C', units='gC13/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C13 fine root to litter due to landcover change', &
            ptr_col=clm3%g%l%c%cc13f%dwt_frootc_to_litr_met_c, default='inactive')
       
       call hist_addfld_decomp (fname='C13_DWT_FROOTC_TO_LITR_CEL_C', units='gC13/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C13 fine root to litter due to landcover change', &
            ptr_col=clm3%g%l%c%cc13f%dwt_frootc_to_litr_cel_c, default='inactive')
       
       call hist_addfld_decomp (fname='C13_DWT_FROOTC_TO_LITR_LIG_C', units='gC13/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C13 fine root to litter due to landcover change', &
            ptr_col=clm3%g%l%c%cc13f%dwt_frootc_to_litr_lig_c, default='inactive')
       
       call hist_addfld_decomp (fname='C13_DWT_LIVECROOTC_TO_CWDC', units='gC13/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C13 live coarse root to CWD due to landcover change', &
            ptr_col=clm3%g%l%c%cc13f%dwt_livecrootc_to_cwdc, default='inactive')
       
       call hist_addfld_decomp (fname='C13_DWT_DEADCROOTC_TO_CWDC', units='gC13/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C13 dead coarse root to CWD due to landcover change', &
            ptr_col=clm3%g%l%c%cc13f%dwt_deadcrootc_to_cwdc, default='inactive')
       
       call hist_addfld1d (fname='C13_DWT_CLOSS', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total carbon loss from land cover conversion', &
            ptr_col=clm3%g%l%c%cc13f%dwt_closs)
       
       call hist_addfld1d (fname='C13_PRODUCT_CLOSS', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total carbon loss from wood product pools', &
            ptr_col=clm3%g%l%c%cc13f%product_closs)
    endif

    if ( use_c14 ) then
       !-------------------------------
       ! C14 flux variables - native to column 
       !-------------------------------
       

       
       do k = 1, ndecomp_pools
          if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
             data1dptr => clm3%g%l%c%cc14f%m_decomp_cpools_to_fire(:,k)
             fieldname = 'C14_M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'
             longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
             call hist_addfld1d (fname=fieldname, units='gC14/m^2',  &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr, default='inactive')
             
             if ( nlevdecomp_full .gt. 1 ) then
                data2dptr => clm3%g%l%c%cc14f%m_decomp_cpools_to_fire_vr(:,:,k)
                fieldname = 'C14_M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'//trim(vr_suffix)
                longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
                call hist_addfld_decomp (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                     avgflag='A', long_name=longname, &
                     ptr_col=data2dptr, default='inactive')
             end if
          endif
          
       end do
       
       do l = 1, ndecomp_cascade_transitions
          !-- HR fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             data2dptr => clm3%g%l%c%cc14f%decomp_cascade_hr_vr(:,:,l)
             ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
             ii = 0
             do jj = 1, ndecomp_cascade_transitions
                if ( decomp_cascade_con%cascade_donor_pool(jj) .eq. decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
             end do
             if ( ii .eq. 1 ) then
                fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR'//trim(vr_suffix)
             else
                fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_'//&
                     trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))//trim(vr_suffix)
             endif
             longname =  'C14 Het. Resp. from '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
             call hist_addfld_decomp (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
          !-- transfer fluxes (none from terminal pool, if present)
          if ( decomp_cascade_con%cascade_receiver_pool(l) .ne. 0 ) then
             data2dptr => clm3%g%l%c%cc14f%decomp_cascade_ctransfer_vr(:,:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'C_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'C'//trim(vr_suffix)
             longname =  'C14 decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
             call hist_addfld_decomp (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       end do
       
       
       call hist_addfld1d (fname='C14_LITHR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C litterfall to litter 3 C', &
            ptr_col=clm3%g%l%c%cc14f%lithr)
       
       call hist_addfld1d (fname='C14_SOMHR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 soil organic matter heterotrophic respiration', &
            ptr_col=clm3%g%l%c%cc14f%somhr)
       
       call hist_addfld1d (fname='C14_HR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total heterotrophic respiration', &
            ptr_col=clm3%g%l%c%cc14f%hr)
       
       call hist_addfld1d (fname='C14_SR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total soil respiration (HR + root resp)', &
            ptr_col=clm3%g%l%c%cc14f%sr)
       
       call hist_addfld1d (fname='C14_ER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total ecosystem respiration, autotrophic + heterotrophic', &
            ptr_col=clm3%g%l%c%cc14f%er)
       
       call hist_addfld1d (fname='C14_LITFIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 litter fire losses', &
            ptr_col=clm3%g%l%c%cc14f%litfire, default='inactive')
       
       call hist_addfld1d (fname='C14_SOMFIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 soil organic matter fire losses', &
            ptr_col=clm3%g%l%c%cc14f%somfire, default='inactive')
       
       call hist_addfld1d (fname='C14_TOTFIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total ecosystem fire losses', &
            ptr_col=clm3%g%l%c%cc14f%totfire, default='inactive')
       
       call hist_addfld1d (fname='C14_NEP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 net ecosystem production, excludes fire flux, positive for sink', &
            ptr_col=clm3%g%l%c%cc14f%nep)
       
       call hist_addfld1d (fname='C14_NEE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 net ecosystem exchange of carbon, includes fire flux, positive for source', &
            ptr_col=clm3%g%l%c%cc14f%nee)
       
       call hist_addfld1d (fname='C14_COL_FIRE_CLOSS', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total column-level fire C loss', &
            ptr_col=clm3%g%l%c%cc14f%col_fire_closs)
       
       call hist_addfld1d (fname='C14_DWT_SEEDC_TO_LEAF', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 seed source to PFT-level leaf', &
            ptr_col=clm3%g%l%c%cc14f%dwt_seedc_to_leaf)
       
       call hist_addfld1d (fname='C14_DWT_SEEDC_TO_DEADSTEM', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 seed source to PFT-level deadstem', &
            ptr_col=clm3%g%l%c%cc14f%dwt_seedc_to_deadstem)
       
       call hist_addfld1d (fname='C14_DWT_CONV_CFLUX', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 conversion C flux (immediate loss to atm)', &
            ptr_col=clm3%g%l%c%cc14f%dwt_conv_cflux)
       
       call hist_addfld1d (fname='C14_DWT_PROD10C_GAIN', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 addition to 10-yr wood product pool', &
            ptr_col=clm3%g%l%c%cc14f%dwt_prod10c_gain)
       
       call hist_addfld1d (fname='C14_PROD10C_LOSS', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 loss from 10-yr wood product pool', &
            ptr_col=clm3%g%l%c%cc14f%prod10c_loss)
       
       call hist_addfld1d (fname='C14_DWT_PROD100C_GAIN', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 addition to 100-yr wood product pool', &
            ptr_col=clm3%g%l%c%cc14f%dwt_prod100c_gain)
       
       call hist_addfld1d (fname='C14_PROD100C_LOSS', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 loss from 100-yr wood product pool', &
            ptr_col=clm3%g%l%c%cc14f%prod100c_loss)
       
       call hist_addfld_decomp (fname='C14_DWT_FROOTC_TO_LITR_MET_C', units='gC14/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C14 fine root to litter due to landcover change', &
            ptr_col=clm3%g%l%c%cc14f%dwt_frootc_to_litr_met_c, default='inactive')
       
       call hist_addfld_decomp (fname='C14_DWT_FROOTC_TO_LITR_CEL_C', units='gC14/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C14 fine root to litter due to landcover change', &
            ptr_col=clm3%g%l%c%cc14f%dwt_frootc_to_litr_cel_c, default='inactive')
       
       call hist_addfld_decomp (fname='C14_DWT_FROOTC_TO_LITR_LIG_C', units='gC14/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C14 fine root to litter due to landcover change', &
            ptr_col=clm3%g%l%c%cc14f%dwt_frootc_to_litr_lig_c, default='inactive')
       
       call hist_addfld_decomp (fname='C14_DWT_LIVECROOTC_TO_CWDC', units='gC14/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C14 live coarse root to CWD due to landcover change', &
            ptr_col=clm3%g%l%c%cc14f%dwt_livecrootc_to_cwdc, default='inactive')
       
       call hist_addfld_decomp (fname='C14_DWT_DEADCROOTC_TO_CWDC', units='gC14/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name='C14 dead coarse root to CWD due to landcover change', &
            ptr_col=clm3%g%l%c%cc14f%dwt_deadcrootc_to_cwdc, default='inactive')
       
       call hist_addfld1d (fname='C14_DWT_CLOSS', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total carbon loss from land cover conversion', &
            ptr_col=clm3%g%l%c%cc14f%dwt_closs)
       
       call hist_addfld1d (fname='C14_PRODUCT_CLOSS', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total carbon loss from wood product pools', &
            ptr_col=clm3%g%l%c%cc14f%product_closs)
    endif

    !-------------------------------
    ! N flux variables - native to PFT
    !-------------------------------

    call hist_addfld1d (fname='M_LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_RETRANSN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_retransn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N fire loss ', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_RETRANSN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_retransn_to_fire, default='inactive')

    call hist_addfld1d (fname='LEAFN_XFER_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%leafn_xfer_to_leafn, default='inactive')

    call hist_addfld1d (fname='FROOTN_XFER_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%frootn_xfer_to_frootn, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_XFER_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%livestemn_xfer_to_livestemn, default='inactive')

    call hist_addfld1d (fname='DEADSTEMN_XFER_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%deadstemn_xfer_to_deadstemn, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_XFER_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%livecrootn_xfer_to_livecrootn, default='inactive')

    call hist_addfld1d (fname='DEADCROOTN_XFER_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%deadcrootn_xfer_to_deadcrootn, default='inactive')

    call hist_addfld1d (fname='LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N litterfall', &
         ptr_pft=clm3%g%l%c%p%pnf%leafn_to_litter, default='inactive')

    call hist_addfld1d (fname='LEAFN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N to retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pnf%leafn_to_retransn, default='inactive')

    call hist_addfld1d (fname='FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N litterfall', &
         ptr_pft=clm3%g%l%c%p%pnf%frootn_to_litter, default='inactive')

    call hist_addfld1d (fname='RETRANSN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of retranslocated N', &
         ptr_pft=clm3%g%l%c%p%pnf%retransn_to_npool)

    call hist_addfld1d (fname='SMINN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of soil mineral N uptake', &
         ptr_pft=clm3%g%l%c%p%pnf%sminn_to_npool)

    call hist_addfld1d (fname='NPOOL_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_leafn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LEAFN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_leafn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_frootn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_FROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_frootn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_livestemn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_livestemn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_deadstemn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_deadstemn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_livecrootn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_livecrootn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_deadcrootn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_deadcrootn_storage, default='inactive')

    call hist_addfld1d (fname='LEAFN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%leafn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='FROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%frootn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%livestemn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='DEADSTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%deadstemn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%livecrootn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='DEADCROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%deadcrootn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N turnover', &
         ptr_pft=clm3%g%l%c%p%pnf%livestemn_to_deadstemn, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N to retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pnf%livestemn_to_retransn, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N turnover', &
         ptr_pft=clm3%g%l%c%p%pnf%livecrootn_to_deadcrootn, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N to retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pnf%livecrootn_to_retransn, default='inactive')

    call hist_addfld1d (fname='NDEPLOY', units='gN/m^2/s', &
         avgflag='A', long_name='total N deployed in new growth', &
         ptr_pft=clm3%g%l%c%p%pnf%ndeploy)

    call hist_addfld1d (fname='WOOD_HARVESTN', units='gN/m^2/s', &
         avgflag='A', long_name='wood harvest N (to product pools)', &
         ptr_pft=clm3%g%l%c%p%pnf%wood_harvestn)

    call hist_addfld1d (fname='PFT_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total pft-level fire N loss', &
         ptr_pft=clm3%g%l%c%p%pnf%pft_fire_nloss)

    if (crop_prog) then
       call hist_addfld1d (fname='FERT', units='gN/m^2/s', &
            avgflag='A', long_name='fertilizer added', &
            ptr_pft=clm3%g%l%c%p%pnf%fert)

       call hist_addfld1d (fname='SOYFIXN', units='gN/m^2/s', &
            avgflag='A', long_name='soybean fixation', &
            ptr_pft=clm3%g%l%c%p%pnf%soyfixn)
    end if

    !-------------------------------
    ! N flux variables - native to column
    !-------------------------------

    call hist_addfld1d (fname='NDEP_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='atmospheric N deposition to soil mineral N', &
         ptr_col=clm3%g%l%c%cnf%ndep_to_sminn)

    call hist_addfld1d (fname='NFIX_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='symbiotic/asymbiotic N fixation to soil mineral N', &
         ptr_col=clm3%g%l%c%cnf%nfix_to_sminn)

    do k = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
          data1dptr => clm3%g%l%c%cnf%m_decomp_npools_to_fire(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TO_FIRE'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N fire loss'
          call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')
       
          if ( nlevdecomp_full .gt. 1 ) then
             data2dptr => clm3%g%l%c%cnf%m_decomp_npools_to_fire_vr(:,:,k)
             fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TO_FIRE'//trim(vr_suffix)
             longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N fire loss'
             call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

       endif
    end do
    
    do l = 1, ndecomp_cascade_transitions
       ! vertically integrated fluxes
       !-- mineralization/immobilization fluxes (none from CWD)
       if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
          data1dptr => clm3%g%l%c%cnf%decomp_cascade_sminn_flux(:,l)
          if ( decomp_cascade_con%cascade_receiver_pool(l) .ne. 0 ) then
             fieldname = 'SMINN_TO_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N_'//&
                  trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))
             longname =  'mineral N flux for decomp. of '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                  'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
          else
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'N_TO_SMINN'
             longname =  'mineral N flux for decomp. of '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
          endif
             call hist_addfld1d (fname=fieldname, units='gN/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       endif
       !-- transfer fluxes (none from terminal pool, if present)
       if ( decomp_cascade_con%cascade_receiver_pool(l) .ne. 0 ) then
          data1dptr => clm3%g%l%c%cnf%decomp_cascade_ntransfer(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'N_TO_'//&
               trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N'
          longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
               ' N to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' N'
          call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       endif

       ! vertically resolved fluxes
       if ( nlevdecomp_full .gt. 1 ) then
          !-- mineralization/immobilization fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             data2dptr => clm3%g%l%c%cnf%decomp_cascade_sminn_flux_vr(:,:,l)
             if ( decomp_cascade_con%cascade_receiver_pool(l) .ne. 0 ) then
                fieldname = 'SMINN_TO_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N_'//&
                     trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))//trim(vr_suffix)
                longname =  'mineral N flux for decomp. of '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                     'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
             else
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'N_TO_SMINN'//trim(vr_suffix)
                longname =  'mineral N flux for decomp. of '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
             endif
             call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
          !-- transfer fluxes (none from terminal pool, if present)
          if ( decomp_cascade_con%cascade_receiver_pool(l) .ne. 0 ) then
             data2dptr => clm3%g%l%c%cnf%decomp_cascade_ntransfer_vr(:,:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'N_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N'//trim(vr_suffix)
             longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  ' N to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' N'
             call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
          
       endif

    end do

    call hist_addfld1d (fname='DENIT', units='gN/m^2/s', &
         avgflag='A', long_name='total rate of denitrification', &
         ptr_col=clm3%g%l%c%cnf%denit)

    call hist_addfld1d (fname='SOM_N_LEACHED', units='gN/m^2/s', &
         avgflag='A', long_name='total flux of N from SOM pools due to leaching', &
         ptr_col=clm3%g%l%c%cnf%som_n_leached, default='inactive')
    
    do k = 1, ndecomp_pools
       if ( .not. decomp_cascade_con%is_cwd(k) ) then
          data1dptr => clm3%g%l%c%cnf%decomp_npools_leached(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TO_LEACHING'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N leaching loss'
          call hist_addfld1d (fname=fieldname, units='gN/m^2/s', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')
          
          data2dptr => clm3%g%l%c%cnf%decomp_npools_transport_tendency(:,:,k)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TNDNCY_VERT_TRANSPORT'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N tendency due to vertical transport'
          call hist_addfld_decomp (fname=fieldname, units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif
    end do


#ifndef NITRIF_DENITRIF
    do l = 1, ndecomp_cascade_transitions
       !-- denitrification fluxes (none from CWD)
       if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
          data1dptr => clm3%g%l%c%cnf%sminn_to_denit_decomp_cascade(:,l)
          fieldname = 'SMINN_TO_DENIT_'//trim(decomp_cascade_con%cascade_step_name(l))
          longname =  'denitrification for decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
               'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
          call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       endif

       if ( nlevdecomp_full .gt. 1 ) then       
          !-- denitrification fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             data2dptr => clm3%g%l%c%cnf%sminn_to_denit_decomp_cascade_vr(:,:,l)
             fieldname = 'SMINN_TO_DENIT_'//trim(decomp_cascade_con%cascade_step_name(l))//trim(vr_suffix)
             longname =  'denitrification for decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
             call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       endif

    end do

    call hist_addfld1d (fname='SMINN_TO_DENIT_EXCESS', units='gN/m^2/s',  &
         avgflag='A', long_name='denitrification from excess mineral N pool', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_excess, default='inactive')

    call hist_addfld1d (fname='SMINN_LEACHED', units='gN/m^2/s',   &
         avgflag='A', long_name='soil mineral N pool loss to leaching', &
         ptr_col=clm3%g%l%c%cnf%sminn_leached)

    if ( nlevdecomp_full .gt. 1 ) then  
       call hist_addfld_decomp (fname='SMINN_TO_DENIT_EXCESS'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='denitrification from excess mineral N pool', &
            ptr_col=clm3%g%l%c%cnf%sminn_to_denit_excess_vr, default='inactive')   

       call hist_addfld_decomp (fname='SMINN_LEACHED'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral N pool loss to leaching', &
            ptr_col=clm3%g%l%c%cnf%sminn_leached_vr, default='inactive')
    endif
    
#else

    call hist_addfld1d (fname='F_NIT', units='gN/m^2/s',  &
         avgflag='A', long_name='nitrification flux', &
         ptr_col=clm3%g%l%c%cnf%f_nit)
    
    call hist_addfld1d (fname='F_DENIT', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification flux', &
         ptr_col=clm3%g%l%c%cnf%f_denit)
    
    call hist_addfld1d (fname='POT_F_NIT', units='gN/m^2/s', &
         avgflag='A', long_name='potential nitrification flux', &
         ptr_col=clm3%g%l%c%cnf%pot_f_nit)
    
    call hist_addfld1d (fname='POT_F_DENIT', units='gN/m^2/s', &
         avgflag='A', long_name='potential denitrification flux', &
         ptr_col=clm3%g%l%c%cnf%pot_f_denit)
    
    call hist_addfld1d (fname='SMIN_NO3_LEACHED', units='gN/m^2/s', &
         avgflag='A', long_name='soil NO3 pool loss to leaching', &
         ptr_col=clm3%g%l%c%cnf%smin_no3_leached)

    call hist_addfld1d (fname='SMIN_NO3_RUNOFF', units='gN/m^2/s', &
         avgflag='A', long_name='soil NO3 pool loss to runoff', &
         ptr_col=clm3%g%l%c%cnf%smin_no3_runoff)

    if ( nlevdecomp_full .gt. 1 ) then 
       call hist_addfld_decomp (fname='F_NIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='nitrification flux', &
            ptr_col=clm3%g%l%c%cnf%f_nit_vr)
       
       call hist_addfld_decomp (fname='F_DENIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='denitrification flux', &
            ptr_col=clm3%g%l%c%cnf%f_denit_vr)
       
       call hist_addfld_decomp (fname='POT_F_NIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='potential nitrification flux', &
            ptr_col=clm3%g%l%c%cnf%pot_f_nit_vr, default='inactive')
       
       call hist_addfld_decomp (fname='POT_F_DENIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='potential denitrification flux', &
            ptr_col=clm3%g%l%c%cnf%pot_f_denit_vr, default='inactive')

       call hist_addfld_decomp (fname='SMIN_NO3_LEACHED'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='soil NO3 pool loss to leaching', &
            ptr_col=clm3%g%l%c%cnf%smin_no3_leached_vr, default='inactive')

       call hist_addfld_decomp (fname='SMIN_NO3_RUNOFF'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='soil NO3 pool loss to runoff', &
            ptr_col=clm3%g%l%c%cnf%smin_no3_runoff_vr, default='inactive')
    endif

    call hist_addfld_decomp (fname='n2_n2o_ratio_denit', units='gN/gN', type2d='levdcmp', &
         avgflag='A', long_name='n2_n2o_ratio_denit', &
         ptr_col=clm3%g%l%c%cnf%n2_n2o_ratio_denit_vr, default='inactive')


    call hist_addfld_decomp (fname='ACTUAL_IMMOB_NO3', units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='immobilization of NO3', &
         ptr_col=clm3%g%l%c%cnf%actual_immob_no3_vr, default='inactive')

    call hist_addfld_decomp (fname='ACTUAL_IMMOB_NH4', units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='immobilization of NH4', &
         ptr_col=clm3%g%l%c%cnf%actual_immob_nh4_vr, default='inactive')

    call hist_addfld_decomp (fname='SMIN_NO3_TO_PLANT', units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='plant uptake of NO3', &
         ptr_col=clm3%g%l%c%cnf%smin_no3_to_plant_vr, default='inactive')

    call hist_addfld_decomp (fname='SMIN_NH4_TO_PLANT', units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='plant uptake of NH4', &
         ptr_col=clm3%g%l%c%cnf%smin_nh4_to_plant_vr, default='inactive')

    call hist_addfld2d (fname='watfc', units='m^3/m^3', type2d='levgrnd', &
         avgflag='A', long_name='water field capacity', &
         ptr_pft=clm3%g%l%c%cps%watfc, default='inactive')

    call hist_addfld2d (fname='watsat', units='m^3/m^3', type2d='levgrnd', &
         avgflag='A', long_name='water saturated', &
         ptr_pft=clm3%g%l%c%cps%watsat, default='inactive')

    call hist_addfld2d (fname='bsw', units='unitless', type2d='levgrnd', &
         avgflag='A', long_name='clapp and hornberger B', &
         ptr_pft=clm3%g%l%c%cps%bsw, default='inactive')

    call hist_addfld_decomp (fname='SMIN_NO3_MASSDENS', units='ugN/cm^3 soil', type2d='levdcmp', &
         avgflag='A', long_name='SMIN_NO3_MASSDENS', &
         ptr_col=clm3%g%l%c%cnf%smin_no3_massdens_vr, default='inactive')

    call hist_addfld_decomp (fname='K_NITR_T', units='unitless', type2d='levdcmp', &
         avgflag='A', long_name='K_NITR_T', &
         ptr_col=clm3%g%l%c%cnf%k_nitr_t_vr, default='inactive')

    call hist_addfld_decomp (fname='K_NITR_PH', units='unitless', type2d='levdcmp', &
         avgflag='A', long_name='K_NITR_PH', &
         ptr_col=clm3%g%l%c%cnf%k_nitr_ph_vr, default='inactive')

    call hist_addfld_decomp (fname='K_NITR_H2O', units='unitless', type2d='levdcmp', &
         avgflag='A', long_name='K_NITR_H2O', &
         ptr_col=clm3%g%l%c%cnf%k_nitr_h2o_vr, default='inactive')

    call hist_addfld_decomp (fname='K_NITR', units='1/s', type2d='levdcmp', &
         avgflag='A', long_name='K_NITR', &
         ptr_col=clm3%g%l%c%cnf%k_nitr_vr, default='inactive')

    call hist_addfld_decomp (fname='WFPS', units='percent', type2d='levdcmp', &
         avgflag='A', long_name='WFPS', &
         ptr_col=clm3%g%l%c%cnf%wfps_vr, default='inactive')

    call hist_addfld_decomp (fname='FMAX_DENIT_CARBONSUBSTRATE', units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='FMAX_DENIT_CARBONSUBSTRATE', &
         ptr_col=clm3%g%l%c%cnf%fmax_denit_carbonsubstrate_vr, default='inactive')

    call hist_addfld_decomp (fname='FMAX_DENIT_NITRATE', units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='FMAX_DENIT_NITRATE', &
         ptr_col=clm3%g%l%c%cnf%fmax_denit_nitrate_vr, default='inactive')

    call hist_addfld_decomp (fname='F_DENIT_BASE', units='gN/m^3/s', type2d='levdcmp', &
         avgflag='A', long_name='F_DENIT_BASE', &
         ptr_col=clm3%g%l%c%cnf%f_denit_base_vr, default='inactive')

    call hist_addfld_decomp (fname='diffus', units='m^2/s', type2d='levdcmp', &
         avgflag='A', long_name='diffusivity', &
         ptr_col=clm3%g%l%c%cnf%diffus, default='inactive')

    call hist_addfld_decomp (fname='ratio_k1', units='none', type2d='levdcmp', &
         avgflag='A', long_name='ratio_k1', &
         ptr_col=clm3%g%l%c%cnf%ratio_k1, default='inactive')

    call hist_addfld_decomp (fname='ratio_no3_co2', units='ratio', type2d='levdcmp', &
         avgflag='A', long_name='ratio_no3_co2', &
         ptr_col=clm3%g%l%c%cnf%ratio_no3_co2, default='inactive')

    call hist_addfld_decomp (fname='soil_co2_prod', units='ug C / g soil / day', type2d='levdcmp', &
         avgflag='A', long_name='soil_co2_prod', &
         ptr_col=clm3%g%l%c%cnf%soil_co2_prod, default='inactive')

    call hist_addfld_decomp (fname='fr_WFPS', units='fraction', type2d='levdcmp', &
         avgflag='A', long_name='fr_WFPS', &
         ptr_col=clm3%g%l%c%cnf%fr_WFPS, default='inactive')

    call hist_addfld_decomp (fname='soil_bulkdensity', units='kg/m3', type2d='levdcmp', &
         avgflag='A', long_name='soil_bulkdensity', &
         ptr_col=clm3%g%l%c%cnf%soil_bulkdensity, default='inactive')

    call hist_addfld_decomp (fname='anaerobic_frac', units='m3/m3', type2d='levdcmp', &
         avgflag='A', long_name='anaerobic_frac', &
         ptr_col=clm3%g%l%c%cnf%anaerobic_frac, default='inactive')

    call hist_addfld_decomp (fname='r_psi', units='m', type2d='levdcmp', &
         avgflag='A', long_name='r_psi', &
         ptr_col=clm3%g%l%c%cnf%r_psi, default='inactive')

#ifdef LCH4
    if (.not. hist_wrtch4diag) then
       call hist_addfld2d (fname='o2_decomp_depth_unsat', units='mol/m3/2', type2d='levgrnd', &
            avgflag='A', long_name='o2_decomp_depth_unsat', &
            ptr_col=clm3%g%l%c%cch4%o2_decomp_depth_unsat)
    end if
#endif
#endif

    if (nlevdecomp_full .gt. 1 ) then
       call hist_addfld_decomp (fname='POTENTIAL_IMMOB'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='potential N immobilization', &
            ptr_col=clm3%g%l%c%cnf%potential_immob_vr, default='inactive')
       
       call hist_addfld_decomp (fname='ACTUAL_IMMOB'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='actual N immobilization', &
            ptr_col=clm3%g%l%c%cnf%actual_immob_vr, default='inactive')
       
       call hist_addfld_decomp (fname='SMINN_TO_PLANT'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of soil mineral N', &
            ptr_col=clm3%g%l%c%cnf%sminn_to_plant_vr, default='inactive')
       
       call hist_addfld_decomp (fname='SUPPLEMENT_TO_SMINN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='supplemental N supply', &
            ptr_col=clm3%g%l%c%cnf%supplement_to_sminn_vr, default='inactive')
       
       call hist_addfld_decomp (fname='GROSS_NMIN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='gross rate of N mineralization', &
            ptr_col=clm3%g%l%c%cnf%gross_nmin_vr, default='inactive')
       
       call hist_addfld_decomp (fname='NET_NMIN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='net rate of N mineralization', &
            ptr_col=clm3%g%l%c%cnf%net_nmin_vr, default='inactive')
    end if
    
    call hist_addfld1d (fname='POTENTIAL_IMMOB', units='gN/m^2/s', &
         avgflag='A', long_name='potential N immobilization', &
         ptr_col=clm3%g%l%c%cnf%potential_immob)
    
    call hist_addfld1d (fname='ACTUAL_IMMOB', units='gN/m^2/s', &
         avgflag='A', long_name='actual N immobilization', &
         ptr_col=clm3%g%l%c%cnf%actual_immob)
    
    call hist_addfld1d (fname='SMINN_TO_PLANT', units='gN/m^2/s', &
         avgflag='A', long_name='plant uptake of soil mineral N', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_plant)
    
    call hist_addfld1d (fname='SUPPLEMENT_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='supplemental N supply', &
         ptr_col=clm3%g%l%c%cnf%supplement_to_sminn)
    
    call hist_addfld1d (fname='GROSS_NMIN', units='gN/m^2/s', &
         avgflag='A', long_name='gross rate of N mineralization', &
         ptr_col=clm3%g%l%c%cnf%gross_nmin)
    
    call hist_addfld1d (fname='NET_NMIN', units='gN/m^2/s', &
         avgflag='A', long_name='net rate of N mineralization', &
         ptr_col=clm3%g%l%c%cnf%net_nmin)
    

#ifdef NITRIF_DENITRIF
    call hist_addfld1d (fname='F_N2O_NIT', units='gN/m^2/s', &
         avgflag='A', long_name='nitrification N2O flux', &
         ptr_col=clm3%g%l%c%cnf%f_n2o_nit)

    call hist_addfld1d (fname='F_N2O_DENIT', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification N2O flux', &
         ptr_col=clm3%g%l%c%cnf%f_n2o_denit)

#endif


    call hist_addfld1d (fname='COL_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total column-level fire N loss', &
         ptr_col=clm3%g%l%c%cnf%col_fire_nloss)

    call hist_addfld1d (fname='DWT_SEEDN_TO_LEAF', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level leaf', &
         ptr_col=clm3%g%l%c%cnf%dwt_seedn_to_leaf)

    call hist_addfld1d (fname='DWT_SEEDN_TO_DEADSTEM', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level deadstem', &
         ptr_col=clm3%g%l%c%cnf%dwt_seedn_to_deadstem)

    call hist_addfld1d (fname='DWT_CONV_NFLUX', units='gN/m^2/s', &
         avgflag='A', long_name='conversion N flux (immediate loss to atm)', &
         ptr_col=clm3%g%l%c%cnf%dwt_conv_nflux)

    call hist_addfld1d (fname='DWT_PROD10N_GAIN', units='gN/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%cnf%dwt_prod10n_gain)

    call hist_addfld1d (fname='PROD10N_LOSS', units='gN/m^2/s', &
         avgflag='A', long_name='loss from 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%cnf%prod10n_loss)

    call hist_addfld1d (fname='DWT_PROD100N_GAIN', units='gN/m^2/s', &
         avgflag='A', long_name='addition to 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%cnf%dwt_prod100n_gain)

    call hist_addfld1d (fname='PROD100N_LOSS', units='gN/m^2/s', &
         avgflag='A', long_name='loss from 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%cnf%prod100n_loss)

    call hist_addfld1d (fname='PRODUCT_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total N loss from wood product pools', &
         ptr_col=clm3%g%l%c%cnf%product_nloss)

    call hist_addfld_decomp (fname='DWT_FROOTN_TO_LITR_MET_N', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_frootn_to_litr_met_n, default='inactive')

    call hist_addfld_decomp (fname='DWT_FROOTN_TO_LITR_CEL_N', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_frootn_to_litr_cel_n, default='inactive')

    call hist_addfld_decomp (fname='DWT_FROOTN_TO_LITR_LIG_N', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_frootn_to_litr_lig_n, default='inactive')

    call hist_addfld_decomp (fname='DWT_LIVECROOTN_TO_CWDN', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_livecrootn_to_cwdn, default='inactive')

    call hist_addfld_decomp (fname='DWT_DEADCROOTN_TO_CWDN', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_deadcrootn_to_cwdn, default='inactive')

    call hist_addfld1d (fname='DWT_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total nitrogen loss from landcover conversion', &
         ptr_col=clm3%g%l%c%cnf%dwt_nloss)

    if (crop_prog) then
       call hist_addfld1d (fname='FERT_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='fertilizer to soil mineral N', &
            ptr_col=clm3%g%l%c%cnf%fert_to_sminn)

       call hist_addfld1d (fname='SOYFIXN_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='Soybean fixation to soil mineral N', &
            ptr_col=clm3%g%l%c%cnf%soyfixn_to_sminn)
    end if

    !-------------------------------
    ! PFT ecophysiological variables (pepv) 
    !-------------------------------

    call hist_addfld1d (fname='DORMANT_FLAG', units='none', &
         avgflag='A', long_name='dormancy flag', &
         ptr_pft=clm3%g%l%c%p%pepv%dormant_flag, default='inactive')

    call hist_addfld1d (fname='DAYS_ACTIVE', units='days', &
         avgflag='A', long_name='number of days since last dormancy', &
         ptr_pft=clm3%g%l%c%p%pepv%days_active, default='inactive')

    call hist_addfld1d (fname='ONSET_FLAG', units='none', &
         avgflag='A', long_name='onset flag', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_flag, default='inactive')

    call hist_addfld1d (fname='ONSET_COUNTER', units='days', &
         avgflag='A', long_name='onset days counter', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_counter, default='inactive')

    call hist_addfld1d (fname='ONSET_GDDFLAG', units='none', &
         avgflag='A', long_name='onset flag for growing degree day sum', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_gddflag, default='inactive')

    call hist_addfld1d (fname='ONSET_FDD', units='C degree-days', &
         avgflag='A', long_name='onset freezing degree days counter', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_fdd, default='inactive')

    call hist_addfld1d (fname='ONSET_GDD', units='C degree-days', &
         avgflag='A', long_name='onset growing degree days', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_gdd, default='inactive')

    call hist_addfld1d (fname='ONSET_SWI', units='none', &
         avgflag='A', long_name='onset soil water index', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_swi, default='inactive')

    call hist_addfld1d (fname='OFFSET_FLAG', units='none', &
         avgflag='A', long_name='offset flag', &
         ptr_pft=clm3%g%l%c%p%pepv%offset_flag, default='inactive')

    call hist_addfld1d (fname='OFFSET_COUNTER', units='days', &
         avgflag='A', long_name='offset days counter', &
         ptr_pft=clm3%g%l%c%p%pepv%offset_counter, default='inactive')

    call hist_addfld1d (fname='OFFSET_FDD', units='C degree-days', &
         avgflag='A', long_name='offset freezing degree days counter', &
         ptr_pft=clm3%g%l%c%p%pepv%offset_fdd, default='inactive')

    call hist_addfld1d (fname='OFFSET_SWI', units='none', &
         avgflag='A', long_name='offset soil water index', &
         ptr_pft=clm3%g%l%c%p%pepv%offset_swi, default='inactive')

    if (crop_prog) then
       call hist_addfld1d (fname='FERT_COUNTER', units='seconds', &
            avgflag='A', long_name='time left to fertilize', &
            ptr_pft=clm3%g%l%c%p%pepv%fert_counter)
    end if

    call hist_addfld1d (fname='LGSF', units='proportion', &
         avgflag='A', long_name='long growing season factor', &
         ptr_pft=clm3%g%l%c%p%pepv%lgsf, default='inactive')

    call hist_addfld1d (fname='BGLFR', units='1/s', &
         avgflag='A', long_name='background litterfall rate', &
         ptr_pft=clm3%g%l%c%p%pepv%bglfr, default='inactive')

    call hist_addfld1d (fname='BGTR', units='1/s', &
         avgflag='A', long_name='background transfer growth rate', &
         ptr_pft=clm3%g%l%c%p%pepv%bgtr, default='inactive')

    call hist_addfld1d (fname='DAYL',  units='s', &
         avgflag='A', long_name='daylength', &
         ptr_pft=clm3%g%l%c%p%pepv%dayl, default='inactive')

    call hist_addfld1d (fname='PREV_DAYL', units='s', &
         avgflag='A', long_name='daylength from previous timestep', &
         ptr_pft=clm3%g%l%c%p%pepv%prev_dayl, default='inactive')

    call hist_addfld1d (fname='ANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average 2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pepv%annavg_t2m, default='inactive')

    call hist_addfld1d (fname='TEMPAVG_T2M', units='K', &
         avgflag='A', long_name='temporary average 2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pepv%tempavg_t2m, default='inactive')

    call hist_addfld1d (fname='INIT_GPP', units='gC/m^2/s', &
         avgflag='A', long_name='GPP flux before downregulation', &
         ptr_pft=clm3%g%l%c%p%pepv%gpp, default='inactive')

    call hist_addfld1d (fname='AVAILC', units='gC/m^2/s', &
         avgflag='A', long_name='C flux available for allocation', &
         ptr_pft=clm3%g%l%c%p%pepv%availc, default='inactive')

    call hist_addfld1d (fname='XSMRPOOL_RECOVER', units='gC/m^2/s', &
         avgflag='A', long_name='C flux assigned to recovery of negative xsmrpool', &
         ptr_pft=clm3%g%l%c%p%pepv%xsmrpool_recover)

    if ( use_c13 ) then
       call hist_addfld1d (fname='XSMRPOOL_C13RATIO', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) ratio for xsmrpool', &
            ptr_pft=clm3%g%l%c%p%pepv%xsmrpool_c13ratio, default='inactive')
    endif

    call hist_addfld1d (fname='ALLOC_PNOW', units='proportion', &
         avgflag='A', long_name='fraction of current allocation to display as new growth', &
         ptr_pft=clm3%g%l%c%p%pepv%alloc_pnow, default='inactive')

    call hist_addfld1d (fname='C_ALLOMETRY', units='none', &
         avgflag='A', long_name='C allocation index', &
         ptr_pft=clm3%g%l%c%p%pepv%c_allometry, default='inactive')

    call hist_addfld1d (fname='N_ALLOMETRY', units='none', &
         avgflag='A', long_name='N allocation index', &
         ptr_pft=clm3%g%l%c%p%pepv%n_allometry, default='inactive')

    call hist_addfld1d (fname='PLANT_NDEMAND', units='gN/m^2/s', &
         avgflag='A', long_name='N flux required to support initial GPP', &
         ptr_pft=clm3%g%l%c%p%pepv%plant_ndemand)

    call hist_addfld1d (fname='TEMPSUM_POTENTIAL_GPP', units='gC/m^2/yr', &
         avgflag='A', long_name='temporary annual sum of potential GPP', &
         ptr_pft=clm3%g%l%c%p%pepv%tempsum_potential_gpp, default='inactive')

    call hist_addfld1d (fname='ANNSUM_POTENTIAL_GPP', units='gN/m^2/yr', &
         avgflag='A', long_name='annual sum of potential GPP', &
         ptr_pft=clm3%g%l%c%p%pepv%annsum_potential_gpp, default='inactive')

    call hist_addfld1d (fname='TEMPMAX_RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='temporary annual max of retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pepv%tempmax_retransn, default='inactive')

    call hist_addfld1d (fname='ANNMAX_RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='annual max of retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pepv%annmax_retransn, default='inactive')

    call hist_addfld1d (fname='AVAIL_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='N flux available from retranslocation pool', &
         ptr_pft=clm3%g%l%c%p%pepv%avail_retransn, default='inactive')

    call hist_addfld1d (fname='PLANT_NALLOC', units='gN/m^2/s', &
         avgflag='A', long_name='total allocated N flux', &
         ptr_pft=clm3%g%l%c%p%pepv%plant_nalloc, default='inactive')

    call hist_addfld1d (fname='PLANT_CALLOC', units='gC/m^2/s', &
         avgflag='A', long_name='total allocated C flux', &
         ptr_pft=clm3%g%l%c%p%pepv%plant_calloc, default='inactive')

    call hist_addfld1d (fname='EXCESS_CFLUX', units='gC/m^2/s', &
         avgflag='A', long_name='C flux not allocated due to downregulation', &
         ptr_pft=clm3%g%l%c%p%pepv%excess_cflux, default='inactive')

    call hist_addfld1d (fname='DOWNREG', units='proportion', &
         avgflag='A', long_name='fractional reduction in GPP due to N limitation', &
         ptr_pft=clm3%g%l%c%p%pepv%downreg, default='inactive')

    call hist_addfld1d (fname='PREV_LEAFC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='previous timestep leaf C litterfall flux', &
         ptr_pft=clm3%g%l%c%p%pepv%prev_leafc_to_litter, default='inactive')

    call hist_addfld1d (fname='PREV_FROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='previous timestep froot C litterfall flux', &
         ptr_pft=clm3%g%l%c%p%pepv%prev_frootc_to_litter, default='inactive')

    call hist_addfld1d (fname='ANNSUM_NPP', units='gC/m^2/yr', &
         avgflag='A', long_name='annual sum of NPP', &
         ptr_pft=clm3%g%l%c%p%pepv%annsum_npp, default='inactive')

    if ( use_c13 ) then
       call hist_addfld1d (fname='RC13_CANAIR', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for canopy air', &
            ptr_pft=clm3%g%l%c%p%pepv%rc13_canair)
       
       call hist_addfld1d (fname='RC13_PSNSUN', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for sunlit photosynthesis', &
            ptr_pft=clm3%g%l%c%p%pepv%rc13_psnsun)
       
       call hist_addfld1d (fname='RC13_PSNSHA', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for shaded photosynthesis', &
            ptr_pft=clm3%g%l%c%p%pepv%rc13_psnsha)
    endif
    
    !-------------------------------
    ! PFT physical state variables not already defined by default
    !-------------------------------
    
    call hist_addfld1d (fname='EMV', units='proportion', &
         avgflag='A', long_name='vegetation emissivity', &
         ptr_pft=clm3%g%l%c%p%pps%emv, default='inactive')

    call hist_addfld1d (fname='Z0MV', units='m', &
         avgflag='A', long_name='roughness length over vegetation, momentum', &
         ptr_pft=clm3%g%l%c%p%pps%z0mv, default='inactive')

    call hist_addfld1d (fname='Z0HV', units='m', &
         avgflag='A', long_name='roughness length over vegetation, sensible heat', &
         ptr_pft=clm3%g%l%c%p%pps%z0hv, default='inactive')

    call hist_addfld1d (fname='Z0QV', units='m', &
         avgflag='A', long_name='roughness length over vegetation, latent heat', &
         ptr_pft=clm3%g%l%c%p%pps%z0qv, default='inactive')

    call hist_addfld1d (fname='DEWMX', units='mm', &
         avgflag='A', long_name='Maximum allowed dew', &
         ptr_pft=clm3%g%l%c%p%pps%dewmx, default='inactive')

    call hist_addfld1d (fname='FSUN', units='proportion', &
         avgflag='A', long_name='sunlit fraction of canopy', &
         ptr_pft=clm3%g%l%c%p%pps%fsun, default='inactive')

    if ( use_c13 ) then
       call hist_addfld1d (fname='ALPHAPSNSUN', units='proportion', &
            avgflag='A', long_name='sunlit c13 fractionation', &
            ptr_pft=clm3%g%l%c%p%pps%alphapsnsun, default='inactive')
       
       call hist_addfld1d (fname='ALPHAPSNSHA', units='proportion', &
            avgflag='A', long_name='shaded c13 fractionation', &
            ptr_pft=clm3%g%l%c%p%pps%alphapsnsha, default='inactive')
    endif
    
    call hist_addfld1d (fname='FWET', units='proportion', &
         avgflag='A', long_name='fraction of canopy that is wet', &
         ptr_pft=clm3%g%l%c%p%pps%fwet, default='inactive')
    
    call hist_addfld1d (fname='FDRY', units='proportion', &
         avgflag='A', long_name='fraction of foliage that is green and dry', &
         ptr_pft=clm3%g%l%c%p%pps%fdry, default='inactive')

    call hist_addfld1d (fname='DT_VEG', units='K', &
         avgflag='A', long_name='change in t_veg, last iteration', &
         ptr_pft=clm3%g%l%c%p%pps%dt_veg, default='inactive')

    call hist_addfld1d (fname='HTOP', units='m', &
         avgflag='A', long_name='canopy top', &
         ptr_pft=clm3%g%l%c%p%pps%htop)

    call hist_addfld1d (fname='HBOT', units='m', &
         avgflag='A', long_name='canopy bottom', &
         ptr_pft=clm3%g%l%c%p%pps%hbot, default='inactive')

    call hist_addfld1d (fname='Z0M', units='m', &
         avgflag='A', long_name='momentum roughness length', &
         ptr_pft=clm3%g%l%c%p%pps%z0m, default='inactive')

    call hist_addfld1d (fname='DISPLA', units='m', &
         avgflag='A', long_name='displacement height', &
         ptr_pft=clm3%g%l%c%p%pps%displa, default='inactive')

    call hist_addfld1d (fname='U10_DUST', units='m/s', &
         avgflag='A', long_name='10-m wind for dust model', &
         ptr_pft=clm3%g%l%c%p%pps%u10, default='inactive')

    call hist_addfld1d (fname='RAM1', units='s/m', &
         avgflag='A', long_name='aerodynamical resistance ', &
         ptr_pft=clm3%g%l%c%p%pps%ram1, default='inactive')

    call hist_addfld1d (fname='FV', units='m/s', &
         avgflag='A', long_name='friction velocity for dust model', &
         ptr_pft=clm3%g%l%c%p%pps%fv, default='inactive')

    call hist_addfld2d (fname='ROOTFR', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='fraction of roots in each soil layer', &
         ptr_pft=clm3%g%l%c%p%pps%rootfr, default='inactive')
                                                                       
    call hist_addfld2d (fname='ROOTR', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='effective fraction of roots in each soil layer', &
         ptr_pft=clm3%g%l%c%p%pps%rootr, default='inactive')
                                                                       
    call hist_addfld2d (fname='RRESIS', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='root resistance in each soil layer', &
         ptr_pft=clm3%g%l%c%p%pps%rresis, default='inactive')
                                                                       
    call hist_addfld2d (fname='ALBD', units='proportion', type2d='numrad', &
         avgflag='A', long_name='surface albedo (direct)', &
         ptr_pft=clm3%g%l%c%p%pps%albd, default='inactive', c2l_scale_type='urbanf')
                                                                        
    call hist_addfld2d (fname='ALBI', units='proportion', type2d='numrad', &
         avgflag='A', long_name='surface albedo (indirect)', &
         ptr_pft=clm3%g%l%c%p%pps%albi, default='inactive', c2l_scale_type='urbanf')
                                                                       
    if ( crop_prog )then

        call hist_addfld1d (fname='GDD0', units='ddays', &
             avgflag='A', long_name='Growing degree days base  0C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd0, default='inactive')

        call hist_addfld1d (fname='GDD8', units='ddays', &
             avgflag='A', long_name='Growing degree days base  8C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd8, default='inactive')

        call hist_addfld1d (fname='GDD10', units='ddays', &
             avgflag='A', long_name='Growing degree days base 10C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd10, default='inactive')

        call hist_addfld1d (fname='GDD020', units='ddays', &
             avgflag='A', long_name='Twenty year average of growing degree days base  0C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd020, default='inactive')

        call hist_addfld1d (fname='GDD820', units='ddays', &
             avgflag='A', long_name='Twenty year average of growing degree days base  8C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd820, default='inactive')

        call hist_addfld1d (fname='GDD1020', units='ddays', &
             avgflag='A', long_name='Twenty year average of growing degree days base 10C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd1020, default='inactive')

        call hist_addfld1d (fname='GDDPLANT', units='ddays', &
             avgflag='A', long_name='Accumulated growing degree days past planting date for crop', &
             ptr_pft=clm3%g%l%c%p%pps%gddplant, default='inactive')

        call hist_addfld1d (fname='GDDHARV', units='ddays', &
             avgflag='A', long_name='Growing degree days (gdd) needed to harvest', &
             ptr_pft=clm3%g%l%c%p%pps%gddmaturity, default='inactive')

        call hist_addfld1d (fname='GDDTSOI', units='ddays', &
             avgflag='A', long_name='Growing degree-days from planting (top two soil layers)', &
             ptr_pft=clm3%g%l%c%p%pps%gddtsoi, default='inactive')

    end if

    call hist_addfld_decomp (fname='CROOT_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from coarse roots', &
         ptr_pft=clm3%g%l%c%p%pps%croot_prof, default='inactive')

    call hist_addfld_decomp (fname='FROOT_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from fine roots', &
         ptr_pft=clm3%g%l%c%p%pps%froot_prof, default='inactive')

    call hist_addfld_decomp (fname='LEAF_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from leaves', &
         ptr_pft=clm3%g%l%c%p%pps%leaf_prof, default='inactive')

    call hist_addfld_decomp (fname='STEM_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from stems', &
         ptr_pft=clm3%g%l%c%p%pps%stem_prof, default='inactive')

    !-------------------------------
    ! Column physical state variables not already defined by default
    !-------------------------------

    call hist_addfld1d (fname='EMG', units='proportion', &
         avgflag='A', long_name='ground emissivity', &
         ptr_col=clm3%g%l%c%cps%emg, default='inactive')

    call hist_addfld1d (fname='BETA', units='none', &
         avgflag='A', long_name='coefficient of convective velocity', &
         ptr_col=clm3%g%l%c%cps%beta, default='inactive')

    call hist_addfld1d (fname='ZII', units='m', &
         avgflag='A', long_name='convective boundary height', &
         ptr_col=clm3%g%l%c%cps%zii, default='inactive')

    call hist_addfld1d (fname='WF', units='proportion', &
         avgflag='A', long_name='soil water as frac. of whc for top 0.05 m', &
         ptr_col=clm3%g%l%c%cps%wf)
  
   ! add by F. Li and S. Levis  
     call hist_addfld1d (fname='LFC2', units='per timestep', &
         avgflag='A', long_name='conversion area fraction of BET and BDT that burned in this timestep', &
         ptr_col=clm3%g%l%c%cps%lfc2)

    if ( nlevdecomp_full .gt. 1 ) then
       call hist_addfld1d (fname='FPI', units='proportion', &
            avgflag='A', long_name='fraction of potential immobilization', &
            ptr_col=clm3%g%l%c%cps%fpi)
    endif

    call hist_addfld_decomp (fname='FPI'//trim(vr_suffix), units='proportion', type2d='levdcmp', & 
         avgflag='A', long_name='fraction of potential immobilization', &
         ptr_col=clm3%g%l%c%cps%fpi_vr)

    call hist_addfld1d (fname='FPG', units='proportion', &
         avgflag='A', long_name='fraction of potential gpp', &
         ptr_col=clm3%g%l%c%cps%fpg)

    call hist_addfld1d (fname='ANNSUM_COUNTER', units='s', &
         avgflag='A', long_name='seconds since last annual accumulator turnover', &
         ptr_col=clm3%g%l%c%cps%annsum_counter, default='inactive')

    call hist_addfld1d (fname='CANNSUM_NPP', units='gC/m^2/s', &
         avgflag='A', long_name='annual sum of column-level NPP', &
         ptr_col=clm3%g%l%c%cps%cannsum_npp, default='inactive')

    call hist_addfld1d (fname='CANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average of 2m air temperature', &
         ptr_col=clm3%g%l%c%cps%cannavg_t2m, default='inactive')

    call hist_addfld2d (fname='FRAC_ICEOLD', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='fraction of ice relative to the tot water', &
         ptr_col=clm3%g%l%c%cps%frac_iceold, default='inactive')

    call hist_addfld2d (fname='EFF_POROSITY', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='effective porosity = porosity - vol_ice', &
         ptr_col=clm3%g%l%c%cps%eff_porosity, default='inactive')

    call hist_addfld2d (fname='ROOTR_COLUMN', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='effective fraction of roots in each soil layer', &
         ptr_col=clm3%g%l%c%cps%rootr_column, default='inactive')

    call hist_addfld2d (fname='ALBGRD', units='proportion', type2d='numrad', &
         avgflag='A', long_name='ground albedo (direct)', &
         ptr_col=clm3%g%l%c%cps%albgrd, default='inactive')

    call hist_addfld2d (fname='ALBGRI', units='proportion', type2d='numrad', &
         avgflag='A', long_name='ground albedo (indirect)', &
         ptr_col=clm3%g%l%c%cps%albgri, default='inactive')

   call hist_addfld1d (fname='NFIRE',  units='counts/km2/timestep', &
         avgflag='A', long_name='timestep fire counts valid only in Reg.C', &
         ptr_col=clm3%g%l%c%cps%nfire)
    
   call hist_addfld1d (fname='FAREA_BURNED',  units='proportion', &
         avgflag='A', long_name='timestep fractional area burned', &
         ptr_col=clm3%g%l%c%cps%farea_burned)

   call hist_addfld1d (fname='BAF_CROP',  units='proportion', &
         avgflag='A', long_name='timestep fractional area burned for crop', &
         ptr_col=clm3%g%l%c%cps%baf_crop)
         
    call hist_addfld1d (fname='BAF_PEATF',  units='proportion', &
         avgflag='A', long_name='timestep fractional area burned in peatland', &
         ptr_col=clm3%g%l%c%cps%baf_peatf)

    !-------------------------------
    ! Energy flux variables not already defined by default - native PFT 
    !-------------------------------

    call hist_addfld1d (fname='DLRAD', units='W/m^2', &
         avgflag='A', long_name='downward longwave radiation below the canopy', &
         ptr_pft=clm3%g%l%c%p%pef%dlrad, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='ULRAD', units='W/m^2', &
         avgflag='A', long_name='upward longwave radiation above the canopy', &
         ptr_pft=clm3%g%l%c%p%pef%ulrad, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='EFLX_SOIL_GRND', units='W/m^2', &
         avgflag='A', long_name='soil heat flux [+ into soil]', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='CGRND', units='W/m^2/K', &
         avgflag='A', long_name='deriv. of soil energy flux wrt to soil temp', &
         ptr_pft=clm3%g%l%c%p%pef%cgrnd, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='CGRNDL', units='W/m^2/K', &
         avgflag='A', long_name='deriv. of soil latent heat flux wrt soil temp', &
         ptr_pft=clm3%g%l%c%p%pef%cgrndl, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='CGRNDS', units='W/m^2/K', &
         avgflag='A', long_name='deriv. of soil sensible heat flux wrt soil temp', &
         ptr_pft=clm3%g%l%c%p%pef%cgrnds, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='EFLX_GNET', units='W/m^2', &
         avgflag='A', long_name='net heat flux into ground', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_gnet, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='DGNETDT', units='W/m^2/K', &
         avgflag='A', long_name='derivative of net ground heat flux wrt soil temp', &
         ptr_pft=clm3%g%l%c%p%pef%dgnetdT, default='inactive', c2l_scale_type='urbanf')

    !-------------------------------
    ! Water flux variables not already defined by default  - native PFT
    !-------------------------------

    call hist_addfld1d (fname='QFLX_RAIN_GRND', units='mm H2O/s', &
         avgflag='A', long_name='rain on ground after interception', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_rain_grnd, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_SNOW_GRND', units='mm H2O/s', &
         avgflag='A', long_name='snow on ground after interception', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_snow_grnd, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_EVAP_GRND', units='mm H2O/s', &
         avgflag='A', long_name='ground surface evaporation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_grnd, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_EVAP_VEG', units='mm H2O/s', &
         avgflag='A', long_name='vegetation evaporation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_veg, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_EVAP_TOT', units='mm H2O/s', &
         avgflag='A', long_name='qflx_evap_soi + qflx_evap_can + qflx_tran_veg', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_tot, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_DEW_GRND', units='mm H2O/s', &
         avgflag='A', long_name='ground surface dew formation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_dew_grnd, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_SUB_SNOW', units='mm H2O/s', &
         avgflag='A', long_name='sublimation rate from snow pack', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_sub_snow, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_DEW_SNOW', units='mm H2O/s', &
         avgflag='A', long_name='surface dew added to snow pacK', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_dew_snow, default='inactive', c2l_scale_type='urbanf')

#endif

    call hist_addfld1d (fname='SNORDSL', units='m^-6', &
         avgflag='A', long_name='top snow layer effective grain radius', &
         ptr_col=clm3%g%l%c%cps%snw_rds_top, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOTTOPL', units='K/m', &
         avgflag='A', long_name='snow temperature (top layer)', &
         ptr_col=clm3%g%l%c%cps%snot_top, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOdTdzL', units='K/m', &
         avgflag='A', long_name='top snow layer temperature gradient (land)', &
         ptr_col=clm3%g%l%c%cps%dTdz_top, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOLIQFL', units='fraction', &
         avgflag='A', long_name='top snow layer liquid water fraction (land)', &
         ptr_col=clm3%g%l%c%cps%sno_liq_top, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation from snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_sno_vd, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation from snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_sno_nd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation from snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_sno_vi, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation from snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_sno_ni, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation on snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_sno_vd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation on snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_sno_nd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation on snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_sno_vi, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation on snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_sno_ni, &
         default='inactive')

! New lake code
    call hist_addfld1d (fname='H2OSNO_TOP', units='kg/m2', &
         avgflag='A', long_name='mass of snow in top snow layer', &
         ptr_col=clm3%g%l%c%cps%h2osno_top, set_urb=spval)

    call hist_addfld1d (fname='SNOBCMCL', units='kg/m2', &
         avgflag='A', long_name='mass of BC in snow column', &
         ptr_col=clm3%g%l%c%cps%mss_bc_col, set_urb=spval)

    call hist_addfld1d (fname='SNOBCMSL', units='kg/m2', &
         avgflag='A', long_name='mass of BC in top snow layer', &
         ptr_col=clm3%g%l%c%cps%mss_bc_top, set_urb=spval)

    call hist_addfld1d (fname='SNOOCMCL', units='kg/m2', &
         avgflag='A', long_name='mass of OC in snow column', &
         ptr_col=clm3%g%l%c%cps%mss_oc_col, set_urb=spval)
   
    call hist_addfld1d (fname='SNOOCMSL', units='kg/m2', &
         avgflag='A', long_name='mass of OC in top snow layer', &
         ptr_col=clm3%g%l%c%cps%mss_oc_top, set_urb=spval)

    call hist_addfld1d (fname='SNODSTMCL', units='kg/m2', &
         avgflag='A', long_name='mass of dust in snow column', &
         ptr_col=clm3%g%l%c%cps%mss_dst_col, set_urb=spval)
    
    call hist_addfld1d (fname='SNODSTMSL', units='kg/m2', &
         avgflag='A', long_name='mass of dust in top snow layer', &
         ptr_col=clm3%g%l%c%cps%mss_dst_top, set_urb=spval)

    call hist_addfld1d (fname='DSTDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total dust deposition (dry+wet) from atmosphere', &
         ptr_col=clm3%g%l%c%cwf%flx_dst_dep, set_urb=spval)

    call hist_addfld1d (fname='BCDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total BC deposition (dry+wet) from atmosphere', &
         ptr_col=clm3%g%l%c%cwf%flx_bc_dep, set_urb=spval)
   
    call hist_addfld1d (fname='OCDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total OC deposition (dry+wet) from atmosphere', &
         ptr_col=clm3%g%l%c%cwf%flx_oc_dep, set_urb=spval)

#if (defined SNICAR_FRC)
    call hist_addfld1d (fname='SNOAERFRCL', units='W/m^2', &
         avgflag='A', long_name='surface forcing of all aerosols in snow (land) ', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_aer, set_urb=spval)
   
    call hist_addfld1d (fname='SNOAERFRC2L', units='W/m^2', &
         avgflag='A', long_name='surface forcing of all aerosols in snow, averaged only when snow is present (land)', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_aer_sno, set_urb=spval)

    call hist_addfld1d (fname='SNOBCFRCL', units='W/m^2', &
         avgflag='A', long_name='surface forcing of BC in snow (land) ', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_bc, set_urb=spval)
   
    call hist_addfld1d (fname='SNOBCFRC2L', units='W/m^2', &
         avgflag='A', long_name='surface forcing of BC in snow, averaged only when snow is present (land)', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_bc_sno, set_urb=spval)

    call hist_addfld1d (fname='SNOOCFRCL', units='W/m^2', &
         avgflag='A', long_name='surface forcing of OC in snow (land) ', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_oc, set_urb=spval)
   
    call hist_addfld1d (fname='SNOOCFRC2L', units='W/m^2', &
         avgflag='A', long_name='surface forcing of OC in snow, averaged only when snow is present (land)', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_oc_sno, set_urb=spval)

    call hist_addfld1d (fname='SNODSTFRCL', units='W/m^2', &
         avgflag='A', long_name='surface forcing of dust in snow (land) ', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_dst, set_urb=spval)
    
    call hist_addfld1d (fname='SNODSTFRC2L', units='W/m^2', &
         avgflag='A', long_name='surface forcing of dust in snow, averaged only when snow is present (land)', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_dst_sno, set_urb=spval)

#endif

    !-------------------------------
    ! Forcings sent to GLC
    !-------------------------------

    if (maxpatch_glcmec > 0) then

       call hist_addfld2d (fname='QICE_FORC', units='mm/s', type2d='glc_nec', &
            avgflag='A', long_name='qice forcing sent to GLC', &
            ptr_lnd=clm_s2x%qice, default='inactive')

       call hist_addfld2d (fname='TSRF_FORC', units='K', type2d='glc_nec', &
            avgflag='A', long_name='surface temperature sent to GLC', &
            ptr_lnd=clm_s2x%tsrf, default='inactive')

       call hist_addfld2d (fname='TOPO_FORC', units='m', type2d='glc_nec', &
            avgflag='A', long_name='topographic height sent to GLC', &
            ptr_lnd=clm_s2x%topo, default='inactive')

    end if

    ! Print masterlist of history fields

    call hist_printflds()

  end subroutine hist_initFlds


  subroutine hist_addfld_decomp (fname, type2d, units, avgflag, long_name, ptr_col, ptr_pft, default)

! !USES:
    use clm_varpar, only : nlevdecomp_full

    use histFileMod, only : hist_addfld1d, hist_addfld2d
    use abortutils, only: endrun
    use clm_varctl, only : iulog
    use surfrdMod , only : crop_prog
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname                    ! field name
    character(len=*), intent(in) :: type2d                   ! 2d output type
    character(len=*), intent(in) :: units                    ! units of field
    character(len=1), intent(in) :: avgflag                  ! time averaging flag
    character(len=*), intent(in) :: long_name                ! long name of field
    real(r8)        , optional, pointer    :: ptr_col(:,:)   ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_pft(:,:)   ! pointer to pft array
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape

    real(r8)        , pointer  :: ptr_1d(:)

    if (present(ptr_col)) then

       ! column-level data
       if (present(default)) then
          if ( nlevdecomp_full .gt. 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_col, default=default)
          else
             ptr_1d => ptr_col(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_1d, default=default)
          endif
       else
          if ( nlevdecomp_full .gt. 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_col)
          else
             ptr_1d => ptr_col(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_1d)
          endif
       endif
       
    else if (present(ptr_pft)) then
       
       ! pft-level data
       if (present(default)) then
          if ( nlevdecomp_full .gt. 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_pft=ptr_pft, default=default)
          else
             ptr_1d => ptr_pft(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_pft=ptr_1d, default=default)
          endif
       else
          if ( nlevdecomp_full .gt. 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_pft=ptr_pft)
          else
             ptr_1d => ptr_pft(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_pft=ptr_1d)
          endif
       endif
       
    else
       write(iulog, *) ' error: hist_addfld_decomp needs either pft or column level pointer'
       write(iulog, *) fname
       call endrun()
    endif
  end subroutine hist_addfld_decomp
  

end module histFldsMod
