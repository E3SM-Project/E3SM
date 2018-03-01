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
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public hist_initFlds ! Build master field list of all possible history
                       ! file fields
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 03/2003
! heald (11/28/06)
!
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
    use clm_glclnd , only : clm_s2x
    use clm_varctl , only : create_glacier_mec_landunit, use_c13, use_cndv, use_cn, use_crop, &
                            use_snicar_frc
    use histFileMod, only : hist_add_subscript, hist_addfld1d, hist_addfld2d, &
                            hist_printflds
    use surfrdMod  , only : crop_prog
    use shr_megan_mod  , only : shr_megan_linkedlist, shr_megan_megcomp_t, shr_megan_megcomps_n

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

    ! Determine what subscripts to add
    ! (uncomment the following call and modify it appropriately)

    ! call hist_add_subscript(subname='subscript_name', subdim=subscript_dim)

    ! NOTE: make a field not appear on the primary history tape by default -
    ! add the keyword to default='inactive' to the call to addfld_1d or addfld_2d

    ! Snow properties
    ! These will be vertically averaged over the snow profile

    call hist_addfld1d (fname='SNOWDP',  units='m',  &
         avgflag='A', long_name='snow height', &
         ptr_col=cps%snowdp, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSNO',  units='unitless',  &
         avgflag='A', long_name='fraction of ground covered by snow', &
         ptr_col=cps%frac_sno, c2l_scale_type='urbanf')

    ! Temperatures

    call hist_addfld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_pft=pes%t_ref2m)

    call hist_addfld1d (fname='TSA_U', units='K',  &
         avgflag='A', long_name='Urban 2m air temperature', &
         ptr_pft=pes%t_ref2m_u, set_nourb=spval)

    call hist_addfld1d (fname='TSA_R', units='K',  &
         avgflag='A', long_name='Rural 2m air temperature', &
         ptr_pft=pes%t_ref2m_r, set_spec=spval)

    call hist_addfld1d(fname='TBUILD', units='K',  &
         avgflag='A', long_name='internal urban building temperature', &
         ptr_lunit=lps%t_building, set_nourb=spval, l2g_scale_type='unity')

    call hist_addfld1d (fname='TREFMNAV', units='K',  &
         avgflag='A', long_name='daily minimum of average 2-m temperature', &
         ptr_pft=pes%t_ref2m_min)

    call hist_addfld1d (fname='TREFMXAV', units='K',  &
         avgflag='A', long_name='daily maximum of average 2-m temperature', &
         ptr_pft=pes%t_ref2m_max)

    call hist_addfld1d (fname='TREFMNAV_U', units='K',  &
         avgflag='A', long_name='Urban daily minimum of average 2-m temperature', &
         ptr_pft=pes%t_ref2m_min_u, set_nourb=spval)

    call hist_addfld1d (fname='TREFMXAV_U', units='K',  &
         avgflag='A', long_name='Urban daily maximum of average 2-m temperature', &
         ptr_pft=pes%t_ref2m_max_u, set_nourb=spval)

    call hist_addfld1d (fname='TREFMNAV_R', units='K',  &
         avgflag='A', long_name='Rural daily minimum of average 2-m temperature', &
         ptr_pft=pes%t_ref2m_min_r, set_spec=spval)

    call hist_addfld1d (fname='TREFMXAV_R', units='K',  &
         avgflag='A', long_name='Rural daily maximum of average 2-m temperature', &
         ptr_pft=pes%t_ref2m_max_r, set_spec=spval)

    call hist_addfld1d (fname='TV', units='K',  &
         avgflag='A', long_name='vegetation temperature', &
         ptr_pft=pes%t_veg)

    call hist_addfld1d (fname='TV24', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 24hrs)', &
         ptr_pft=pvs%t_veg24, default='inactive')

    call hist_addfld1d (fname='TV240', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 240hrs)', &
         ptr_pft=pvs%t_veg240, default='inactive')

    call hist_addfld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=ces%t_grnd, c2l_scale_type='urbans')

    call hist_addfld1d (fname='TG_U', units='K',  &
         avgflag='A', long_name='Urban ground temperature', &
         ptr_col=ces%t_grnd_u, set_nourb=spval, c2l_scale_type='urbans')

    call hist_addfld1d (fname='TG_R', units='K',  &
         avgflag='A', long_name='Rural ground temperature', &
         ptr_col=ces%t_grnd_r, set_spec=spval)

    call hist_addfld1d (fname='HCSOI',  units='MJ/m2',  &
         avgflag='A', long_name='soil heat content', &
         ptr_col=ces%hc_soi, set_lake=spval, set_urb=spval, l2g_scale_type='veg')

    call hist_addfld1d (fname='HC',  units='MJ/m2',  &
         avgflag='A', long_name='heat content of soil/snow/lake', &
         ptr_col=ces%hc_soisno, set_urb=spval)

    call hist_addfld2d (fname='TSOI',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (vegetated landunits only)', &
         ptr_col=ces%t_soisno, l2g_scale_type='veg')

    call hist_addfld2d (fname='TSOI_ICE',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (ice landunits only)', &
         ptr_col=ces%t_soisno, l2g_scale_type='ice')

    call hist_addfld1d (fname='TSOI_10CM',  units='K', &
         avgflag='A', long_name='soil temperature in top 10cm of soil', &
         ptr_col=ces%t_soi_10cm, set_urb=spval)

    call hist_addfld2d (fname='TLAKE',  units='K', type2d='levlak', &
         avgflag='A', long_name='lake temperature', &
         ptr_col=ces%t_lake)

    ! Specific humidity

    call hist_addfld1d (fname='Q2M', units='kg/kg',  &
         avgflag='A', long_name='2m specific humidity', &
         ptr_pft=pes%q_ref2m)

    ! Relative humidity

    call hist_addfld1d (fname='RH2M', units='%',  &
         avgflag='A', long_name='2m relative humidity', &
         ptr_pft=pes%rh_ref2m)

    call hist_addfld1d (fname='RH2M_U', units='%',  &
         avgflag='A', long_name='Urban 2m relative humidity', &
         ptr_pft=pes%rh_ref2m_u, set_nourb=spval)

    call hist_addfld1d (fname='RH2M_R', units='%',  &
         avgflag='A', long_name='Rural 2m specific humidity', &
         ptr_pft=pes%rh_ref2m_r, set_spec=spval)

    ! Wind

    call hist_addfld1d (fname='U10', units='m/s', &
         avgflag='A', long_name='10-m wind', &
         ptr_pft=pps%u10_clm)
    call hist_addfld1d (fname='VA', units='m/s', &
         avgflag='A', long_name='atmospheric wind speed plus convective velocity', &
         ptr_pft=pps%va, default='inactive')

    ! Surface radiation

    call hist_addfld1d (fname='SABV', units='W/m^2',  &
         avgflag='A', long_name='solar rad absorbed by veg', &
         ptr_pft=pef%sabv, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='SABG', units='W/m^2',  &
         avgflag='A', long_name='solar rad absorbed by ground', &
         ptr_pft=pef%sabg, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation', &
         ptr_pft=pef%fsds_vis_d)

    call hist_addfld1d (fname='FSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation', &
         ptr_pft=pef%fsds_nir_d)

    call hist_addfld1d (fname='FSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation', &
         ptr_pft=pef%fsds_vis_i)

    call hist_addfld1d (fname='FSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation', &
         ptr_pft=pef%fsds_nir_i)

    call hist_addfld1d (fname='FSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation', &
         ptr_pft=pef%fsr_vis_d, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation', &
         ptr_pft=pef%fsr_nir_d, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation', &
         ptr_pft=pef%fsr_vis_i, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation', &
         ptr_pft=pef%fsr_nir_i, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSDSVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation at local noon', &
         ptr_pft=pef%fsds_vis_d_ln)

    call hist_addfld1d (fname='FSDSNDLN', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation at local noon', &
         ptr_pft=pef%fsds_nir_d_ln)

    call hist_addfld1d (fname='FSRVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation at local noon', &
         ptr_pft=pef%fsr_vis_d_ln, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRNDLN', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation at local noon', &
         ptr_pft=pef%fsr_nir_d_ln, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSA', units='W/m^2',  &
         avgflag='A', long_name='absorbed solar radiation', &
         ptr_pft=pef%fsa, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban absorbed solar radiation', &
         ptr_pft=pef%fsa_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural absorbed solar radiation', &
         ptr_pft=pef%fsa_r, set_spec=spval)

    call hist_addfld1d (fname='FSR', units='W/m^2',  &
         avgflag='A', long_name='reflected solar radiation', &
         ptr_pft=pef%fsr, c2l_scale_type='urbanf')

    ! Rename of FSR for Urban intercomparision project
    call hist_addfld1d (fname='SWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling shortwave radiation', &
         ptr_pft=pef%fsr, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FIRA', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_pft=pef%eflx_lwrad_net, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FIRA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban net infrared (longwave) radiation', &
         ptr_pft=pef%eflx_lwrad_net_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FIRA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural net infrared (longwave) radiation', &
         ptr_pft=pef%eflx_lwrad_net_r, set_spec=spval)

    call hist_addfld1d (fname='FIRE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_pft=pef%eflx_lwrad_out, c2l_scale_type='urbanf')

    ! Rename of FIRE for Urban intercomparision project
    call hist_addfld1d (fname='LWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling longwave radiation', &
         ptr_pft=pef%eflx_lwrad_out, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='BUILDHEAT', units='W/m^2',  &
         avgflag='A', long_name='heat flux from urban building interior to walls and roof', &
         ptr_col=cef%eflx_building_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
         avgflag='A', long_name='urban air conditioning flux', &
         ptr_col=cef%eflx_urban_ac, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
         avgflag='A', long_name='urban heating flux', &
         ptr_col=cef%eflx_urban_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='TRAFFICFLUX', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from urban traffic', &
         ptr_pft=pef%eflx_traffic_pft, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    call hist_addfld1d (fname='WASTEHEAT', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from heating/cooling sources of urban waste heat', &
         ptr_pft=pef%eflx_wasteheat_pft, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='HEAT_FROM_AC', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux put into canyon due to heat removed from air conditioning', &
         ptr_pft=pef%eflx_heat_from_ac_pft, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='Qanth', units='W/m^2',  &
         avgflag='A', long_name='anthropogenic heat flux', &
         ptr_pft=pef%eflx_anthro, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    call hist_addfld1d (fname='Rnet', units='W/m^2',  &
         avgflag='A', long_name='net radiation', &
         ptr_pft=pef%netrad, c2l_scale_type='urbanf', &
         default='inactive')

    ! Solar zenith angle and solar declination angle

    call hist_addfld1d (fname='COSZEN', units='none', &
         avgflag='A', long_name='cosine of solar zenith angle', &
         ptr_col=cps%coszen, default='inactive')

    call hist_addfld1d (fname='DECL', units='radians', &
         avgflag='A', long_name='solar declination angle', &
         ptr_col=cps%decl, default='inactive')

    ! Surface energy fluxes

    call hist_addfld1d (fname='FCTR', units='W/m^2',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_pft=pef%eflx_lh_vegt, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FCEV', units='W/m^2',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_pft=pef%eflx_lh_vege, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGEV', units='W/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_pft=pef%eflx_lh_grnd, c2l_scale_type='urbanf') 

    call hist_addfld1d (fname='FSH_NODYNLNDUSE', units='W/m^2',  &
         avgflag='A', long_name='sensible heat not including correction for land use change', &
         ptr_pft=pef%eflx_sh_tot, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSH', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_lnd=gef%eflx_sh_totg)

    call hist_addfld1d (fname='FSH_U', units='W/m^2',  &
         avgflag='A', long_name='Urban sensible heat', &
         ptr_pft=pef%eflx_sh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSH_R', units='W/m^2',  &
         avgflag='A', long_name='Rural sensible heat', &
         ptr_pft=pef%eflx_sh_tot_r, set_spec=spval)

    call hist_addfld1d (fname='GC_HEAT1',  units='J/m^2',  &
         avgflag='A', long_name='initial gridcell total heat content', &
         ptr_lnd=ges%gc_heat1)

    call hist_addfld1d (fname='GC_HEAT2',  units='J/m^2',  &
         avgflag='A', long_name='post land cover change total heat content', &
         ptr_lnd=ges%gc_heat2, default='inactive')

    call hist_addfld1d (fname='EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_lnd=gef%eflx_dynbal)

    call hist_addfld1d (fname='Qh', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_pft=pef%eflx_sh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='Qle', units='W/m^2',  &
         avgflag='A', long_name='total evaporation', &
         ptr_pft=pef%eflx_lh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='EFLX_LH_TOT_U', units='W/m^2',  &
         avgflag='A', long_name='Urban total evaporation', &
         ptr_pft=pef%eflx_lh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='EFLX_LH_TOT_R', units='W/m^2',  &
         avgflag='A', long_name='Rural total evaporation', &
         ptr_pft=pef%eflx_lh_tot_r, set_spec=spval)

    call hist_addfld1d (fname='Qstor', units='W/m^2',  &
         avgflag='A', long_name='storage heat flux (includes snowmelt)', &
         ptr_pft=pef%eflx_soil_grnd, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='FSH_V', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from veg', &
         ptr_pft=pef%eflx_sh_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSH_G', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_pft=pef%eflx_sh_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGR', units='W/m^2',  &
         avgflag='A', long_name='heat flux into soil/snow including snow melt', &
         ptr_pft=pef%eflx_soil_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGR_U', units='W/m^2',  &
         avgflag='A', long_name='Urban heat flux into soil/snow including snow melt', &
         ptr_pft=pef%eflx_soil_grnd_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FGR_R', units='W/m^2',  &
         avgflag='A', long_name='Rural heat flux into soil/snow including snow melt', &
         ptr_pft=pef%eflx_soil_grnd_r, set_spec=spval)

    call hist_addfld1d (fname='FSM',  units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=cef%eflx_snomelt, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSM_U',  units='W/m^2',  &
         avgflag='A', long_name='Urban snow melt heat flux', &
         ptr_col=cef%eflx_snomelt_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSM_R',  units='W/m^2',  &
         avgflag='A', long_name='Rural snow melt heat flux', &
         ptr_col=cef%eflx_snomelt_r, set_spec=spval)

    call hist_addfld1d (fname='FGR12',  units='W/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=cef%eflx_fgr12, set_lake=spval)

    call hist_addfld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_pft=pmf%taux)

    ! Rename of TAUX for Urban intercomparision project (when U=V)
    call hist_addfld1d (fname='Qtau', units='kg/m/s^2',  &
         avgflag='A', long_name='momentum flux', &
         ptr_pft=pmf%taux, default='inactive')

    call hist_addfld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_pft=pmf%tauy)

    ! Vegetation phenology

    call hist_addfld1d (fname='ELAI', units='m^2/m^2', &
          avgflag='A', long_name='exposed one-sided leaf area index', &
         ptr_pft=pps%elai)

    call hist_addfld1d (fname='ESAI', units='m^2/m^2', &
          avgflag='A', long_name='exposed one-sided stem area index', &
         ptr_pft=pps%esai)

    call hist_addfld1d (fname='LAISUN', units='none', &
         avgflag='A', long_name='sunlit projected leaf area index', &
         ptr_pft=pps%laisun, set_urb=0._r8)

    call hist_addfld1d (fname='LAISHA', units='none', &
         avgflag='A', long_name='shaded projected leaf area index', &
         ptr_pft=pps%laisha, set_urb=0._r8)

    call hist_addfld1d (fname='TLAI', units='none', &
         avgflag='A', long_name='total projected leaf area index', &
         ptr_pft=pps%tlai)

    call hist_addfld1d (fname='TSAI', units='none', &
         avgflag='A', long_name='total projected stem area index', &
         ptr_pft=pps%tsai)

    call hist_addfld1d (fname='SLASUN', units='m^2/gC', &
         avgflag='A', long_name='specific leaf area for sunlit canopy, projected area basis', &
         ptr_pft=pps%slasun, set_urb=0._r8, default='inactive')

    call hist_addfld1d (fname='SLASHA', units='m^2/gC', &
         avgflag='A', long_name='specific leaf area for shaded canopy, projected area basis', &
         ptr_pft=pps%slasha, set_urb=0._r8, default='inactive')

    ! Canopy physiology

    call hist_addfld1d (fname='RSSUN', units='s/m',  &
         avgflag='M', long_name='sunlit leaf stomatal resistance', &
         ptr_pft=pps%rssun, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='RSSHA', units='s/m',  &
         avgflag='M', long_name='shaded leaf stomatal resistance', &
         ptr_pft=pps%rssha, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='BTRAN', units='unitless',  &
         avgflag='A', long_name='transpiration beta factor', &
         ptr_pft=pps%btran, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='FPSN', units='umol/m2s',  &
         avgflag='A', long_name='photosynthesis', &
         ptr_pft=pcf%fpsn, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='DSTFLXT', units='kg/m2/s',  &
         avgflag='A', long_name='total surface dust emission', &
         ptr_pft=pdf%flx_mss_vrt_dst_tot, set_lake=0._r8, set_urb=0._r8)
    call hist_addfld1d (fname='DPVLTRB1', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 1', &
         ptr_pft=pdf%vlc_trb_1, default='inactive')
    call hist_addfld1d (fname='DPVLTRB2', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 2', &
         ptr_pft=pdf%vlc_trb_2, default='inactive')
    call hist_addfld1d (fname='DPVLTRB3', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 3', &
         ptr_pft=pdf%vlc_trb_3, default='inactive')
    call hist_addfld1d (fname='DPVLTRB4', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 4', &
         ptr_pft=pdf%vlc_trb_4, default='inactive')

    ! for MEGAN emissions diagnositics
    if (shr_megan_megcomps_n>0) then
       
       ! loop over megan compounds
       meg_cmp => shr_megan_linkedlist
       do while(associated(meg_cmp))
          imeg = meg_cmp%index

          call hist_addfld1d ( fname='MEG_'//trim(meg_cmp%name), units='kg/m2/sec',  &
               avgflag='A', long_name='MEGAN flux', &
               ptr_pft=pvf%meg(imeg)%flux_out, set_lake=0._r8, set_urb=0._r8 )

          meg_cmp => meg_cmp%next_megcomp
       enddo
       
       call hist_addfld1d (fname='VOCFLXT', units='moles/m2/sec',  &
            avgflag='A', long_name='total VOC flux into atmosphere', &
            ptr_pft=pvf%vocflx_tot, set_lake=0._r8, set_urb=0._r8)

       call hist_addfld1d (fname='GAMMA', units='non',  &
            avgflag='A', long_name='total gamma for VOC calc', &
            ptr_pft=pvf%gamma_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAL', units='non',  &
            avgflag='A', long_name='gamma L for VOC calc', &
            ptr_pft=pvf%gammaL_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAT', units='non',  &
            avgflag='A', long_name='gamma T for VOC calc', &
            ptr_pft=pvf%gammaT_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAP', units='non',  &
            avgflag='A', long_name='gamma P for VOC calc', &
            ptr_pft=pvf%gammaP_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAA', units='non',  &
            avgflag='A', long_name='gamma A for VOC calc', &
            ptr_pft=pvf%gammaA_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAS', units='non',  &
            avgflag='A', long_name='gamma S for VOC calc', &
            ptr_pft=pvf%gammaS_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='GAMMAC', units='non',  &
            avgflag='A', long_name='gamma C for VOC calc', &
            ptr_pft=pvf%gammaC_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='EOPT', units='non',  &
            avgflag='A', long_name='Eopt coefficient for VOC calc', &
            ptr_pft=pvf%Eopt_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='TOPT', units='non',  &
            avgflag='A', long_name='topt coefficient for VOC calc', &
            ptr_pft=pvf%topt_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='ALPHA', units='non',  &
            avgflag='A', long_name='alpha coefficient for VOC calc', &
            ptr_pft=pvf%alpha_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='CP', units='non',  &
            avgflag='A', long_name='cp coefficient for VOC calc', &
            ptr_pft=pvf%cp_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR', &
            ptr_pft=pvf%paru_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR24_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR (24 hrs)', &
            ptr_pft=pvf%par24u_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR240_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR (240 hrs)', &
            ptr_pft=pvf%par240u_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR', &
            ptr_pft=pvf%para_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR24_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR (24 hrs)', &
            ptr_pft=pvf%par24a_out, set_lake=0._r8, default='inactive')

       call hist_addfld1d (fname='PAR240_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR (240 hrs)', &
            ptr_pft=pvf%par240a_out, set_lake=0._r8, default='inactive')

    endif

    call hist_addfld1d (fname='FSUN24', units='K',  &
         avgflag='A', long_name='fraction sunlit (last 24hrs)', &
         ptr_pft=pvs%fsun24, default='inactive')

    call hist_addfld1d (fname='FSUN240', units='K',  &
         avgflag='A', long_name='fraction sunlit (last 240hrs)', &
         ptr_pft=pvs%fsun240, default='inactive')

    call hist_addfld1d (fname='FSI24', units='K',  &
         avgflag='A', long_name='indirect radiation (last 24hrs)', &
         ptr_pft=pvs%fsi24, default='inactive')

    call hist_addfld1d (fname='FSI240', units='K',  &
         avgflag='A', long_name='indirect radiation (last 240hrs)', &
         ptr_pft=pvs%fsi240, default='inactive')

    call hist_addfld1d (fname='FSD24', units='K',  &
         avgflag='A', long_name='direct radiation (last 24hrs)', &
         ptr_pft=pvs%fsd24, default='inactive')

    call hist_addfld1d (fname='FSD240', units='K',  &
         avgflag='A', long_name='direct radiation (last 240hrs)', &
         ptr_pft=pvs%fsd240, default='inactive')

    ! Hydrology

    call hist_addfld1d (fname='SoilAlpha',  units='unitless',  &
         avgflag='A', long_name='factor limiting ground evap', &
         ptr_col=cws%soilalpha, set_urb=spval)

    call hist_addfld1d (fname='SoilAlpha_U',  units='unitless',  &
         avgflag='A', long_name='urban factor limiting ground evap', &
         ptr_col=cws%soilalpha_u, set_nourb=spval)

    call hist_addfld1d (fname='FCOV',  units='unitless',  &
         avgflag='A', long_name='fractional impermeable area', &
         ptr_col=cws%fcov, l2g_scale_type='veg')
    call hist_addfld1d (fname='FSAT',  units='unitless',  &
         avgflag='A', long_name='fractional area with water table at surface', &
         ptr_col=cws%fsat, l2g_scale_type='veg')
    call hist_addfld1d (fname='ZWT',  units='m',  &
         avgflag='A', long_name='water table depth (vegetated landunits only)', &
         ptr_col=cws%zwt, l2g_scale_type='veg')

    call hist_addfld1d (fname='WA',  units='mm',  &
         avgflag='A', long_name='water in the unconfined aquifer (vegetated landunits only)', &
         ptr_col=cws%wa, l2g_scale_type='veg')

    call hist_addfld1d (fname='WT',  units='mm',  &
         avgflag='A', long_name='total water storage (unsaturated soil water + groundwater, veg landunits)', &
         ptr_col=cws%wt, l2g_scale_type='veg')

    call hist_addfld1d (fname='QCHARGE',  units='mm/s',  &
         avgflag='A', long_name='aquifer recharge rate (vegetated landunits only)', &
         ptr_col=cws%qcharge, l2g_scale_type='veg')

    call hist_addfld2d (fname='SMP',  units='mm', type2d='levgrnd',  &
         avgflag='A', long_name='soil matric potential (vegetated landunits only)', &
         ptr_col=cws%smp_l, set_spec=spval, l2g_scale_type='veg', default='inactive')

    call hist_addfld2d (fname='HK',  units='mm/s', type2d='levgrnd',  &
         avgflag='A', long_name='hydraulic conductivity (vegetated landunits only)', &
         ptr_col=cws%hk_l, set_spec=spval, l2g_scale_type='veg', default='inactive')

    call hist_addfld1d (fname='H2OSNO',  units='mm',  &
         avgflag='A', long_name='snow depth (liquid water)', &
         ptr_col=cws%h2osno, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='ERRH2OSNO',  units='mm',  &
         avgflag='A', long_name='imbalance in snow depth (liquid water)', &
         ptr_col=cws%errh2osno, c2l_scale_type='urbanf')

    ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
    ! any given time step but only if there is at least one snow layer (for all landunits 
    ! except lakes).  h2osno also includes snow that is part of the soil column (an 
    ! initial snow layer is only created if h2osno > 10mm). Also note that monthly average
    ! files of snow_sources and snow sinks must be weighted by number of days in the month to 
    ! diagnose, for example, an annual value of the change in h2osno. 

    call hist_addfld1d (fname='SNOW_SOURCES',  units='mm/s',  &
         avgflag='A', long_name='snow sources (liquid water)', &
         ptr_col=cws%snow_sources, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='SNOW_SINKS',  units='mm/s',  &
         avgflag='A', long_name='snow sinks (liquid water)', &
         ptr_col=cws%snow_sinks, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='H2OCAN', units='mm',  &
         avgflag='A', long_name='intercepted water', &
         ptr_pft=pws%h2ocan, set_lake=0._r8)

    call hist_addfld2d (fname='H2OSOI',  units='mm3/mm3', type2d='levgrnd', &
         avgflag='A', long_name='volumetric soil water (vegetated landunits only)', &
         ptr_col=cws%h2osoi_vol, l2g_scale_type='veg')

    call hist_addfld2d (fname='SOILLIQ',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='soil liquid water (vegetated landunits only)', &
         ptr_col=cws%h2osoi_liq, l2g_scale_type='veg')

    call hist_addfld2d (fname='SOILICE',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='soil ice (vegetated landunits only)', &
         ptr_col=cws%h2osoi_ice, l2g_scale_type='veg')

    call hist_addfld1d (fname='SOILWATER_10CM',  units='kg/m2', &
         avgflag='A', long_name='soil liquid water + ice in top 10cm of soil (veg landunits only)', &
         ptr_col=cws%h2osoi_liqice_10cm, set_urb=spval, l2g_scale_type='veg')

    call hist_addfld1d (fname='SNOWLIQ',  units='kg/m2',  &
         avgflag='A', long_name='snow liquid water', &
         ptr_col=cws%snowliq)

    call hist_addfld1d (fname='SNOWICE',  units='kg/m2',  &
         avgflag='A', long_name='snow ice', &
         ptr_col=cws%snowice)

    call hist_addfld1d (fname='QTOPSOIL',  units='mm/s',  &
         avgflag='A', long_name='water input to surface', &
         ptr_col=cwf%qflx_top_soil, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='QINFL',  units='mm/s',  &
         avgflag='A', long_name='infiltration', &
         ptr_col=cwf%qflx_infl, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QOVER',  units='mm/s',  &
         avgflag='A', long_name='surface runoff', &
         ptr_col=cwf%qflx_surf, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRGWL',  units='mm/s',  &
         avgflag='A', long_name='surface runoff at glaciers (liquid only), wetlands, lakes', &
         ptr_col=cwf%qflx_qrgwl, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNWCPLIQ', units='mm H2O/s', &
         avgflag='A', long_name='excess rainfall due to snow capping', &
         ptr_pft=pwf%qflx_snwcp_liq, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='QSNWCPICE_NODYNLNDUSE', units='mm H2O/s', &
         avgflag='A', &
         long_name='excess snowfall due to snow capping not including correction for land use change', &
         ptr_pft=pwf%qflx_snwcp_ice, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNWCPICE',  units='mm/s',  &
         avgflag='A', long_name='excess snowfall due to snow capping', &
         ptr_lnd=gwf%qflx_snwcp_iceg)

    call hist_addfld1d (fname='QDRAI',  units='mm/s',  &
         avgflag='A', long_name='sub-surface drainage', &
         ptr_col=cwf%qflx_drain, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF_NODYNLNDUSE',  units='mm/s',  &
         avgflag='A', &
         long_name='total liquid runoff (does not include QSNWCPICE) not including correction for land use change', &
         ptr_col=cwf%qflx_runoff, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF',  units='mm/s',  &
         avgflag='A', long_name='total liquid runoff (does not include QSNWCPICE)', &
         ptr_lnd=gwf%qflx_runoffg)

    call hist_addfld1d (fname='GC_LIQ1',  units='mm',  &
         avgflag='A', long_name='initial gridcell total liq content', &
         ptr_lnd=gws%gc_liq1)

    call hist_addfld1d (fname='GC_LIQ2',  units='mm',  &  
         avgflag='A', long_name='post landuse change gridcell total liq content', &              
         ptr_lnd=gws%gc_liq2, default='inactive')     

    call hist_addfld1d (fname='QFLX_LIQ_DYNBAL',  units='mm/s',  &  
         avgflag='A', long_name='liq dynamic land cover change conversion runoff flux', &              
         ptr_lnd=gwf%qflx_liq_dynbal)     

    call hist_addfld1d (fname='GC_ICE1',  units='mm',  &  
         avgflag='A', long_name='initial gridcell total ice content', &              
         ptr_lnd=gws%gc_ice1)     

    call hist_addfld1d (fname='GC_ICE2',  units='mm',  &  
         avgflag='A', long_name='post land cover change total ice content', &              
         ptr_lnd=gws%gc_ice2, default='inactive')

    call hist_addfld1d (fname='QFLX_ICE_DYNBAL',  units='mm/s',  &
         avgflag='A', long_name='ice dynamic land cover change conversion runoff flux', &                                   
         ptr_lnd=gwf%qflx_ice_dynbal)

    call hist_addfld1d (fname='QRUNOFF_U', units='mm/s',  &
         avgflag='A', long_name='Urban total runoff', &
         ptr_col=cwf%qflx_runoff_u, set_nourb=spval, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF_R', units='mm/s',  &
         avgflag='A', long_name='Rural total runoff', &
         ptr_col=cwf%qflx_runoff_r, set_spec=spval)

    call hist_addfld1d (fname='QINTR', units='mm/s',  &
         avgflag='A', long_name='interception', &
         ptr_pft=pwf%qflx_prec_intr, set_lake=0._r8)

    call hist_addfld1d (fname='QDRIP', units='mm/s',  &
         avgflag='A', long_name='throughfall', &
         ptr_pft=pwf%qflx_prec_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNOMELT',  units='mm/s',  &
         avgflag='A', long_name='snow melt', &
         ptr_col=cwf%qflx_snomelt, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNOFRZ', units='kg/m2/s', &
         avgflag='A', long_name='column-integrated snow freezing rate', &
         ptr_col=cwf%qflx_snofrz_col, default='inactive', &
         set_lake=spval, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSOIL', units='mm/s',  &
         avgflag='A', long_name= &
         'Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)', &
         ptr_pft=pwf%qflx_evap_soi, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QVEGE', units='mm/s',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_pft=pwf%qflx_evap_can, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QVEGT', units='mm/s',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_pft=pwf%qflx_tran_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QIRRIG', units='mm/s', &
         avgflag='A', long_name='water added through irrigation', &
         ptr_col=cwf%qflx_irrig, set_lake=0._r8)

    if (create_glacier_mec_landunit) then

       call hist_addfld1d (fname='QICE',  units='mm/s',  &
            avgflag='A', long_name='ice growth/melt', &
            ptr_col=cwf%qflx_glcice, set_noglcmec=spval)

       call hist_addfld1d (fname='QICE_FRZ',  units='mm/s',  &
            avgflag='A', long_name='ice growth', &
            ptr_col=cwf%qflx_glcice_frz, set_noglcmec=spval)

       call hist_addfld1d (fname='QICE_MELT',  units='mm/s',  &
            avgflag='A', long_name='ice melt', &
            ptr_col=cwf%qflx_glcice_melt, set_noglcmec=spval)

       call hist_addfld1d (fname='gris_mask',  units='unitless',  &
            avgflag='A', long_name='Greenland mask', &
            ptr_gcell=grc%gris_mask)

       call hist_addfld1d (fname='gris_area',  units='km^2',  &
            avgflag='A', long_name='Greenland ice area', &
            ptr_gcell=grc%gris_area)

       call hist_addfld1d (fname='aais_mask',  units='unitless',  &
            avgflag='A', long_name='Antarctic mask', &
            ptr_gcell=grc%aais_mask)

       call hist_addfld1d (fname='aais_area',  units='km^2',  &
            avgflag='A', long_name='Antarctic ice area', &
            ptr_gcell=grc%aais_area)

   endif

    ! Water and energy balance checks

    call hist_addfld1d (fname='ERRSOI',  units='W/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=cebal%errsoi)

    call hist_addfld1d (fname='ERRSEB',  units='W/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_pft=pebal%errseb)

    call hist_addfld1d (fname='ERRSOL',  units='W/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_pft=pebal%errsol, set_urb=spval)

    call hist_addfld1d (fname='ERRH2O', units='mm',  &
         avgflag='A', long_name='total water conservation error', &
         ptr_col=cwbal%errh2o)

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

    if (use_cndv .or. use_crop) then
       call hist_addfld1d (fname='T10', units='K',  &
            avgflag='A', long_name='10-day running mean of 2-m temperature', &
            ptr_pft=pes%t10)
    end if

    if (use_cndv) then
       call hist_addfld1d (fname='TDA', units='K',  &
            avgflag='A', long_name='daily average 2-m temperature', &
            ptr_pft=pdgvs%t_mo)
       
       call hist_addfld1d (fname='AGDD', units='K',  &
            avgflag='A', long_name='growing degree-days base 5C', &
            ptr_pft=pdgvs%agdd)
    end if

    if (use_cn) then
       call hist_addfld2d (fname='SOILPSI', units='MPa', type2d='levgrnd', &
            avgflag='A', long_name='soil water potential in each soil layer', &
            ptr_col=cps%soilpsi)
    end if

    if (use_cn) then
       ! add history fields for all CN variables, always set as default='inactive'
       if ( crop_prog )then
          
          call hist_addfld1d (fname='A5TMIN', units='K',  &
               avgflag='A', long_name='5-day running mean of min 2-m temperature', &
               ptr_pft=pes%a5tmin, default='inactive')
          
          call hist_addfld1d (fname='A10TMIN', units='K',  &
               avgflag='A', long_name='10-day running mean of min 2-m temperature', &
               ptr_pft=pes%a10tmin, default='inactive')
          
       end if
    
       !-------------------------------
       ! C state variables - native to PFT 
       !-------------------------------
       ! add history fields for all CLAMP CN variables
       
       call hist_addfld1d (fname='WOODC', units='gC/m^2', &
            avgflag='A', long_name='wood C', &
            ptr_pft=pcs%woodc)
       
       call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
            avgflag='A', long_name='leaf C', &
            ptr_pft=pcs%leafc)
       
       call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='leaf C storage', &
            ptr_pft=pcs%leafc_storage, default='inactive')

       call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
            avgflag='A', long_name='leaf C transfer', &
            ptr_pft=pcs%leafc_xfer, default='inactive')

       call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
            avgflag='A', long_name='fine root C', &
            ptr_pft=pcs%frootc)

       call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='fine root C storage', &
            ptr_pft=pcs%frootc_storage, default='inactive')

       call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='fine root C transfer', &
            ptr_pft=pcs%frootc_xfer, default='inactive')

       call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
            avgflag='A', long_name='live stem C', &
            ptr_pft=pcs%livestemc)

       call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live stem C storage', &
            ptr_pft=pcs%livestemc_storage, default='inactive')

       call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live stem C transfer', &
            ptr_pft=pcs%livestemc_xfer, default='inactive')

       call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
            avgflag='A', long_name='dead stem C', &
            ptr_pft=pcs%deadstemc)

       call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead stem C storage', &
            ptr_pft=pcs%deadstemc_storage, default='inactive')

       call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead stem C transfer', &
            ptr_pft=pcs%deadstemc_xfer, default='inactive')

       call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C', &
            ptr_pft=pcs%livecrootc)

       call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C storage', &
            ptr_pft=pcs%livecrootc_storage, default='inactive')

       call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C transfer', &
            ptr_pft=pcs%livecrootc_xfer, default='inactive')

       call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C', &
            ptr_pft=pcs%deadcrootc)

       call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C storage', &
            ptr_pft=pcs%deadcrootc_storage,  default='inactive')

       call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C transfer', &
            ptr_pft=pcs%deadcrootc_xfer, default='inactive')

       call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='growth respiration storage', &
            ptr_pft=pcs%gresp_storage, default='inactive')

       call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
            avgflag='A', long_name='growth respiration transfer', &
            ptr_pft=pcs%gresp_xfer, default='inactive')

       call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_pft=pcs%cpool)

       call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_pft=pcs%xsmrpool)

       call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
            avgflag='A', long_name='pft-level sink for C truncation', &
            ptr_pft=pcs%pft_ctrunc)

       call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
            avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
            ptr_pft=pcs%dispvegc)

       call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
            avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
            ptr_pft=pcs%storvegc)

       call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
            avgflag='A', long_name='total vegetation carbon, excluding cpool', &
            ptr_pft=pcs%totvegc)

       call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
            avgflag='A', long_name='total pft-level carbon, including cpool', &
            ptr_pft=pcs%totpftc)

       if (use_c13) then
          !-------------------------------
          ! C13 state variables - native to PFT 
          !-------------------------------

          call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
               avgflag='A', long_name='C13 leaf C', &
               ptr_pft=pc13s%leafc)

          call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
               avgflag='A', long_name='C13 leaf C storage', &
               ptr_pft=pc13s%leafc_storage, default='inactive')

          call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
               avgflag='A', long_name='C13 leaf C transfer', &
               ptr_pft=pc13s%leafc_xfer, default='inactive')

          call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
               avgflag='A', long_name='C13 fine root C', &
               ptr_pft=pc13s%frootc)

          call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
               avgflag='A', long_name='C13 fine root C storage', &
               ptr_pft=pc13s%frootc_storage, default='inactive')

          call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
               avgflag='A', long_name='C13 fine root C transfer', &
               ptr_pft=pc13s%frootc_xfer, default='inactive')

          call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
               avgflag='A', long_name='C13 live stem C', &
               ptr_pft=pc13s%livestemc)

          call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
               avgflag='A', long_name='C13 live stem C storage', &
               ptr_pft=pc13s%livestemc_storage, default='inactive')

          call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
               avgflag='A', long_name='C13 live stem C transfer', &
               ptr_pft=pc13s%livestemc_xfer, default='inactive')

          call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
               avgflag='A', long_name='C13 dead stem C', &
               ptr_pft=pc13s%deadstemc)

          call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
               avgflag='A', long_name='C13 dead stem C storage', &
               ptr_pft=pc13s%deadstemc_storage, default='inactive')

          call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
               avgflag='A', long_name='C13 dead stem C transfer', &
               ptr_pft=pc13s%deadstemc_xfer, default='inactive')

          call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
               avgflag='A', long_name='C13 live coarse root C', &
               ptr_pft=pc13s%livecrootc)

          call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
               avgflag='A', long_name='C13 live coarse root C storage', &
               ptr_pft=pc13s%livecrootc_storage, default='inactive')

          call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
               avgflag='A', long_name='C13 live coarse root C transfer', &
               ptr_pft=pc13s%livecrootc_xfer, default='inactive')

          call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
               avgflag='A', long_name='C13 dead coarse root C', &
               ptr_pft=pc13s%deadcrootc)

          call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
               avgflag='A', long_name='C13 dead coarse root C storage', &
               ptr_pft=pc13s%deadcrootc_storage,  default='inactive')

          call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
               avgflag='A', long_name='C13 dead coarse root C transfer', &
               ptr_pft=pc13s%deadcrootc_xfer, default='inactive')

          call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
               avgflag='A', long_name='C13 growth respiration storage', &
               ptr_pft=pc13s%gresp_storage, default='inactive')

          call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
               avgflag='A', long_name='C13 growth respiration transfer', &
               ptr_pft=pc13s%gresp_xfer, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
               avgflag='A', long_name='C13 temporary photosynthate C pool', &
               ptr_pft=pc13s%cpool)

          call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
               avgflag='A', long_name='C13 temporary photosynthate C pool', &
               ptr_pft=pc13s%xsmrpool)

          call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
               avgflag='A', long_name='C13 pft-level sink for C truncation', &
               ptr_pft=pc13s%pft_ctrunc)

          call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
               avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
               ptr_pft=pc13s%dispvegc)

          call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
               avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
               ptr_pft=pc13s%storvegc)

          call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
               avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
               ptr_pft=pc13s%totvegc)

          call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
               avgflag='A', long_name='C13 total pft-level carbon, including cpool', &
               ptr_pft=pc13s%totpftc)
       endif

       !-------------------------------
       ! C state variables - native to column
       !-------------------------------
       ! add history fields for all CLAMP CN variables
       call hist_addfld1d (fname='SOILC', units='gC/m^2', &
            avgflag='A', long_name='soil C', &
            ptr_col=ccs%totsomc)

       call hist_addfld1d (fname='LITTERC', units='gC/m^2', &
            avgflag='A', long_name='litter C', &
            ptr_col=ccs%totlitc)

       call hist_addfld1d (fname='CWDC', units='gC/m^2', &
            avgflag='A', long_name='coarse woody debris C', &
            ptr_col=ccs%cwdc)

       call hist_addfld1d (fname='LITR1C', units='gC/m^2', &
            avgflag='A', long_name='litter labile C', &
            ptr_col=ccs%litr1c)

       call hist_addfld1d (fname='LITR2C', units='gC/m^2', &
            avgflag='A', long_name='litter cellulose C', &
            ptr_col=ccs%litr2c)

       call hist_addfld1d (fname='LITR3C', units='gC/m^2', &
            avgflag='A', long_name='litter lignin C', &
            ptr_col=ccs%litr3c)

       call hist_addfld1d (fname='SOIL1C', units='gC/m^2', &
            avgflag='A', long_name='soil organic matter C (fast pool)', &
            ptr_col=ccs%soil1c)

       call hist_addfld1d (fname='SOIL2C', units='gC/m^2', &
            avgflag='A', long_name='soil organic matter C (medium pool)', &
            ptr_col=ccs%soil2c)

       call hist_addfld1d (fname='SOIL3C', units='gC/m^2', &
            avgflag='A', long_name='soil organic matter C (slow pool)', &
            ptr_col=ccs%soil3c)

       call hist_addfld1d (fname='SOIL4C', units='gC/m^2', &
            avgflag='A', long_name='soil organic matter C (slowest pool)', &
            ptr_col=ccs%soil4c)

       call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
            avgflag='A', long_name='pool for seeding new PFTs', &
            ptr_col=ccs%seedc)

       call hist_addfld1d (fname='COL_CTRUNC', units='gC/m^2', &
            avgflag='A', long_name='column-level sink for C truncation', &
            ptr_col=ccs%col_ctrunc)

       call hist_addfld1d (fname='TOTLITC', units='gC/m^2', &
            avgflag='A', long_name='total litter carbon', &
            ptr_col=ccs%totlitc)

       call hist_addfld1d (fname='TOTSOMC', units='gC/m^2', &
            avgflag='A', long_name='total soil organic matter carbon', &
            ptr_col=ccs%totsomc)

       call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
            avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=ccs%totecosysc)

       call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
            avgflag='A', long_name='total column carbon, incl veg and cpool', &
            ptr_col=ccs%totcolc)

       call hist_addfld1d (fname='PROD10C', units='gC/m^2', &
            avgflag='A', long_name='10-yr wood product C', &
            ptr_col=ccs%prod10c)

       call hist_addfld1d (fname='PROD100C', units='gC/m^2', &
            avgflag='A', long_name='100-yr wood product C', &
            ptr_col=ccs%prod100c)

       call hist_addfld1d (fname='TOTPRODC', units='gC/m^2', &
            avgflag='A', long_name='total wood product C', &
            ptr_col=ccs%totprodc)


       if (use_c13) then
          !-------------------------------
          ! C13 state variables - native to column
          !-------------------------------

          call hist_addfld1d (fname='C13_CWDC', units='gC13/m^2', &
               avgflag='A', long_name='C13 coarse woody debris C', &
               ptr_col=cc13s%cwdc)

          call hist_addfld1d (fname='C13_LITR1C', units='gC13/m^2', &
               avgflag='A', long_name='C13 litter labile C', &
               ptr_col=cc13s%litr1c)

          call hist_addfld1d (fname='C13_LITR2C', units='gC13/m^2', &
               avgflag='A', long_name='C13 litter cellulose C', &
               ptr_col=cc13s%litr2c)

          call hist_addfld1d (fname='C13_LITR3C', units='gC13/m^2', &
               avgflag='A', long_name='C13 litter lignin C', &
               ptr_col=cc13s%litr3c)

          call hist_addfld1d (fname='C13_SOIL1C', units='gC13/m^2', &
               avgflag='A', long_name='C13 soil organic matter C (fast pool)', &
               ptr_col=cc13s%soil1c)

          call hist_addfld1d (fname='C13_SOIL2C', units='gC13/m^2', &
               avgflag='A', long_name='C13 soil organic matter C (medium pool)', &
               ptr_col=cc13s%soil2c)

          call hist_addfld1d (fname='C13_SOIL3C', units='gC13/m^2', &
               avgflag='A', long_name='C13 soil organic matter C (slow pool)', &
               ptr_col=cc13s%soil3c)

          call hist_addfld1d (fname='C13_SOIL4C', units='gC13/m^2', &
               avgflag='A', long_name='C13 soil organic matter C (slowest pool)', &
               ptr_col=cc13s%soil4c)

          call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
               avgflag='A', long_name='C13 pool for seeding new PFTs', &
               ptr_col=ccs%seedc)

          call hist_addfld1d (fname='C13_COL_CTRUNC', units='gC13/m^2', &
               avgflag='A', long_name='C13 column-level sink for C truncation', &
               ptr_col=cc13s%col_ctrunc)

          call hist_addfld1d (fname='C13_TOTLITC', units='gC13/m^2', &
               avgflag='A', long_name='C13 total litter carbon', &
               ptr_col=cc13s%totlitc)

          call hist_addfld1d (fname='C13_TOTSOMC', units='gC13/m^2', &
               avgflag='A', long_name='C13 total soil organic matter carbon', &
               ptr_col=cc13s%totsomc)

          call hist_addfld1d (fname='C13_TOTECOSYSC', units='gC13/m^2', &
               avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool', &
               ptr_col=cc13s%totecosysc)

          call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
               avgflag='A', long_name='C13 total column carbon, incl veg and cpool', &
               ptr_col=cc13s%totcolc)

          call hist_addfld1d (fname='C13_PROD10C', units='gC13/m^2', &
               avgflag='A', long_name='C13 10-yr wood product C', &
               ptr_col=cc13s%prod10c)

          call hist_addfld1d (fname='C13_PROD100C', units='gC13/m^2', &
               avgflag='A', long_name='C13 100-yr wood product C', &
               ptr_col=cc13s%prod100c)

          call hist_addfld1d (fname='C13_TOTPRODC', units='gC13/m^2', &
               avgflag='A', long_name='C13 total wood product C', &
               ptr_col=cc13s%totprodc)
       endif

       !-------------------------------
       ! N state variables - native to PFT
       !-------------------------------

       call hist_addfld1d (fname='LEAFN', units='gN/m^2', &
            avgflag='A', long_name='leaf N', &
            ptr_pft=pns%leafn)

       call hist_addfld1d (fname='LEAFN_STORAGE', units='gN/m^2', &
            avgflag='A', long_name='leaf N storage', &
            ptr_pft=pns%leafn_storage, default='inactive')

       call hist_addfld1d (fname='LEAFN_XFER', units='gN/m^2', &
            avgflag='A', long_name='leaf N transfer', &
            ptr_pft=pns%leafn_xfer, default='inactive')

       call hist_addfld1d (fname='FROOTN', units='gN/m^2', &
            avgflag='A', long_name='fine root N', &
            ptr_pft=pns%frootn)

       call hist_addfld1d (fname='FROOTN_STORAGE', units='gN/m^2', &
            avgflag='A', long_name='fine root N storage', &
            ptr_pft=pns%frootn_storage, default='inactive')

       call hist_addfld1d (fname='FROOTN_XFER', units='gN/m^2', &
            avgflag='A', long_name='fine root N transfer', &
            ptr_pft=pns%frootn_xfer, default='inactive')

       call hist_addfld1d (fname='LIVESTEMN', units='gN/m^2', &
            avgflag='A', long_name='live stem N', &
            ptr_pft=pns%livestemn)

       call hist_addfld1d (fname='LIVESTEMN_STORAGE', units='gN/m^2', &
            avgflag='A', long_name='live stem N storage', &
            ptr_pft=pns%livestemn_storage, default='inactive')

       call hist_addfld1d (fname='LIVESTEMN_XFER', units='gN/m^2', &
            avgflag='A', long_name='live stem N transfer', &
            ptr_pft=pns%livestemn_xfer, default='inactive')

       call hist_addfld1d (fname='DEADSTEMN', units='gN/m^2', &
            avgflag='A', long_name='dead stem N', &
            ptr_pft=pns%deadstemn)

       call hist_addfld1d (fname='DEADSTEMN_STORAGE', units='gN/m^2', &
            avgflag='A', long_name='dead stem N storage', &
            ptr_pft=pns%deadstemn_storage, default='inactive')

       call hist_addfld1d (fname='DEADSTEMN_XFER', units='gN/m^2', &
            avgflag='A', long_name='dead stem N transfer', &
            ptr_pft=pns%deadstemn_xfer, default='inactive')

       call hist_addfld1d (fname='LIVECROOTN', units='gN/m^2', &
            avgflag='A', long_name='live coarse root N', &
            ptr_pft=pns%livecrootn)

       call hist_addfld1d (fname='LIVECROOTN_STORAGE', units='gN/m^2', &
            avgflag='A', long_name='live coarse root N storage', &
            ptr_pft=pns%livecrootn_storage, default='inactive')

       call hist_addfld1d (fname='LIVECROOTN_XFER', units='gN/m^2', &
            avgflag='A', long_name='live coarse root N transfer', &
            ptr_pft=pns%livecrootn_xfer, default='inactive')

       call hist_addfld1d (fname='DEADCROOTN', units='gN/m^2', &
            avgflag='A', long_name='dead coarse root N', &
            ptr_pft=pns%deadcrootn)

       call hist_addfld1d (fname='DEADCROOTN_STORAGE', units='gN/m^2', &
            avgflag='A', long_name='dead coarse root N storage', &
            ptr_pft=pns%deadcrootn_storage, default='inactive')

       call hist_addfld1d (fname='DEADCROOTN_XFER', units='gN/m^2', &
            avgflag='A', long_name='dead coarse root N transfer', &
            ptr_pft=pns%deadcrootn_xfer, default='inactive')

       call hist_addfld1d (fname='RETRANSN', units='gN/m^2', &
            avgflag='A', long_name='plant pool of retranslocated N', &
            ptr_pft=pns%retransn)

       call hist_addfld1d (fname='NPOOL', units='gN/m^2', &
            avgflag='A', long_name='temporary plant N pool', &
            ptr_pft=pns%npool, default='inactive')

       call hist_addfld1d (fname='PFT_NTRUNC', units='gN/m^2', &
            avgflag='A', long_name='pft-level sink for N truncation', &
            ptr_pft=pns%pft_ntrunc)

       call hist_addfld1d (fname='DISPVEGN', units='gN/m^2', &
            avgflag='A', long_name='displayed vegetation nitrogen', &
            ptr_pft=pns%dispvegn)

       call hist_addfld1d (fname='STORVEGN', units='gN/m^2', &
            avgflag='A', long_name='stored vegetation nitrogen', &
            ptr_pft=pns%storvegn)

       call hist_addfld1d (fname='TOTVEGN', units='gN/m^2', &
            avgflag='A', long_name='total vegetation nitrogen', &
            ptr_pft=pns%totvegn)

       call hist_addfld1d (fname='TOTPFTN', units='gN/m^2', &
            avgflag='A', long_name='total PFT-level nitrogen', &
            ptr_pft=pns%totpftn)

       !-------------------------------
       ! N state variables - native to column
       !-------------------------------

       call hist_addfld1d (fname='CWDN', units='gN/m^2', &
            avgflag='A', long_name='coarse woody debris N', &
            ptr_col=cns%cwdn)

       call hist_addfld1d (fname='LITR1N', units='gN/m^2', &
            avgflag='A', long_name='litter labile N', &
            ptr_col=cns%litr1n)

       call hist_addfld1d (fname='LITR2N', units='gN/m^2', &
            avgflag='A', long_name='litter cellulose N', &
            ptr_col=cns%litr2n)

       call hist_addfld1d (fname='LITR3N', units='gN/m^2', &
            avgflag='A', long_name='litter lignin N', &
            ptr_col=cns%litr3n)

       call hist_addfld1d (fname='SOIL1N', units='gN/m^2', &
            avgflag='A', long_name='soil organic matter N (fast pool)', &
            ptr_col=cns%soil1n)

       call hist_addfld1d (fname='SOIL2N', units='gN/m^2', &
            avgflag='A', long_name='soil organic matter N (medium pool)', &
            ptr_col=cns%soil2n)

       call hist_addfld1d (fname='SOIL3N', units='gN/m^2', &
            avgflag='A', long_name='soil orgainc matter N (slow pool)', &
            ptr_col=cns%soil3n)

       call hist_addfld1d (fname='SOIL4N', units='gN/m^2', &
            avgflag='A', long_name='soil orgainc matter N (slowest pool)', &
            ptr_col=cns%soil4n)

       call hist_addfld1d (fname='SMINN', units='gN/m^2', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=cns%sminn)

       call hist_addfld1d (fname='COL_NTRUNC', units='gN/m^2', &
            avgflag='A', long_name='column-level sink for N truncation', &
            ptr_col=cns%col_ntrunc)

       call hist_addfld1d (fname='TOTLITN', units='gN/m^2', &
            avgflag='A', long_name='total litter N', &
            ptr_col=cns%totlitn)

       call hist_addfld1d (fname='TOTSOMN', units='gN/m^2', &
            avgflag='A', long_name='total soil organic matter N', &
            ptr_col=cns%totsomn)

       call hist_addfld1d (fname='TOTECOSYSN', units='gN/m^2', &
            avgflag='A', long_name='total ecosystem N', &
            ptr_col=cns%totecosysn)

       call hist_addfld1d (fname='TOTCOLN', units='gN/m^2', &
            avgflag='A', long_name='total column-level N', &
            ptr_col=cns%totcoln)

       call hist_addfld1d (fname='SEEDN', units='gN/m^2', &
            avgflag='A', long_name='pool for seeding new PFTs ', &
            ptr_col=cns%seedn)

       call hist_addfld1d (fname='PROD10N', units='gN/m^2', &
            avgflag='A', long_name='10-yr wood product N', &
            ptr_col=cns%prod10n)

       call hist_addfld1d (fname='PROD100N', units='gN/m^2', &
            avgflag='A', long_name='100-yr wood product N', &
            ptr_col=cns%prod100n)

       call hist_addfld1d (fname='TOTPRODN', units='gN/m^2', &
            avgflag='A', long_name='total wood product N', &
            ptr_col=cns%totprodn)

       !-------------------------------
       ! C flux variables - native to PFT
       !-------------------------------

       ! add history fields for all CLAMP CN variables

       call hist_addfld1d (fname='WOODC_ALLOC', units='gC/m^2/s', &
            avgflag='A', long_name='wood C allocation', &
            ptr_pft=pcf%woodc_alloc)

       call hist_addfld1d (fname='WOODC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='wood C loss', &
            ptr_pft=pcf%woodc_loss)

       call hist_addfld1d (fname='LEAFC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C loss', &
            ptr_pft=pcf%leafc_loss)

       call hist_addfld1d (fname='LEAFC_ALLOC', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C allocation', &
            ptr_pft=pcf%leafc_alloc)

       call hist_addfld1d (fname='FROOTC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C loss', &
            ptr_pft=pcf%frootc_loss)

       call hist_addfld1d (fname='FROOTC_ALLOC', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C allocation', &
            ptr_pft=pcf%frootc_alloc)

       call hist_addfld1d (fname='PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='sunlit leaf photosynthesis', &
            ptr_pft=pcf%psnsun)

       call hist_addfld1d (fname='PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='shaded leaf photosynthesis', &
            ptr_pft=pcf%psnsha)

       call hist_addfld1d (fname='M_LEAFC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C mortality', &
            ptr_pft=pcf%m_leafc_to_litter, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C mortality', &
            ptr_pft=pcf%m_frootc_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C storage mortality', &
            ptr_pft=pcf%m_leafc_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C storage mortality', &
            ptr_pft=pcf%m_frootc_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C storage mortality', &
            ptr_pft=pcf%m_livestemc_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C storage mortality', &
            ptr_pft=pcf%m_deadstemc_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C storage mortality', &
            ptr_pft=pcf%m_livecrootc_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C storage mortality', &
            ptr_pft=pcf%m_deadcrootc_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C transfer mortality', &
            ptr_pft=pcf%m_leafc_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C transfer mortality', &
            ptr_pft=pcf%m_frootc_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C transfer mortality', &
            ptr_pft=pcf%m_livestemc_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C transfer mortality', &
            ptr_pft=pcf%m_deadstemc_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C transfer mortality', &
            ptr_pft=pcf%m_livecrootc_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C transfer mortality', &
            ptr_pft=pcf%m_deadcrootc_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C mortality', &
            ptr_pft=pcf%m_livestemc_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C mortality', &
            ptr_pft=pcf%m_deadstemc_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C mortality', &
            ptr_pft=pcf%m_livecrootc_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C mortality', &
            ptr_pft=pcf%m_deadcrootc_to_litter, default='inactive')

       call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration storage mortality', &
            ptr_pft=pcf%m_gresp_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_GRESP_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration transfer mortality', &
            ptr_pft=pcf%m_gresp_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LEAFC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C fire loss', &
            ptr_pft=pcf%m_leafc_to_fire, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C fire loss', &
            ptr_pft=pcf%m_frootc_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C storage fire loss', &
            ptr_pft=pcf%m_leafc_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C storage fire loss', &
            ptr_pft=pcf%m_frootc_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C storage fire loss', &
            ptr_pft=pcf%m_livestemc_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C storage fire loss', &
            ptr_pft=pcf%m_deadstemc_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C storage fire loss', &
            ptr_pft=pcf%m_livecrootc_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C storage fire loss', &
            ptr_pft=pcf%m_deadcrootc_storage_to_fire,  default='inactive')

       call hist_addfld1d (fname='M_LEAFC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C transfer fire loss', &
            ptr_pft=pcf%m_leafc_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C transfer fire loss', &
            ptr_pft=pcf%m_frootc_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C transfer fire loss', &
            ptr_pft=pcf%m_livestemc_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C transfer fire loss', &
            ptr_pft=pcf%m_deadstemc_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C transfer fire loss', &
            ptr_pft=pcf%m_livecrootc_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C transfer fire loss', &
            ptr_pft=pcf%m_deadcrootc_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C fire loss', &
            ptr_pft=pcf%m_livestemc_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C fire loss', &
            ptr_pft=pcf%m_deadstemc_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C fire mortality to litter', &
            ptr_pft=pcf%m_deadstemc_to_litter_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C fire loss', &
            ptr_pft=pcf%m_livecrootc_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C fire loss', &
            ptr_pft=pcf%m_deadcrootc_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C fire mortality to litter', &
            ptr_pft=pcf%m_deadcrootc_to_litter_fire, default='inactive')

       call hist_addfld1d (fname='M_GRESP_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration storage fire loss', &
            ptr_pft=pcf%m_gresp_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_GRESP_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration transfer fire loss', &
            ptr_pft=pcf%m_gresp_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='LEAFC_XFER_TO_LEAFC', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C growth from storage', &
            ptr_pft=pcf%leafc_xfer_to_leafc, default='inactive')

       call hist_addfld1d (fname='FROOTC_XFER_TO_FROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C growth from storage', &
            ptr_pft=pcf%frootc_xfer_to_frootc, default='inactive')

       call hist_addfld1d (fname='LIVESTEMC_XFER_TO_LIVESTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C growth from storage', &
            ptr_pft=pcf%livestemc_xfer_to_livestemc, default='inactive')

       call hist_addfld1d (fname='DEADSTEMC_XFER_TO_DEADSTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C growth from storage', &
            ptr_pft=pcf%deadstemc_xfer_to_deadstemc, default='inactive')

       call hist_addfld1d (fname='LIVECROOTC_XFER_TO_LIVECROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C growth from storage', &
            ptr_pft=pcf%livecrootc_xfer_to_livecrootc, default='inactive')

       call hist_addfld1d (fname='DEADCROOTC_XFER_TO_DEADCROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C growth from storage', &
            ptr_pft=pcf%deadcrootc_xfer_to_deadcrootc, default='inactive')

       call hist_addfld1d (fname='LEAFC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C litterfall', &
            ptr_pft=pcf%leafc_to_litter, default='inactive')

       call hist_addfld1d (fname='FROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C litterfall', &
            ptr_pft=pcf%frootc_to_litter, default='inactive')

       call hist_addfld1d (fname='LEAF_MR', units='gC/m^2/s', &
            avgflag='A', long_name='leaf maintenance respiration', &
            ptr_pft=pcf%leaf_mr, default='inactive')

       call hist_addfld1d (fname='FROOT_MR', units='gC/m^2/s', &
            avgflag='A', long_name='fine root maintenance respiration', &
            ptr_pft=pcf%froot_mr, default='inactive')

       call hist_addfld1d (fname='LIVESTEM_MR', units='gC/m^2/s', &
            avgflag='A', long_name='live stem maintenance respiration', &
            ptr_pft=pcf%livestem_mr, default='inactive')

       call hist_addfld1d (fname='LIVECROOT_MR', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root maintenance respiration', &
            ptr_pft=pcf%livecroot_mr, default='inactive')

       call hist_addfld1d (fname='PSNSUN_TO_CPOOL', units='gC/m^2/s', &
            avgflag='A', long_name='C fixation from sunlit canopy', &
            ptr_pft=pcf%psnsun_to_cpool)

       call hist_addfld1d (fname='PSNSHADE_TO_CPOOL', units='gC/m^2/s', &
            avgflag='A', long_name='C fixation from shaded canopy', &
            ptr_pft=pcf%psnshade_to_cpool)

       call hist_addfld1d (fname='CPOOL_TO_LEAFC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to leaf C', &
            ptr_pft=pcf%cpool_to_leafc, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_LEAFC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to leaf C storage', &
            ptr_pft=pcf%cpool_to_leafc_storage, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_FROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to fine root C', &
            ptr_pft=pcf%cpool_to_frootc, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_FROOTC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to fine root C storage', &
            ptr_pft=pcf%cpool_to_frootc_storage, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to live stem C', &
            ptr_pft=pcf%cpool_to_livestemc, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to live stem C storage', &
            ptr_pft=pcf%cpool_to_livestemc_storage, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to dead stem C', &
            ptr_pft=pcf%cpool_to_deadstemc, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to dead stem C storage', &
            ptr_pft=pcf%cpool_to_deadstemc_storage, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to live coarse root C', &
            ptr_pft=pcf%cpool_to_livecrootc, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to live coarse root C storage', &
            ptr_pft=pcf%cpool_to_livecrootc_storage, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to dead coarse root C', &
            ptr_pft=pcf%cpool_to_deadcrootc, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to dead coarse root C storage', &
            ptr_pft=pcf%cpool_to_deadcrootc_storage, default='inactive')

       call hist_addfld1d (fname='CPOOL_TO_GRESP_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to growth respiration storage', &
            ptr_pft=pcf%cpool_to_gresp_storage, default='inactive')

       call hist_addfld1d (fname='CPOOL_LEAF_GR', units='gC/m^2/s', &
            avgflag='A', long_name='leaf growth respiration', &
            ptr_pft=pcf%cpool_leaf_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_LEAF_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='leaf growth respiration to storage', &
            ptr_pft=pcf%cpool_leaf_storage_gr, default='inactive')

       call hist_addfld1d (fname='TRANSFER_LEAF_GR', units='gC/m^2/s', &
            avgflag='A', long_name='leaf growth respiration from storage', &
            ptr_pft=pcf%transfer_leaf_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_FROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='fine root growth respiration', &
            ptr_pft=pcf%cpool_froot_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_FROOT_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='fine root  growth respiration to storage', &
            ptr_pft=pcf%cpool_froot_storage_gr, default='inactive')

       call hist_addfld1d (fname='TRANSFER_FROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='fine root  growth respiration from storage', &
            ptr_pft=pcf%transfer_froot_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_LIVESTEM_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live stem growth respiration', &
            ptr_pft=pcf%cpool_livestem_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_LIVESTEM_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live stem growth respiration to storage', &
            ptr_pft=pcf%cpool_livestem_storage_gr, default='inactive')

       call hist_addfld1d (fname='TRANSFER_LIVESTEM_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live stem growth respiration from storage', &
            ptr_pft=pcf%transfer_livestem_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_DEADSTEM_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem growth respiration', &
            ptr_pft=pcf%cpool_deadstem_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_DEADSTEM_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem growth respiration to storage', &
            ptr_pft=pcf%cpool_deadstem_storage_gr, default='inactive')

       call hist_addfld1d (fname='TRANSFER_DEADSTEM_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem growth respiration from storage', &
            ptr_pft=pcf%transfer_deadstem_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_LIVECROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root growth respiration', &
            ptr_pft=pcf%cpool_livecroot_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_LIVECROOT_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root growth respiration to storage', &
            ptr_pft=pcf%cpool_livecroot_storage_gr, default='inactive')

       call hist_addfld1d (fname='TRANSFER_LIVECROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root growth respiration from storage', &
            ptr_pft=pcf%transfer_livecroot_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_DEADCROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root growth respiration', &
            ptr_pft=pcf%cpool_deadcroot_gr, default='inactive')

       call hist_addfld1d (fname='CPOOL_DEADCROOT_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root growth respiration to storage', &
            ptr_pft=pcf%cpool_deadcroot_storage_gr, default='inactive')

       call hist_addfld1d (fname='TRANSFER_DEADCROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root growth respiration from storage', &
            ptr_pft=pcf%transfer_deadcroot_gr, default='inactive')

       call hist_addfld1d (fname='LEAFC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C shift storage to transfer', &
            ptr_pft=pcf%leafc_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='FROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C shift storage to transfer', &
            ptr_pft=pcf%frootc_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='LIVESTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C shift storage to transfer', &
            ptr_pft=pcf%livestemc_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='DEADSTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C shift storage to transfer', &
            ptr_pft=pcf%deadstemc_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='LIVECROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C shift storage to transfer', &
            ptr_pft=pcf%livecrootc_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='DEADCROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C shift storage to transfer', &
            ptr_pft=pcf%deadcrootc_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='GRESP_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration shift storage to transfer', &
            ptr_pft=pcf%gresp_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='LIVESTEMC_TO_DEADSTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C turnover', &
            ptr_pft=pcf%livestemc_to_deadstemc, default='inactive')

       call hist_addfld1d (fname='LIVECROOTC_TO_DEADCROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C turnover', &
            ptr_pft=pcf%livecrootc_to_deadcrootc, default='inactive')

       call hist_addfld1d (fname='GPP', units='gC/m^2/s', &
            avgflag='A', long_name='gross primary production', &
            ptr_pft=pcf%gpp)

       call hist_addfld1d (fname='MR', units='gC/m^2/s', &
            avgflag='A', long_name='maintenance respiration', &
            ptr_pft=pcf%mr)

       call hist_addfld1d (fname='CURRENT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='growth resp for new growth displayed in this timestep', &
            ptr_pft=pcf%current_gr, default='inactive')

       call hist_addfld1d (fname='TRANSFER_GR', units='gC/m^2/s', &
            avgflag='A', long_name='growth resp for transfer growth displayed in this timestep', &
            ptr_pft=pcf%transfer_gr, default='inactive')

       call hist_addfld1d (fname='STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='growth resp for growth sent to storage for later display', &
            ptr_pft=pcf%storage_gr, default='inactive')

       call hist_addfld1d (fname='GR', units='gC/m^2/s', &
            avgflag='A', long_name='total growth respiration', &
            ptr_pft=pcf%gr)

       call hist_addfld1d (fname='AR', units='gC/m^2/s', &
            avgflag='A', long_name='autotrophic respiration (MR + GR)', &
            ptr_pft=pcf%ar)

       call hist_addfld1d (fname='RR', units='gC/m^2/s', &
            avgflag='A', long_name='root respiration (fine root MR + total root GR)', &
            ptr_pft=pcf%rr)

       call hist_addfld1d (fname='NPP', units='gC/m^2/s', &
            avgflag='A', long_name='net primary production', &
            ptr_pft=pcf%npp)

       call hist_addfld1d (fname='AGNPP', units='gC/m^2/s', &
            avgflag='A', long_name='aboveground NPP', &
            ptr_pft=pcf%agnpp)

       call hist_addfld1d (fname='BGNPP', units='gC/m^2/s', &
            avgflag='A', long_name='belowground NPP', &
            ptr_pft=pcf%bgnpp)

       call hist_addfld1d (fname='LITFALL', units='gC/m^2/s', &
            avgflag='A', long_name='litterfall (leaves and fine roots)', &
            ptr_pft=pcf%litfall)

       call hist_addfld1d (fname='VEGFIRE', units='gC/m^2/s', &
            avgflag='A', long_name='pft-level fire loss', &
            ptr_pft=pcf%vegfire, default='inactive')

       call hist_addfld1d (fname='WOOD_HARVESTC', units='gC/m^2/s', &
            avgflag='A', long_name='wood harvest carbon (to product pools)', &
            ptr_pft=pcf%wood_harvestc)

       call hist_addfld1d (fname='PFT_FIRE_CLOSS', units='gC/m^2/s', &
            avgflag='A', long_name='total pft-level fire C loss', &
            ptr_pft=pcf%pft_fire_closs)

       if (use_c13) then
          !-------------------------------
          ! C13 flux variables - native to PFT
          !-------------------------------

          call hist_addfld1d (fname='C13_PSNSUN', units='umolCO2/m^2/s', &
               avgflag='A', long_name='C13 sunlit leaf photosynthesis', &
               ptr_pft=pc13f%psnsun)

          call hist_addfld1d (fname='C13_PSNSHA', units='umolCO2/m^2/s', &
               avgflag='A', long_name='C13 shaded leaf photosynthesis', &
               ptr_pft=pc13f%psnsha)

          call hist_addfld1d (fname='C13_M_LEAFC_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C mortality', &
               ptr_pft=pc13f%m_leafc_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C mortality', &
               ptr_pft=pc13f%m_frootc_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C storage mortality', &
               ptr_pft=pc13f%m_leafc_storage_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C storage mortality', &
               ptr_pft=pc13f%m_frootc_storage_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C storage mortality', &
               ptr_pft=pc13f%m_livestemc_storage_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C storage mortality', &
               ptr_pft=pc13f%m_deadstemc_storage_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C storage mortality', &
               ptr_pft=pc13f%m_livecrootc_storage_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C storage mortality', &
               ptr_pft=pc13f%m_deadcrootc_storage_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C transfer mortality', &
               ptr_pft=pc13f%m_leafc_xfer_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C transfer mortality', &
               ptr_pft=pc13f%m_frootc_xfer_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C transfer mortality', &
               ptr_pft=pc13f%m_livestemc_xfer_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C transfer mortality', &
               ptr_pft=pc13f%m_deadstemc_xfer_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C transfer mortality', &
               ptr_pft=pc13f%m_livecrootc_xfer_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C transfer mortality', &
               ptr_pft=pc13f%m_deadcrootc_xfer_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C mortality', &
               ptr_pft=pc13f%m_livestemc_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C mortality', &
               ptr_pft=pc13f%m_deadstemc_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C mortality', &
               ptr_pft=pc13f%m_livecrootc_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C mortality', &
               ptr_pft=pc13f%m_deadcrootc_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth respiration storage mortality', &
               ptr_pft=pc13f%m_gresp_storage_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth respiration transfer mortality', &
               ptr_pft=pc13f%m_gresp_xfer_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C fire loss', &
               ptr_pft=pc13f%m_leafc_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C fire loss', &
               ptr_pft=pc13f%m_frootc_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C storage fire loss', &
               ptr_pft=pc13f%m_leafc_storage_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C storage fire loss', &
               ptr_pft=pc13f%m_frootc_storage_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C storage fire loss', &
               ptr_pft=pc13f%m_livestemc_storage_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C storage fire loss', &
               ptr_pft=pc13f%m_deadstemc_storage_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C storage fire loss', &
               ptr_pft=pc13f%m_livecrootc_storage_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C storage fire loss', &
               ptr_pft=pc13f%m_deadcrootc_storage_to_fire,  default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C transfer fire loss', &
               ptr_pft=pc13f%m_leafc_xfer_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C transfer fire loss', &
               ptr_pft=pc13f%m_frootc_xfer_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C transfer fire loss', &
               ptr_pft=pc13f%m_livestemc_xfer_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C transfer fire loss', &
               ptr_pft=pc13f%m_deadstemc_xfer_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C transfer fire loss', &
               ptr_pft=pc13f%m_livecrootc_xfer_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C transfer fire loss', &
               ptr_pft=pc13f%m_deadcrootc_xfer_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C fire loss', &
               ptr_pft=pc13f%m_livestemc_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C fire loss', &
               ptr_pft=pc13f%m_deadstemc_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C fire mortality to litter', &
               ptr_pft=pc13f%m_deadstemc_to_litter_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C fire loss', &
               ptr_pft=pc13f%m_livecrootc_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C fire loss', &
               ptr_pft=pc13f%m_deadcrootc_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C fire mortality to litter', &
               ptr_pft=pc13f%m_deadcrootc_to_litter_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth respiration storage fire loss', &
               ptr_pft=pc13f%m_gresp_storage_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth respiration transfer fire loss', &
               ptr_pft=pc13f%m_gresp_xfer_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_LEAFC_XFER_TO_LEAFC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C growth from storage', &
               ptr_pft=pc13f%leafc_xfer_to_leafc, default='inactive')

          call hist_addfld1d (fname='C13_FROOTC_XFER_TO_FROOTC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C growth from storage', &
               ptr_pft=pc13f%frootc_xfer_to_frootc, default='inactive')

          call hist_addfld1d (fname='C13_LIVESTEMC_XFER_TO_LIVESTEMC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C growth from storage', &
               ptr_pft=pc13f%livestemc_xfer_to_livestemc, default='inactive')

          call hist_addfld1d (fname='C13_DEADSTEMC_XFER_TO_DEADSTEMC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C growth from storage', &
               ptr_pft=pc13f%deadstemc_xfer_to_deadstemc, default='inactive')

          call hist_addfld1d (fname='C13_LIVECROOTC_XFER_TO_LIVECROOTC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C growth from storage', &
               ptr_pft=pc13f%livecrootc_xfer_to_livecrootc, default='inactive')

          call hist_addfld1d (fname='C13_DEADCROOTC_XFER_TO_DEADCROOTC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C growth from storage', &
               ptr_pft=pc13f%deadcrootc_xfer_to_deadcrootc, default='inactive')

          call hist_addfld1d (fname='C13_LEAFC_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C litterfall', &
               ptr_pft=pc13f%leafc_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_FROOTC_TO_LITTER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C litterfall', &
               ptr_pft=pc13f%frootc_to_litter, default='inactive')

          call hist_addfld1d (fname='C13_LEAF_MR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf maintenance respiration', &
               ptr_pft=pc13f%leaf_mr, default='inactive')

          call hist_addfld1d (fname='C13_FROOT_MR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root maintenance respiration', &
               ptr_pft=pc13f%froot_mr, default='inactive')

          call hist_addfld1d (fname='C13_LIVESTEM_MR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem maintenance respiration', &
               ptr_pft=pc13f%livestem_mr, default='inactive')

          call hist_addfld1d (fname='C13_LIVECROOT_MR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root maintenance respiration', &
               ptr_pft=pc13f%livecroot_mr, default='inactive')

          call hist_addfld1d (fname='C13_PSNSUN_TO_CPOOL', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 C fixation from sunlit canopy', &
               ptr_pft=pc13f%psnsun_to_cpool)

          call hist_addfld1d (fname='C13_PSNSHADE_TO_CPOOL', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 C fixation from shaded canopy', &
               ptr_pft=pc13f%psnshade_to_cpool)

          call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to leaf C', &
               ptr_pft=pc13f%cpool_to_leafc, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC_STORAGE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to leaf C storage', &
               ptr_pft=pc13f%cpool_to_leafc_storage, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to fine root C', &
               ptr_pft=pc13f%cpool_to_frootc, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC_STORAGE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to fine root C storage', &
               ptr_pft=pc13f%cpool_to_frootc_storage, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to live stem C', &
               ptr_pft=pc13f%cpool_to_livestemc, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC_STORAGE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to live stem C storage', &
               ptr_pft=pc13f%cpool_to_livestemc_storage, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to dead stem C', &
               ptr_pft=pc13f%cpool_to_deadstemc, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC_STORAGE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to dead stem C storage', &
               ptr_pft=pc13f%cpool_to_deadstemc_storage, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to live coarse root C', &
               ptr_pft=pc13f%cpool_to_livecrootc, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC_STORAGE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to live coarse root C storage', &
               ptr_pft=pc13f%cpool_to_livecrootc_storage, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to dead coarse root C', &
               ptr_pft=pc13f%cpool_to_deadcrootc, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC_STORAGE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to dead coarse root C storage', &
               ptr_pft=pc13f%cpool_to_deadcrootc_storage, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_TO_GRESP_STORAGE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 allocation to growth respiration storage', &
               ptr_pft=pc13f%cpool_to_gresp_storage, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_LEAF_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf growth respiration', &
               ptr_pft=pc13f%cpool_leaf_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_LEAF_STORAGE_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf growth respiration to storage', &
               ptr_pft=pc13f%cpool_leaf_storage_gr, default='inactive')

          call hist_addfld1d (fname='C13_TRANSFER_LEAF_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf growth respiration from storage', &
               ptr_pft=pc13f%transfer_leaf_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_FROOT_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root growth respiration', &
               ptr_pft=pc13f%cpool_froot_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_FROOT_STORAGE_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root  growth respiration to storage', &
               ptr_pft=pc13f%cpool_froot_storage_gr, default='inactive')

          call hist_addfld1d (fname='C13_TRANSFER_FROOT_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root  growth respiration from storage', &
               ptr_pft=pc13f%transfer_froot_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem growth respiration', &
               ptr_pft=pc13f%cpool_livestem_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_STORAGE_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem growth respiration to storage', &
               ptr_pft=pc13f%cpool_livestem_storage_gr, default='inactive')

          call hist_addfld1d (fname='C13_TRANSFER_LIVESTEM_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem growth respiration from storage', &
               ptr_pft=pc13f%transfer_livestem_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem growth respiration', &
               ptr_pft=pc13f%cpool_deadstem_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_STORAGE_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem growth respiration to storage', &
               ptr_pft=pc13f%cpool_deadstem_storage_gr, default='inactive')

          call hist_addfld1d (fname='C13_TRANSFER_DEADSTEM_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem growth respiration from storage', &
               ptr_pft=pc13f%transfer_deadstem_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root growth respiration', &
               ptr_pft=pc13f%cpool_livecroot_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_STORAGE_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root growth respiration to storage', &
               ptr_pft=pc13f%cpool_livecroot_storage_gr, default='inactive')

          call hist_addfld1d (fname='C13_TRANSFER_LIVECROOT_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root growth respiration from storage', &
               ptr_pft=pc13f%transfer_livecroot_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root growth respiration', &
               ptr_pft=pc13f%cpool_deadcroot_gr, default='inactive')

          call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_STORAGE_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root growth respiration to storage', &
               ptr_pft=pc13f%cpool_deadcroot_storage_gr, default='inactive')

          call hist_addfld1d (fname='C13_TRANSFER_DEADCROOT_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root growth respiration from storage', &
               ptr_pft=pc13f%transfer_deadcroot_gr, default='inactive')

          call hist_addfld1d (fname='C13_LEAFC_STORAGE_TO_XFER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C shift storage to transfer', &
               ptr_pft=pc13f%leafc_storage_to_xfer, default='inactive')

          call hist_addfld1d (fname='C13_FROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C shift storage to transfer', &
               ptr_pft=pc13f%frootc_storage_to_xfer, default='inactive')

          call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C shift storage to transfer', &
               ptr_pft=pc13f%livestemc_storage_to_xfer, default='inactive')

          call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C shift storage to transfer', &
               ptr_pft=pc13f%deadstemc_storage_to_xfer, default='inactive')

          call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C shift storage to transfer', &
               ptr_pft=pc13f%livecrootc_storage_to_xfer, default='inactive')

          call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C shift storage to transfer', &
               ptr_pft=pc13f%deadcrootc_storage_to_xfer, default='inactive')

          call hist_addfld1d (fname='C13_GRESP_STORAGE_TO_XFER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth respiration shift storage to transfer', &
               ptr_pft=pc13f%gresp_storage_to_xfer, default='inactive')

          call hist_addfld1d (fname='C13_LIVESTEMC_TO_DEADSTEMC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C turnover', &
               ptr_pft=pc13f%livestemc_to_deadstemc, default='inactive')

          call hist_addfld1d (fname='C13_LIVECROOTC_TO_DEADCROOTC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C turnover', &
               ptr_pft=pc13f%livecrootc_to_deadcrootc, default='inactive')

          call hist_addfld1d (fname='C13_GPP', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 gross primary production', &
               ptr_pft=pc13f%gpp)

          call hist_addfld1d (fname='C13_MR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 maintenance respiration', &
               ptr_pft=pc13f%mr)

          call hist_addfld1d (fname='C13_CURRENT_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth resp for new growth displayed in this timestep', &
               ptr_pft=pc13f%current_gr, default='inactive')

          call hist_addfld1d (fname='C13_TRANSFER_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth resp for transfer growth displayed in this timestep', &
               ptr_pft=pc13f%transfer_gr, default='inactive')

          call hist_addfld1d (fname='C13_STORAGE_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth resp for growth sent to storage for later display', &
               ptr_pft=pc13f%storage_gr, default='inactive')

          call hist_addfld1d (fname='C13_GR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total growth respiration', &
               ptr_pft=pc13f%gr)

          call hist_addfld1d (fname='C13_AR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 autotrophic respiration (MR + GR)', &
               ptr_pft=pc13f%ar)

          call hist_addfld1d (fname='C13_RR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 root respiration (fine root MR + total root GR)', &
               ptr_pft=pc13f%rr)

          call hist_addfld1d (fname='C13_NPP', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 net primary production', &
               ptr_pft=pc13f%npp)

          call hist_addfld1d (fname='C13_AGNPP', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 aboveground NPP', &
               ptr_pft=pc13f%agnpp)

          call hist_addfld1d (fname='C13_BGNPP', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 belowground NPP', &
               ptr_pft=pc13f%bgnpp)

          call hist_addfld1d (fname='C13_LITFALL', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 litterfall (leaves and fine roots)', &
               ptr_pft=pc13f%litfall, default='inactive')

          call hist_addfld1d (fname='C13_VEGFIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 pft-level fire loss', &
               ptr_pft=pc13f%vegfire, default='inactive')

          call hist_addfld1d (fname='C13_PFT_FIRE_CLOSS', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total pft-level fire C loss', &
               ptr_pft=pc13f%pft_fire_closs)
       endif

       !-------------------------------
       ! C flux variables - native to column 
       !-------------------------------
       ! add history fields for all CLAMP CN variables

       call hist_addfld1d (fname='CWDC_HR', units='gC/m^2/s', &
            avgflag='A', long_name='coarse woody debris C heterotrophic respiration', &
            ptr_col=ccf%cwdc_hr)

       call hist_addfld1d (fname='CWDC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='coarse woody debris C loss', &
            ptr_col=ccf%cwdc_loss)

       call hist_addfld1d (fname='LITTERC_HR', units='gC/m^2/s', &
            avgflag='A', long_name='litter C heterotrophic respiration', &
            ptr_col=ccf%lithr)

       call hist_addfld1d (fname='LITTERC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='litter C loss', &
            ptr_col=ccf%litterc_loss)

       call hist_addfld1d (fname='SOILC_HR', units='gC/m^2/s', &
            avgflag='A', long_name='soil C heterotrophic respiration', &
            ptr_col=ccf%somhr)

       call hist_addfld1d (fname='SOILC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='soil C loss', &
            ptr_col=ccf%somhr)

       call hist_addfld1d (fname='M_LEAFC_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C mortality to litter 1 C', &
            ptr_col=ccf%m_leafc_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_LEAFC_TO_LITR2C', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C mortality to litter 2 C', &
            ptr_col=ccf%m_leafc_to_litr2c, default='inactive')

       call hist_addfld1d (fname='M_LEAFC_TO_LITR3C', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C mortality to litter 3 C', &
            ptr_col=ccf%m_leafc_to_litr3c, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C mortality to litter 1 C', &
            ptr_col=ccf%m_frootc_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_TO_LITR2C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C mortality to litter 2 C', &
            ptr_col=ccf%m_frootc_to_litr2c, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_TO_LITR3C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C mortality to litter 3 C', &
            ptr_col=ccf%m_frootc_to_litr3c, default='inactive')

       call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C storage mortality to litter 1 C', &
            ptr_col=ccf%m_leafc_storage_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C storage mortality to litter 1 C', &
            ptr_col=ccf%m_frootc_storage_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C storage mortality to litter 1 C', &
            ptr_col=ccf%m_livestemc_storage_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C storage mortality to litter 1 C', &
            ptr_col=ccf%m_deadstemc_storage_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C storage mortality to litter 1 C', &
            ptr_col=ccf%m_livecrootc_storage_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C storage mortality to litter 1 C', &
            ptr_col=ccf%m_deadcrootc_storage_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C transfer mortality to litter 1 C', &
            ptr_col=ccf%m_leafc_xfer_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C transfer mortality to litter 1 C', &
            ptr_col=ccf%m_frootc_xfer_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C transfer mortality to litter 1 C', &
            ptr_col=ccf%m_livestemc_xfer_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C transfer mortality to litter 1 C', &
            ptr_col=ccf%m_deadstemc_xfer_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C transfer mortality to litter 1 C', &
            ptr_col=ccf%m_livecrootc_xfer_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C transfer mortality to litter 1 C', &
            ptr_col=ccf%m_deadcrootc_xfer_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMC_TO_CWDC', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C mortality to coarse woody debris C', &
            ptr_col=ccf%m_livestemc_to_cwdc, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_TO_CWDC', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C mortality to coarse woody debris C', &
            ptr_col=ccf%m_deadstemc_to_cwdc, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTC_TO_CWDC', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C mortality to coarse woody debris C', &
            ptr_col=ccf%m_livecrootc_to_cwdc, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_TO_CWDC', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C mortality to coarse woody debris C', &
            ptr_col=ccf%m_deadcrootc_to_cwdc, default='inactive')

       call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration storage mortality to litter 1 C', &
            ptr_col=ccf%m_gresp_storage_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_GRESP_XFER_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration transfer mortality to litter 1 C', &
            ptr_col=ccf%m_gresp_xfer_to_litr1c, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMC_TO_CWDC_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C to coarse woody debris C by fire', &
            ptr_col=ccf%m_deadstemc_to_cwdc_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTC_TO_CWDC_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C to to woody debris C by fire', &
            ptr_col=ccf%m_deadcrootc_to_cwdc_fire, default='inactive')

       call hist_addfld1d (fname='M_LITR1C_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='litter 1 C fire loss', &
            ptr_col=ccf%m_litr1c_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LITR2C_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='litter 2 C fire loss', &
            ptr_col=ccf%m_litr2c_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LITR3C_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='litter 3 C fire loss', &
            ptr_col=ccf%m_litr3c_to_fire, default='inactive')

       call hist_addfld1d (fname='M_CWDC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='coarse woody debris C fire loss', &
            ptr_col=ccf%m_cwdc_to_fire, default='inactive')

       call hist_addfld1d (fname='LEAFC_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C litterfall to litter 1 C', &
            ptr_col=ccf%leafc_to_litr1c, default='inactive')

       call hist_addfld1d (fname='LEAFC_TO_LITR2C', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C litterfall to litter 2 C', &
            ptr_col=ccf%leafc_to_litr2c, default='inactive')

       call hist_addfld1d (fname='LEAFC_TO_LITR3C', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C litterfall to litter 3 C', &
            ptr_col=ccf%leafc_to_litr3c, default='inactive')

       call hist_addfld1d (fname='FROOTC_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C litterfall to litter 1 C', &
            ptr_col=ccf%frootc_to_litr1c, default='inactive')

       call hist_addfld1d (fname='FROOTC_TO_LITR2C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C litterfall to litter 2 C', &
            ptr_col=ccf%frootc_to_litr2c, default='inactive')

       call hist_addfld1d (fname='FROOTC_TO_LITR3C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C litterfall to litter 3 C', &
            ptr_col=ccf%frootc_to_litr3c, default='inactive')

       call hist_addfld1d (fname='CWDC_TO_LITR2C', units='gC/m^2/s', &
            avgflag='A', long_name='decomp. of coarse woody debris C to litter 2 C', &
            ptr_col=ccf%cwdc_to_litr2c, default='inactive')

       call hist_addfld1d (fname='CWDC_TO_LITR3C', units='gC/m^2/s', &
            avgflag='A', long_name='decomp. of coarse woody debris C to litter 3 C', &
            ptr_col=ccf%cwdc_to_litr3c, default='inactive')

       call hist_addfld1d (fname='LITR1_HR', units='gC/m^2/s', &
            avgflag='A', long_name='het. resp. from litter 1 C', &
            ptr_col=ccf%litr1_hr, default='inactive')

       call hist_addfld1d (fname='LITR1C_TO_SOIL1C', units='gC/m^2/s', &
            avgflag='A', long_name='decomp. of litter 1 C to SOM 1 C', &
            ptr_col=ccf%litr1c_to_soil1c)

       call hist_addfld1d (fname='LITR2_HR', units='gC/m^2/s', &
            avgflag='A', long_name='het. resp. from litter 2 C', &
            ptr_col=ccf%litr2_hr, default='inactive')

       call hist_addfld1d (fname='LITR2C_TO_SOIL2C', units='gC/m^2/s', &
            avgflag='A', long_name='decomp. of litter 2 C to SOM 2 C', &
            ptr_col=ccf%litr2c_to_soil2c)

       call hist_addfld1d (fname='LITR3_HR', units='gC/m^2/s', &
            avgflag='A', long_name='het. resp. from litter 3 C', &
            ptr_col=ccf%litr3_hr, default='inactive')

       call hist_addfld1d (fname='LITR3C_TO_SOIL3C', units='gC/m^2/s', &
            avgflag='A', long_name='decomp. of litter 3 C to SOM 3 C', &
            ptr_col=ccf%litr3c_to_soil3c)

       call hist_addfld1d (fname='SOIL1_HR', units='gC/m^2/s', &
            avgflag='A', long_name='het. resp. from SOM 1 C', &
            ptr_col=ccf%soil1_hr, default='inactive')

       call hist_addfld1d (fname='SOIL1C_TO_SOIL2C', units='gC/m^2/s', &
            avgflag='A', long_name='decomp. of SOM 1 C to SOM 2 C', &
            ptr_col=ccf%soil1c_to_soil2c, default='inactive')

       call hist_addfld1d (fname='SOIL2_HR', units='gC/m^2/s', &
            avgflag='A', long_name='het. resp. from SOM 2 C', &
            ptr_col=ccf%soil2_hr, default='inactive')

       call hist_addfld1d (fname='SOIL2C_TO_SOIL3C', units='gC/m^2/s', &
            avgflag='A', long_name='decomp. of SOM 2 C to SOM 3 C', &
            ptr_col=ccf%soil2c_to_soil3c, default='inactive')

       call hist_addfld1d (fname='SOIL3_HR', units='gC/m^2/s', &
            avgflag='A', long_name='het. resp. from SOM 3 C', &
            ptr_col=ccf%soil3_hr, default='inactive')

       call hist_addfld1d (fname='SOIL3C_TO_SOIL4C', units='gC/m^2/s', &
            avgflag='A', long_name='decomp. of SOM 3 C to SOM 4 C', &
            ptr_col=ccf%soil3c_to_soil4c, default='inactive')

       call hist_addfld1d (fname='SOIL4_HR', units='gC/m^2/s', &
            avgflag='A', long_name='het. resp. from SOM 4 C', &
            ptr_col=ccf%soil4_hr, default='inactive')

       call hist_addfld1d (fname='LITHR', units='gC/m^2/s', &
            avgflag='A', long_name='litter heterotrophic respiration', &
            ptr_col=ccf%lithr)

       call hist_addfld1d (fname='SOMHR', units='gC/m^2/s', &
            avgflag='A', long_name='soil organic matter heterotrophic respiration', &
            ptr_col=ccf%somhr)

       call hist_addfld1d (fname='HR', units='gC/m^2/s', &
            avgflag='A', long_name='total heterotrophic respiration', &
            ptr_col=ccf%hr)

       call hist_addfld1d (fname='SR', units='gC/m^2/s', &
            avgflag='A', long_name='total soil respiration (HR + root resp)', &
            ptr_col=ccf%sr)

       call hist_addfld1d (fname='ER', units='gC/m^2/s', &
            avgflag='A', long_name='total ecosystem respiration, autotrophic + heterotrophic', &
            ptr_col=ccf%er)

       call hist_addfld1d (fname='LITFIRE', units='gC/m^2/s', &
            avgflag='A', long_name='litter fire losses', &
            ptr_col=ccf%litfire, default='inactive')

       call hist_addfld1d (fname='SOMFIRE', units='gC/m^2/s', &
            avgflag='A', long_name='soil organic matter fire losses', &
            ptr_col=ccf%somfire, default='inactive')

       call hist_addfld1d (fname='TOTFIRE', units='gC/m^2/s', &
            avgflag='A', long_name='total ecosystem fire losses', &
            ptr_col=ccf%totfire, default='inactive')

       call hist_addfld1d (fname='NEP', units='gC/m^2/s', &
            avgflag='A', long_name='net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink', &
            ptr_col=ccf%nep)

       call hist_addfld1d (fname='NBP', units='gC/m^2/s', &
            avgflag='A', long_name='net biome production, includes fire, landuse, and harvest flux, positive for sink', &
            ptr_col=ccf%nbp)

       call hist_addfld1d (fname='NEE', units='gC/m^2/s', &
            avgflag='A', long_name=&
            'net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source', &
            ptr_col=ccf%nee)

       call hist_addfld1d (fname='COL_FIRE_CLOSS', units='gC/m^2/s', &
            avgflag='A', long_name='total column-level fire C loss', &
            ptr_col=ccf%col_fire_closs)

       call hist_addfld1d (fname='DWT_SEEDC_TO_LEAF', units='gC/m^2/s', &
            avgflag='A', long_name='seed source to PFT-level leaf', &
            ptr_col=ccf%dwt_seedc_to_leaf)

       call hist_addfld1d (fname='DWT_SEEDC_TO_DEADSTEM', units='gC/m^2/s', &
            avgflag='A', long_name='seed source to PFT-level deadstem', &
            ptr_col=ccf%dwt_seedc_to_deadstem)

       call hist_addfld1d (fname='DWT_CONV_CFLUX', units='gC/m^2/s', &
            avgflag='A', long_name='conversion C flux (immediate loss to atm)', &
            ptr_col=ccf%dwt_conv_cflux)

       call hist_addfld1d (fname='DWT_PROD10C_GAIN', units='gC/m^2/s', &
            avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
            ptr_col=ccf%dwt_prod10c_gain)

       call hist_addfld1d (fname='PROD10C_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='loss from 10-yr wood product pool', &
            ptr_col=ccf%prod10c_loss)

       call hist_addfld1d (fname='DWT_PROD100C_GAIN', units='gC/m^2/s', &
            avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
            ptr_col=ccf%dwt_prod100c_gain)

       call hist_addfld1d (fname='PROD100C_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='loss from 100-yr wood product pool', &
            ptr_col=ccf%prod100c_loss)

       call hist_addfld1d (fname='DWT_FROOTC_TO_LITR1C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root to litter due to landcover change', &
            ptr_col=ccf%dwt_frootc_to_litr1c, default='inactive')

       call hist_addfld1d (fname='DWT_FROOTC_TO_LITR2C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root to litter due to landcover change', &
            ptr_col=ccf%dwt_frootc_to_litr2c, default='inactive')

       call hist_addfld1d (fname='DWT_FROOTC_TO_LITR3C', units='gC/m^2/s', &
            avgflag='A', long_name='fine root to litter due to landcover change', &
            ptr_col=ccf%dwt_frootc_to_litr3c, default='inactive')

       call hist_addfld1d (fname='DWT_LIVECROOTC_TO_CWDC', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root to CWD due to landcover change', &
            ptr_col=ccf%dwt_livecrootc_to_cwdc, default='inactive')

       call hist_addfld1d (fname='DWT_DEADCROOTC_TO_CWDC', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
            ptr_col=ccf%dwt_deadcrootc_to_cwdc, default='inactive')

       call hist_addfld1d (fname='DWT_CLOSS', units='gC/m^2/s', &
            avgflag='A', long_name='total carbon loss from land cover conversion', &
            ptr_col=ccf%dwt_closs)

       call hist_addfld1d (fname='PRODUCT_CLOSS', units='gC/m^2/s', &
            avgflag='A', long_name='total carbon loss from wood product pools', &
            ptr_col=ccf%product_closs)

       call hist_addfld1d (fname='LAND_USE_FLUX', units='gC/m^2/s', &
            avgflag='A', long_name='total C emitted from land cover conversion and wood product pools', &
            ptr_col=ccf%landuseflux)

       call hist_addfld1d (fname='LAND_UPTAKE', units='gC/m^2/s', &
            avgflag='A', long_name='NEE minus LAND_USE_FLUX, negative for update', &
            ptr_col=ccf%landuptake)

       if (use_c13) then
          !-------------------------------
          ! C13 flux variables - native to column 
          !-------------------------------

          call hist_addfld1d (fname='C13_M_LEAFC_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C mortality to litter 1 C', &
               ptr_col=cc13f%m_leafc_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_TO_LITR2C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C mortality to litter 2 C', &
               ptr_col=cc13f%m_leafc_to_litr2c, default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_TO_LITR3C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C mortality to litter 3 C', &
               ptr_col=cc13f%m_leafc_to_litr3c, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C mortality to litter 1 C', &
               ptr_col=cc13f%m_frootc_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_TO_LITR2C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C mortality to litter 2 C', &
               ptr_col=cc13f%m_frootc_to_litr2c, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_TO_LITR3C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C mortality to litter 3 C', &
               ptr_col=cc13f%m_frootc_to_litr3c, default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C storage mortality to litter 1 C', &
               ptr_col=cc13f%m_leafc_storage_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C storage mortality to litter 1 C', &
               ptr_col=cc13f%m_frootc_storage_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C storage mortality to litter 1 C', &
               ptr_col=cc13f%m_livestemc_storage_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C storage mortality to litter 1 C', &
               ptr_col=cc13f%m_deadstemc_storage_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C storage mortality to litter 1 C', &
               ptr_col=cc13f%m_livecrootc_storage_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C storage mortality to litter 1 C', &
               ptr_col=cc13f%m_deadcrootc_storage_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C transfer mortality to litter 1 C', &
               ptr_col=cc13f%m_leafc_xfer_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C transfer mortality to litter 1 C', &
               ptr_col=cc13f%m_frootc_xfer_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C transfer mortality to litter 1 C', &
               ptr_col=cc13f%m_livestemc_xfer_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C transfer mortality to litter 1 C', &
               ptr_col=cc13f%m_deadstemc_xfer_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C transfer mortality to litter 1 C', &
               ptr_col=cc13f%m_livecrootc_xfer_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C transfer mortality to litter 1 C', &
               ptr_col=cc13f%m_deadcrootc_xfer_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_CWDC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live stem C mortality to coarse woody debris C', &
               ptr_col=cc13f%m_livestemc_to_cwdc, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_CWDC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C mortality to coarse woody debris C', &
               ptr_col=cc13f%m_deadstemc_to_cwdc, default='inactive')

          call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_CWDC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root C mortality to coarse woody debris C', &
               ptr_col=cc13f%m_livecrootc_to_cwdc, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_CWDC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C mortality to coarse woody debris C', &
               ptr_col=cc13f%m_deadcrootc_to_cwdc, default='inactive')

          call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth respiration storage mortality to litter 1 C', &
               ptr_col=cc13f%m_gresp_storage_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 growth respiration transfer mortality to litter 1 C', &
               ptr_col=cc13f%m_gresp_xfer_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_CWDC_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead stem C to coarse woody debris C by fire', &
               ptr_col=cc13f%m_deadstemc_to_cwdc_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_CWDC_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root C to to woody debris C by fire', &
               ptr_col=cc13f%m_deadcrootc_to_cwdc_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LITR1C_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 litter 1 C fire loss', &
               ptr_col=cc13f%m_litr1c_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LITR2C_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 litter 2 C fire loss', &
               ptr_col=cc13f%m_litr2c_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_LITR3C_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 litter 3 C fire loss', &
               ptr_col=cc13f%m_litr3c_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_M_CWDC_TO_FIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 coarse woody debris C fire loss', &
               ptr_col=cc13f%m_cwdc_to_fire, default='inactive')

          call hist_addfld1d (fname='C13_LEAFC_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C litterfall to litter 1 C', &
               ptr_col=cc13f%leafc_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_LEAFC_TO_LITR2C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C litterfall to litter 2 C', &
               ptr_col=cc13f%leafc_to_litr2c, default='inactive')

          call hist_addfld1d (fname='C13_LEAFC_TO_LITR3C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 leaf C litterfall to litter 3 C', &
               ptr_col=cc13f%leafc_to_litr3c, default='inactive')

          call hist_addfld1d (fname='C13_FROOTC_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C litterfall to litter 1 C', &
               ptr_col=cc13f%frootc_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_FROOTC_TO_LITR2C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C litterfall to litter 2 C', &
               ptr_col=cc13f%frootc_to_litr2c, default='inactive')

          call hist_addfld1d (fname='C13_FROOTC_TO_LITR3C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root C litterfall to litter 3 C', &
               ptr_col=cc13f%frootc_to_litr3c, default='inactive')

          call hist_addfld1d (fname='C13_CWDC_TO_LITR2C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 decomp. of coarse woody debris C to litter 2 C', &
               ptr_col=cc13f%cwdc_to_litr2c, default='inactive')

          call hist_addfld1d (fname='C13_CWDC_TO_LITR3C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 decomp. of coarse woody debris C to litter 3 C', &
               ptr_col=cc13f%cwdc_to_litr3c, default='inactive')

          call hist_addfld1d (fname='C13_LITR1_HR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 het. resp. from litter 1 C', &
               ptr_col=cc13f%litr1_hr, default='inactive')

          call hist_addfld1d (fname='C13_LITR1C_TO_SOIL1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 decomp. of litter 1 C to SOM 1 C', &
               ptr_col=cc13f%litr1c_to_soil1c, default='inactive')

          call hist_addfld1d (fname='C13_LITR2_HR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 het. resp. from litter 2 C', &
               ptr_col=cc13f%litr2_hr, default='inactive')

          call hist_addfld1d (fname='C13_LITR2C_TO_SOIL2C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 decomp. of litter 2 C to SOM 2 C', &
               ptr_col=cc13f%litr2c_to_soil2c, default='inactive')

          call hist_addfld1d (fname='C13_LITR3_HR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 het. resp. from litter 3 C', &
               ptr_col=cc13f%litr3_hr, default='inactive')

          call hist_addfld1d (fname='C13_LITR3C_TO_SOIL3C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 decomp. of litter 3 C to SOM 3 C', &
               ptr_col=cc13f%litr3c_to_soil3c, default='inactive')

          call hist_addfld1d (fname='C13_SOIL1_HR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 het. resp. from SOM 1 C', &
               ptr_col=cc13f%soil1_hr, default='inactive')

          call hist_addfld1d (fname='C13_SOIL1C_TO_SOIL2C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 decomp. of SOM 1 C to SOM 2 C', &
               ptr_col=cc13f%soil1c_to_soil2c, default='inactive')

          call hist_addfld1d (fname='C13_SOIL2_HR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 het. resp. from SOM 2 C', &
               ptr_col=cc13f%soil2_hr, default='inactive')

          call hist_addfld1d (fname='C13_SOIL2C_TO_SOIL3C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 decomp. of SOM 2 C to SOM 3 C', &
               ptr_col=cc13f%soil2c_to_soil3c, default='inactive')

          call hist_addfld1d (fname='C13_SOIL3_HR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 het. resp. from SOM 3 C', &
               ptr_col=cc13f%soil3_hr, default='inactive')

          call hist_addfld1d (fname='C13_SOIL3C_TO_SOIL4C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 decomp. of SOM 3 C to SOM 4 C', &
               ptr_col=cc13f%soil3c_to_soil4c, default='inactive')

          call hist_addfld1d (fname='C13_SOIL4_HR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 het. resp. from SOM 4 C', &
               ptr_col=cc13f%soil4_hr, default='inactive')

          call hist_addfld1d (fname='C13_LITHR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 litter heterotrophic respiration', &
               ptr_col=cc13f%lithr)

          call hist_addfld1d (fname='C13_SOMHR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 soil organic matter heterotrophic respiration', &
               ptr_col=cc13f%somhr)

          call hist_addfld1d (fname='C13_HR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total heterotrophic respiration', &
               ptr_col=cc13f%hr)

          call hist_addfld1d (fname='C13_SR', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total soil respiration (HR + root resp)', &
               ptr_col=cc13f%sr)

          call hist_addfld1d (fname='C13_ER', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total ecosystem respiration, autotrophic + heterotrophic', &
               ptr_col=cc13f%er)

          call hist_addfld1d (fname='C13_LITFIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 litter fire losses', &
               ptr_col=cc13f%litfire, default='inactive')

          call hist_addfld1d (fname='C13_SOMFIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 soil organic matter fire losses', &
               ptr_col=cc13f%somfire, default='inactive')

          call hist_addfld1d (fname='C13_TOTFIRE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total ecosystem fire losses', &
               ptr_col=cc13f%totfire, default='inactive')

          call hist_addfld1d (fname='C13_NEP', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 net ecosystem production, excludes fire flux, positive for sink', &
               ptr_col=cc13f%nep)

          call hist_addfld1d (fname='C13_NEE', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 net ecosystem exchange of carbon, includes fire flux, positive for source', &
               ptr_col=cc13f%nee)

          call hist_addfld1d (fname='C13_COL_FIRE_CLOSS', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total column-level fire C loss', &
               ptr_col=cc13f%col_fire_closs)

          call hist_addfld1d (fname='C13_DWT_SEEDC_TO_LEAF', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 seed source to PFT-level leaf', &
               ptr_col=cc13f%dwt_seedc_to_leaf)

          call hist_addfld1d (fname='C13_DWT_SEEDC_TO_DEADSTEM', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 seed source to PFT-level deadstem', &
               ptr_col=cc13f%dwt_seedc_to_deadstem)

          call hist_addfld1d (fname='C13_DWT_CONV_CFLUX', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 conversion C flux (immediate loss to atm)', &
               ptr_col=cc13f%dwt_conv_cflux)

          call hist_addfld1d (fname='C13_DWT_PROD10C_GAIN', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 addition to 10-yr wood product pool', &
               ptr_col=cc13f%dwt_prod10c_gain)

          call hist_addfld1d (fname='C13_PROD10C_LOSS', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 loss from 10-yr wood product pool', &
               ptr_col=cc13f%prod10c_loss)

          call hist_addfld1d (fname='C13_DWT_PROD100C_GAIN', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 addition to 100-yr wood product pool', &
               ptr_col=cc13f%dwt_prod100c_gain)

          call hist_addfld1d (fname='C13_PROD100C_LOSS', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 loss from 100-yr wood product pool', &
               ptr_col=cc13f%prod100c_loss)

          call hist_addfld1d (fname='C13_DWT_FROOTC_TO_LITR1C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root to litter due to landcover change', &
               ptr_col=cc13f%dwt_frootc_to_litr1c, default='inactive')

          call hist_addfld1d (fname='C13_DWT_FROOTC_TO_LITR2C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root to litter due to landcover change', &
               ptr_col=cc13f%dwt_frootc_to_litr2c, default='inactive')

          call hist_addfld1d (fname='C13_DWT_FROOTC_TO_LITR3C', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 fine root to litter due to landcover change', &
               ptr_col=cc13f%dwt_frootc_to_litr3c, default='inactive')

          call hist_addfld1d (fname='C13_DWT_LIVECROOTC_TO_CWDC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 live coarse root to CWD due to landcover change', &
               ptr_col=cc13f%dwt_livecrootc_to_cwdc, default='inactive')

          call hist_addfld1d (fname='C13_DWT_DEADCROOTC_TO_CWDC', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 dead coarse root to CWD due to landcover change', &
               ptr_col=cc13f%dwt_deadcrootc_to_cwdc, default='inactive')

          call hist_addfld1d (fname='C13_DWT_CLOSS', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total carbon loss from land cover conversion', &
               ptr_col=cc13f%dwt_closs)

          call hist_addfld1d (fname='C13_PRODUCT_CLOSS', units='gC13/m^2/s', &
               avgflag='A', long_name='C13 total carbon loss from wood product pools', &
               ptr_col=cc13f%product_closs)
       endif

       !-------------------------------
       ! N flux variables - native to PFT
       !-------------------------------

       call hist_addfld1d (fname='M_LEAFN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N mortality', &
            ptr_pft=pnf%m_leafn_to_litter, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N mortality', &
            ptr_pft=pnf%m_frootn_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N storage mortality', &
            ptr_pft=pnf%m_leafn_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N storage mortality', &
            ptr_pft=pnf%m_frootn_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N storage mortality', &
            ptr_pft=pnf%m_livestemn_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N storage mortality', &
            ptr_pft=pnf%m_deadstemn_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N storage mortality', &
            ptr_pft=pnf%m_livecrootn_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N storage mortality', &
            ptr_pft=pnf%m_deadcrootn_storage_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_XFER_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N transfer mortality', &
            ptr_pft=pnf%m_leafn_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N transfer mortality', &
            ptr_pft=pnf%m_frootn_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N transfer mortality', &
            ptr_pft=pnf%m_livestemn_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N transfer mortality', &
            ptr_pft=pnf%m_deadstemn_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N transfer mortality', &
            ptr_pft=pnf%m_livecrootn_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N transfer mortality', &
            ptr_pft=pnf%m_deadcrootn_xfer_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N mortality', &
            ptr_pft=pnf%m_livestemn_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N mortality', &
            ptr_pft=pnf%m_deadstemn_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N mortality', &
            ptr_pft=pnf%m_livecrootn_to_litter, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N mortality', &
            ptr_pft=pnf%m_deadcrootn_to_litter, default='inactive')

       call hist_addfld1d (fname='M_RETRANSN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='retranslocated N pool mortality', &
            ptr_pft=pnf%m_retransn_to_litter, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N fire loss', &
            ptr_pft=pnf%m_leafn_to_fire, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N fire loss ', &
            ptr_pft=pnf%m_frootn_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N storage fire loss', &
            ptr_pft=pnf%m_leafn_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N storage fire loss', &
            ptr_pft=pnf%m_frootn_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N storage fire loss', &
            ptr_pft=pnf%m_livestemn_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N storage fire loss', &
            ptr_pft=pnf%m_deadstemn_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N storage fire loss', &
            ptr_pft=pnf%m_livecrootn_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N storage fire loss', &
            ptr_pft=pnf%m_deadcrootn_storage_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_XFER_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N transfer fire loss', &
            ptr_pft=pnf%m_leafn_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N transfer fire loss', &
            ptr_pft=pnf%m_frootn_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N transfer fire loss', &
            ptr_pft=pnf%m_livestemn_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N transfer fire loss', &
            ptr_pft=pnf%m_deadstemn_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N transfer fire loss', &
            ptr_pft=pnf%m_livecrootn_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N transfer fire loss', &
            ptr_pft=pnf%m_deadcrootn_xfer_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N fire loss', &
            ptr_pft=pnf%m_livestemn_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N fire loss', &
            ptr_pft=pnf%m_deadstemn_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N fire mortality to litter', &
            ptr_pft=pnf%m_deadstemn_to_litter_fire, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N fire loss', &
            ptr_pft=pnf%m_livecrootn_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N fire loss', &
            ptr_pft=pnf%m_deadcrootn_to_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N fire mortality to litter', &
            ptr_pft=pnf%m_deadcrootn_to_litter_fire, default='inactive')

       call hist_addfld1d (fname='M_RETRANSN_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='retranslocated N pool fire loss', &
            ptr_pft=pnf%m_retransn_to_fire, default='inactive')

       call hist_addfld1d (fname='LEAFN_XFER_TO_LEAFN', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N growth from storage', &
            ptr_pft=pnf%leafn_xfer_to_leafn, default='inactive')

       call hist_addfld1d (fname='FROOTN_XFER_TO_FROOTN', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N growth from storage', &
            ptr_pft=pnf%frootn_xfer_to_frootn, default='inactive')

       call hist_addfld1d (fname='LIVESTEMN_XFER_TO_LIVESTEMN', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N growth from storage', &
            ptr_pft=pnf%livestemn_xfer_to_livestemn, default='inactive')

       call hist_addfld1d (fname='DEADSTEMN_XFER_TO_DEADSTEMN', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N growth from storage', &
            ptr_pft=pnf%deadstemn_xfer_to_deadstemn, default='inactive')

       call hist_addfld1d (fname='LIVECROOTN_XFER_TO_LIVECROOTN', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N growth from storage', &
            ptr_pft=pnf%livecrootn_xfer_to_livecrootn, default='inactive')

       call hist_addfld1d (fname='DEADCROOTN_XFER_TO_DEADCROOTN', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N growth from storage', &
            ptr_pft=pnf%deadcrootn_xfer_to_deadcrootn, default='inactive')

       call hist_addfld1d (fname='LEAFN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N litterfall', &
            ptr_pft=pnf%leafn_to_litter, default='inactive')

       call hist_addfld1d (fname='LEAFN_TO_RETRANSN', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N to retranslocated N pool', &
            ptr_pft=pnf%leafn_to_retransn, default='inactive')

       call hist_addfld1d (fname='FROOTN_TO_LITTER', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N litterfall', &
            ptr_pft=pnf%frootn_to_litter, default='inactive')

       call hist_addfld1d (fname='RETRANSN_TO_NPOOL', units='gN/m^2/s', &
            avgflag='A', long_name='deployment of retranslocated N', &
            ptr_pft=pnf%retransn_to_npool)

       call hist_addfld1d (fname='SMINN_TO_NPOOL', units='gN/m^2/s', &
            avgflag='A', long_name='deployment of soil mineral N uptake', &
            ptr_pft=pnf%sminn_to_npool)

       call hist_addfld1d (fname='NPOOL_TO_LEAFN', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to leaf N', &
            ptr_pft=pnf%npool_to_leafn, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_LEAFN_STORAGE', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to leaf N storage', &
            ptr_pft=pnf%npool_to_leafn_storage, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_FROOTN', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to fine root N', &
            ptr_pft=pnf%npool_to_frootn, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_FROOTN_STORAGE', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to fine root N storage', &
            ptr_pft=pnf%npool_to_frootn_storage, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to live stem N', &
            ptr_pft=pnf%npool_to_livestemn, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN_STORAGE', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to live stem N storage', &
            ptr_pft=pnf%npool_to_livestemn_storage, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to dead stem N', &
            ptr_pft=pnf%npool_to_deadstemn, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN_STORAGE', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to dead stem N storage', &
            ptr_pft=pnf%npool_to_deadstemn_storage, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to live coarse root N', &
            ptr_pft=pnf%npool_to_livecrootn, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN_STORAGE', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to live coarse root N storage', &
            ptr_pft=pnf%npool_to_livecrootn_storage, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to dead coarse root N', &
            ptr_pft=pnf%npool_to_deadcrootn, default='inactive')

       call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN_STORAGE', units='gN/m^2/s', &
            avgflag='A', long_name='allocation to dead coarse root N storage', &
            ptr_pft=pnf%npool_to_deadcrootn_storage, default='inactive')

       call hist_addfld1d (fname='LEAFN_STORAGE_TO_XFER', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N shift storage to transfer', &
            ptr_pft=pnf%leafn_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='FROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N shift storage to transfer', &
            ptr_pft=pnf%frootn_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='LIVESTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N shift storage to transfer', &
            ptr_pft=pnf%livestemn_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='DEADSTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N shift storage to transfer', &
            ptr_pft=pnf%deadstemn_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='LIVECROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N shift storage to transfer', &
            ptr_pft=pnf%livecrootn_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='DEADCROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N shift storage to transfer', &
            ptr_pft=pnf%deadcrootn_storage_to_xfer, default='inactive')

       call hist_addfld1d (fname='LIVESTEMN_TO_DEADSTEMN', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N turnover', &
            ptr_pft=pnf%livestemn_to_deadstemn, default='inactive')

       call hist_addfld1d (fname='LIVESTEMN_TO_RETRANSN', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N to retranslocated N pool', &
            ptr_pft=pnf%livestemn_to_retransn, default='inactive')

       call hist_addfld1d (fname='LIVECROOTN_TO_DEADCROOTN', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N turnover', &
            ptr_pft=pnf%livecrootn_to_deadcrootn, default='inactive')

       call hist_addfld1d (fname='LIVECROOTN_TO_RETRANSN', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N to retranslocated N pool', &
            ptr_pft=pnf%livecrootn_to_retransn, default='inactive')

       call hist_addfld1d (fname='NDEPLOY', units='gN/m^2/s', &
            avgflag='A', long_name='total N deployed in new growth', &
            ptr_pft=pnf%ndeploy)

       call hist_addfld1d (fname='WOOD_HARVESTN', units='gN/m^2/s', &
            avgflag='A', long_name='wood harvest N (to product pools)', &
            ptr_pft=pnf%wood_harvestn)

       call hist_addfld1d (fname='PFT_FIRE_NLOSS', units='gN/m^2/s', &
            avgflag='A', long_name='total pft-level fire N loss', &
            ptr_pft=pnf%pft_fire_nloss)

       !-------------------------------
       ! N flux variables - native to column
       !-------------------------------

       call hist_addfld1d (fname='NDEP_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='atmospheric N deposition to soil mineral N', &
            ptr_col=cnf%ndep_to_sminn)

       call hist_addfld1d (fname='NFIX_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='symbiotic/asymbiotic N fixation to soil mineral N', &
            ptr_col=cnf%nfix_to_sminn)

       call hist_addfld1d (fname='M_LEAFN_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N mortality to litter 1 N', &
            ptr_col=cnf%m_leafn_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_TO_LITR2N', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N mortality to litter 2 N', &
            ptr_col=cnf%m_leafn_to_litr2n, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_TO_LITR3N', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N mortality to litter 3 N', &
            ptr_col=cnf%m_leafn_to_litr3n, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N mortality to litter 1 N', &
            ptr_col=cnf%m_frootn_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_TO_LITR2N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N mortality to litter 2 N', &
            ptr_col=cnf%m_frootn_to_litr2n, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_TO_LITR3N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N mortality to litter 3 N', &
            ptr_col=cnf%m_frootn_to_litr3n, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N storage mortality to litter 1 N', &
            ptr_col=cnf%m_leafn_storage_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N storage mortality to litter 1 N', &
            ptr_col=cnf%m_frootn_storage_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N storage mortality to litter 1 N', &
            ptr_col=cnf%m_livestemn_storage_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N storage mortality to litter 1 N', &
            ptr_col=cnf%m_deadstemn_storage_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N storage mortality to litter 1 N', &
            ptr_col=cnf%m_livecrootn_storage_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N storage mortality to litter 1 N', &
            ptr_col=cnf%m_deadcrootn_storage_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_LEAFN_XFER_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N transfer mortality to litter 1 N', &
            ptr_col=cnf%m_leafn_xfer_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_FROOTN_XFER_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N transfer mortality to litter 1 N', &
            ptr_col=cnf%m_frootn_xfer_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N transfer mortality to litter 1 N', &
            ptr_col=cnf%m_livestemn_xfer_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N transfer mortality to litter 1 N', &
            ptr_col=cnf%m_deadstemn_xfer_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N transfer mortality to litter 1 N', &
            ptr_col=cnf%m_livecrootn_xfer_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N transfer mortality to litter 1 N', &
            ptr_col=cnf%m_deadcrootn_xfer_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_LIVESTEMN_TO_CWDN', units='gN/m^2/s', &
            avgflag='A', long_name='live stem N mortality to coarse woody debris N', &
            ptr_col=cnf%m_livestemn_to_cwdn, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_TO_CWDN', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N mortality to coarse woody debris N', &
            ptr_col=cnf%m_deadstemn_to_cwdn, default='inactive')

       call hist_addfld1d (fname='M_LIVECROOTN_TO_CWDN', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root N mortality to coarse woody debris N', &
            ptr_col=cnf%m_livecrootn_to_cwdn, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_TO_CWDN', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N mortality to coarse woody debris N', &
            ptr_col=cnf%m_deadcrootn_to_cwdn, default='inactive')

       call hist_addfld1d (fname='M_RETRANSN_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='retranslocated N pool mortality to litter 1 N', &
            ptr_col=cnf%m_retransn_to_litr1n, default='inactive')

       call hist_addfld1d (fname='M_DEADSTEMN_TO_CWDN_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead stem N to coarse woody debris N by fire', &
            ptr_col=cnf%m_deadstemn_to_cwdn_fire, default='inactive')

       call hist_addfld1d (fname='M_DEADCROOTN_TO_CWDN_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root N to to woody debris N by fire', &
            ptr_col=cnf%m_deadcrootn_to_cwdn_fire, default='inactive')

       call hist_addfld1d (fname='M_LITR1N_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='litter 1 N fire loss', &
            ptr_col=cnf%m_litr1n_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LITR2N_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='litter 2 N fire loss', &
            ptr_col=cnf%m_litr2n_to_fire, default='inactive')

       call hist_addfld1d (fname='M_LITR3N_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='litter 3 N fire loss', &
            ptr_col=cnf%m_litr3n_to_fire, default='inactive')

       call hist_addfld1d (fname='M_CWDN_TO_FIRE', units='gN/m^2/s', &
            avgflag='A', long_name='coarse woody debris N fire loss', &
            ptr_col=cnf%m_cwdn_to_fire, default='inactive')

       call hist_addfld1d (fname='LEAFN_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N litterfall to litter 1 N', &
            ptr_col=cnf%leafn_to_litr1n, default='inactive')

       call hist_addfld1d (fname='LEAFN_TO_LITR2N', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N litterfall to litter 2 N', &
            ptr_col=cnf%leafn_to_litr2n, default='inactive')

       call hist_addfld1d (fname='LEAFN_TO_LITR3N', units='gN/m^2/s', &
            avgflag='A', long_name='leaf N litterfall to litter 3 N', &
            ptr_col=cnf%leafn_to_litr3n, default='inactive')

       call hist_addfld1d (fname='FROOTN_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N litterfall to litter 1 N', &
            ptr_col=cnf%frootn_to_litr1n, default='inactive')

       call hist_addfld1d (fname='FROOTN_TO_LITR2N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N litterfall to litter 2 N', &
            ptr_col=cnf%frootn_to_litr2n, default='inactive')

       call hist_addfld1d (fname='FROOTN_TO_LITR3N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root N litterfall to litter 3 N ', &
            ptr_col=cnf%frootn_to_litr3n, default='inactive')

       call hist_addfld1d (fname='CWDN_TO_LITR2N', units='gN/m^2/s', &
            avgflag='A', long_name='decomp. of coarse woody debris N to litter 2 N', &
            ptr_col=cnf%cwdn_to_litr2n, default='inactive')

       call hist_addfld1d (fname='CWDN_TO_LITR3N', units='gN/m^2/s', &
            avgflag='A', long_name='decomp. of coarse woody debris N to litter 3 N', &
            ptr_col=cnf%cwdn_to_litr3n, default='inactive')

       call hist_addfld1d (fname='LITR1N_TO_SOIL1N', units='gN/m^2/s', &
            avgflag='A', long_name='decomp. of litter 1 N to SOM 1 N', &
            ptr_col=cnf%litr1n_to_soil1n, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_SOIL1N_L1', units='gN/m^2/s', &
            avgflag='A', long_name='mineral N flux for decomp. of litter 1 to SOM 1', &
            ptr_col=cnf%sminn_to_soil1n_l1, default='inactive')

       call hist_addfld1d (fname='LITR2N_TO_SOIL2N', units='gN/m^2/s', &
            avgflag='A', long_name='decomp. of litter 2 N to SOM 2 N', &
            ptr_col=cnf%litr2n_to_soil2n, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_SOIL2N_L2', units='gN/m^2/s', &
            avgflag='A', long_name='mineral N flux for decomp. of litter 2 to SOM 2', &
            ptr_col=cnf%sminn_to_soil2n_l2, default='inactive')

       call hist_addfld1d (fname='LITR3N_TO_SOIL3N', units='gN/m^2/s', &
            avgflag='A', long_name='decomp. of litter 3 N to SOM 3 N', &
            ptr_col=cnf%litr3n_to_soil3n, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_SOIL3N_L3', units='gN/m^2/s', &
            avgflag='A', long_name='mineral N flux for decomp. of litter 3 to SOM 3', &
            ptr_col=cnf%sminn_to_soil3n_l3, default='inactive')

       call hist_addfld1d (fname='SOIL1N_TO_SOIL2n', units='gN/m^2/s', &
            avgflag='A', long_name='decomp. of SOM 1 N to SOM 2 N', &
            ptr_col=cnf%soil1n_to_soil2n, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_SOIL2N_S1', units='gN/m^2/s', &
            avgflag='A', long_name='mineral N flux for decomp. of SOM 1 to SOM 2', &
            ptr_col=cnf%sminn_to_soil2n_s1, default='inactive')

       call hist_addfld1d (fname='SOIL2N_TO_SOIL3N', units='gN/m^2/s', &
            avgflag='A', long_name='decomp. of SOM 2 N to SOM 3 N', &
            ptr_col=cnf%soil2n_to_soil3n, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_SOIL3N_S2', units='gN/m^2/s', &
            avgflag='A', long_name='mineral N flux for decomp. of SOM 2 to SOM 3', &
            ptr_col=cnf%sminn_to_soil3n_s2, default='inactive')

       call hist_addfld1d (fname='SOIL3N_TO_SOIL4N', units='gN/m^2/s', &
            avgflag='A', long_name='decomp. of SOM 3 N to SOM 4 N', &
            ptr_col=cnf%soil3n_to_soil4n, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_SOIL4N_S3', units='gN/m^2/s', &
            avgflag='A', long_name='mineral N flux for decomp. of SOM 3 to SOM 4', &
            ptr_col=cnf%sminn_to_soil4n_s3, default='inactive')

       call hist_addfld1d (fname='SOIL4N_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='N mineralization for decomp. of SOM 4', &
            ptr_col=cnf%soil4n_to_sminn, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_DENIT_L1S1', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification for decomp. of litter 1 to SOM 1', &
            ptr_col=cnf%sminn_to_denit_l1s1, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_DENIT_L2S2', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification for decomp. of litter 2 to SOM 2', &
            ptr_col=cnf%sminn_to_denit_l2s2, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_DENIT_L3S3', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification for decomp. of litter 3 to SOM 3', &
            ptr_col=cnf%sminn_to_denit_l3s3, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_DENIT_S1S2', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification for decomp. of SOM 1 to SOM 2', &
            ptr_col=cnf%sminn_to_denit_s1s2, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_DENIT_S2S3', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification for decomp. of SOM 2 to SOM 3', &
            ptr_col=cnf%sminn_to_denit_s2s3, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_DENIT_S3S4', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification for decomp. of SOM 3 to SOM 4', &
            ptr_col=cnf%sminn_to_denit_s3s4, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_DENIT_S4', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification for decomp. of SOM 4', &
            ptr_col=cnf%sminn_to_denit_s4, default='inactive')

       call hist_addfld1d (fname='SMINN_TO_DENIT_EXCESS', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification from excess mineral N pool', &
            ptr_col=cnf%sminn_to_denit_excess, default='inactive')

       call hist_addfld1d (fname='SMINN_LEACHED', units='gN/m^2/s', &
            avgflag='A', long_name='soil mineral N pool loss to leaching', &
            ptr_col=cnf%sminn_leached)

       call hist_addfld1d (fname='POTENTIAL_IMMOB', units='gN/m^2/s', &
            avgflag='A', long_name='potential N immobilization', &
            ptr_col=cnf%potential_immob)

       call hist_addfld1d (fname='ACTUAL_IMMOB', units='gN/m^2/s', &
            avgflag='A', long_name='actual N immobilization', &
            ptr_col=cnf%actual_immob)

       call hist_addfld1d (fname='SMINN_TO_PLANT', units='gN/m^2/s', &
            avgflag='A', long_name='plant uptake of soil mineral N', &
            ptr_col=cnf%sminn_to_plant)

       call hist_addfld1d (fname='SUPPLEMENT_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='supplemental N supply', &
            ptr_col=cnf%supplement_to_sminn)

       call hist_addfld1d (fname='GROSS_NMIN', units='gN/m^2/s', &
            avgflag='A', long_name='gross rate of N mineralization', &
            ptr_col=cnf%gross_nmin)

       call hist_addfld1d (fname='NET_NMIN', units='gN/m^2/s', &
            avgflag='A', long_name='net rate of N mineralization', &
            ptr_col=cnf%net_nmin)

       call hist_addfld1d (fname='DENIT', units='gN/m^2/s', &
            avgflag='A', long_name='total rate of denitrification', &
            ptr_col=cnf%denit)

       call hist_addfld1d (fname='COL_FIRE_NLOSS', units='gN/m^2/s', &
            avgflag='A', long_name='total column-level fire N loss', &
            ptr_col=cnf%col_fire_nloss)

       call hist_addfld1d (fname='DWT_SEEDN_TO_LEAF', units='gN/m^2/s', &
            avgflag='A', long_name='seed source to PFT-level leaf', &
            ptr_col=cnf%dwt_seedn_to_leaf)

       call hist_addfld1d (fname='DWT_SEEDN_TO_DEADSTEM', units='gN/m^2/s', &
            avgflag='A', long_name='seed source to PFT-level deadstem', &
            ptr_col=cnf%dwt_seedn_to_deadstem)

       call hist_addfld1d (fname='DWT_CONV_NFLUX', units='gN/m^2/s', &
            avgflag='A', long_name='conversion N flux (immediate loss to atm)', &
            ptr_col=cnf%dwt_conv_nflux)

       call hist_addfld1d (fname='DWT_PROD10N_GAIN', units='gN/m^2/s', &
            avgflag='A', long_name='addition to 10-yr wood product pool', &
            ptr_col=cnf%dwt_prod10n_gain)

       call hist_addfld1d (fname='PROD10N_LOSS', units='gN/m^2/s', &
            avgflag='A', long_name='loss from 10-yr wood product pool', &
            ptr_col=cnf%prod10n_loss)

       call hist_addfld1d (fname='DWT_PROD100N_GAIN', units='gN/m^2/s', &
            avgflag='A', long_name='addition to 100-yr wood product pool', &
            ptr_col=cnf%dwt_prod100n_gain)

       call hist_addfld1d (fname='PROD100N_LOSS', units='gN/m^2/s', &
            avgflag='A', long_name='loss from 100-yr wood product pool', &
            ptr_col=cnf%prod100n_loss)

       call hist_addfld1d (fname='PRODUCT_NLOSS', units='gN/m^2/s', &
            avgflag='A', long_name='total N loss from wood product pools', &
            ptr_col=cnf%product_nloss)

       call hist_addfld1d (fname='DWT_FROOTN_TO_LITR1N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root to litter due to landcover change', &
            ptr_col=cnf%dwt_frootn_to_litr1n, default='inactive')

       call hist_addfld1d (fname='DWT_FROOTN_TO_LITR2N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root to litter due to landcover change', &
            ptr_col=cnf%dwt_frootn_to_litr2n, default='inactive')

       call hist_addfld1d (fname='DWT_FROOTN_TO_LITR3N', units='gN/m^2/s', &
            avgflag='A', long_name='fine root to litter due to landcover change', &
            ptr_col=cnf%dwt_frootn_to_litr3n, default='inactive')

       call hist_addfld1d (fname='DWT_LIVECROOTN_TO_CWDN', units='gN/m^2/s', &
            avgflag='A', long_name='live coarse root to CWD due to landcover change', &
            ptr_col=cnf%dwt_livecrootn_to_cwdn, default='inactive')

       call hist_addfld1d (fname='DWT_DEADCROOTN_TO_CWDN', units='gN/m^2/s', &
            avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
            ptr_col=cnf%dwt_deadcrootn_to_cwdn, default='inactive')

       call hist_addfld1d (fname='DWT_NLOSS', units='gN/m^2/s', &
            avgflag='A', long_name='total nitrogen loss from landcover conversion', &
            ptr_col=cnf%dwt_nloss)

       !-------------------------------
       ! PFT ecophysiological variables (pepv) 
       !-------------------------------

       call hist_addfld1d (fname='DORMANT_FLAG', units='none', &
            avgflag='A', long_name='dormancy flag', &
            ptr_pft=pepv%dormant_flag, default='inactive')

       call hist_addfld1d (fname='DAYS_ACTIVE', units='days', &
            avgflag='A', long_name='number of days since last dormancy', &
            ptr_pft=pepv%days_active, default='inactive')

       call hist_addfld1d (fname='ONSET_FLAG', units='none', &
            avgflag='A', long_name='onset flag', &
            ptr_pft=pepv%onset_flag, default='inactive')

       call hist_addfld1d (fname='ONSET_COUNTER', units='days', &
            avgflag='A', long_name='onset days counter', &
            ptr_pft=pepv%onset_counter, default='inactive')

       call hist_addfld1d (fname='ONSET_GDDFLAG', units='none', &
            avgflag='A', long_name='onset flag for growing degree day sum', &
            ptr_pft=pepv%onset_gddflag, default='inactive')

       call hist_addfld1d (fname='ONSET_FDD', units='C degree-days', &
            avgflag='A', long_name='onset freezing degree days counter', &
            ptr_pft=pepv%onset_fdd, default='inactive')

       call hist_addfld1d (fname='ONSET_GDD', units='C degree-days', &
            avgflag='A', long_name='onset growing degree days', &
            ptr_pft=pepv%onset_gdd, default='inactive')

       call hist_addfld1d (fname='ONSET_SWI', units='none', &
            avgflag='A', long_name='onset soil water index', &
            ptr_pft=pepv%onset_swi, default='inactive')

       call hist_addfld1d (fname='OFFSET_FLAG', units='none', &
            avgflag='A', long_name='offset flag', &
            ptr_pft=pepv%offset_flag, default='inactive')

       call hist_addfld1d (fname='OFFSET_COUNTER', units='days', &
            avgflag='A', long_name='offset days counter', &
            ptr_pft=pepv%offset_counter, default='inactive')

       call hist_addfld1d (fname='OFFSET_FDD', units='C degree-days', &
            avgflag='A', long_name='offset freezing degree days counter', &
            ptr_pft=pepv%offset_fdd, default='inactive')

       call hist_addfld1d (fname='OFFSET_SWI', units='none', &
            avgflag='A', long_name='offset soil water index', &
            ptr_pft=pepv%offset_swi, default='inactive')

       call hist_addfld1d (fname='LGSF', units='proportion', &
            avgflag='A', long_name='long growing season factor', &
            ptr_pft=pepv%lgsf, default='inactive')

       call hist_addfld1d (fname='BGLFR', units='1/s', &
            avgflag='A', long_name='background litterfall rate', &
            ptr_pft=pepv%bglfr, default='inactive')

       call hist_addfld1d (fname='BGTR', units='1/s', &
            avgflag='A', long_name='background transfer growth rate', &
            ptr_pft=pepv%bgtr, default='inactive')

       call hist_addfld1d (fname='DAYL',  units='s', &
            avgflag='A', long_name='daylength', &
            ptr_pft=pepv%dayl, default='inactive')

       call hist_addfld1d (fname='PREV_DAYL', units='s', &
            avgflag='A', long_name='daylength from previous timestep', &
            ptr_pft=pepv%prev_dayl, default='inactive')

       call hist_addfld1d (fname='ANNAVG_T2M', units='K', &
            avgflag='A', long_name='annual average 2m air temperature', &
            ptr_pft=pepv%annavg_t2m, default='inactive')

       call hist_addfld1d (fname='TEMPAVG_T2M', units='K', &
            avgflag='A', long_name='temporary average 2m air temperature', &
            ptr_pft=pepv%tempavg_t2m, default='inactive')

       call hist_addfld1d (fname='INIT_GPP', units='gC/m^2/s', &
            avgflag='A', long_name='GPP flux before downregulation', &
            ptr_pft=pepv%gpp, default='inactive')

       call hist_addfld1d (fname='AVAILC', units='gC/m^2/s', &
            avgflag='A', long_name='C flux available for allocation', &
            ptr_pft=pepv%availc, default='inactive')

       call hist_addfld1d (fname='XSMRPOOL_RECOVER', units='gC/m^2/s', &
            avgflag='A', long_name='C flux assigned to recovery of negative xsmrpool', &
            ptr_pft=pepv%xsmrpool_recover)

       if (use_c13) then
          call hist_addfld1d (fname='XSMRPOOL_C13RATIO', units='proportion', &
               avgflag='A', long_name='C13/C(12+13) ratio for xsmrpool', &
               ptr_pft=pepv%xsmrpool_c13ratio, default='inactive')
       endif

       call hist_addfld1d (fname='ALLOC_PNOW', units='proportion', &
            avgflag='A', long_name='fraction of current allocation to display as new growth', &
            ptr_pft=pepv%alloc_pnow, default='inactive')

       call hist_addfld1d (fname='C_ALLOMETRY', units='none', &
            avgflag='A', long_name='C allocation index', &
            ptr_pft=pepv%c_allometry, default='inactive')

       call hist_addfld1d (fname='N_ALLOMETRY', units='none', &
            avgflag='A', long_name='N allocation index', &
            ptr_pft=pepv%n_allometry, default='inactive')

       call hist_addfld1d (fname='PLANT_NDEMAND', units='gN/m^2/s', &
            avgflag='A', long_name='N flux required to support initial GPP', &
            ptr_pft=pepv%plant_ndemand)

       call hist_addfld1d (fname='TEMPSUM_POTENTIAL_GPP', units='gC/m^2/yr', &
            avgflag='A', long_name='temporary annual sum of potential GPP', &
            ptr_pft=pepv%tempsum_potential_gpp, default='inactive')

       call hist_addfld1d (fname='ANNSUM_POTENTIAL_GPP', units='gN/m^2/yr', &
            avgflag='A', long_name='annual sum of potential GPP', &
            ptr_pft=pepv%annsum_potential_gpp, default='inactive')

       call hist_addfld1d (fname='TEMPMAX_RETRANSN', units='gN/m^2', &
            avgflag='A', long_name='temporary annual max of retranslocated N pool', &
            ptr_pft=pepv%tempmax_retransn, default='inactive')

       call hist_addfld1d (fname='ANNMAX_RETRANSN', units='gN/m^2', &
            avgflag='A', long_name='annual max of retranslocated N pool', &
            ptr_pft=pepv%annmax_retransn, default='inactive')

       call hist_addfld1d (fname='AVAIL_RETRANSN', units='gN/m^2/s', &
            avgflag='A', long_name='N flux available from retranslocation pool', &
            ptr_pft=pepv%avail_retransn, default='inactive')

       call hist_addfld1d (fname='PLANT_NALLOC', units='gN/m^2/s', &
            avgflag='A', long_name='total allocated N flux', &
            ptr_pft=pepv%plant_nalloc, default='inactive')

       call hist_addfld1d (fname='PLANT_CALLOC', units='gC/m^2/s', &
            avgflag='A', long_name='total allocated C flux', &
            ptr_pft=pepv%plant_calloc, default='inactive')

       call hist_addfld1d (fname='EXCESS_CFLUX', units='gC/m^2/s', &
            avgflag='A', long_name='C flux not allocated due to downregulation', &
            ptr_pft=pepv%excess_cflux, default='inactive')

       call hist_addfld1d (fname='DOWNREG', units='proportion', &
            avgflag='A', long_name='fractional reduction in GPP due to N limitation', &
            ptr_pft=pepv%downreg, default='inactive')

       call hist_addfld1d (fname='PREV_LEAFC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='previous timestep leaf C litterfall flux', &
            ptr_pft=pepv%prev_leafc_to_litter, default='inactive')

       call hist_addfld1d (fname='PREV_FROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='previous timestep froot C litterfall flux', &
            ptr_pft=pepv%prev_frootc_to_litter, default='inactive')

       call hist_addfld1d (fname='ANNSUM_NPP', units='gC/m^2/yr', &
            avgflag='A', long_name='annual sum of NPP', &
            ptr_pft=pepv%annsum_npp, default='inactive')

       if (use_c13) then
          call hist_addfld1d (fname='RC13_CANAIR', units='proportion', &
               avgflag='A', long_name='C13/C(12+13) for canopy air', &
               ptr_pft=pepv%rc13_canair, default='inactive')

          call hist_addfld1d (fname='RC13_PSNSUN', units='proportion', &
               avgflag='A', long_name='C13/C(12+13) for sunlit photosynthesis', &
               ptr_pft=pepv%rc13_psnsun, default='inactive')

          call hist_addfld1d (fname='RC13_PSNSHA', units='proportion', &
               avgflag='A', long_name='C13/C(12+13) for shaded photosynthesis', &
               ptr_pft=pepv%rc13_psnsha, default='inactive')
       endif

       !-------------------------------
       ! PFT physical state variables not already defined by default
       !-------------------------------

       call hist_addfld1d (fname='EMV', units='proportion', &
            avgflag='A', long_name='vegetation emissivity', &
            ptr_pft=pps%emv, default='inactive')

       call hist_addfld1d (fname='Z0MV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, momentum', &
            ptr_pft=pps%z0mv, default='inactive')

       call hist_addfld1d (fname='Z0HV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, sensible heat', &
            ptr_pft=pps%z0hv, default='inactive')

       call hist_addfld1d (fname='Z0QV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, latent heat', &
            ptr_pft=pps%z0qv, default='inactive')

       call hist_addfld1d (fname='DEWMX', units='mm', &
            avgflag='A', long_name='Maximum allowed dew', &
            ptr_pft=pps%dewmx, default='inactive')

       call hist_addfld1d (fname='LNCSUN', units='gN/m^2', &
            avgflag='A', long_name='leaf N concentration per unit projected LAI', &
            ptr_pft=pps%lncsun, default='inactive')

       call hist_addfld1d (fname='LNCSHA', units='gN/m^2', &
            avgflag='A', long_name='leaf N concentration per unit projected LAI', &
            ptr_pft=pps%lncsha, default='inactive')

       call hist_addfld1d (fname='VCMXSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='sunlit leaf Vcmax', &
            ptr_pft=pps%vcmxsun, default='inactive')

       call hist_addfld1d (fname='VCMXSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='shaded leaf Vcmax', &
            ptr_pft=pps%vcmxsha, default='inactive')

       call hist_addfld1d (fname='FSUN', units='proportion', &
            avgflag='A', long_name='sunlit fraction of canopy', &
            ptr_pft=pps%fsun, default='inactive')

       call hist_addfld1d (fname='GDIR', units='proportion', &
            avgflag='A', long_name='leaf projection in solar direction', &
            ptr_pft=pps%gdir, default='inactive')

       call hist_addfld1d (fname='CISUN', units='Pa', &
            avgflag='A', long_name='sunlit intracellular CO2', &
            ptr_pft=pps%cisun, default='inactive')

       call hist_addfld1d (fname='CISHA', units='Pa', &
            avgflag='A', long_name='shaded intracellular CO2', &
            ptr_pft=pps%cisha, default='inactive')

       if (use_c13) then
          call hist_addfld1d (fname='ALPHAPSNSUN', units='proportion', &
               avgflag='A', long_name='sunlit c13 fractionation', &
               ptr_pft=pps%alphapsnsun, default='inactive')

          call hist_addfld1d (fname='ALPHAPSNSHA', units='proportion', &
               avgflag='A', long_name='shaded c13 fractionation', &
               ptr_pft=pps%alphapsnsha, default='inactive')
       endif

       call hist_addfld1d (fname='FWET', units='proportion', &
            avgflag='A', long_name='fraction of canopy that is wet', &
            ptr_pft=pps%fwet, default='inactive')

       call hist_addfld1d (fname='FDRY', units='proportion', &
            avgflag='A', long_name='fraction of foliage that is green and dry', &
            ptr_pft=pps%fdry, default='inactive')

       call hist_addfld1d (fname='DT_VEG', units='K', &
            avgflag='A', long_name='change in t_veg, last iteration', &
            ptr_pft=pps%dt_veg, default='inactive')

       call hist_addfld1d (fname='HTOP', units='m', &
            avgflag='A', long_name='canopy top', &
            ptr_pft=pps%htop)

       call hist_addfld1d (fname='HBOT', units='m', &
            avgflag='A', long_name='canopy bottom', &
            ptr_pft=pps%hbot, default='inactive')

       call hist_addfld1d (fname='Z0M', units='m', &
            avgflag='A', long_name='momentum roughness length', &
            ptr_pft=pps%z0m, default='inactive')

       call hist_addfld1d (fname='DISPLA', units='m', &
            avgflag='A', long_name='displacement height', &
            ptr_pft=pps%displa, default='inactive')

       call hist_addfld1d (fname='U10_DUST', units='m/s', &
            avgflag='A', long_name='10-m wind for dust model', &
            ptr_pft=pps%u10, default='inactive')

       call hist_addfld1d (fname='RAM1', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_pft=pps%ram1, default='inactive')

       call hist_addfld1d (fname='FV', units='m/s', &
            avgflag='A', long_name='friction velocity for dust model', &
            ptr_pft=pps%fv, default='inactive')

       call hist_addfld2d (fname='ROOTFR', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='fraction of roots in each soil layer', &
            ptr_pft=pps%rootfr, default='inactive')

       call hist_addfld2d (fname='ROOTR', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective fraction of roots in each soil layer', &
            ptr_pft=pps%rootr, default='inactive')

       call hist_addfld2d (fname='RRESIS', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='root resistance in each soil layer', &
            ptr_pft=pps%rresis, default='inactive')

       call hist_addfld2d (fname='ALBD', units='proportion', type2d='numrad', &
            avgflag='A', long_name='surface albedo (direct)', &
            ptr_pft=pps%albd, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld2d (fname='ALBI', units='proportion', type2d='numrad', &
            avgflag='A', long_name='surface albedo (indirect)', &
            ptr_pft=pps%albi, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld2d (fname='FABD', units='proportion', type2d='numrad', &
            avgflag='A', long_name='flux absorbed by veg per unit direct flux', &
            ptr_pft=pps%fabd, default='inactive')

       call hist_addfld2d (fname='FABI', units='proportion', type2d='numrad', &
            avgflag='A', long_name='flux absorbed by veg per unit indirect flux', &
            ptr_pft=pps%fabi, default='inactive')

       call hist_addfld2d (fname='FTDD', units='proportion', type2d='numrad', &
            avgflag='A', long_name='down direct flux below veg per unit dir flx', &
            ptr_pft=pps%ftdd, default='inactive')

       call hist_addfld2d (fname='FTID', units='proportion', type2d='numrad', &
            avgflag='A', long_name='down indirect flux below veg per unit dir flx', &
            ptr_pft=pps%ftid, default='inactive')

       call hist_addfld2d (fname='FTII', units='proportion', type2d='numrad', &
            avgflag='A', long_name='down indirect flux below veg per unit indirect flx', &
            ptr_pft=pps%ftii, default='inactive')

       call hist_addfld2d (fname='OMEGA', units='proportion', type2d='numrad', &
            avgflag='A', long_name='fraction of intercepted radiation that is scattered', &
            ptr_pft=pps%omega, default='inactive')

       call hist_addfld2d (fname='EFF_KID', units='none', type2d='numrad', &
            avgflag='A', long_name='effective extinction coefficient for indirect from direct', &
            ptr_pft=pps%eff_kid, default='inactive')

       call hist_addfld2d (fname='EFF_KII', units='none', type2d='numrad', &
            avgflag='A', long_name='effective extinction coefficient for indirect from indirect', &
            ptr_pft=pps%eff_kii, default='inactive')

       call hist_addfld2d (fname='SUN_FAID', units='proportion', type2d='numrad', &
            avgflag='A', long_name='fraction sun canopy absorbed indirect from direct', &
            ptr_pft=pps%sun_faid, default='inactive')

       call hist_addfld2d (fname='SUN_FAII', units='proportion', type2d='numrad', &
            avgflag='A', long_name='fraction sun canopy absorbed indirect from indirect', &
            ptr_pft=pps%sun_faii, default='inactive')

       call hist_addfld2d (fname='SHA_FAID', units='proportion', type2d='numrad', &
            avgflag='A', long_name='fraction shade canopy absorbed indirect from direct', &
            ptr_pft=pps%sha_faid, default='inactive')

       call hist_addfld2d (fname='SHA_FAII', units='proportion', type2d='numrad', &
            avgflag='A', long_name='fraction shade canopy absorbed indirect from indirect', &
            ptr_pft=pps%sha_faii, default='inactive')

       if ( crop_prog )then

          call hist_addfld1d (fname='GDD0', units='ddays', &
               avgflag='A', long_name='Growing degree days base  0C from planting', &
               ptr_pft=pps%gdd0, default='inactive')

          call hist_addfld1d (fname='GDD8', units='ddays', &
               avgflag='A', long_name='Growing degree days base  8C from planting', &
               ptr_pft=pps%gdd8, default='inactive')

          call hist_addfld1d (fname='GDD10', units='ddays', &
               avgflag='A', long_name='Growing degree days base 10C from planting', &
               ptr_pft=pps%gdd10, default='inactive')

          call hist_addfld1d (fname='GDD020', units='ddays', &
               avgflag='A', long_name='Twenty year average of growing degree days base  0C from planting', &
               ptr_pft=pps%gdd020, default='inactive')

          call hist_addfld1d (fname='GDD820', units='ddays', &
               avgflag='A', long_name='Twenty year average of growing degree days base  8C from planting', &
               ptr_pft=pps%gdd820, default='inactive')

          call hist_addfld1d (fname='GDD1020', units='ddays', &
               avgflag='A', long_name='Twenty year average of growing degree days base 10C from planting', &
               ptr_pft=pps%gdd1020, default='inactive')

          call hist_addfld1d (fname='GDDPLANT', units='ddays', &
               avgflag='A', long_name='Accumulated growing degree days past planting date for crop', &
               ptr_pft=pps%gddplant, default='inactive')

          call hist_addfld1d (fname='GDDHARV', units='ddays', &
               avgflag='A', long_name='Growing degree days (gdd) needed to harvest', &
               ptr_pft=pps%gddmaturity, default='inactive')

          call hist_addfld1d (fname='GDDTSOI', units='ddays', &
               avgflag='A', long_name='Growing degree-days from planting (top two soil layers)', &
               ptr_pft=pps%gddtsoi, default='inactive')

       end if

       !-------------------------------
       ! Column physical state variables not already defined by default
       !-------------------------------
       
       call hist_addfld1d (fname='EMG', units='proportion', &
            avgflag='A', long_name='ground emissivity', &
            ptr_col=cps%emg, default='inactive')

       call hist_addfld1d (fname='Z0MG', units='m', &
            avgflag='A', long_name='roughness length over ground, momentum', &
            ptr_col=cps%z0mg, default='inactive')

       call hist_addfld1d (fname='Z0HG', units='m', &
            avgflag='A', long_name='roughness length over ground, sensible heat', &
            ptr_col=cps%z0hg, default='inactive')

       call hist_addfld1d (fname='Z0QG', units='m', &
            avgflag='A', long_name='roughness length over ground, latent heat', &
            ptr_col=cps%z0qg, default='inactive')

       call hist_addfld1d (fname='BETA', units='none', &
            avgflag='A', long_name='coefficient of convective velocity', &
            ptr_col=cps%beta, default='inactive')

       call hist_addfld1d (fname='ZII', units='m', &
            avgflag='A', long_name='convective boundary height', &
            ptr_col=cps%zii, default='inactive')

       call hist_addfld1d (fname='WF', units='proportion', &
            avgflag='A', long_name='soil water as frac. of whc for top 0.5 m', &
            ptr_col=cps%wf, default='inactive')

       call hist_addfld1d (fname='FPI', units='proportion', &
            avgflag='A', long_name='fraction of potential immobilization', &
            ptr_col=cps%fpi)

       call hist_addfld1d (fname='FPG', units='proportion', &
            avgflag='A', long_name='fraction of potential gpp', &
            ptr_col=cps%fpg)

       call hist_addfld1d (fname='ANNSUM_COUNTER', units='s', &
            avgflag='A', long_name='seconds since last annual accumulator turnover', &
            ptr_col=cps%annsum_counter, default='inactive')

       call hist_addfld1d (fname='CANNSUM_NPP', units='gC/m^2/s', &
            avgflag='A', long_name='annual sum of column-level NPP', &
            ptr_col=cps%cannsum_npp, default='inactive')

       call hist_addfld1d (fname='CANNAVG_T2M', units='K', &
            avgflag='A', long_name='annual average of 2m air temperature', &
            ptr_col=cps%cannavg_t2m, default='inactive')

       call hist_addfld2d (fname='FRAC_ICEOLD', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='fraction of ice relative to the tot water', &
            ptr_col=cps%frac_iceold, default='inactive')

       call hist_addfld2d (fname='EFF_POROSITY', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective porosity = porosity - vol_ice', &
            ptr_col=cps%eff_porosity, default='inactive')

       call hist_addfld2d (fname='ROOTR_COLUMN', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective fraction of roots in each soil layer', &
            ptr_col=cps%rootr_column, default='inactive')

       call hist_addfld2d (fname='ALBGRD', units='proportion', type2d='numrad', &
            avgflag='A', long_name='ground albedo (direct)', &
            ptr_col=cps%albgrd, default='inactive')

       call hist_addfld2d (fname='ALBGRI', units='proportion', type2d='numrad', &
            avgflag='A', long_name='ground albedo (indirect)', &
            ptr_col=cps%albgri, default='inactive')

       call hist_addfld1d (fname='ME',  units='proportion', &
            avgflag='A', long_name='moisture of extinction', &
            ptr_col=cps%me, default='inactive')

       call hist_addfld1d (fname='FIRE_PROB',  units='0-1', &
            avgflag='A', long_name='daily fire probability', &
            ptr_col=cps%fire_prob, default='inactive')

       call hist_addfld1d (fname='MEAN_FIRE_PROB',  units='0-1', &
            avgflag='A', long_name='e-folding mean of daily fire probability', &
            ptr_col=cps%mean_fire_prob)

       call hist_addfld1d (fname='FIRESEASONL',  units='days', &
            avgflag='A', long_name='annual fire season length', &
            ptr_col=cps%fireseasonl)

       call hist_addfld1d (fname='FAREA_BURNED',  units='proportion', &
            avgflag='A', long_name='timestep fractional area burned', &
            ptr_col=cps%farea_burned, default='inactive')

       call hist_addfld1d (fname='ANN_FAREA_BURNED',  units='proportion', &
            avgflag='A', long_name='annual total fractional area burned', &
            ptr_col=cps%ann_farea_burned)

       !-------------------------------
       ! Energy flux variables not already defined by default - native PFT 
       !-------------------------------

       call hist_addfld1d (fname='PARSUN', units='W/m^2', &
            avgflag='A', long_name='average absorbed PAR for sunlit leaves', &
            ptr_pft=pef%parsun, default='inactive')

       call hist_addfld1d (fname='PARSHA', units='W/m^2', &
            avgflag='A', long_name='average absorbed PAR for shaded leaves', &
            ptr_pft=pef%parsha, default='inactive')

       call hist_addfld1d (fname='DLRAD', units='W/m^2', &
            avgflag='A', long_name='downward longwave radiation below the canopy', &
            ptr_pft=pef%dlrad, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='ULRAD', units='W/m^2', &
            avgflag='A', long_name='upward longwave radiation above the canopy', &
            ptr_pft=pef%ulrad, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='EFLX_LH_TOT', units='W/m^2', &
            avgflag='A', long_name='total latent heat flux [+ to atm]', &
            ptr_pft=pef%eflx_lh_tot, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='EFLX_SOIL_GRND', units='W/m^2', &
            avgflag='A', long_name='soil heat flux [+ into soil]', &
            ptr_pft=pef%eflx_soil_grnd, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='CGRND', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil energy flux wrt to soil temp', &
            ptr_pft=pef%cgrnd, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='CGRNDL', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil latent heat flux wrt soil temp', &
            ptr_pft=pef%cgrndl, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='CGRNDS', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil sensible heat flux wrt soil temp', &
            ptr_pft=pef%cgrnds, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='EFLX_GNET', units='W/m^2', &
            avgflag='A', long_name='net heat flux into ground', &
            ptr_pft=pef%eflx_gnet, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='DGNETDT', units='W/m^2/K', &
            avgflag='A', long_name='derivative of net ground heat flux wrt soil temp', &
            ptr_pft=pef%dgnetdT, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld2d (fname='SUN_ADD', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='sun canopy absorbed direct from direct', &
            ptr_pft=pef%sun_add, default='inactive')

       call hist_addfld2d (fname='TOT_AID', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='total canopy absorbed indirect from direct', &
            ptr_pft=pef%tot_aid, default='inactive')

       call hist_addfld2d (fname='SUN_AID', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='sun canopy absorbed indirect from direct', &
            ptr_pft=pef%sun_aid, default='inactive')

       call hist_addfld2d (fname='SUN_AII', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='sun canopy absorbed indirect from indirect', &
            ptr_pft=pef%sun_aii, default='inactive')

       call hist_addfld2d (fname='SHA_AID', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='shade canopy absorbed indirect from direct', &
            ptr_pft=pef%sha_aid, default='inactive')

       call hist_addfld2d (fname='SHA_AII', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='shade canopy absorbed indirect from indirect', &
            ptr_pft=pef%sha_aii, default='inactive')

       call hist_addfld2d (fname='SUN_ATOT', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='sun canopy total absorbed', &
            ptr_pft=pef%sun_atot, default='inactive')

       call hist_addfld2d (fname='SHA_ATOT', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='shade canopy total absorbed', &
            ptr_pft=pef%sha_atot, default='inactive')

       call hist_addfld2d (fname='SUN_ALF', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='sun canopy total absorbed by leaves', &
            ptr_pft=pef%sun_alf, default='inactive')

       call hist_addfld2d (fname='SHA_ALF', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='shade canopy total absored by leaves', &
            ptr_pft=pef%sha_alf, default='inactive')

       call hist_addfld2d (fname='SUN_APERLAI', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='sun canopy total absorbed per unit LAI', &
            ptr_pft=pef%sun_aperlai, default='inactive')

       call hist_addfld2d (fname='SHA_APERLAI', units='W/m^2', type2d='numrad', &
            avgflag='A', long_name='shade canopy total absorbed per unit LAI', &
            ptr_pft=pef%sha_aperlai, default='inactive')

       !-------------------------------
       ! Water flux variables not already defined by default  - native PFT
       !-------------------------------

       call hist_addfld1d (fname='QFLX_RAIN_GRND', units='mm H2O/s', &
            avgflag='A', long_name='rain on ground after interception', &
            ptr_pft=pwf%qflx_rain_grnd, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='QFLX_SNOW_GRND', units='mm H2O/s', &
            avgflag='A', long_name='snow on ground after interception', &
            ptr_pft=pwf%qflx_snow_grnd, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='QFLX_EVAP_GRND', units='mm H2O/s', &
            avgflag='A', long_name='ground surface evaporation', &
            ptr_pft=pwf%qflx_evap_grnd, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='QFLX_EVAP_VEG', units='mm H2O/s', &
            avgflag='A', long_name='vegetation evaporation', &
            ptr_pft=pwf%qflx_evap_veg, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='QFLX_EVAP_TOT', units='mm H2O/s', &
            avgflag='A', long_name='qflx_evap_soi + qflx_evap_can + qflx_tran_veg', &
            ptr_pft=pwf%qflx_evap_tot, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='QFLX_DEW_GRND', units='mm H2O/s', &
            avgflag='A', long_name='ground surface dew formation', &
            ptr_pft=pwf%qflx_dew_grnd, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='QFLX_SUB_SNOW', units='mm H2O/s', &
            avgflag='A', long_name='sublimation rate from snow pack', &
            ptr_pft=pwf%qflx_sub_snow, default='inactive', c2l_scale_type='urbanf')

       call hist_addfld1d (fname='QFLX_DEW_SNOW', units='mm H2O/s', &
            avgflag='A', long_name='surface dew added to snow pacK', &
            ptr_pft=pwf%qflx_dew_snow, default='inactive', c2l_scale_type='urbanf')

    end if

    call hist_addfld1d (fname='SNORDSL', units='m^-6', &
         avgflag='A', long_name='top snow layer effective grain radius', &
         ptr_col=cps%snw_rds_top, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOTTOPL', units='K/m', &
         avgflag='A', long_name='snow temperature (top layer)', &
         ptr_col=cps%snot_top, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOdTdzL', units='K/m', &
         avgflag='A', long_name='top snow layer temperature gradient (land)', &
         ptr_col=cps%dTdz_top, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOLIQFL', units='fraction', &
         avgflag='A', long_name='top snow layer liquid water fraction (land)', &
         ptr_col=cps%sno_liq_top, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation from snow', &
         ptr_pft=pef%fsr_sno_vd, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation from snow', &
         ptr_pft=pef%fsr_sno_nd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation from snow', &
         ptr_pft=pef%fsr_sno_vi, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation from snow', &
         ptr_pft=pef%fsr_sno_ni, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation on snow', &
         ptr_pft=pef%fsds_sno_vd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation on snow', &
         ptr_pft=pef%fsds_sno_nd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation on snow', &
         ptr_pft=pef%fsds_sno_vi, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation on snow', &
         ptr_pft=pef%fsds_sno_ni, &
         default='inactive')

    call hist_addfld1d (fname='H2OSNO_TOP', units='kg/m2', &
         avgflag='A', long_name='mass of snow in top snow layer', &
         ptr_col=cps%h2osno_top, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNOBCMCL', units='kg/m2', &
         avgflag='A', long_name='mass of BC in snow column', &
         ptr_col=cps%mss_bc_col, set_lake=spval, set_urb=spval)
    
    call hist_addfld1d (fname='SNOBCMSL', units='kg/m2', &
         avgflag='A', long_name='mass of BC in top snow layer', &
         ptr_col=cps%mss_bc_top, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNOOCMCL', units='kg/m2', &
         avgflag='A', long_name='mass of OC in snow column', &
         ptr_col=cps%mss_oc_col, set_lake=spval, set_urb=spval)
   
    call hist_addfld1d (fname='SNOOCMSL', units='kg/m2', &
         avgflag='A', long_name='mass of OC in top snow layer', &
         ptr_col=cps%mss_oc_top, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNODSTMCL', units='kg/m2', &
         avgflag='A', long_name='mass of dust in snow column', &
         ptr_col=cps%mss_dst_col, set_lake=spval, set_urb=spval)
    
    call hist_addfld1d (fname='SNODSTMSL', units='kg/m2', &
         avgflag='A', long_name='mass of dust in top snow layer', &
         ptr_col=cps%mss_dst_top, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='DSTDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total dust deposition (dry+wet) from atmosphere', &
         ptr_col=cwf%flx_dst_dep, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='BCDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total BC deposition (dry+wet) from atmosphere', &
         ptr_col=cwf%flx_bc_dep, set_lake=spval, set_urb=spval)
   
    call hist_addfld1d (fname='OCDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total OC deposition (dry+wet) from atmosphere', &
         ptr_col=cwf%flx_oc_dep, set_lake=spval, set_urb=spval)

    if (use_snicar_frc) then
       call hist_addfld1d (fname='SNOAERFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of all aerosols in snow (land) ', &
            ptr_pft=pef%sfc_frc_aer, set_lake=spval, set_urb=spval)
   
       call hist_addfld1d (fname='SNOAERFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of all aerosols in snow, averaged only when snow is present (land)', &
            ptr_pft=pef%sfc_frc_aer_sno, set_lake=spval, set_urb=spval)

       call hist_addfld1d (fname='SNOBCFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of BC in snow (land) ', &
            ptr_pft=pef%sfc_frc_bc, set_lake=spval, set_urb=spval)

       call hist_addfld1d (fname='SNOBCFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of BC in snow, averaged only when snow is present (land)', &
            ptr_pft=pef%sfc_frc_bc_sno, set_lake=spval, set_urb=spval)

       call hist_addfld1d (fname='SNOOCFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of OC in snow (land) ', &
            ptr_pft=pef%sfc_frc_oc, set_lake=spval, set_urb=spval)

       call hist_addfld1d (fname='SNOOCFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of OC in snow, averaged only when snow is present (land)', &
            ptr_pft=pef%sfc_frc_oc_sno, set_lake=spval, set_urb=spval)

       call hist_addfld1d (fname='SNODSTFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of dust in snow (land) ', &
            ptr_pft=pef%sfc_frc_dst, set_lake=spval, set_urb=spval)

       call hist_addfld1d (fname='SNODSTFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of dust in snow, averaged only when snow is present (land)', &
            ptr_pft=pef%sfc_frc_dst_sno, set_lake=spval, set_urb=spval)
    end if

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

end module histFldsMod
