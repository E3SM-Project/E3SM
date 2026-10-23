module lnd2atmMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use abortutils             , only : endrun
  use shr_megan_mod        , only : shr_megan_mechcomps_n
  use shr_fan_mod          , only : shr_fan_to_atm
  use elm_varpar           , only : numrad, ndst, nlevgrnd, nlevsno, nlevsoi !ndst = number of dust bins.
  use elm_varcon           , only : rair, grav, cpair, hfus, tfrz, spval
  use elm_varctl           , only : iulog, use_c13, use_cn, use_lch4, use_voc, use_fates, use_atm_downscaling_to_topunit, use_fan
  use elm_varctl           , only : use_lnd_rof_two_way, use_finetop_rad
  use seq_drydep_mod   , only : n_drydep, drydep_method, DD_XLND
  use decompMod            , only : bounds_type
  use subgridAveMod        , only : p2t
  use lnd2atmType          , only : lnd2atm_type
  use atm2lndType          , only : atm2lnd_type
  use CH4Mod               , only : ch4_type
  use DUSTMod              , only : dust_type
  use DryDepVelocity       , only : drydepvel_type
  use VocEmissionMod       , only : vocemis_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityType , only : frictionvel_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use GridcellType         , only : grc_pp
  use TopounitDataType     , only : top_es, top_af                 ! To calculate t_rad at topounit level needed in downscaling
  use GridcellDataType     , only : grc_ef, grc_ws, grc_wf
  use ColumnDataType       , only : col_ws, col_wf, col_cf, col_es, col_nf
  use VegetationDataType   , only : veg_es, veg_ef, veg_ws, veg_wf
  use SoilHydrologyType    , only : soilhydrology_type 
  use SedFluxType          , only : sedflux_type
  use spmdmod          , only: masterproc
  use elm_varctl     , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm
  public :: lnd2atm_minimal
  private :: l2a_map_init
  private :: l2a_p2g_1d
  private :: l2a_p2g_2d
  private :: l2a_c2g_1d
  private :: l2a_c2g_2d

  integer, parameter :: unity = 0, urbanf = 1, urbans = 2
  integer, parameter :: natveg = 3, veg =4, ice=5, nonurb=6, lake=7
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine lnd2atm_minimal(bounds, &
      surfalb_vars, solarabs_vars, energyflux_vars, atm2lnd_vars, lnd2atm_vars)
    !
    ! !DESCRIPTION:
    ! Compute elm_l2a_vars component of gridcell derived type. This routine computes
    ! the bare minimum of components necessary to get the first step of a
    ! run started.
    !
    ! !USES:
      !$acc routine seq
    use elm_varcon, only : sb

    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    type(surfalb_type)    , intent(in)    :: surfalb_vars
    type(solarabs_type)   , intent(in)    :: solarabs_vars
    type(energyflux_type) , intent(in)    :: energyflux_vars
    type(atm2lnd_type)    , intent(in)    :: atm2lnd_vars
    type(lnd2atm_type)    , intent(inout) :: lnd2atm_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g, t, ib                                  ! index
    integer  :: npts_p                                          ! number of entries in pts_p
    integer  :: npts_c                                          ! number of entries in pts_c
    integer  :: pts_p(bounds%endp-bounds%begp+1)                ! active patches with nonzero wtgcell
    integer  :: pts_c(bounds%endc-bounds%begc+1)                ! active columns with nonzero wtgcell
    real(r8) :: scale_c2l(bounds%begc:bounds%endc,unity:urbans) ! c2l scale factors [col, scale type]
    real(r8) :: scale_l2g(bounds%begl:bounds%endl)              ! l2g scale factors for unity [landunit]
    !------------------------------------------------------------------------
    associate( &
      h2osno => col_ws%h2osno  , &
      h2osno_grc => lnd2atm_vars%h2osno_grc , &
      h2osoi_vol => col_ws%h2osoi_vol , &
      h2osoi_vol_grc => lnd2atm_vars%h2osoi_vol_grc , &
      albd_patch => surfalb_vars%albd_patch , &
      albd_grc   => lnd2atm_vars%albd_grc   , &
      apparent_albd_grc   => lnd2atm_vars%apparent_albd_grc   , &
      albi_patch => surfalb_vars%albi_patch , &
      albi_grc   => lnd2atm_vars%albi_grc   , &
      apparent_albi_grc   => lnd2atm_vars%apparent_albi_grc   , &
      eflx_lwrad_out => veg_ef%eflx_lwrad_out , &
      eflx_lwrad_out_grc => lnd2atm_vars%eflx_lwrad_out_grc   , &
      forc_solad_pp_grc  => atm2lnd_vars%forc_solad_pp_grc    , &
      forc_solai_pp_grc  => atm2lnd_vars%forc_solai_pp_grc    , &
      fsr_vis_d_patch    => solarabs_vars%fsr_vis_d_patch     , &
      fsr_vis_i_patch    => solarabs_vars%fsr_vis_i_patch     , &
      fsr_nir_d_patch    => solarabs_vars%fsr_nir_d_patch     , &
      fsr_nir_i_patch    => solarabs_vars%fsr_nir_i_patch     , &
      fsr_vis_d_grc      => lnd2atm_vars%fsr_vis_d_grc        , &
      fsr_vis_i_grc      => lnd2atm_vars%fsr_vis_i_grc        , &
      fsr_nir_d_grc      => lnd2atm_vars%fsr_nir_d_grc        , &
      fsr_nir_i_grc      => lnd2atm_vars%fsr_nir_i_grc        , &
      sky_view_factor    => grc_pp%sky_view_factor              &
      )

    ! Map to gridcells over the active points with nonzero gridcell weight only
    call l2a_map_init(bounds, npts_p, pts_p, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,unity:urbans), scale_l2g(bounds%begl:bounds%endl))

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         h2osno (bounds%begc:bounds%endc), &
         h2osno_grc(bounds%begg:bounds%endg))

    do g = bounds%begg,bounds%endg
       h2osno_grc(g) = h2osno_grc(g)/1000._r8
    end do

    call l2a_c2g_2d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         nlevgrnd, &
         h2osoi_vol (bounds%begc:bounds%endc,:), &
         h2osoi_vol_grc(bounds%begg:bounds%endg,:))

    call l2a_p2g_2d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         numrad, &
         albd_patch(bounds%begp:bounds%endp,:), &
         albd_grc (bounds%begg:bounds%endg,:))

    call l2a_p2g_2d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         numrad, &
         albi_patch(bounds%begp:bounds%endp,:), &
         albi_grc (bounds%begg:bounds%endg,:))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         eflx_lwrad_out (bounds%begp:bounds%endp), &
         eflx_lwrad_out_grc (bounds%begg:bounds%endg))

    if (use_finetop_rad) then
         do g = bounds%begg,bounds%endg
            eflx_lwrad_out_grc(g) = eflx_lwrad_out_grc(g)*sky_view_factor(g)
         end do
    end if

    do g = bounds%begg,bounds%endg
       lnd2atm_vars%t_rad_grc(g) = sqrt(sqrt(eflx_lwrad_out_grc(g)/sb))
    end do

    ! Calculate topounit level eflx_lwrad_out_topo for downscaling purpose
    if (use_atm_downscaling_to_topunit) then
       call p2t(bounds, &
            eflx_lwrad_out (bounds%begp:bounds%endp), &
            top_es%eflx_lwrad_out_topo      (bounds%begt:bounds%endt), &
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2t_scale_type='unity')
    
       do t = bounds%begt,bounds%endt
          top_es%t_rad(t) = sqrt(sqrt(top_es%eflx_lwrad_out_topo(t)/sb))   
       end do
    end if
    
    if (use_finetop_rad) then
       ! calculate gridcell level reflected radiation
       call l2a_p2g_1d(bounds, npts_p, pts_p, &
            scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
            fsr_vis_d_patch(bounds%begp:bounds%endp), &
            fsr_vis_d_grc (bounds%begg:bounds%endg))

       call l2a_p2g_1d(bounds, npts_p, pts_p, &
            scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
            fsr_vis_i_patch(bounds%begp:bounds%endp), &
            fsr_vis_i_grc (bounds%begg:bounds%endg))

       call l2a_p2g_1d(bounds, npts_p, pts_p, &
            scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
            fsr_nir_d_patch(bounds%begp:bounds%endp), &
            fsr_nir_d_grc (bounds%begg:bounds%endg))

       call l2a_p2g_1d(bounds, npts_p, pts_p, &
            scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
            fsr_nir_i_patch(bounds%begp:bounds%endp), &
            fsr_nir_i_grc (bounds%begg:bounds%endg))

       ! calculate gridcell level apparent albedo 
       do g = bounds%begg,bounds%endg 
           do ib = 1, numrad
               apparent_albd_grc(g,ib) = 1._r8
               apparent_albi_grc(g,ib) = 1._r8
           end do

           ! calculate direct albedo
           ib = 1
           if (forc_solad_pp_grc(g,ib) > 0._r8) then
               apparent_albd_grc(g,ib) = fsr_vis_d_grc(g)/forc_solad_pp_grc(g,ib)*sky_view_factor(g)
           end if
           ib = 2
           if (forc_solad_pp_grc(g,ib) > 0._r8) then
               apparent_albd_grc(g,ib) = fsr_nir_d_grc(g)/forc_solad_pp_grc(g,ib)*sky_view_factor(g)
           end if

           ! calculate diffuse albedo
           ib = 1
           if (forc_solai_pp_grc(g,ib) > 0._r8) then
               apparent_albi_grc(g,ib) = fsr_vis_i_grc(g)/forc_solai_pp_grc(g,ib)*sky_view_factor(g)
           end if
           ib = 2
           if (forc_solai_pp_grc(g,ib) > 0._r8) then
               apparent_albi_grc(g,ib) = fsr_nir_i_grc(g)/forc_solai_pp_grc(g,ib)*sky_view_factor(g)
           end if

           ! limit albedo to be <= 1
           do ib = 1, numrad
               if (apparent_albd_grc(g,ib) > 1._r8) then
                   apparent_albd_grc(g,ib) = 1._r8
               end if

               if (apparent_albi_grc(g,ib) > 1._r8) then
                   apparent_albi_grc(g,ib) = 1._r8
               end if
           end do
        end do
     end if

    end associate

  end subroutine lnd2atm_minimal

  !------------------------------------------------------------------------
  subroutine lnd2atm(bounds, &
       atm2lnd_vars, surfalb_vars, frictionvel_vars, &
       energyflux_vars, &
       solarabs_vars, drydepvel_vars, &
       vocemis_vars, dust_vars, ch4_vars, soilhydrology_vars, &
       sedflux_vars, lnd2atm_vars)
    !
    ! !DESCRIPTION:
    ! Compute lnd2atm_vars component of gridcell derived type
    !
    ! !USES:
      !$acc routine seq
    use CH4varcon  , only : ch4offline
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)     :: bounds
    type(atm2lnd_type)     , intent(in)     :: atm2lnd_vars
    type(surfalb_type)     , intent(in)     :: surfalb_vars
    type(frictionvel_type) , intent(in)     :: frictionvel_vars
    type(energyflux_type)  , intent(in)     :: energyflux_vars
    type(solarabs_type)    , intent(in)     :: solarabs_vars
    type(drydepvel_type)   , intent(in)     :: drydepvel_vars
    type(vocemis_type)     , intent(in)     :: vocemis_vars
    type(dust_type)        , intent(in)     :: dust_vars
    type(ch4_type)         , intent(in)     :: ch4_vars
    type(soilhydrology_type), intent(in)    :: soilhydrology_vars
    type(sedflux_type)     , intent(in)     :: sedflux_vars
    type(lnd2atm_type)     , intent(inout)  :: lnd2atm_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g, lvl             ! index
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    integer  :: npts_p                                          ! number of entries in pts_p
    integer  :: npts_c                                          ! number of entries in pts_c
    integer  :: pts_p(bounds%endp-bounds%begp+1)                ! active patches with nonzero wtgcell
    integer  :: pts_c(bounds%endc-bounds%begc+1)                ! active columns with nonzero wtgcell
    real(r8) :: scale_c2l(bounds%begc:bounds%endc,unity:urbans) ! c2l scale factors [col, scale type]
    real(r8) :: scale_l2g(bounds%begl:bounds%endl)              ! l2g scale factors for unity [landunit]
    !------------------------------------------------------------------------
    associate( &
      t_ref2m     => veg_es%t_ref2m , &
      t_ref2m_grc => lnd2atm_vars%t_ref2m_grc       , &
      q_ref2m     => veg_ws%q_ref2m , &
      q_ref2m_grc => lnd2atm_vars%q_ref2m_grc      , &
      u10_elm_patch => frictionvel_vars%u10_elm_patch , &
      u10_with_gusts_elm_patch => frictionvel_vars%u10_with_gusts_elm_patch, &
      u_ref10m_grc => lnd2atm_vars%u_ref10m_grc      , &
      u_ref10m_with_gusts_grc => lnd2atm_vars%u_ref10m_with_gusts_grc      , &
      taux     => veg_ef%taux , &
      taux_grc => lnd2atm_vars%taux_grc      , &
      tauy     => veg_ef%tauy , &
      tauy_grc => lnd2atm_vars%tauy_grc      , &
      qflx_evap_tot     => veg_wf%qflx_evap_tot , &
      qflx_evap_tot_grc => lnd2atm_vars%qflx_evap_tot_grc     , &
      fsa_patch => solarabs_vars%fsa_patch , &
      fsa_grc   =>lnd2atm_vars%fsa_grc    , &
      fv_patch  => frictionvel_vars%fv_patch , &
      fv_grc    => lnd2atm_vars%fv_grc       , &
      ram1_patch  => frictionvel_vars%ram1_patch , &
      ram1_grc    => lnd2atm_vars%ram1_grc       , &
      eflx_sh_tot => veg_ef%eflx_sh_tot , &
      eflx_sh_tot_grc => lnd2atm_vars%eflx_sh_tot_grc      , &
      eflx_lh_tot => veg_ef%eflx_lh_tot , &
      eflx_lh_tot_grc => lnd2atm_vars%eflx_lh_tot_grc      , &
      nee             => col_cf%nee, &
      nee_grc         => lnd2atm_vars%nee_grc   , &
      velocity_patch  => drydepvel_vars%velocity_patch , &
      ddvel_grc       => lnd2atm_vars%ddvel_grc        , &
      flx_mss_vrt_dst_patch => dust_vars%flx_mss_vrt_dst_patch, &
      flxdst_grc      => lnd2atm_vars%flxdst_grc        , &
      ch4_surf_flux_tot_col => ch4_vars%ch4_surf_flux_tot_col , &
      flux_ch4_grc  => lnd2atm_vars%flux_ch4_grc      , &
      qflx_runoff   => col_wf%qflx_runoff , &
      qflx_rofliq_grc => lnd2atm_vars%qflx_rofliq_grc   , &
      qflx_surf      => col_wf%qflx_surf , &
      qflx_rofliq_qsur_grc  => lnd2atm_vars%qflx_rofliq_qsur_grc   , &
      qflx_h2osfc_surf      => col_wf%qflx_h2osfc_surf , &
      qflx_rofliq_qsurp_grc => lnd2atm_vars%qflx_rofliq_qsurp_grc  , &
      qflx_irr_demand       => col_wf%qflx_irr_demand , &
      qflx_irr_demand_grc   => lnd2atm_vars%qflx_irr_demand_grc   , &
      qflx_drain            => col_wf%qflx_drain , &
      qflx_rofliq_qsub_grc  => lnd2atm_vars%qflx_rofliq_qsub_grc   , &
      qflx_drain_perched    => col_wf%qflx_drain_perched , &
      qflx_rofliq_qsubp_grc => lnd2atm_vars%qflx_rofliq_qsubp_grc  , &
      qflx_qrgwl            => col_wf%qflx_qrgwl                   , &
      qflx_rofliq_qgwl_grc  => lnd2atm_vars%qflx_rofliq_qgwl_grc   , &
      qflx_snwcp_ice        => col_wf%qflx_snwcp_ice,  &
      qflx_ice_runoff_xs    => col_wf%qflx_ice_runoff_xs,  &
      qflx_rofice_grc       => lnd2atm_vars%qflx_rofice_grc     ,  &
      endwb                 => col_ws%endwb , &
      tws                   => grc_ws%tws   , &
      t_grnd                => col_es%t_grnd , &
      t_grnd_grc            => lnd2atm_vars%t_grnd_grc   , &
      t_soisno              =>  col_es%t_soisno, &
      t_soisno_grc          =>  lnd2atm_vars%t_soisno_grc, &
      zwt_col          =>   soilhydrology_vars%zwt_col , &
      zwt_grc          =>   lnd2atm_vars%zwt_grc,    &
      coszen_col       => surfalb_vars%coszen_col , &
      coszen_str       => lnd2atm_vars%coszen_str , &
      qflx_h2orof_drain     => col_wf%qflx_h2orof_drain , &
      qflx_h2orof_drain_grc => lnd2atm_vars%qflx_h2orof_drain_grc, &
      nh3_total        => col_nf%nh3_total &
      )
    !----------------------------------------------------
    ! lnd -> atm
    !----------------------------------------------------

    ! First, compute the "minimal" set of fields.
    call lnd2atm_minimal(bounds, &
         surfalb_vars, solarabs_vars, energyflux_vars, atm2lnd_vars, lnd2atm_vars)

    ! Map to gridcells over the active points with nonzero gridcell weight only
    call l2a_map_init(bounds, npts_p, pts_p, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,unity:urbans), scale_l2g(bounds%begl:bounds%endl))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         t_ref2m (bounds%begp:bounds%endp), &
         t_ref2m_grc(bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         q_ref2m (bounds%begp:bounds%endp), &
         q_ref2m_grc(bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         u10_elm_patch(bounds%begp:bounds%endp), &
         u_ref10m_grc (bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         u10_with_gusts_elm_patch(bounds%begp:bounds%endp), &
         u_ref10m_with_gusts_grc (bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         taux (bounds%begp:bounds%endp), &
         taux_grc (bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         tauy (bounds%begp:bounds%endp), &
         tauy_grc (bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_evap_tot (bounds%begp:bounds%endp), &
         qflx_evap_tot_grc(bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         fsa_patch(bounds%begp:bounds%endp), &
         fsa_grc (bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         fv_patch(bounds%begp:bounds%endp), &
         fv_grc (bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         ram1_patch(bounds%begp:bounds%endp), &
         ram1_grc (bounds%begg:bounds%endg))

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         eflx_sh_tot (bounds%begp:bounds%endp), &
         eflx_sh_tot_grc (bounds%begg:bounds%endg))

    do g = bounds%begg, bounds%endg
       eflx_sh_tot_grc(g) =  eflx_sh_tot_grc(g) - grc_ef%eflx_dynbal(g)
    enddo

    call l2a_p2g_1d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         eflx_lh_tot (bounds%begp:bounds%endp), &
         eflx_lh_tot_grc(bounds%begg:bounds%endg))

    if (use_cn .or. use_fates) then
       call l2a_c2g_1d(bounds, npts_c, pts_c, &
            scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
            nee (bounds%begc:bounds%endc), &
            nee_grc(bounds%begg:bounds%endg))

       if (use_lch4) then
          if (.not. ch4offline) then
             ! Adjust flux of CO2 by the net conversion of mineralizing C to CH4
             do g = bounds%begg,bounds%endg
                ! nem is in g C/m2/sec
                nee_grc(g) = nee_grc(g) + lnd2atm_vars%nem_grc(g)
             end do
          end if
       end if

       ! Convert from gC/m2/s to kgCO2/m2/s
       do g = bounds%begg,bounds%endg
          nee_grc(g) = nee_grc(g)*convertgC2kgCO2
       end do
    else
       do g = bounds%begg,bounds%endg
          nee_grc(g) = 0._r8
       end do
    end if

    ! drydepvel
    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
       call l2a_p2g_2d(bounds, npts_p, pts_p, &
            scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
            n_drydep, &
            velocity_patch(bounds%begp:bounds%endp,:), &
            ddvel_grc (bounds%begg:bounds%endg,:))
    endif

    ! voc emission flux
    if (use_voc .and. shr_megan_mechcomps_n>0) then
       call l2a_p2g_2d(bounds, npts_p, pts_p, &
            scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
            shr_megan_mechcomps_n, &
            vocemis_vars%vocflx_patch, &
            lnd2atm_vars%flxvoc_grc)
    end if

    ! dust emission flux
    call l2a_p2g_2d(bounds, npts_p, pts_p, &
         scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
         ndst, &
         flx_mss_vrt_dst_patch(bounds%begp:bounds%endp,:), &
         flxdst_grc (bounds%begg:bounds%endg,:))


    ! ch4 flux
    if (use_lch4) then
       call l2a_c2g_1d(bounds, npts_c, pts_c, &
            scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
            ch4_surf_flux_tot_col(bounds%begc:bounds%endc), &
            flux_ch4_grc (bounds%begg:bounds%endg))
    end if

    ! nh3 flux
    if (shr_fan_to_atm) then
       call l2a_c2g_1d(bounds, npts_c, pts_c, &
            scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
            nh3_total (bounds%begc:bounds%endc), &
            lnd2atm_vars%flux_nh3_grc (bounds%begg:bounds%endg))
    end if

    !----------------------------------------------------
    ! lnd -> rof
    !----------------------------------------------------

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_runoff (bounds%begc:bounds%endc), &
         qflx_rofliq_grc(bounds%begg:bounds%endg))

    do g = bounds%begg, bounds%endg
       qflx_rofliq_grc(g) = qflx_rofliq_grc(g) - grc_wf%qflx_liq_dynbal(g)
    enddo

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_surf (bounds%begc:bounds%endc), &
         qflx_rofliq_qsur_grc(bounds%begg:bounds%endg))

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_h2osfc_surf (bounds%begc:bounds%endc), &
         qflx_rofliq_qsurp_grc(bounds%begg:bounds%endg))

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_irr_demand (bounds%begc:bounds%endc), &
         qflx_irr_demand_grc(bounds%begg:bounds%endg))

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_drain (bounds%begc:bounds%endc), &
         qflx_rofliq_qsub_grc(bounds%begg:bounds%endg))

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_drain_perched (bounds%begc:bounds%endc), &
         qflx_rofliq_qsubp_grc(bounds%begg:bounds%endg))

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_qrgwl (bounds%begc:bounds%endc), &
         qflx_rofliq_qgwl_grc(bounds%begg:bounds%endg))
    
    qflx_snwcp_ice(bounds%begc:bounds%endc) = qflx_snwcp_ice(bounds%begc:bounds%endc) + qflx_ice_runoff_xs(bounds%begc:bounds%endc) 
    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         qflx_snwcp_ice (bounds%begc:bounds%endc), &
         qflx_rofice_grc(bounds%begg:bounds%endg))

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         coszen_col (bounds%begc:bounds%endc), &
         coszen_str (bounds%begg:bounds%endg))

    do g = bounds%begg, bounds%endg
       qflx_rofice_grc(g) = qflx_rofice_grc(g) - grc_wf%qflx_ice_dynbal(g)
    enddo

    ! land river two-way coupling
    ! Average up  to gridcell for the inundation drainage
    if (use_lnd_rof_two_way) then
          call l2a_c2g_1d(bounds, npts_c, pts_c, &
               scale_c2l(bounds%begc:bounds%endc,unity), scale_l2g(bounds%begl:bounds%endl), &
               qflx_h2orof_drain(bounds%begc:bounds%endc), &
               qflx_h2orof_drain_grc(bounds%begg:bounds%endg))
    endif

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         col_ws%wslake_col(bounds%begc:bounds%endc), &
         lnd2atm_vars%wslake_grc(bounds%begg:bounds%endg))

    ! calculate total water storage for history files
    ! first set tws to gridcell total endwb
    ! second add river storage as gridcell average depth (1.e-3 converts [m3/km2] to [mm])
    ! TODO - this was in BalanceCheckMod - not sure where it belongs?

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         endwb(bounds%begc:bounds%endc), &
         tws (bounds%begg:bounds%endg))
    do g = bounds%begg, bounds%endg
       tws(g) = tws(g) + atm2lnd_vars%volr_grc(g) / grc_pp%area(g) * 1.e-3_r8
    enddo


    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbans), scale_l2g(bounds%begl:bounds%endl), &
         t_grnd (bounds%begc:bounds%endc), &
         t_grnd_grc(bounds%begg:bounds%endg))

    !do lvl = -nlevsno+1, nlevgrnd
        call l2a_c2g_2d(bounds, npts_c, pts_c, &
             scale_c2l(bounds%begc:bounds%endc,urbans), scale_l2g(bounds%begl:bounds%endl), &
             nlevgrnd+nlevsno, &
             t_soisno (bounds%begc:bounds%endc,:), &
             t_soisno_grc(bounds%begg:bounds%endg,:))
    !enddo

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbans), scale_l2g(bounds%begl:bounds%endl), &
         zwt_col(bounds%begc:bounds%endc), &
         zwt_grc(bounds%begg:bounds%endg))

    do g = bounds%begg,bounds%endg
       ! TODO temperary treatment in case weird values after c2g
       if(lnd2atm_vars%t_soisno_grc(g, 1) > 400._r8) then
             write(iulog,*)'lnd2atm_vars%t_soisno_grc(g, 1) is',lnd2atm_vars%t_soisno_grc(g, 1)
             call endrun( msg=' lnd2atm ERROR: lnd2atm_vars%t_soisno_grc >  400 Kelvin degree.'//errMsg(__FILE__, __LINE__))
       end if
       lnd2atm_vars%Tqsur_grc(g) = avg_tsoil_surf(t_soisno_grc(g,:))
       lnd2atm_vars%Tqsub_grc(g) = avg_tsoil(zwt_grc(g),t_soisno_grc(g,:))

    end do

    call l2a_c2g_1d(bounds, npts_c, pts_c, &
         scale_c2l(bounds%begc:bounds%endc,urbanf), scale_l2g(bounds%begl:bounds%endl), &
         sedflux_vars%sed_yld_col(bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofmud_grc(bounds%begg:bounds%endg))
    
    end associate
  end subroutine lnd2atm


    function avg_tsoil_surf(Tsoil_) result(avgT_)
      !$acc routine seq
    ! Function for estimating average soil temperature within the top few layers (which closely interacts with surface runoff)
        implicit none
        real(r8), intent(in) :: Tsoil_(-nlevsno+1:nlevgrnd)       ! water table depth, soil temperature
        real(r8) :: avgT_             ! average soil temperature within the saturated layers

        integer :: ilvl   !local index
        real(r8) :: depth_(nlevsoi), h_layer(nlevsoi), sum_h, sum_ht

        ! calculate the thickness of each 15 soil layer, refer to Eqn. (6.5) and (6.6) in CLM4.0 tech note
        do ilvl = 1, nlevsoi
            depth_(ilvl) = 0.025_r8*(EXP(0.5_r8*(REAL(ilvl)-0.5_r8))-1._r8)
        enddo
        h_layer(1) = 0.5_r8*(depth_(1)+depth_(2))
        do ilvl = 2, nlevsoi-1
            h_layer(ilvl) = 0.5_r8*(depth_(ilvl+1)-depth_(ilvl-1))
        end do
        h_layer(nlevsoi) = depth_(nlevsoi) - depth_(nlevsoi-1)

        sum_h = 0._r8
        sum_ht = 0._r8
        do ilvl = 1, 3
            sum_h = sum_h + h_layer(ilvl)
            sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
        enddo
        avgT_ = sum_ht/sum_h

        return
    end function avg_tsoil_surf

    function avg_tsoil(zwt_, Tsoil_) result(avgT_)
      !$acc routine seq
    ! Function for estimating average soil temperature within the saturated soil zone (which produces subsurface runoff)
        implicit none
        real(r8), intent(in) :: zwt_, Tsoil_(-nlevsno+1:nlevgrnd)       ! water table depth, soil temperature
        real(r8) :: avgT_             ! average soil temperature within the saturated layers

        integer :: ilvl,izwt   !local index
        real(r8) :: depth_(nlevsoi), h_layer(nlevsoi), sum_h, sum_ht

        ! calculate the thickness of each 15 soil layer, refer to Eqn. (6.5) and (6.6) in CLM4.0 tech note
        do ilvl = 1, nlevsoi
            depth_(ilvl) = 0.025_r8*(EXP(0.5_r8*(REAL(ilvl)-0.5_r8))-1._r8)
        enddo
        h_layer(1) = 0.5_r8*(depth_(1)+depth_(2))
        do ilvl = 2, nlevsoi-1
            h_layer(ilvl) = 0.5_r8*(depth_(ilvl+1)-depth_(ilvl-1))
        end do
        h_layer(nlevsoi) = depth_(nlevsoi) - depth_(nlevsoi-1)

        if(zwt_ <= 0._r8) then ! water table close to ground surface, average over the whole soil column
            sum_h = 0._r8
            sum_ht = 0._r8
            do ilvl = 1, nlevsoi
                sum_h = sum_h + h_layer(ilvl)
                sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
            enddo
            avgT_ = sum_ht/sum_h
        else if(zwt_ >= depth_(nlevsoi)) then ! water table deeper than the total soil depth, taking the temperature of the deepes
            avgT_ = Tsoil_(nlevsoi)
        else
            sum_h = 0._r8
            sum_ht = 0._r8
            ! find out which soil layer the water table is located
            do ilvl = 1, nlevsoi
                if(zwt_ <= depth_(ilvl)) then
                    izwt = ilvl
                    sum_h = depth_(ilvl) - zwt_
                    sum_ht = (depth_(ilvl) - zwt_)*Tsoil_(ilvl)
                    exit
                end if
            enddo
            ! calculate mean soil temperature of the total saturated soil zone
            do ilvl = izwt + 1, nlevsoi
                sum_h = sum_h + h_layer(ilvl)
                sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
            enddo
            avgT_ = sum_ht/sum_h
        end if

        return
    end function avg_tsoil

  !------------------------------------------------------------------------
  subroutine l2a_map_init(bounds, npts_p, pts_p, npts_c, pts_c, scale_c2l, scale_l2g)
    !
    ! !DESCRIPTION:
    ! Set up the active patches and columns with nonzero gridcell weight, in
    ! index order, and the subgrid averaging scale factors (c2l: unity, urbanf
    ! and urbans; l2g: unity) used by lnd2atm_minimal and lnd2atm, so that
    ! l2a_p2g_1d/2d and l2a_c2g_1d/2d map fields to gridcells over these points
    ! only. p2g and c2g (used before) loop over all allocated points, skip the
    ! same points and set up the same scale factors on every call, so the
    ! results are the same.
    !
    ! !USES:
      !$acc routine seq
    use subgridAveMod , only : create_scale_c2l, build_scale_l2g_gpu
    use VegetationType, only : veg_pp
    use ColumnType    , only : col_pp
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds
    integer          , intent(out)   :: npts_p                        ! number of entries in pts_p
    integer          , intent(out)   :: pts_p(:)                      ! active patches with nonzero wtgcell
    integer          , intent(out)   :: npts_c                        ! number of entries in pts_c
    integer          , intent(out)   :: pts_c(:)                      ! active columns with nonzero wtgcell
    real(r8)         , intent(inout) :: scale_c2l(bounds%begc:,unity:) ! c2l scale factors [col, scale type]
    real(r8)         , intent(inout) :: scale_l2g(bounds%begl:)        ! l2g scale factors for unity [landunit]
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, ic
    !------------------------------------------------------------------------

    npts_p = 0
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .and. veg_pp%wtgcell(p) /= 0._r8) then
          npts_p = npts_p + 1
          pts_p(npts_p) = p
       end if
    end do

    npts_c = 0
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .and. col_pp%wtgcell(c) /= 0._r8) then
          npts_c = npts_c + 1
          pts_c(npts_c) = c
       end if
    end do

    do ic = unity, urbans
       call create_scale_c2l(bounds, ic, scale_c2l(bounds%begc:bounds%endc,ic))
    end do

    call build_scale_l2g_gpu(bounds, unity, scale_l2g(bounds%begl:bounds%endl))

  end subroutine l2a_map_init

  !------------------------------------------------------------------------
  subroutine l2a_p2g_1d(bounds, npts_p, pts_p, scale_c2l, scale_l2g, parr, garr)
    !
    ! !DESCRIPTION:
    ! Same as p2g_1d_gpu with p2c_scale_type=unity and l2g_scale_type=unity,
    ! over the patches pts_p and with the scale factors set up in l2a_map_init
    !
    ! !USES:
      !$acc routine seq
    use VegetationType, only : veg_pp
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)  :: bounds
    integer , intent(in)  :: npts_p                ! number of entries in pts_p
    integer , intent(in)  :: pts_p(:)              ! active patches with nonzero wtgcell
    real(r8), intent(in)  :: scale_c2l(bounds%begc:) ! c2l scale factors of the averaging type
    real(r8), intent(in)  :: scale_l2g(bounds%begl:) ! l2g scale factors
    real(r8), intent(in)  :: parr( bounds%begp: )  ! input pft array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,p,c,l,g,index                   ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_p2c                          ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------

    scale_p2c = 1.0_r8

    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do fp = 1,npts_p
       p = pts_p(fp)
       c = veg_pp%column(p)
       l = veg_pp%landunit(p)
       if (parr(p) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
          g = veg_pp%gridcell(p)
          if (sumwt(g) == 0._r8) garr(g) = 0._r8
          garr(g) = garr(g) + parr(p) * scale_p2c * scale_c2l(c) * scale_l2g(l) * veg_pp%wtgcell(p)
          sumwt(g) = sumwt(g) + veg_pp%wtgcell(p)
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       ! as in p2g_1d_gpu, report and continue
       print *,'p2g_1d error: sumwt is greater than 1.0 at g= ',index
    end if

  end subroutine l2a_p2g_1d

  !------------------------------------------------------------------------
  subroutine l2a_p2g_2d(bounds, npts_p, pts_p, scale_c2l, scale_l2g, num2d, parr, garr)
    !
    ! !DESCRIPTION:
    ! Same as p2g_2d_gpu with p2c_scale_type=unity and l2g_scale_type=unity,
    ! over the patches pts_p and with the scale factors set up in l2a_map_init
    !
    ! !USES:
      !$acc routine seq
    use VegetationType, only : veg_pp
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)  :: bounds
    integer , intent(in)  :: npts_p                    ! number of entries in pts_p
    integer , intent(in)  :: pts_p(:)                  ! active patches with nonzero wtgcell
    real(r8), intent(in)  :: scale_c2l(bounds%begc:)   ! c2l scale factors of the averaging type
    real(r8), intent(in)  :: scale_l2g(bounds%begl:)   ! l2g scale factors
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: parr( bounds%begp: , 1: ) ! input pft array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    !
    ! !LOCAL VARIABLES:
    integer  :: j,fp,p,c,l,g,index                     ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_p2c                              ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------

    scale_p2c = 1.0_r8

    garr(bounds%begg : bounds%endg, :) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do fp = 1,npts_p
          p = pts_p(fp)
          c = veg_pp%column(p)
          l = veg_pp%landunit(p)
          if (parr(p,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
             g = veg_pp%gridcell(p)
             if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
             garr(g,j) = garr(g,j) + parr(p,j) * scale_p2c * scale_c2l(c) * scale_l2g(l) * veg_pp%wtgcell(p)
             sumwt(g) = sumwt(g) + veg_pp%wtgcell(p)
          end if
       end do
       found = .false.
       do g = bounds%begg, bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          ! as in p2g_2d_gpu, report and stop
          print *,'p2g_2d error: sumwt gt 1.0 at g/sumwt = ',index,sumwt(index)
          stop
       end if
    end do

  end subroutine l2a_p2g_2d

  !------------------------------------------------------------------------
  subroutine l2a_c2g_1d(bounds, npts_c, pts_c, scale_c2l, scale_l2g, carr, garr)
    !
    ! !DESCRIPTION:
    ! Same as c2g_1d_gpu (and c2g_1d) with l2g_scale_type=unity, over the
    ! columns pts_c and with the scale factors set up in l2a_map_init
    !
    ! !USES:
      !$acc routine seq
    use ColumnType  , only : col_pp
    use elm_varcon  , only : nameg
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)  :: bounds
    integer , intent(in)  :: npts_c                ! number of entries in pts_c
    integer , intent(in)  :: pts_c(:)              ! active columns with nonzero wtgcell
    real(r8), intent(in)  :: scale_c2l(bounds%begc:) ! c2l scale factors of the averaging type
    real(r8), intent(in)  :: scale_l2g(bounds%begl:) ! l2g scale factors
    real(r8), intent(in)  :: carr( bounds%begc: )  ! input column array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,l,g,index                     ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------

    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do fc = 1,npts_c
       c = pts_c(fc)
       l = col_pp%landunit(c)
       if (carr(c) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
          g = col_pp%gridcell(c)
          if (sumwt(g) == 0._r8) garr(g) = 0._r8
          garr(g) = garr(g) + carr(c) * scale_c2l(c) * scale_l2g(l) * col_pp%wtgcell(c)
          sumwt(g) = sumwt(g) + col_pp%wtgcell(c)
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
#ifndef _OPENACC
       write(iulog,*) 'c2g_1d error: sumwt is greater than 1.0 at g= ',index
       call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
#endif
    end if

  end subroutine l2a_c2g_1d

  !------------------------------------------------------------------------
  subroutine l2a_c2g_2d(bounds, npts_c, pts_c, scale_c2l, scale_l2g, num2d, carr, garr)
    !
    ! !DESCRIPTION:
    ! Same as c2g_2d_gpu with l2g_scale_type=unity, over the columns pts_c
    ! and with the scale factors set up in l2a_map_init
    !
    ! !USES:
      !$acc routine seq
    use ColumnType  , only : col_pp
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)  :: bounds
    integer , intent(in)  :: npts_c                    ! number of entries in pts_c
    integer , intent(in)  :: pts_c(:)                  ! active columns with nonzero wtgcell
    real(r8), intent(in)  :: scale_c2l(bounds%begc:)   ! c2l scale factors of the averaging type
    real(r8), intent(in)  :: scale_l2g(bounds%begl:)   ! l2g scale factors
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: carr( bounds%begc: , 1: ) ! input column array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    !
    ! !LOCAL VARIABLES:
    integer  :: j,fc,c,l,g,index                       ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------

    garr(bounds%begg : bounds%endg,:) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do fc = 1,npts_c
          c = pts_c(fc)
          l = col_pp%landunit(c)
          if (carr(c,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
             g = col_pp%gridcell(c)
             if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
             garr(g,j) = garr(g,j) + carr(c,j) * scale_c2l(c) * scale_l2g(l) * col_pp%wtgcell(c)
             sumwt(g) = sumwt(g) + col_pp%wtgcell(c)
          end if
       end do
       found = .false.
       do g = bounds%begg, bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          ! as in c2g_2d_gpu, report and continue
          print *, 'c2g_2d error: sumwt is greater than 1.0 at g= ',index
       end if
    end do

  end subroutine l2a_c2g_2d

end module lnd2atmMod
