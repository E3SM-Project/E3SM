
module mo_chemini

  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmd_utils,   only : masterproc
  use cam_logfile,  only : iulog

  implicit none

  private
  public :: chemini

contains

  subroutine chemini &
       ( solar_parms_file &
       , euvac_file &
       , euvacdat_file &
       , photon_file &
       , electron_file &
       , airpl_emis_file &
       , sulf_file &
       , sad_file &
       , sad_timing &
       , depvel_file &
       , depvel_lnd_file &
       , clim_soilw_file &
       , season_wes_file &
       , xs_coef_file &
       , xs_short_file &
       , xs_long_file &
       , rsf_file &
       , fstrat_file &
       , fstrat_list &
       , srf_emis_specifier &
       , srf_emis_type &
       , srf_emis_cycle_yr &
       , srf_emis_fixed_ymd &
       , srf_emis_fixed_tod &
       , ext_frc_specifier &
       , ext_frc_type &
       , ext_frc_cycle_yr &
       , ext_frc_fixed_ymd &
       , ext_frc_fixed_tod &
       , xactive_prates &
       , exo_coldens_file &
       , tuv_xsect_file &
       , o2_xsect_file &
       , lght_no_prd_factor &
       , chem_name &
       , pbuf2d &
       )

    !-----------------------------------------------------------------------
    ! 	... Chemistry module intialization
    !-----------------------------------------------------------------------

    use mo_airplane,       only : airpl_src
    use mo_srf_emissions,  only : srf_emissions_inti
    use mo_sulf,           only : sulf_inti
    use mo_photo,          only : photo_inti
    use mo_lightning,      only : lightning_inti
    use mo_drydep,         only : drydep_inti
    use seq_drydep_mod,    only : DD_XLND, DD_XATM, drydep_method
    use mo_imp_sol,        only : imp_slv_inti
    use mo_exp_sol,        only : exp_sol_inti
    use spmd_utils,        only : iam
    use mo_fstrat,         only : fstrat_inti
    use m_types,           only : time_ramp
    use abortutils,        only : endrun
    use pmgrid,            only : plev           
    use mo_sethet,         only : sethet_inti
    use mo_usrrxt,         only : usrrxt_inti
    use mo_extfrc,         only : extfrc_inti
    use mo_setext,         only : setext_inti
    use mo_setinv,         only : setinv_inti
    use mo_gas_phase_chemdr,only: gas_phase_chemdr_inti
    
    use tracer_cnst,       only : tracer_cnst_init
    use tracer_srcs,       only : tracer_srcs_init
    use mo_synoz,          only : synoz_inti
    use mo_chem_utls,      only : get_spc_ndx
    use mo_airglow,        only : init_airglow
    use mo_mean_mass,      only : init_mean_mass
    use mo_mass_xforms,    only : init_mass_xforms
    use mo_strato_rates,   only : init_strato_rates
    use mo_cph,            only : init_cph
    use mo_strato_sad,     only : strato_sad_inti
    use mo_sad,            only : sad_inti
    use mo_solarproton,    only : spe_init
    use mo_solar_parms,    only : solar_parms_init, get_solar_parms
    use euvac,             only : euvac_init, euvac_set_etf
    use mo_heatnirco2,     only : heatnirco2_init
    use mo_waccm_hrates,   only : init_hrates
    use mo_aurora,         only : aurora_inti
    use clybry_fam,        only : clybry_fam_init
    use mo_neu_wetdep,     only : neu_wetdep_init 
    use physics_buffer,    only : physics_buffer_desc

    implicit none

    character(len=*), intent(in) :: solar_parms_file
    character(len=*), intent(in) :: euvac_file
    character(len=*), intent(in) :: euvacdat_file
    character(len=*), intent(in) :: photon_file
    character(len=*), intent(in) :: electron_file

    character(len=*), intent(in) :: airpl_emis_file
    character(len=*), intent(in) :: sulf_file
    character(len=*), intent(in) :: sad_file
    type(time_ramp),  intent(in) :: sad_timing 
    character(len=*), intent(in) :: depvel_file
    character(len=*), intent(in) :: depvel_lnd_file
    character(len=*), intent(in) :: clim_soilw_file
    character(len=*), intent(in) :: season_wes_file
    character(len=*), intent(in) :: xs_coef_file
    character(len=*), intent(in) :: xs_short_file
    character(len=*), intent(in) :: xs_long_file
    character(len=*), intent(in) :: rsf_file
    character(len=*), intent(in) :: fstrat_file
    character(len=*), intent(in) :: fstrat_list(:)
    character(len=*), dimension(:), intent(in) :: srf_emis_specifier
    character(len=*), dimension(:), intent(in) :: ext_frc_specifier
    logical, intent(in)          :: xactive_prates
    character(len=*), intent(in) :: exo_coldens_file
    character(len=*), intent(in) :: tuv_xsect_file
    character(len=*), intent(in) :: o2_xsect_file
    real(r8),         intent(in) :: lght_no_prd_factor
    character(len=*), intent(in) :: ext_frc_type
    integer,          intent(in) :: ext_frc_cycle_yr
    integer,          intent(in) :: ext_frc_fixed_ymd
    integer,          intent(in) :: ext_frc_fixed_tod
    character(len=*), intent(in) :: srf_emis_type
    integer,          intent(in) :: srf_emis_cycle_yr
    integer,          intent(in) :: srf_emis_fixed_ymd
    integer,          intent(in) :: srf_emis_fixed_tod

    character(len=*), intent(in) :: chem_name

    integer :: ndx
    logical :: is_waccm
    real(r8)          ::   f107
    real(r8)          ::   f107a
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    call gas_phase_chemdr_inti()

    call init_mean_mass
    call init_mass_xforms

    call setinv_inti()
    call sethet_inti()
    call usrrxt_inti()
    call init_hrates
    call init_airglow

    call init_strato_rates
    call init_cph

    !-----------------------------------------------------------------------
    ! 	... initialize tracer modules
    !-----------------------------------------------------------------------
    call tracer_cnst_init()
    call tracer_srcs_init()

    !-----------------------------------------------------------------------
    ! 	... read time-independent airplane emissions
    !-----------------------------------------------------------------------
    call airpl_src(airpl_emis_file)
    if (masterproc) write(iulog,*) 'chemini: after airpl_src on node ',iam

    !-----------------------------------------------------------------------
    ! 	... read time-dependent surface flux dataset
    !-----------------------------------------------------------------------
    call srf_emissions_inti ( srf_emis_specifier, srf_emis_type, srf_emis_cycle_yr, srf_emis_fixed_ymd, srf_emis_fixed_tod)

    if (masterproc) write(iulog,*) 'chemini: after srf_emissions_inti on node ',iam

    !-----------------------------------------------------------------------
    ! 	... initialize external forcings module
    !-----------------------------------------------------------------------
    call setext_inti()
    call extfrc_inti(ext_frc_specifier, ext_frc_type, ext_frc_cycle_yr, ext_frc_fixed_ymd, ext_frc_fixed_tod)
    if (masterproc) write(iulog,*) 'chemini: after extfrc_inti on node ',iam

    !-----------------------------------------------------------------------
    ! 	... read the stratospheric sad dataset
    !-----------------------------------------------------------------------
    call strato_sad_inti(sad_file, sad_timing)
    if (masterproc) write(iulog,*) 'chemini: after strato_sad_inti on node ',iam


    call sulf_inti(sulf_file)
    if (masterproc) write(iulog,*) 'chemini: after sulf_inti on node ',iam

    !-----------------------------------------------------------------------
    !	... initialize the sad module
    !-----------------------------------------------------------------------
    call sad_inti(pbuf2d)
    if (masterproc) write(iulog,*) 'chemini: after sad_inti on node ',iam

    !-----------------------------------------------------------------------
    !	... initialize the lightning module
    !-----------------------------------------------------------------------
    call lightning_inti(lght_no_prd_factor)
    if (masterproc) write(iulog,*) 'chemini: after lightning_inti on node ',iam

    !-----------------------------------------------------------------------
    !	... initialize the dry deposition module
    !-----------------------------------------------------------------------
    if ( drydep_method == DD_XATM .or. drydep_method == DD_XLND ) then
       call drydep_inti(depvel_lnd_file, clim_soilw_file, season_wes_file )
    else
       call drydep_inti( depvel_file )
    endif

    if (masterproc) write(iulog,*) 'chemini: after drydep_inti on node ',iam

    !-----------------------------------------------------------------------
    !	... Initialize the upper boundary module
    !-----------------------------------------------------------------------
    call fstrat_inti( fstrat_file, fstrat_list )
    if (masterproc) write(iulog,*) 'chemini: after fstrat_inti on node ',iam

    !-----------------------------------------------------------------------
    ! 	... initialize the co2 nir heating module
    !-----------------------------------------------------------------------
    call heatnirco2_init

    !-----------------------------------------------------------------------
    ! 	... initialize photorate module
    !-----------------------------------------------------------------------

    is_waccm = ((trim(chem_name) == 'waccm_mozart') .or. &
      (trim(chem_name) == 'waccm_mozart_mam3')) 
    
    if (len_trim(solar_parms_file)>0) then
       call solar_parms_init (solar_parms_file)
    endif

    if ((len_trim(solar_parms_file)>0) .and. (.not.trim(chem_name) == 'waccm_ghg') ) then       
       !-----------------------------------------------------------------------
       ! 	... initialize the solar parameters module
       !-----------------------------------------------------------------------
       call get_solar_parms( f107_s = f107, f107a_s = f107a )
       if (masterproc) write(iulog,*) 'chemini: f107,f107a = ',f107,f107a

    endif
    
    if (is_waccm) then
       !-----------------------------------------------------------------------
       ! 	... initialize the euvac etf module
       !-----------------------------------------------------------------------
       call euvac_init (euvac_file)
       call euvac_set_etf( f107, f107a )
       !-----------------------------------------------------------------------
       ! 	... initialize the solar proton event module
       !-----------------------------------------------------------------------
       call spe_init()

    endif

    call photo_inti( xs_coef_file, xs_short_file, xs_long_file, rsf_file, &
         euvacdat_file, photon_file, electron_file, &
         exo_coldens_file, tuv_xsect_file, o2_xsect_file, xactive_prates, is_waccm=is_waccm )

    if (masterproc) write(iulog,*) 'chemini: after waccm_prate_inti on node ',iam

    !-----------------------------------------------------------------------
    !	... initialize the implicit solver
    !-----------------------------------------------------------------------
    call imp_slv_inti()
    call exp_sol_inti()

    !-----------------------------------------------------------------------
    !       ... initialize the stratospheric ozone source
    !-----------------------------------------------------------------------
    if( get_spc_ndx( 'SYNOZ' ) > 0 ) then
       call synoz_inti( )
       ! over ride the ozone constituent used for radiation feedbacks
    end if

    !-----------------------------------------------------------------------
    !	... initialize ion production
    !-----------------------------------------------------------------------
    call aurora_inti
    if (masterproc) write(iulog,*) 'chemini: after aurora_inti'

    call neu_wetdep_init()
    if (masterproc) write(iulog,*) 'chemini: after wetdep_init'

    call clybry_fam_init()

    if (masterproc) write(iulog,*) 'chemini: finished on node ',iam

  end subroutine chemini

end module mo_chemini
