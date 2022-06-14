module mo_gas_phase_chemdr

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use shr_const_mod,    only : pi => shr_const_pi
  use constituents,     only : pcnst
  use cam_history,      only : fieldname_len
  use chem_mods,        only : phtcnt, rxntot, gas_pcnst
  use chem_mods,        only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map, extcnt
  use dust_model,       only : dust_names, ndust => dust_nbin
  use ppgrid,           only : pcols, pver
  use spmd_utils,       only : iam, masterproc
  use phys_control,     only : phys_getopts
  use cam_logfile,      only : iulog
  use physics_buffer,   only : pbuf_get_index

  implicit none
  save

  private
  public :: gas_phase_chemdr, gas_phase_chemdr_inti 
  public :: map2chm
  public :: gas_ac_name, gas_ac_name_2D

  integer :: map2chm(pcnst) = 0           ! index map to/from chemistry/constituents list
  character(len=fieldname_len) :: gas_ac_name(gas_pcnst)
  character(len=fieldname_len) :: gas_ac_name_2D(gas_pcnst)

  integer :: o3_ndx, synoz_ndx, so4_ndx, h2o_ndx, o2_ndx, o_ndx, hno3_ndx, dst_ndx, cldice_ndx, e90_ndx
  !integer :: o3lnz_ndx, n2olnz_ndx, noylnz_ndx, ch4lnz_ndx
  integer :: o3lnz_ndx, ch4lnz_ndx
  integer :: uci1_ndx
  integer :: het1_ndx
  integer :: ndx_cldfr, ndx_cmfdqr, ndx_nevapr, ndx_cldtop, ndx_prain, ndx_sadsulf
  integer :: ndx_h2so4
  integer :: inv_ndx_cnst_o3, inv_ndx_m
! for ozone budget 
  integer :: jo2_b_ndx, lch3o2_no_ndx, lno_ho2_ndx, lc2h5o2_no_ndx, lch3co3_no_ndx, lroho2_no_ndx, lisopo2_no_ndx, lmvko2_no_ndx
  integer :: jo1dU_ndx, jno2_ndx, jno3_a_ndx, jn2o5_b_ndx, po3_oh_ndx
  integer :: uci2_ndx, uci3_ndx, lc2h4_o3_ndx, lisop_o3_ndx, lo3_oh_ndx, lo3_ho2_ndx, lmvkmacr_o3_ndx  
  integer :: ch3o2_ndx, no_ndx, ho2_ndx, c2h5o2_ndx, ch3co3_ndx, roho2_ndx, isopo2_ndx, mvko2_ndx
  integer :: h2_ndx, ch4_ndx, c2h4_ndx, isop_ndx, oh_ndx, mvkmacr_ndx
  integer :: no2_ndx, no3_ndx, n2o5_ndx
  integer :: lo3_no_ndx, lo3_no2_ndx 

  character(len=fieldname_len),dimension(rxntot-phtcnt) :: rxn_names
  character(len=fieldname_len),dimension(phtcnt)        :: pht_names
  character(len=fieldname_len),dimension(rxt_tag_cnt)   :: tag_names
  character(len=fieldname_len),dimension(extcnt)        :: extfrc_name
  logical :: convproc_do_aer 

  real(r8), parameter :: low_limit = 0.1_r8  ! lower limit factor for dry deposition
  integer,  parameter :: srf_emit_nlayer = 3 ! number of layers adding surface emission

contains

  subroutine gas_phase_chemdr_inti(chem_name)

    use mo_chem_utls,      only : get_spc_ndx, get_extfrc_ndx, get_inv_ndx, get_rxt_ndx
    use cam_history,       only : addfld, horiz_only
    use cam_abortutils,        only : endrun
    use cam_history,       only : add_default
    use mo_chm_diags,      only : chm_diags_inti
    use mo_tracname,       only : solsym
    use constituents,      only : cnst_get_ind
    use rate_diags,        only : rate_diags_init
    use rad_constituents,  only : rad_cnst_get_info

    implicit none

    character(len=*), intent(in) :: chem_name

    character(len=3)  :: string
    integer           :: n, m
    logical           :: history_aerosol      ! Output the MAM aerosol tendencies

    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out = history_aerosol, &
         convproc_do_aer_out = convproc_do_aer ) 
   
    ndx_h2so4 = get_spc_ndx('H2SO4')

    uci1_ndx= get_rxt_ndx('uci1')
    het1_ndx= get_rxt_ndx('het1')
    o3_ndx  = get_spc_ndx('O3')
    o_ndx   = get_spc_ndx('O')
    o2_ndx  = get_spc_ndx('O2')
    so4_ndx = get_spc_ndx('SO4')
    h2o_ndx = get_spc_ndx('H2O')
    hno3_ndx = get_spc_ndx('HNO3')
    dst_ndx = get_spc_ndx( dust_names(1) )
    synoz_ndx = get_extfrc_ndx( 'SYNOZ' )
    e90_ndx  = get_spc_ndx('E90')
    o3lnz_ndx =   get_spc_ndx('O3LNZ')
!    n2olnz_ndx =  get_spc_ndx('N2OLNZ')
!    noylnz_ndx =  get_spc_ndx('NOYLNZ')
    ch4lnz_ndx =  get_spc_ndx('CH4LNZ')
    ! for ozone budget 
    jo2_b_ndx      = get_rxt_ndx('jo2_b')
    lch3o2_no_ndx  = get_rxt_ndx('lch3o2_no')
    lno_ho2_ndx    = get_rxt_ndx('lno_ho2')
    lc2h5o2_no_ndx = get_rxt_ndx('lc2h5o2_no')
    lch3co3_no_ndx = get_rxt_ndx('lch3co3_no')
    lroho2_no_ndx  = get_rxt_ndx('lroho2_no')
    lisopo2_no_ndx = get_rxt_ndx('lisopo2_no')
    lmvko2_no_ndx  = get_rxt_ndx('lmvko2_no')
    ch3o2_ndx      = get_spc_ndx('CH3O2')
    no_ndx         = get_spc_ndx('NO')
    ho2_ndx        = get_spc_ndx('HO2')
    c2h5o2_ndx     = get_spc_ndx('C2H5O2')
    ch3co3_ndx     = get_spc_ndx('CH3CO3')
    roho2_ndx      = get_spc_ndx('ROHO2')
    isopo2_ndx     = get_spc_ndx('ISOPO2')
    mvko2_ndx      = get_spc_ndx('MVKO2')
    no2_ndx        = get_spc_ndx('NO2')
    no3_ndx        = get_spc_ndx('NO3')
    n2o5_ndx       = get_spc_ndx('N2O5')
 
    uci2_ndx       = get_rxt_ndx('uci2')
    uci3_ndx       = get_rxt_ndx('uci3')
    lc2h4_o3_ndx   = get_rxt_ndx('lc2h4_o3')
    lisop_o3_ndx   = get_rxt_ndx('lisop_o3')
    lo3_oh_ndx     = get_rxt_ndx('lo3_oh')
    lo3_ho2_ndx    = get_rxt_ndx('lo3_ho2')
    lmvkmacr_o3_ndx = get_rxt_ndx('lmvkmacr_o3')
    h2_ndx         = get_spc_ndx('H2')
    ch4_ndx        = get_spc_ndx('CH4')
    c2h4_ndx       = get_spc_ndx('C2H4')
    isop_ndx       = get_spc_ndx('ISOP')
    oh_ndx         = get_spc_ndx('OH')
    mvkmacr_ndx    = get_spc_ndx('MVKMACR')
    lo3_no_ndx     = get_rxt_ndx('lo3_no')
    lo3_no2_ndx    = get_rxt_ndx('lo3_no2')
    jo1dU_ndx      = get_rxt_ndx('jo1dU')
    jno2_ndx       = get_rxt_ndx('jno2')
    jno3_a_ndx     = get_rxt_ndx('jno3_a')
    jn2o5_b_ndx    = get_rxt_ndx('jn2o5_b')
    po3_oh_ndx     = get_rxt_ndx('po3_oh')

    call cnst_get_ind( 'CLDICE', cldice_ndx )

    call addfld( 'Photolysis_CLOUD',       (/ 'lev' /), 'I', 'fraction', 'Cloud fraction in photolysis')
    call addfld( 'Photolysis_CLOUD_water', (/ 'lev' /), 'I', 'g/g',      'Cloud liquid water in photolysis')
    call addfld( 'Photolysis_CLOUD_ice',   (/ 'lev' /), 'I', 'g/g',      'Cloud ice in photolysis')
    
    do m = 1,extcnt
       WRITE(UNIT=string, FMT='(I2.2)') m
       extfrc_name(m) = 'extfrc_'// trim(string)
       call addfld( extfrc_name(m), (/ 'lev' /), 'I', ' ', 'ext frcing' )
       !call add_default( extfrc_name(m), 3, ' ' )
    end do

    do n = 1,rxt_tag_cnt
       tag_names(n) = trim(rxt_tag_lst(n))
       if (n<=phtcnt) then
          call addfld(      tag_names(n),              (/ 'lev' /), 'I', '/s', 'photolysis rate' )
          call addfld( trim(tag_names(n))//"_nocloud", (/ 'lev' /), 'I', '/s', 'photolysis rate (cloud free)' )
          call addfld( trim(tag_names(n))//"_aerosol", (/ 'lev' /), 'I', '/s', 'photolysis rate (including aerosols)' )
          call addfld( trim(tag_names(n))//"_lookup", (/ 'lev' /), 'I', '/s', 'photolysis rate (lookup table: diagnostic)' )
       else
          call addfld( tag_names(n),             (/ 'lev' /), 'I', '/cm3/s', 'reaction rate' )
       endif
       !call add_default( tag_names(n), 1, ' ' )
    enddo

    do n = 1,phtcnt
       WRITE(UNIT=string, FMT='(I3.3)') n
       pht_names(n) = 'J_' // trim(string)
       call addfld(      pht_names(n),              (/ 'lev' /), 'I', '/s', 'photolysis rate' )
       call addfld( trim(pht_names(n))//"_nocloud", (/ 'lev' /), 'I', '/s', 'photolysis rate (cloud free)' )
       call addfld( trim(pht_names(n))//"_aerosol", (/ 'lev' /), 'I', '/s', 'photolysis rate (including aerosols)' )
       call addfld( trim(pht_names(n))//"_lookup", (/ 'lev' /), 'I', '/s', 'photolysis rate (lookup table: diagnostic)' )
       !call add_default( pht_names(n), 3, ' ' )
    enddo

    do n = 1,rxntot-phtcnt
       WRITE(UNIT=string, FMT='(I3.3)') n
       rxn_names(n) = 'R_' // trim(string)
       call addfld( rxn_names(n), (/ 'lev' /), 'I', '/cm3/s', 'reaction rate' )
       !call add_default( rxn_names(n), 1, ' ' )
    enddo

    call addfld( 'DTCBS',  horiz_only, 'I',   ' ','photolysis diagnostic black carbon OD' )
    call addfld( 'DTOCS',  horiz_only, 'I',   ' ','photolysis diagnostic organic carbon OD' )
    call addfld( 'DTSO4',  horiz_only, 'I',   ' ','photolysis diagnostic SO4 OD' )
    call addfld( 'DTSOA',  horiz_only, 'I',   ' ','photolysis diagnostic SOA OD' )
    call addfld( 'DTANT',  horiz_only, 'I',   ' ','photolysis diagnostic NH4SO4 OD' )
    call addfld( 'DTSAL',  horiz_only, 'I',   ' ','photolysis diagnostic salt OD' )
    call addfld( 'DTDUST',  horiz_only, 'I',  ' ','photolysis diagnostic dust OD' )
    call addfld( 'DTTOTAL',  horiz_only, 'I', ' ','photolysis diagnostic total aerosol OD' )   
    call addfld( 'FRACDAY',  horiz_only, 'I', ' ','photolysis diagnostic fraction of day' )

    call addfld( 'QDSAD', (/ 'lev' /), 'I', '/s', 'water vapor sad delta' )
    call addfld( 'SAD', (/ 'lev' /), 'I', 'cm2/cm3', 'sulfate aerosol SAD' )
    call addfld( 'SAD_SULFC', (/ 'lev' /), 'I', 'cm2/cm3', 'chemical sulfate aerosol SAD' )
    call addfld( 'SAD_SAGE', (/ 'lev' /), 'I', 'cm2/cm3', 'SAGE sulfate aerosol SAD' )
    call addfld( 'SAD_LNAT', (/ 'lev' /), 'I', 'cm2/cm3', 'large-mode NAT aerosol SAD' )
    call addfld( 'SAD_ICE', (/ 'lev' /), 'I', 'cm2/cm3', 'water-ice aerosol SAD' )
    call addfld( 'RAD_SULFC', (/ 'lev' /), 'I', 'cm', 'chemical sad sulfate' )
    call addfld( 'RAD_LNAT', (/ 'lev' /), 'I', 'cm', 'large nat radius' )
    call addfld( 'RAD_ICE', (/ 'lev' /), 'I', 'cm', 'sad ice' )
    call addfld( 'SAD_TROP', (/ 'lev' /), 'I', 'cm2/cm3', 'tropospheric aerosol SAD' )
    call addfld( 'HNO3_STS', (/ 'lev' /), 'I', 'mol/mol', 'STS condensed HNO3' )
    call addfld( 'HNO3_NAT', (/ 'lev' /), 'I', 'mol/mol', 'NAT condensed HNO3' )
    call addfld( 'QDSETT', (/ 'lev' /), 'I', '/s', 'water vapor settling delta' )
    call addfld( 'QDCHEM', (/ 'lev' /), 'I', '/s', 'water vapor chemistry delta' )
    call addfld( 'HNO3_GAS', (/ 'lev' /), 'I', 'mol/mol', 'gas-phase hno3' )
    call addfld( 'H2O_GAS', (/ 'lev' /), 'I', 'mol/mol', 'gas-phase h2o' )
    if (het1_ndx>0) then
       call addfld( 'het1_total', (/ 'lev' /), 'I', '/s', 'total N2O5 + H2O het rate constant' )
    endif
    call addfld( 'SZA', horiz_only, 'I', 'degrees', 'solar zenith angle' )

    call chm_diags_inti()
    call rate_diags_init()

!-----------------------------------------------------------------------
! get pbuf indicies
!-----------------------------------------------------------------------
    ndx_cldfr  = pbuf_get_index('CLD')
    ndx_cmfdqr = pbuf_get_index('RPRDTOT')
    ndx_nevapr = pbuf_get_index('NEVAPR')
    ndx_prain  = pbuf_get_index('PRAIN')
    ndx_cldtop = pbuf_get_index('CLDTOP')
    
!-----------------------------------------------------------------------
! get fixed oxidant (troposphere) index for Linoz_MAM
!-----------------------------------------------------------------------
   
     inv_ndx_cnst_o3 = get_inv_ndx( 'cnst_O3' ) ! prescribed O3 oxidant field
     inv_ndx_m       = get_inv_ndx( 'M' )        ! airmass.  Elsewhere this variable is known as m_ndx
     
     if ( chem_name == 'linoz_mam3'.or.chem_name == 'linoz_mam4_resus'.or.chem_name == 'linoz_mam4_resus_mom' &
         .or.chem_name == 'linoz_mam4_resus_soag'.or.chem_name == 'linoz_mam4_resus_mom_soag') then
       if ( inv_ndx_cnst_o3 < 1 ) then
          call endrun('ERROR: chem_name = '//trim(chem_name)//&
          ' requies cnst_O3 fixed oxidant field. Use cnst_O3:O3 in namelist tracer_cnst_specifier')
       end if
     end if


  end subroutine gas_phase_chemdr_inti


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine gas_phase_chemdr(lchnk, ncol, imozart, q, &
                              phis, zm, zi, calday, &
                              tfld, pmid, pdel, pdeldry, pint,  &
                              cldw, troplev, &
                              ncldwtr, ufld, vfld,  &
                              delt, ps, xactive_prates, &
                              do_cloudj_photolysis, do_cloudj_clouds, do_cloudj_aerosols, &
                              do_cloudj_lookup_diag, do_cloudj_aerosol_diag, do_cloudj_nocloud_diag, &
                              fsds, ts, asdir, ocnfrac, icefrac, &
                              precc, precl, snowhland, ghg_chem, latmapback, &
                              chem_name, drydepflx, cflx, qtend, pbuf, ixcldliq, ixcldice, tropFlag, tropFlagInt)

    !-----------------------------------------------------------------------
    !     ... Chem_solver advances the volumetric mixing ratio
    !         forward one time step via a combination of explicit,
    !         ebi, hov, fully implicit, and/or rodas algorithms.
    !-----------------------------------------------------------------------

    use chem_mods,         only : nabscol, nfs, indexm
    use physconst,         only : rga
    use mo_photo,          only : set_ub_col, setcol, table_photo, xactive_photo
    use mo_exp_sol,        only : exp_sol
    use mo_imp_sol,        only : imp_sol
    use mo_setrxt,         only : setrxt
    use mo_adjrxt,         only : adjrxt
    use mo_phtadj,         only : phtadj
    use llnl_O1D_to_2OH_adj,only : O1D_to_2OH_adj
    use mo_usrrxt,         only : usrrxt
    use mo_setinv,         only : setinv
    use mo_negtrc,         only : negtrc
    use mo_sulf,           only : sulf_interp
    use mo_lightning,      only : prod_no
    use mo_setext,         only : setext
    use mo_sethet,         only : sethet
    use mo_drydep,         only : drydep, set_soilw
    use seq_drydep_mod,    only : DD_XLND, DD_XATM, DD_TABL, drydep_method
    use mo_fstrat,         only : set_fstrat_vals, set_fstrat_h2o
    use mo_flbc,           only : flbc_set
    use phys_grid,         only : get_rlat_all_p, get_rlon_all_p, get_lat_all_p, get_lon_all_p
    use mo_mean_mass,      only : set_mean_mass
    use cam_history,       only : outfld
    use wv_saturation,     only : qsat
    use constituents,      only : cnst_mw
    use mo_drydep,         only : has_drydep
    use time_manager,      only : get_ref_date
    use mo_ghg_chem,       only : ghg_chem_set_rates, ghg_chem_set_flbc
    use mo_sad,            only : sad_strat_calc
    use charge_neutrality, only : charge_balance
    use mo_strato_sad,     only : strato_sad_set
    use mo_strato_rates,   only : ratecon_sfstrat
    use mo_aero_settling,  only : strat_aer_settling
    use shr_orb_mod,       only : shr_orb_decl
    use cam_control_mod,   only : lambm0, eccen, mvelpp, obliqr
    use mo_strato_rates,   only : has_strato_chem
    use short_lived_species,only: set_short_lived_species,get_short_lived_species
    use mo_chm_diags,      only : chm_diags, het_diags, gaschmmass_diags, aer_species
    use perf_mod,          only : t_startf, t_stopf
    use mo_tracname,       only : solsym
    use physconst,         only : gravit
    use chem_mods,         only : adv_mass
    use mo_neu_wetdep,     only : do_neu_wetdep
    use physics_buffer,    only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use infnan,            only : nan, assignment(=)
    use rate_diags,        only : rate_diags_calc
    use mo_mass_xforms,    only : mmr2vmr, vmr2mmr, h2o_to_vmr, mmr2vmri
    use orbit,             only : zenith
    use UCI_cloudJ_interface, only : cloudJ_interface
    use cam_abortutils,    only : endrun
    use mo_constants, only : pi, rgrav, rearth, avogadro
!
! LINOZ
!
    use lin_strat_chem,    only : do_lin_strat_chem,linoz_v2, linoz_v3, linv2_strat_chem_solve, linv3_strat_chem_solve, lin_strat_sfcsink, fstrat_efold_decay 
    use linoz_data,        only : has_linoz_data
!
! for aqueous chemistry and aerosol growth
!
    use aero_model,        only : aero_model_gasaerexch

    implicit none

    !-----------------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------------
    integer,        intent(in)    :: lchnk                          ! chunk index
    integer,        intent(in)    :: ncol                           ! number columns in chunk
    integer,        intent(in)    :: imozart                        ! gas phase start index in q
    real(r8),       intent(in)    :: delt                           ! timestep (s)
    real(r8),       intent(in)    :: calday                         ! day of year
    real(r8),       intent(in)    :: ps(pcols)                      ! surface pressure
    real(r8),       intent(in)    :: phis(pcols)                    ! surface geopotential
    real(r8),       intent(in)    :: tfld(pcols,pver)               ! midpoint temperature (K)
    real(r8),       intent(in)    :: pmid(pcols,pver)               ! midpoint pressures (Pa)
    real(r8),       intent(in)    :: pdel(pcols,pver)               ! pressure delta about midpoints (Pa)
    real(r8),       intent(in)    :: pdeldry(pcols,pver)            ! dry pressure delta about midpoints (Pa)
    real(r8),       intent(in)    :: ufld(pcols,pver)               ! zonal velocity (m/s)
    real(r8),       intent(in)    :: vfld(pcols,pver)               ! meridional velocity (m/s)
    real(r8),       intent(in)    :: cldw(pcols,pver)               ! cloud water (kg/kg)
    real(r8),       intent(in)    :: ncldwtr(pcols,pver)            ! droplet number concentration (#/kg)
    real(r8),       intent(in)    :: zm(pcols,pver)                 ! midpoint geopotential height above the surface (m)
    real(r8),       intent(in)    :: zi(pcols,pver+1)               ! interface geopotential height above the surface (m)
    real(r8),       intent(in)    :: pint(pcols,pver+1)             ! interface pressures (Pa)
    real(r8),       intent(in)    :: q(pcols,pver,pcnst)            ! species concentrations (kg/kg)
    logical,        intent(in)    :: xactive_prates
    logical,        intent(in)    :: do_cloudj_photolysis           ! Use Cloud-J for photolysis rates.
    logical,        intent(in)    :: do_cloudj_clouds               ! Include Clouds in Cloud-J.
    logical,        intent(in)    :: do_cloudj_aerosols             ! Include Aerosols in Cloud-J photolysis rates.
    logical,        intent(in)    :: do_cloudj_lookup_diag          ! Output lookup table photolysis rates to compare to cloud-J.
    logical,        intent(in)    :: do_cloudj_aerosol_diag         ! Output Cloud-J photolysis rates that include aerosols.
    logical,        intent(in)    :: do_cloudj_nocloud_diag         ! Output Cloud-J photolysis rates with no clouds.
    real(r8),       intent(in)    :: fsds(pcols)                    ! longwave down at sfc
    real(r8),       intent(in)    :: icefrac(pcols)                 ! sea-ice areal fraction
    real(r8),       intent(in)    :: ocnfrac(pcols)                 ! ocean areal fraction
    real(r8),       intent(in)    :: asdir(pcols)                   ! albedo: shortwave, direct
    real(r8),       intent(in)    :: ts(pcols)                      ! sfc temp (merged w/ocean if coupled)
    real(r8),       intent(in)    :: precc(pcols)                   !
    real(r8),       intent(in)    :: precl(pcols)                   !
    real(r8),       intent(in)    :: snowhland(pcols)               !
    logical,        intent(in)    :: ghg_chem 
    integer,        intent(in)    :: latmapback(pcols)
    character(len=*), intent(in)  :: chem_name
    integer,         intent(in)   :: troplev(pcols)

    real(r8),       intent(inout) :: qtend(pcols,pver,pcnst)        ! species tendencies (kg/kg/s)
    real(r8),       intent(inout) :: cflx(pcols,pcnst)              ! constituent surface flux (kg/m^2/s)
    real(r8),       intent(out)   :: drydepflx(pcols,pcnst)         ! dry deposition flux (kg/m^2/s)

    type(physics_buffer_desc), pointer :: pbuf(:)

    integer,        intent(in)    :: ixcldliq, ixcldice                     ! indicies of liquid and ice cloud water
    logical, optional, intent(in) :: tropFlag(pcols,pver)      ! 3D tropospheric level flag
                                                               ! true: tropospheric box
    real(r8), optional, intent(in) :: tropFlagInt(pcols,pver)      ! 3D tropospheric level flag
                                                                   ! 1: tropospheric box

    
    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    real(r8), parameter :: m2km  = 1.e-3_r8
    real(r8), parameter :: Pa2mb = 1.e-2_r8

    real(r8),       pointer    :: prain(:,:)
    real(r8),       pointer    :: nevapr(:,:)
    real(r8),       pointer    :: cmfdqr(:,:)
    real(r8),       pointer    :: cldfr(:,:)
    real(r8),       pointer    :: cldtop(:)

    integer      ::  i, k, m, n
    integer      ::  tim_ndx
    real(r8)     ::  delt_inverse
    real(r8)     ::  esfact
    integer      ::  latndx(pcols)                         ! chunk lat indicies
    integer      ::  lonndx(pcols)                         ! chunk lon indicies
    real(r8)     ::  invariants(ncol,pver,nfs)
    real(r8)     ::  col_dens(ncol,pver,   max(1, nabscol+1) )    ! column densities (molecules/cm^2)
    real(r8)     ::  col_delta(ncol,0:pver,max(1, nabscol+1) )    ! layer column densities (molecules/cm^2)
    real(r8)     ::  extfrc(ncol,pver,max(1,extcnt))
    real(r8)     ::  vmr(ncol,pver,gas_pcnst)         ! xported species (vmr)
    real(r8)     ::  reaction_rates(ncol,pver,max(1,rxntot))      ! reaction rates
    real(r8)     ::  depvel(ncol,gas_pcnst)                ! dry deposition velocity (cm/s)
    real(r8)     ::  het_rates(ncol,pver,max(1,gas_pcnst)) ! washout rate (1/s)
    real(r8), dimension(ncol,pver) :: &
         h2ovmr, &                                         ! water vapor volume mixing ratio
         mbar, &                                           ! mean wet atmospheric mass ( amu )
         zmid, &                                           ! midpoint geopotential in km
         zmidr, &                                          ! midpoint geopotential in km realitive to surf
         sulfate, &                                        ! trop sulfate aerosols
         pmb                                               ! pressure at midpoints ( hPa )
    real(r8), dimension(ncol,pver) :: &
         cwat, &                                           ! cloud water mass mixing ratio (kg/kg)
         wrk
    real(r8), dimension(ncol,pver+1) :: &
         zintr                                              ! interface geopotential in km realitive to surf
    real(r8), dimension(ncol,pver+1) :: &
         zint                                              ! interface geopotential in km
    real(r8), dimension(ncol)  :: &
         zen_angle, &                                      ! solar zenith angles
         zsurf, &                                          ! surface height (m)
         rlats, rlons                                      ! chunk latitudes and longitudes (radians)
    real(r8) :: sza(ncol)                                  ! solar zenith angles (degrees)
    real(r8), parameter :: rad2deg = 180._r8/pi                ! radians to degrees conversion factor
    real(r8) :: relhum(ncol,pver)                          ! relative humidity
    real(r8) :: satv(ncol,pver)                            ! wrk array for relative humidity
    real(r8) :: satq(ncol,pver)                            ! wrk array for relative humidity

    integer                   :: j,wd_index
    real(r8), dimension(ncol) :: noy_wk, sox_wk, nhx_wk
    integer                   ::  ltrop_sol(pcols)                 ! tropopause vertical index used in chem solvers
    real(r8)                  ::  strato_sad(pcols,pver)           ! stratospheric SAD (cm2/cm3)
    real(r8)                  ::  sad_total(pcols,pver)            ! total trop. SAD (cm2/cm3)
    real(r8)                  ::  sad_sage(pcols,pver)             ! SAGE SAD (cm2/cm3)

    real(r8) :: tvs(pcols)
    integer  :: ncdate,yr,mon,day,sec
    real(r8) :: wind_speed(pcols)        ! surface wind speed (m/s)
    logical, parameter :: dyn_soilw = .false.
    logical  :: table_soilw
    real(r8) :: soilw(pcols)
    real(r8) :: prect(pcols)
    real(r8) :: sflx(pcols,gas_pcnst)
    real(r8) :: sflx2(pcols,gas_pcnst)
    real(r8) :: dust_vmr(ncol,pver,ndust)
    real(r8) :: noy(ncol,pver)
    real(r8) :: sox(ncol,pver)
    real(r8) :: nhx(ncol,pver)
    real(r8) :: dt_diag(pcols,8)               ! od diagnostics
    real(r8) :: fracday(pcols)                 ! fraction of day
    real(r8) :: o2mmr(ncol,pver)               ! o2 concentration (kg/kg)
    real(r8) :: ommr(ncol,pver)                ! o concentration (kg/kg)
    real(r8) :: mmr(pcols,pver,gas_pcnst)      ! chem working concentrations (kg/kg)
    real(r8) :: mmr_new(pcols,pver,gas_pcnst)      ! chem working concentrations (kg/kg)
    real(r8) :: hno3_gas(ncol,pver)            ! hno3 gas phase concentration (mol/mol)
    real(r8) :: hno3_cond(ncol,pver,2)         ! hno3 condensed phase concentration (mol/mol)
    real(r8) :: h2o_gas(ncol,pver)             ! h2o gas phase concentration (mol/mol)
    real(r8) :: h2o_cond(ncol,pver)            ! h2o condensed phase concentration (mol/mol)
    real(r8) :: cldice(pcols,pver)             ! cloud water "ice" (kg/kg)
    real(r8) :: radius_strat(ncol,pver,3)      ! radius of sulfate, nat, & ice ( cm )
    real(r8) :: sad_strat(ncol,pver,3)         ! surf area density of sulfate, nat, & ice ( cm^2/cm^3 )
    real(r8) :: mmr_tend(pcols,pver,gas_pcnst) ! chemistry species tendencies (kg/kg/s)
    real(r8) :: qh2o(pcols,pver)               ! specific humidity (kg/kg)
    real(r8) :: qvmr(ncol,pver)                ! water vapor vmr
    real(r8) :: delta
     !linoz v3
    real(r8) :: xsfc(4, ncol)                  ! constant surface concentration for o3/o3lnz, n2o, noylnz, and ch4 if called
    real(r8) :: o3lsfcsink(ncol)               ! linoz o3lnz surface sink from call lin_strat_sfcsink
  !  
  ! for aerosol formation....  
    real(r8) :: del_h2so4_gasprod(ncol,pver)

    real(r8) :: vmr0(ncol,pver,gas_pcnst)
    real(r8) :: vmr_old(ncol,pver,gas_pcnst)
    real(r8) :: vmr_old2(ncol,pver,gas_pcnst)
    real(r8) :: tmp(ncol)
    
    ! for explicit chem prodcution and loss
    real(r8)     ::  chem_prod(ncol,pver,gas_pcnst)         ! xported species (vmr/delt)
    real(r8)     ::  chem_loss(ncol,pver,gas_pcnst)         ! xported species (vmr/delt)
    real(r8)     ::  chemmp_prod(ncol,pver,gas_pcnst)         ! xported species (vmr/delt)
    real(r8)     ::  chemmp_loss(ncol,pver,gas_pcnst)         ! xported species (vmr/delt)
    real(r8)     ::  diags_reaction_rates(ncol,pver,max(1,rxntot))       
   
    ! flags for MMF configuration
    logical :: use_MMF, use_ECPP

    ! output gas chemistry tracer concentrations and tendencies
    logical :: history_gaschmbudget
    logical :: history_gaschmbudget_2D
    logical :: history_gaschmbudget_2D_levels
    logical :: history_UCIgaschmbudget_2D
    logical :: history_UCIgaschmbudget_2D_levels

    integer ::  gas_ac_idx
    real(r8), pointer, dimension(:) :: gas_ac_2D
    real(r8), pointer, dimension(:,:) :: gas_ac
    real(r8) :: ftem(pcols,pver) ! tmp space

    call phys_getopts (use_MMF_out = use_MMF)
    call phys_getopts (use_ECPP_out  = use_ECPP )
    call phys_getopts( history_gaschmbudget_out = history_gaschmbudget, &
                    history_gaschmbudget_2D_out = history_gaschmbudget_2D, &
             history_gaschmbudget_2D_levels_out = history_gaschmbudget_2D_levels,&
                 history_UCIgaschmbudget_2D_out = history_UCIgaschmbudget_2D, &
          history_UCIgaschmbudget_2D_levels_out = history_UCIgaschmbudget_2D_levels)

    call t_startf('chemdr_init')

    ! initialize to NaN to hopefully catch user defined rxts that go unset
    reaction_rates(:,:,:) = nan

    delt_inverse = 1._r8 / delt
    !-----------------------------------------------------------------------      
    !        ... Get chunck latitudes and longitudes
    !-----------------------------------------------------------------------      
    call get_lat_all_p( lchnk, ncol, latndx )
    call get_lon_all_p( lchnk, ncol, lonndx )
    call get_rlat_all_p( lchnk, ncol, rlats )
    call get_rlon_all_p( lchnk, ncol, rlons )
    tim_ndx = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ndx_prain,      prain,  start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cldfr,      cldfr,  start=(/1,1,tim_ndx/), kount=(/ncol,pver,1/) )
    call pbuf_get_field(pbuf, ndx_cmfdqr,     cmfdqr, start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_nevapr,     nevapr, start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cldtop,     cldtop )

    !-----------------------------------------------------------------------      
    !        ... Calculate cosine of zenith angle
    !            then cast back to angle (radians)
    !-----------------------------------------------------------------------      
    call zenith( calday, rlats, rlons, zen_angle, ncol )
    zen_angle(:) = acos( zen_angle(:) )

    sza(:) = zen_angle(:) * rad2deg
    call outfld( 'SZA',   sza,    ncol, lchnk )

    !-----------------------------------------------------------------------      
    !        ... Xform geopotential height from m to km 
    !            and pressure from Pa to mb
    !-----------------------------------------------------------------------      
    zsurf(:ncol) = rga * phis(:ncol)
    do k = 1,pver
       zintr(:ncol,k) = m2km * zi(:ncol,k)
       zmidr(:ncol,k) = m2km * zm(:ncol,k)
       zmid(:ncol,k) = m2km * (zm(:ncol,k) + zsurf(:ncol))
       zint(:ncol,k) = m2km * (zi(:ncol,k) + zsurf(:ncol))
       pmb(:ncol,k)  = Pa2mb * pmid(:ncol,k)
    end do
    zint(:ncol,pver+1) = m2km * (zi(:ncol,pver+1) + zsurf(:ncol))
    zintr(:ncol,pver+1)= m2km *  zi(:ncol,pver+1)

    !-----------------------------------------------------------------------      
    !        ... map incoming concentrations to working array
    !-----------------------------------------------------------------------      
    do m = 1,pcnst
       n = map2chm(m)
       if( n > 0 ) then
          mmr(:ncol,:,n) = q(:ncol,:,m)
       end if
    end do

    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

    !-----------------------------------------------------------------------      
    !        ... Set atmosphere mean mass
    !-----------------------------------------------------------------------      
    call set_mean_mass( ncol, mmr, mbar )

    !-----------------------------------------------------------------------      
    !        ... Xform from mmr to vmr
    !-----------------------------------------------------------------------      
    call mmr2vmr( mmr, vmr, mbar, ncol )

    if (h2o_ndx>0) then
       !-----------------------------------------------------------------------      
       !        ... store water vapor in wrk variable
       !-----------------------------------------------------------------------      
       qh2o(:ncol,:) = mmr(:ncol,:,h2o_ndx)
       h2ovmr(:ncol,:) = vmr(:ncol,:,h2o_ndx)
    else
       qh2o(:ncol,:) = q(:ncol,:,1)
       !-----------------------------------------------------------------------      
       !        ... Xform water vapor from mmr to vmr and set upper bndy values
       !-----------------------------------------------------------------------      
       call h2o_to_vmr( q(:,:,1), h2ovmr, mbar, ncol )

       call set_fstrat_h2o( h2ovmr, pmid, troplev, calday, ncol, lchnk )

    endif

    !-----------------------------------------------------------------------      
    !        ... force ion/electron balance
    !-----------------------------------------------------------------------      
    call charge_balance( ncol, vmr )

    !-----------------------------------------------------------------------      
    !        ... Set the "invariants"
    !-----------------------------------------------------------------------  
    call setinv( invariants, tfld, h2ovmr, vmr, pmid, ncol, lchnk, pbuf )

    !-----------------------------------------------------------------------      
    !        ... interpolate SAGEII data for surface area
    !-----------------------------------------------------------------------
    strato_sad(:,:) = 0.0_r8
    sad_sage(:,:) = 0.0_r8
    call strato_sad_set( pmid, sad_sage, ncol, lchnk)

    !-----------------------------------------------------------------------      
    !        ... set tropospheric ozone for Linoz_MAM  (pjc, 2015)
    !-----------------------------------------------------------------------
!    if ( chem_name == 'linoz_mam3'.or.chem_name == 'linoz_mam4_resus'.or.chem_name == 'linoz_mam4_resus_mom' &
!       .or.chem_name == 'linoz_mam4_resus_soag'.or.chem_name == 'linoz_mam4_resus_mom_soag' ) then
!     write(iulog,*) 'Set tropospheric ozone for linoz_mam: inv_ndx_cnst_o3 =',inv_ndx_cnst_o3
!      do k = 1, pver                !Following loop logic from below.  However, reordering loops can get rid of IF statement.
!         do i = 1, ncol
!            if( k > troplev(i) ) then
!              vmr(i,k,o3_ndx) = invariants(i,k,inv_ndx_cnst_o3) / invariants(i,k,inv_ndx_m)   ! O3 and cnst_o3
!            endif
!         end do
!      end do
!    end if

    !-----------------------------------------------------------------
    ! ... zero out sulfate below tropopause
    !-----------------------------------------------------------------
    do k = 1, pver
       do i = 1, ncol
          if( k < troplev(i) ) then
             strato_sad(i,k) = sad_sage(i,k)
          else
             strato_sad(i,k) = 0.0_r8
          endif
       end do
    end do

    if ( has_strato_chem ) then
       !-----------------------------------------------------------------------      
       !        ... initialize condensed and gas phases; all hno3 to gas
       !-----------------------------------------------------------------------      
       do k = 1,pver
          hno3_gas(:,k)   = vmr(:,k,hno3_ndx)
          h2o_gas(:,k)    = h2ovmr(:,k)
          wrk(:,k)        = h2ovmr(:,k)
          cldice(:ncol,k) = q(:ncol,k,cldice_ndx)
       end do
       do m = 1,2
          do k = 1,pver
             hno3_cond(:,k,m) = 0._r8
          end do
       end do
       call mmr2vmri( cldice, h2o_cond, mbar, cnst_mw(cldice_ndx), ncol )

       !-----------------------------------------------------------------------      
       !        ... call SAD routine
       !-----------------------------------------------------------------------      
       call sad_strat_calc( lchnk, invariants(:ncol,:,indexm), pmb, tfld, hno3_gas, &
            hno3_cond, h2o_gas, h2o_cond, strato_sad(:ncol,:), radius_strat, &
            sad_strat, ncol, pbuf )
       do k = 1,pver
          vmr(:,k,hno3_ndx) = hno3_gas(:,k)
          h2ovmr(:,k)       = h2o_gas(:,k)
          vmr(:,k,h2o_ndx)  = h2o_gas(:,k)
          wrk(:,k)          = (h2ovmr(:,k) - wrk(:,k))*delt_inverse
       end do
       call outfld( 'QDSAD', wrk(:,:), ncol, lchnk )
       call outfld( 'SAD', strato_sad(:ncol,:), ncol, lchnk )
       call outfld( 'SAD_SULFC', sad_strat(:,:,1), ncol, lchnk )
       call outfld( 'SAD_SAGE',   sad_sage(:,:), ncol, lchnk )
       call outfld( 'SAD_LNAT', sad_strat(:,:,2), ncol, lchnk )
       call outfld( 'SAD_ICE', sad_strat(:,:,3), ncol, lchnk )
       call outfld( 'RAD_SULFC', radius_strat(:,:,1), ncol, lchnk )
       call outfld( 'RAD_LNAT', radius_strat(:,:,2), ncol, lchnk )
       call outfld( 'RAD_ICE', radius_strat(:,:,3), ncol, lchnk )

       !-----------------------------------------------------------------------      
       !        ... call aerosol reaction rates
       !-----------------------------------------------------------------------      
       call ratecon_sfstrat( invariants(:,:,indexm), pmid, tfld, &
            radius_strat(:,:,1), sad_strat(:,:,1), sad_strat(:,:,2), &
            sad_strat(:,:,3), h2ovmr, vmr, reaction_rates, ncol )

    endif

    !-----------------------------------------------------------------------      
    !        ... Set the column densities at the upper boundary
    !-----------------------------------------------------------------------      
    call set_ub_col( col_delta, vmr, invariants, pint(:,1), pdel, ncol, lchnk)

    !-----------------------------------------------------------------------      
    !       ...  Set rates for "tabular" and user specified reactions
    !-----------------------------------------------------------------------      
    call setrxt( reaction_rates, tfld, invariants(1,1,indexm), ncol )

    sulfate(:,:) = 0._r8
    if( so4_ndx < 1 ) then ! get offline so4 field if not prognostic
       call sulf_interp( ncol, lchnk, sulfate )
    else
       sulfate(:,:) = vmr(:,:,so4_ndx)
    endif

    !-----------------------------------------------------------------
    ! ... zero out sulfate above tropopause
    !-----------------------------------------------------------------
    do k = 1, pver
       do i = 1, ncol
          if( k < troplev(i) ) then
             sulfate(i,k) = 0.0_r8
          else
             sulfate(i,k) = sulfate(i,k)
          end if
       end do
    end do

    !-----------------------------------------------------------------
    !	... compute the relative humidity
    !-----------------------------------------------------------------
    call qsat(tfld(:ncol,:), pmid(:ncol,:), satv, satq)

    do k = 1,pver
       relhum(:,k) = .622_r8 * h2ovmr(:,k) / satq(:,k)
       relhum(:,k) = max( 0._r8,min( 1._r8,relhum(:,k) ) )
    end do
    
    cwat(:ncol,:pver) = cldw(:ncol,:pver)

    call usrrxt( reaction_rates, tfld, tfld, tfld, invariants, h2ovmr, ps, &
                 pmid, invariants(:,:,indexm), sulfate, mmr, relhum, strato_sad, &
                 troplev, ncol, sad_total, cwat, mbar, pbuf )

    call outfld( 'SAD_TROP', sad_total(:ncol,:), ncol, lchnk )

    if (het1_ndx>0) then
       call outfld( 'het1_total', reaction_rates(:,:,het1_ndx), ncol, lchnk )
    endif

    if (ghg_chem) then
       call ghg_chem_set_rates( reaction_rates, latmapback, zen_angle, ncol, lchnk )
    endif

    do i = phtcnt+1,rxntot
       call outfld( rxn_names(i-phtcnt), reaction_rates(:,:,i), ncol, lchnk )
    enddo

    call adjrxt( reaction_rates, invariants, invariants(1,1,indexm), ncol )

    !-----------------------------------------------------------------------
    !        ... Compute the photolysis rates at time = t(n+1)
    !-----------------------------------------------------------------------      
    !-----------------------------------------------------------------------      
    !     	... Set the column densities
    !-----------------------------------------------------------------------      
    call setcol( col_delta, col_dens, vmr, pdel,  ncol )

    !-----------------------------------------------------------------------      
    !     	... Calculate the photodissociation rates
    !-----------------------------------------------------------------------      

    esfact = 1._r8
    call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr  , &
         delta, esfact )

    if ( do_cloudj_photolysis ) then 

       !------------------------------------------------------------------
       !	... compute the photolysis rates using CLoud-J and Fast-JX
       !------------------------------------------------------------------
       
!!! Photolosys rates for clear sky (Diagnostic sensitivity)
       if ( do_cloudj_nocloud_diag ) then
          
          call cloudJ_interface( reaction_rates, vmr, mmr, invariants, tfld, &
            cldfr, q(:ncol,:,ixcldliq), q(:ncol,:,ixcldice), &
               pmid, pint, zmidr, zint, rlats, rlons, col_dens, zen_angle, asdir, &
               invariants(1,1,indexm), ps, ts, &
               esfact, relhum, dust_vmr, &
               ncol, lchnk, .false., do_cloudj_aerosols ) 
          
          do i = 1,phtcnt
             call outfld( trim(pht_names(i))//"_nocloud", reaction_rates(:ncol,:,i), ncol, lchnk )
             call outfld( trim(tag_names(i))//"_nocloud", reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
          enddo
       endif
       
!!! Photolosys rates with aerosols and clouds (Diagnostic sensitivity)
       if ( do_cloudj_aerosol_diag ) then
       
          call cloudJ_interface( reaction_rates, vmr, mmr, invariants, tfld, &
               cldfr, q(:ncol,:,ixcldliq), q(:ncol,:,ixcldice), &
               pmid, pint, zmidr, zint, rlats, rlons, col_dens, zen_angle, asdir, &
               invariants(1,1,indexm), ps, ts, &
               esfact, relhum, dust_vmr, &
               ncol, lchnk, do_cloudj_clouds, .true. ) 
          
          do i = 1,phtcnt
             call outfld( trim(pht_names(i))//"_aerosol", reaction_rates(:ncol,:,i), ncol, lchnk )
             call outfld( trim(tag_names(i))//"_aerosol", reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
          enddo
       endif

!!! Photolosys rates from lookup table (Diagnostic sensitivity)
       if ( do_cloudj_lookup_diag ) then

          call table_photo( reaction_rates, pmid, pdel, tfld, zmid, zint, &
               col_dens, zen_angle, asdir, cwat, cldfr, &
               esfact, vmr, invariants, ncol, lchnk, pbuf )
          
          do i = 1,phtcnt
             call outfld( trim(pht_names(i))//"_lookup", reaction_rates(:ncol,:,i), ncol, lchnk )
             call outfld( trim(tag_names(i))//"_lookup", reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
          enddo
       endif

!!! Photolosys rates for E3SM
       
       call cloudJ_interface( reaction_rates, vmr, mmr, invariants, tfld, &
            cldfr, q(:ncol,:,ixcldliq), q(:ncol,:,ixcldice), &
            pmid, pint, zmidr, zint, rlats, rlons, col_dens, zen_angle, asdir, &
            invariants(1,1,indexm), ps, ts, &
            esfact, relhum, dust_vmr, &
            ncol, lchnk, do_cloudj_clouds, do_cloudj_aerosols ) 

    else if ( xactive_prates ) then
       if ( dst_ndx > 0 ) then
          dust_vmr(:ncol,:,1:ndust) = vmr(:ncol,:,dst_ndx:dst_ndx+ndust-1)
       else 
          dust_vmr(:ncol,:,:) = 0._r8
       endif

       !-----------------------------------------------------------------
       !	... compute the photolysis rates
       !-----------------------------------------------------------------
       call xactive_photo( reaction_rates, vmr, tfld, cwat, cldfr, &
            pmid, zmidr, col_dens, zen_angle, asdir, &
            invariants(1,1,indexm), ps, ts, &
            esfact, relhum, dust_vmr, dt_diag, fracday, ncol, lchnk )

       call outfld('DTCBS',   dt_diag(:ncol,1), ncol, lchnk )
       call outfld('DTOCS',   dt_diag(:ncol,2), ncol, lchnk )
       call outfld('DTSO4',   dt_diag(:ncol,3), ncol, lchnk )
       call outfld('DTANT',   dt_diag(:ncol,4), ncol, lchnk )
       call outfld('DTSAL',   dt_diag(:ncol,5), ncol, lchnk )
       call outfld('DTDUST',  dt_diag(:ncol,6), ncol, lchnk )
       call outfld('DTSOA',   dt_diag(:ncol,7), ncol, lchnk )
       call outfld('DTTOTAL', dt_diag(:ncol,8), ncol, lchnk )
       call outfld('FRACDAY', fracday(:ncol), ncol, lchnk )

    else

       !-----------------------------------------------------------------
       !	... lookup the photolysis rates from table
       !-----------------------------------------------------------------
       call table_photo( reaction_rates, pmid, pdel, tfld, zmid, zint, &
                         col_dens, zen_angle, asdir, cwat, cldfr, &
                         esfact, vmr, invariants, ncol, lchnk, pbuf )
    endif

    do i = 1,phtcnt
       call outfld( pht_names(i), reaction_rates(:ncol,:,i), ncol, lchnk )
       call outfld( tag_names(i), reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    !! NOTE: the lookup table and xactive_photo use cwat, which is the sum of liquid and ice.
    call outfld( 'Photolysis_CLOUD',       cldfr,               ncol, lchnk )
    call outfld( 'Photolysis_CLOUD_water', q(:ncol,:,ixcldliq), ncol, lchnk )
    call outfld( 'Photolysis_CLOUD_ice',   q(:ncol,:,ixcldice), ncol, lchnk )
     
    !-----------------------------------------------------------------------      
    !     	... Adjust the photodissociation rates
    !-----------------------------------------------------------------------  
    call O1D_to_2OH_adj( reaction_rates, invariants, invariants(:,:,indexm), ncol, tfld )
    call phtadj( reaction_rates, invariants, invariants(:,:,indexm), ncol )

    !-----------------------------------------------------------------------
    !        ... Compute the extraneous frcing at time = t(n+1)
    !-----------------------------------------------------------------------      
    if ( o2_ndx > 0 .and. o_ndx > 0 ) then
       do k = 1,pver
          o2mmr(:ncol,k) = mmr(:ncol,k,o2_ndx)
          ommr(:ncol,k)  = mmr(:ncol,k,o_ndx)
       end do
    endif
    !-----------------------------------------------------------------------
    !        ... Compute the extraneous frcing at time = t(n+1)
    !-----------------------------------------------------------------------      
    call setext( extfrc, zint, zintr, cldtop, &
                 zmid, lchnk, tfld, o2mmr, ommr, &
                 pmid, mbar, rlats, calday, ncol, rlons, pbuf )

    do m = 1,extcnt
       if( m /= synoz_ndx ) then
          do k = 1,pver
             extfrc(:ncol,k,m) = extfrc(:ncol,k,m) / invariants(:ncol,k,indexm)
          end do
       endif
       call outfld( extfrc_name(m), extfrc(:ncol,:,m), ncol, lchnk )
    end do

    !-----------------------------------------------------------------------
    !        ... Form the washout rates
    !-----------------------------------------------------------------------      
    if ( do_neu_wetdep ) then
      het_rates = 0._r8
    else
      if (use_MMF .and. use_ECPP) then
        ! ECPP handles wet removal, so set to zero here
        het_rates(:,:,:) = 0.0
      else
        call sethet( het_rates, pmid, zmid, phis, tfld, &
                     cmfdqr, prain, nevapr, delt, invariants(:,:,indexm), &
                     vmr, ncol, lchnk )
      end if

       if(.not. convproc_do_aer) then 
          call het_diags( het_rates(:ncol,:,:), mmr(:ncol,:,:), pdel(:ncol,:), lchnk, ncol )
       endif
    end if

    do i = phtcnt+1,rxt_tag_cnt
       call outfld( tag_names(i), reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    if ( has_linoz_data .and. .not. &
       (chem_name == 'linoz_mam3'.or.chem_name == 'linoz_mam4_resus'.or.chem_name == 'linoz_mam4_resus_mom' &
       .or.chem_name == 'linoz_mam4_resus_soag'.or.chem_name == 'linoz_mam4_resus_mom_soag' ) ) then
       ltrop_sol(:ncol) = troplev(:ncol)
    else
       ltrop_sol(:ncol) = 0 ! apply solver to all levels
    endif

    ! save h2so4 before gas phase chem (for later new particle nucleation)
    if (ndx_h2so4 > 0) then
       del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4)
    else
       del_h2so4_gasprod(:,:) = 0.0_r8
    endif

    vmr0(:ncol,:,:) = vmr(:ncol,:,:) ! mixing ratios before chemistry changes

    call t_stopf('chemdr_init')

    if ( history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
       call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSW' )
    endif

    !=======================================================================
    !        ... Call the class solution algorithms
    !=======================================================================
    !-----------------------------------------------------------------------
    !	... Solve for "Implicit" species
    !-----------------------------------------------------------------------
    if ( has_strato_chem ) wrk(:,:) = vmr(:,:,h2o_ndx)
    !
    if (uci1_ndx > 0) then
       vmr_old2(:ncol,:,:) = vmr(:ncol,:,:)
    endif
    call t_startf('imp_sol')
    call imp_sol( vmr, reaction_rates, het_rates, extfrc, delt, &
                  invariants(1,1,indexm), ncol, lchnk, ltrop_sol(:ncol) ) 
    call t_stopf('imp_sol')

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TDI' )
       endif
       if ( history_gaschmbudget_2D ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTDI' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSI' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDI_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDI_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    ! ozone production 
    if (uci1_ndx > 0) then
      chem_prod(:,:,:) = 0._r8
      chem_loss(:,:,:) = 0._r8

      level0_loop: do k = 1, pver
      column0_loop: do i = 1, ncol
      IF (k <= ltrop_sol(i)) CYCLE column0_loop
      chem_prod(i,k,o3_ndx ) =  2._r8 * reaction_rates(i,k,jo2_b_ndx) &
         +  reaction_rates(i,k,jno2_ndx)*vmr(i,k,no2_ndx) &
         +  reaction_rates(i,k,jno3_a_ndx)*vmr(i,k,no3_ndx) &
         +  reaction_rates(i,k,jn2o5_b_ndx)*vmr(i,k,n2o5_ndx) &
         +  reaction_rates(i,k,po3_oh_ndx)*vmr(i,k,oh_ndx)*vmr(i,k,oh_ndx)

      chem_loss(i,k,o3_ndx) = (reaction_rates(i,k,uci1_ndx) &
         + reaction_rates(i,k,uci2_ndx) &
         + reaction_rates(i,k,uci3_ndx) * vmr(i,k,ch4lnz_ndx) &
         + reaction_rates(i,k,lc2h4_o3_ndx)*vmr(i,k,c2h4_ndx) &
         + reaction_rates(i,k,lisop_o3_ndx)*vmr(i,k,isop_ndx)  &
         + reaction_rates(i,k,lo3_oh_ndx)*vmr(i,k,oh_ndx)  &
         + reaction_rates(i,k,lo3_ho2_ndx)*vmr(i,k,ho2_ndx)  &
         + reaction_rates(i,k,lo3_no_ndx)*vmr(i,k,no_ndx)  &
         + reaction_rates(i,k,lo3_no2_ndx)*vmr(i,k,no2_ndx)  &
         + reaction_rates(i,k,lmvkmacr_o3_ndx)*vmr(i,k,mvkmacr_ndx)  &
         + het_rates(i,k,o3_ndx)) *vmr(i,k,o3_ndx)  

      end do column0_loop
      end do level0_loop
    endif

    if ( history_UCIgaschmbudget_2D .or. history_UCIgaschmbudget_2D_levels) then
       if ( history_UCIgaschmbudget_2D .and. uci1_ndx > 0) then
          call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTIP' )
          call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTIL' )
       endif
       if ( history_UCIgaschmbudget_2D_levels .and. uci1_ndx > 0) then
         call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTIP_LL' )
         call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTIP_trop', tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTIL_LL' )
         call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTIL_trop', tropFlagInt )
       endif
    endif

    ! reset vmr to pre-imp_sol values for stratospheric boxes
    if (uci1_ndx > 0) then
       diags_reaction_rates(:,:,:) = reaction_rates(:,:,:)
       do i = 1,ncol
          do k = 1,pver
             if ( .not. tropFlag(i,k) ) then
                vmr(i,k,:) = vmr_old2(i,k,:)
                diags_reaction_rates(i,k,:) = 0._r8
             endif
          enddo
       enddo
    endif

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget  .and. uci1_ndx > 0) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TRI' )
       endif
       if ( history_gaschmbudget_2D .and. uci1_ndx > 0) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTRI' )
       endif
       if ( history_gaschmbudget_2D_levels .and. uci1_ndx > 0) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTRI_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTRI_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    if (uci1_ndx > 0) then
      chemmp_prod(:,:,:) = 0._r8
      chemmp_loss(:,:,:) = 0._r8
      chem_prod(:,:,:) = 0._r8
      chem_loss(:,:,:) = 0._r8
      ! ozone production 
      level_loop: do k = 1, pver
      column_loop: do i = 1, ncol
      IF (k <= ltrop_sol(i)) CYCLE column_loop
      chemmp_prod(i,k,o3_ndx) = 2._r8 * diags_reaction_rates(i,k,jo2_b_ndx) & 
         + diags_reaction_rates(i,k,lch3o2_no_ndx)*vmr(i,k,ch3o2_ndx)*vmr(i,k,no_ndx) &
         + diags_reaction_rates(i,k,lno_ho2_ndx)*vmr(i,k,no_ndx)*vmr(i,k,ho2_ndx) &
         + diags_reaction_rates(i,k,lc2h5o2_no_ndx)*vmr(i,k,c2h5o2_ndx)*vmr(i,k,no_ndx) &
         + diags_reaction_rates(i,k,lch3co3_no_ndx)*vmr(i,k,ch3co3_ndx)*vmr(i,k,no_ndx) &
         + diags_reaction_rates(i,k,lroho2_no_ndx)*vmr(i,k,roho2_ndx)*vmr(i,k,no_ndx) &
         + diags_reaction_rates(i,k,lisopo2_no_ndx)*vmr(i,k,isopo2_ndx)*vmr(i,k,no_ndx) &
         + diags_reaction_rates(i,k,lmvko2_no_ndx)*vmr(i,k,mvko2_ndx)*vmr(i,k,no_ndx) 
 
      chemmp_loss(i,k,o3_ndx) = (diags_reaction_rates(i,k,uci1_ndx) &
         + diags_reaction_rates(i,k,uci2_ndx) &
         + diags_reaction_rates(i,k,uci3_ndx) * vmr(i,k,ch4lnz_ndx) &
         + diags_reaction_rates(i,k,lc2h4_o3_ndx)*vmr(i,k,c2h4_ndx) &
         + diags_reaction_rates(i,k,lisop_o3_ndx)*vmr(i,k,isop_ndx)  &
         + diags_reaction_rates(i,k,lo3_oh_ndx)*vmr(i,k,oh_ndx)  &
         + diags_reaction_rates(i,k,lo3_ho2_ndx)*vmr(i,k,ho2_ndx)  &
         + diags_reaction_rates(i,k,lmvkmacr_o3_ndx)*vmr(i,k,mvkmacr_ndx) ) *vmr(i,k,o3_ndx)  

      chem_prod(i,k,o3_ndx ) =  2._r8 * diags_reaction_rates(i,k,jo2_b_ndx) &
         +  diags_reaction_rates(i,k,jno2_ndx)*vmr(i,k,no2_ndx) &
         +  diags_reaction_rates(i,k,jno3_a_ndx)*vmr(i,k,no3_ndx) &
         +  diags_reaction_rates(i,k,jn2o5_b_ndx)*vmr(i,k,n2o5_ndx) &
         +  diags_reaction_rates(i,k,po3_oh_ndx)*vmr(i,k,oh_ndx)*vmr(i,k,oh_ndx)

      chem_loss(i,k,o3_ndx) = (diags_reaction_rates(i,k,uci1_ndx) &
         + diags_reaction_rates(i,k,uci2_ndx) &
         + diags_reaction_rates(i,k,uci3_ndx) * vmr(i,k,ch4lnz_ndx) &
         + diags_reaction_rates(i,k,lc2h4_o3_ndx)*vmr(i,k,c2h4_ndx) &
         + diags_reaction_rates(i,k,lisop_o3_ndx)*vmr(i,k,isop_ndx)  &
         + diags_reaction_rates(i,k,lo3_oh_ndx)*vmr(i,k,oh_ndx)  &
         + diags_reaction_rates(i,k,lo3_ho2_ndx)*vmr(i,k,ho2_ndx)  &
         + diags_reaction_rates(i,k,lo3_no_ndx)*vmr(i,k,no_ndx)  &
         + diags_reaction_rates(i,k,lo3_no2_ndx)*vmr(i,k,no2_ndx)  &
         + diags_reaction_rates(i,k,lmvkmacr_o3_ndx)*vmr(i,k,mvkmacr_ndx)  &
         + het_rates(i,k,o3_ndx)) *vmr(i,k,o3_ndx)  

      end do column_loop
      end do level_loop
    endif !uci1_ndx

    if ( history_UCIgaschmbudget_2D .or. history_UCIgaschmbudget_2D_levels) then
       if ( history_UCIgaschmbudget_2D .and. uci1_ndx > 0) then
          call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DCIP' )
          call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DCIL' )
          call gaschmmass_diags( lchnk, ncol, chemmp_prod(:ncol,:,:), vmr_old(ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMPP' )
          call gaschmmass_diags( lchnk, ncol, chemmp_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMPL' )
       endif
       if ( history_UCIgaschmbudget_2D_levels .and. uci1_ndx > 0) then
         call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCIP_LL' )
         call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCIP_trop', tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCIL_LL' )
         call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCIL_trop', tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, chemmp_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMPP_LL' )
         call gaschmmass_diags( lchnk, ncol, chemmp_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMPP_trop', tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, chemmp_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMPL_LL' )
         call gaschmmass_diags( lchnk, ncol, chemmp_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMPL_trop', tropFlagInt )
       endif
    endif
    !-----------------------------------------------------------------------
    !	... Solve for "Explicit" species
    !-----------------------------------------------------------------------
    if (uci1_ndx > 0) then
       vmr_old2(:ncol,:,:) = vmr(:ncol,:,:)
    endif
    call t_startf('exp_sol')
    call exp_sol( vmr, reaction_rates, het_rates, extfrc, delt, invariants(1,1,indexm), ncol, lchnk, ltrop_sol, &
                  diags_reaction_rates, chem_prod, chem_loss, chemmp_prod, chemmp_loss)
    call t_stopf('exp_sol')

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels .or.&
         history_UCIgaschmbudget_2D .or. history_UCIgaschmbudget_2D_levels ) then
       if ( history_gaschmbudget ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, 'TDE' )
       endif
       if ( history_gaschmbudget_2D ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDE' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMSE' )
       endif
       if ( history_UCIgaschmbudget_2D .and. uci1_ndx > 0) then
         call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTEP' )
         call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTEL' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDE_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDE_trop', tropFlagInt )
       endif
       if ( history_UCIgaschmbudget_2D_levels .and. uci1_ndx > 0) then
         call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTEP_LL' )
         call gaschmmass_diags( lchnk, ncol, chem_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTEP_trop', tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTEL_LL' )
         call gaschmmass_diags( lchnk, ncol, chem_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTEL_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    ! reset vmr to pre-exp_sol values for stratospheric boxes
    if (uci1_ndx > 0) then
       ! exclude E90 from resetting
       vmr_old2(:,:,e90_ndx) = vmr(:,:,e90_ndx)
       do i = 1,ncol
          do k = 1,pver
             if ( .not. tropFlag(i,k) ) then
                vmr(i,k,:) = vmr_old2(i,k,:)
             endif
          enddo
       enddo
    endif

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget  .and. uci1_ndx > 0) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TRE' )
       endif
       if ( history_gaschmbudget_2D .and. uci1_ndx > 0) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTRE' )
       endif
       if ( history_UCIgaschmbudget_2D .and. uci1_ndx > 0) then
         call gaschmmass_diags( lchnk, ncol, chemmp_prod(:ncol,:,:), vmr_old(ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCEP' )
         call gaschmmass_diags( lchnk, ncol, chemmp_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCEL' )
       endif
       if ( history_gaschmbudget_2D_levels .and. uci1_ndx > 0) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTRE_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTRE_trop', tropFlagInt )
       endif
       if ( history_UCIgaschmbudget_2D_levels .and. uci1_ndx > 0) then
         call gaschmmass_diags( lchnk, ncol, chemmp_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCEP_LL' )
         call gaschmmass_diags( lchnk, ncol, chemmp_prod(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCEP_trop', tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, chemmp_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCEL_LL' )
         call gaschmmass_diags( lchnk, ncol, chemmp_loss(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DCEL_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    if( h2o_ndx>0) call outfld( 'H2O_GAS',  vmr(1,1,h2o_ndx),  ncol ,lchnk )
    if(convproc_do_aer) then 
       call vmr2mmr( vmr, mmr_new, mbar, ncol )  !RCE
       mmr_new(:ncol,:,:) = 0.5_r8*( mmr(:ncol,:,:)+mmr_new(:ncol,:,:) )  !RCE
       !RCE - mmr_new = average of mmr values before and after imp_sol
       call het_diags( het_rates(:ncol,:,:), mmr_new(:ncol,:,:), pdel(:ncol,:), lchnk, ncol )  !RCE
    endif

    ! save h2so4 change by gas phase chem (for later new particle nucleation)
    if (ndx_h2so4 > 0) then
       del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - del_h2so4_gasprod(1:ncol,:)
    endif

!
! Aerosol processes ...
!

    if (uci1_ndx > 0) then
       vmr_old2(:ncol,:,:) = vmr(:ncol,:,:)
    endif
    call t_startf('aero_model_gasaerexch')
    call aero_model_gasaerexch( imozart-1, ncol, lchnk, delt, latndx, lonndx, reaction_rates, &
                                tfld, pmid, pdel, mbar, relhum, &
                                zm,  qh2o, cwat, cldfr, ncldwtr, &
                                invariants(:,:,indexm), invariants, del_h2so4_gasprod,  &
                                vmr0, vmr, pbuf )
    call t_stopf('aero_model_gasaerexch')
!
! Remove the impact of aerosol processes on gas chemistry tracers ...
!
    if (uci1_ndx > 0) then
       do m = 1,pcnst
          n = map2chm( m )
          if ( n > 0 ) then
             if ( .not. any( aer_species == n ) ) then
               if (trim(solsym(n))/='DMS' .and. trim(solsym(n))/='SO2' .and. &
                   trim(solsym(n))/='H2SO4' .and. trim(solsym(n))/='SOAG') then
                   !write(iulog,*) 'n=',n,'solsym=',trim(solsym(n))
                   vmr(:ncol,:,n) = vmr_old2(:ncol,:,n)
               endif
             endif
          endif
       enddo
    endif

    if ( has_strato_chem ) then 

       wrk(:ncol,:) = (vmr(:ncol,:,h2o_ndx) - wrk(:ncol,:))*delt_inverse
       call outfld( 'QDCHEM',   wrk(:ncol,:),         ncol, lchnk )
       call outfld( 'HNO3_GAS', vmr(:ncol,:,hno3_ndx), ncol ,lchnk )

       !-----------------------------------------------------------------------      
       !         ... aerosol settling
       !             first settle hno3(2) using radius ice
       !             secnd settle hno3(3) using radius large nat
       !-----------------------------------------------------------------------      
       wrk(:,:) = vmr(:,:,h2o_ndx)
#ifdef ALT_SETTL
       where( h2o_cond(:,:) > 0._r8 )
          settl_rad(:,:) = radius_strat(:,:,3)
       elsewhere
          settl_rad(:,:) = 0._r8
       endwhere
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), settl_rad, ncol, lchnk, 1 )

       where( h2o_cond(:,:) == 0._r8 )
          settl_rad(:,:) = radius_strat(:,:,2)
       elsewhere
          settl_rad(:,:) = 0._r8
       endwhere
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), settl_rad, ncol, lchnk, 2 )
#else
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), radius_strat(1,1,2), ncol, lchnk, 2 )
#endif

       !-----------------------------------------------------------------------      
       !	... reform total hno3 = gas + all condensed
       !-----------------------------------------------------------------------      
       do k = 1,pver
          vmr(:,k,hno3_ndx) = vmr(:,k,hno3_ndx) + hno3_cond(:,k,1) &
               + hno3_cond(:,k,2) 
       end do
       call outfld( 'HNO3_STS', hno3_cond(:ncol,:,1), ncol ,lchnk )
       call outfld( 'HNO3_NAT', hno3_cond(:ncol,:,2), ncol ,lchnk )
       wrk(:,:) = (vmr(:,:,h2o_ndx) - wrk(:,:))*delt_inverse
       call outfld( 'QDSETT', wrk(:,:), ncol, lchnk )

    endif

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TDA' )
       endif
       if ( history_gaschmbudget_2D ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTDA' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSA' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDA_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDA_trop', tropFlagInt)
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

!
    ! LINOZ
    ! col_dens(:,:,1) is for O3; col_dens(:,:,2) is for O2; col_dens(:,:,3) is for O3LNZ
    if (masterproc)write(iulog,*)'do_lin_strat_chem=',do_lin_strat_chem,'linoz_v2=',linoz_v2,'linoz_v3=',linoz_v3
    
    xsfc(:,:)=0    
    if ( do_lin_strat_chem ) then
       if (uci1_ndx <= 0) then
          !if(linoz_v2) call linv2_strat_chem_solve( ncol, lchnk, vmr,              col_dens(:,:,3), tfld, zen_angle, pmid, delt, rlats, troplev )                  
          if(linoz_v2) call linv2_strat_chem_solve( ncol, lchnk, vmr,              col_dens(:,:,1), tfld, zen_angle, pmid, delt, rlats, troplev, pdeldry)                  
          if(linoz_v3) then
                       call h2o_to_vmr(q(:,:,1), qvmr, mbar, ncol )  
                       !pass in column O3 
                       call linv3_strat_chem_solve( ncol, lchnk, vmr, qvmr, xsfc,  col_dens(:,:,1), tfld, zen_angle, pmid, delt, rlats, troplev, pdeldry)
          endif
        else
          !if(linoz_v2) call linv2_strat_chem_solve( ncol, lchnk, vmr,              col_dens(:,:,3), tfld, zen_angle, pmid, delt, rlats, troplev, tropFlag=tropFlag )                  
          if(linoz_v2) call linv2_strat_chem_solve( ncol, lchnk, vmr,              col_dens(:,:,1), tfld, zen_angle, pmid, delt, rlats, troplev, pdeldry, tropFlag=tropFlag )                  
          if(linoz_v3) then
                       call h2o_to_vmr(q(:,:,1), qvmr, mbar, ncol )  
                       !pass in column O3
                       if (o3lnz_ndx > 0) then
                          write(iulog,*) 'mo_gas_phase_chemdr: calling LINOZv3 for chemUCI+linozv3 mechanism, O3LNZ should not exist'
                          call endrun('mo_gas_phase_chemdr: calling LINOZv3 for chemUCI+linozv3 mechanism, O3LNZ should not exist')
                       else 
                          call linv3_strat_chem_solve( ncol, lchnk, vmr, qvmr, xsfc,  col_dens(:,:,1), tfld, zen_angle, pmid, delt, rlats, troplev, pdeldry, tropFlag=tropFlag)
                       endif
          endif
          call fstrat_efold_decay(ncol, vmr, delt, troplev, tropFlag=tropFlag) !if chemuci is on
       endif
       
       call lin_strat_sfcsink(ncol, lchnk, vmr, xsfc, delt,   pdel(:ncol,:) )
    end if

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TDL' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'MSL' )
       endif
       if ( history_gaschmbudget_2D ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTDL' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSL' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDL_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDL_trop', tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMSL_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMSL_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    !-----------------------------------------------------------------------      
    !         ... Check for negative values and reset to zero
    !-----------------------------------------------------------------------      
    call negtrc( 'After chemistry ', vmr, ncol )

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TDN' )
       endif
       if ( history_gaschmbudget_2D ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTDN' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSN' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDN_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDN_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    !-----------------------------------------------------------------------      
    !         ... Set upper boundary mmr values
    !-----------------------------------------------------------------------      
    if(.not. do_lin_strat_chem) call set_fstrat_vals( vmr, pmid, pint, troplev, calday, ncol,lchnk )

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TDU' )
       endif
       if ( history_gaschmbudget_2D ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTDU' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSU' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDU_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDU_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    !-----------------------------------------------------------------------      
    !         ... Set fixed lower boundary mmr values
    !-----------------------------------------------------------------------      
    call flbc_set( vmr, ncol, lchnk, map2chm )

    if ( ghg_chem ) then
       call ghg_chem_set_flbc( vmr, ncol )
    endif

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TDB' )
       endif
       if ( history_gaschmbudget_2D ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTDB' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSC' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDB_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDB_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    !-----------------------------------------------------------------------      
    !         ... Surface emissions
    !----------------------------------------------------------------------- 
    ! if chemUCI is used, apply surface emissions to the bottom (srf_emit_nlayer) layers here
    if (uci1_ndx > 0) then 
       do m = 1,pcnst
          n = map2chm( m )
          if ( n > 0 ) then
            if ( any( cflx(:ncol,m) /= 0._r8 ) ) then
              if ( .not. any( aer_species == n ) .and. adv_mass(n) /= 0._r8 ) then
                do k = pver, pver-srf_emit_nlayer+1, -1 ! loop from bottom to top
                  ! kg/m2, tracer mass
                  wrk(:ncol,k) = adv_mass(n)*vmr(:ncol,k,n)/mbar(:ncol,k) &
                                    *pdeldry(:ncol,k)*rga
                  tmp(:ncol) = wrk(:ncol,k) + cflx(:ncol,m)*delt/dble(srf_emit_nlayer)
                  vmr(:ncol,k,n) = tmp(:ncol)*mbar(:ncol,k)/adv_mass(n)/pdeldry(:ncol,k)/rga
                end do
                cflx(:ncol,m) = 0._r8
              endif
            endif
          endif
       end do
    endif

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TDS' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'MSS' )
       endif
       if ( history_gaschmbudget_2D ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTDS' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSS' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDS_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDS_trop', tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMSS_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMSS_trop', tropFlagInt )
       endif
       vmr_old(:ncol,:,:) = vmr(:ncol,:,:)
    endif

    !-----------------------------------------------------------------------      
    !         ... Xform from vmr to mmr (for dry deposition)
    !-----------------------------------------------------------------------      
    call vmr2mmr( vmr, mmr_tend, mbar, ncol )

    call set_short_lived_species( mmr_tend, lchnk, ncol, pbuf )

    !-----------------------------------------------------------------------      
    !         ... Dry deposition
    !----------------------------------------------------------------------- 
    tvs(:ncol) = tfld(:ncol,pver) * (1._r8 + qh2o(:ncol,pver))

    sflx(:,:) = 0._r8
    call get_ref_date(yr, mon, day, sec)
    ncdate = yr*10000 + mon*100 + day
    wind_speed(:ncol) = sqrt( ufld(:ncol,pver)*ufld(:ncol,pver) + vfld(:ncol,pver)*vfld(:ncol,pver) )
    prect(:ncol) = precc(:ncol) + precl(:ncol)

    call t_startf('drydep')
    if ( drydep_method == DD_XLND ) then
       soilw = -99
       call drydep( ocnfrac, icefrac, ncdate, ts, ps,  &
            wind_speed, qh2o(:,pver), tfld(:,pver), pmid(:,pver), prect, &
            snowhland, fsds, depvel, sflx, mmr, &
            tvs, soilw, relhum(:,pver:pver), ncol, lonndx, latndx, lchnk )
    else if ( drydep_method == DD_XATM ) then
       table_soilw = has_drydep( 'H2' ) .or. has_drydep( 'CO' )
       if( .not. dyn_soilw .and. table_soilw ) then
          call set_soilw( soilw, lchnk, calday )
       end if
       call drydep( ncdate, ts, ps,  &
            wind_speed, qh2o(:,pver), tfld(:,pver), pmid(:,pver), prect, &
            snowhland, fsds, depvel, sflx, mmr, &
            tvs, soilw, relhum(:,pver:pver), ncol, lonndx, latndx, lchnk )
    else if ( drydep_method == DD_TABL ) then
       call drydep( calday, ts, zen_angle, &
            depvel, sflx, mmr, pmid(:,pver), &
            tvs, ncol, icefrac, ocnfrac, lchnk )
    endif

    drydepflx(:,:) = 0._r8
    ! if chemUCI is used, apply dry dep flux here
    if (uci1_ndx <= 0) then 
       do m = 1,pcnst
          n = map2chm( m )
          if ( n > 0 ) then
            cflx(:ncol,m)      = cflx(:ncol,m) - sflx(:ncol,n)
            drydepflx(:ncol,m) = sflx(:ncol,n)
          endif
       end do
    else
       sflx2 = sflx
       do m = 1,pcnst
          n = map2chm( m )
          if ( n > 0 ) then
            drydepflx(:ncol,m) = sflx(:ncol,n)
            if ( any( sflx2(:ncol,n) /= 0._r8 ) ) then
              if ( .not. any( aer_species == n ) .and. adv_mass(n) /= 0._r8 ) then
                do k = pver, 1, -1 ! loop from bottom to top
                  ! kg/m2, tracer mass
                  wrk(:ncol,k) = adv_mass(n)*vmr(:ncol,k,n)/mbar(:ncol,k) &
                                    *pdeldry(:ncol,k)*rga
                  j = 0 ! number of columns w/o enough mass to remove at the current layer
                  do i = 1,ncol
                    if ( sflx2(i,n) /= 0._r8 ) then
                      if ( wrk(i,k)*(1._r8-low_limit) >= sflx2(i,n)*delt ) then
                        tmp(1) = wrk(i,k) - sflx2(i,n)*delt
                        sflx2(i,n) = 0._r8
                      else
                        tmp(1) = wrk(i,k)*low_limit
                        sflx2(i,n) = sflx2(i,n) - wrk(i,k)*(1._r8-low_limit)*delt_inverse
                        j = j + 1
                      endif
                      vmr(i,k,n) = tmp(1)*mbar(i,k)/adv_mass(n)/pdeldry(i,k)/rga
                    endif
                  end do

                  if ( j == 0 ) then
                    exit
                  endif

                end do
              endif
            endif
          endif
       end do
    endif
    call t_stopf('drydep')

    if ( history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels) then
       if ( history_gaschmbudget ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'TDD' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, 'MSD' )
       endif
       if ( history_gaschmbudget_2D ) then
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DTDD' )
          call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                 pdeldry(:ncol,:), mbar, delt_inverse, '2DMSD' )
       endif
       if ( history_gaschmbudget_2D_levels ) then
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDD_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DTDD_trop',tropFlagInt )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMSD_LL' )
         call gaschmmass_diags( lchnk, ncol, vmr(:ncol,:,:), vmr_old(:ncol,:,:), &
                                pdeldry(:ncol,:), mbar, delt_inverse, '2DMSD_trop',tropFlagInt )
       endif
       do m = 1,pcnst
          n = map2chm(m)
          if (n > 0 .and. (.not. any( aer_species == n ))) then
            ! store gas chemistry tracer concentration in pbuf for tendency
            ! calculation
            if (history_gaschmbudget .or. history_gaschmbudget_2D_levels) then
               gas_ac_idx = pbuf_get_index(gas_ac_name(n))
               call pbuf_get_field(pbuf, gas_ac_idx, gas_ac )
               ftem(:ncol,:) = adv_mass(n)*vmr(:ncol,:,n)/mbar(:ncol,:) &
                               *pdeldry(:ncol,:)*rga
               gas_ac(:ncol,:) = ftem(:ncol,:)
            endif
            if (history_gaschmbudget_2D) then
               gas_ac_idx = pbuf_get_index(gas_ac_name_2D(n))
               call pbuf_get_field(pbuf, gas_ac_idx, gas_ac_2D )
               ftem(:ncol,:) = adv_mass(n)*vmr(:ncol,:,n)/mbar(:ncol,:) &
                               *pdeldry(:ncol,:)*rga
               do k=2,pver
                 ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
               end do
               gas_ac_2D(:ncol) = ftem(:ncol,1)
            endif
          endif
       end do

    endif

    !-----------------------------------------------------------------------      
    !         ... Xform from vmr to mmr
    !-----------------------------------------------------------------------      
    call vmr2mmr( vmr, mmr_tend, mbar, ncol )

    call set_short_lived_species( mmr_tend, lchnk, ncol, pbuf )

    !-----------------------------------------------------------------------      
    !         ... Form the tendencies
    !----------------------------------------------------------------------- 
    do m = 1,gas_pcnst 
       mmr_new(:ncol,:,m) = mmr_tend(:ncol,:,m)
       mmr_tend(:ncol,:,m) = (mmr_tend(:ncol,:,m) - mmr(:ncol,:,m))*delt_inverse
    enddo

    do m = 1,pcnst
       n = map2chm(m)
       if( n > 0 ) then
          qtend(:ncol,:,m) = qtend(:ncol,:,m) + mmr_tend(:ncol,:,n) 
       end if
    end do

    call t_startf('chemdr_diags')
    call chm_diags( lchnk, ncol, vmr(:ncol,:,:), mmr_new(:ncol,:,:), &
                    reaction_rates(:ncol,:,:), invariants(:ncol,:,:), depvel(:ncol,:),  sflx(:ncol,:), &
                    mmr_tend(:ncol,:,:), pdel(:ncol,:), pdeldry(:ncol,:), pbuf, troplev(:ncol), tropFlag=tropFlag )

    call rate_diags_calc( reaction_rates(:,:,:), vmr(:,:,:), invariants(:,:,indexm), ncol, lchnk, &
                          pver, pdeldry(:ncol,:), mbar)

    call t_stopf('chemdr_diags')

  end subroutine gas_phase_chemdr

end module mo_gas_phase_chemdr
