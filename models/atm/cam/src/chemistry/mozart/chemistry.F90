module chemistry

!---------------------------------------------------------------------------------
! "Interactive" gas phase module
!---------------------------------------------------------------------------------

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use ppgrid,           only : pcols, pver, pverp, begchunk, endchunk
  use physconst,        only : gravit
  use constituents,     only : pcnst, cnst_add, cnst_name, sflxnam, cnst_fixed_ubc
  use chem_mods,        only : gas_pcnst, nfs
  use cam_history,      only : fieldname_len
  use physics_types,    only : physics_state, physics_ptend, physics_ptend_init
  use m_types,          only : time_ramp
  use dust_intr,        only : dust_names
  use progseasalts_intr,only : progseasalts_names
  use spmd_utils,       only : masterproc
  use cam_logfile,      only : iulog
  use mo_gas_phase_chemdr, only : map2chm
  use shr_megan_mod,    only : shr_megan_mechcomps, shr_megan_mechcomps_n 

  use tracer_data,  only : MAXTRCRS

  implicit none
  private
  save

!---------------------------------------------------------------------------------
! Public interfaces
!---------------------------------------------------------------------------------
  public :: chem_is                        ! identify which chemistry is being used
  public :: chem_register                  ! register consituents
  public :: chem_readnl                    ! read chem namelist 
  public :: chem_is_active                 ! returns true
  public :: chem_implements_cnst           ! returns true if consituent is implemented by this package
  public :: chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public :: chem_init                      ! initialize (history) variables
  public :: chem_timestep_init             ! per timestep initializations
  public :: chem_timestep_tend             ! interface to tendency computation
  public :: chem_final
  public :: chem_write_restart
  public :: chem_read_restart
  public :: chem_init_restart
  public :: chem_reset_fluxes
  
  integer, public :: imozart = -1       ! index of 1st constituent
  
  ! Namelist variables
  
  ! control
  
  integer :: chem_freq = 1 ! time steps

  ! ghg

  character(len=shr_kind_cl) :: bndtvg = ' ' ! pathname for greenhouse gas loss rate
  character(len=shr_kind_cl) :: h2orates = ' ' ! pathname for greenhouse gas (lyman-alpha H2O loss)

  ! lightning

  real(r8)           :: lght_no_prd_factor = 1._r8

  ! photolysis

  logical            :: xactive_prates = .false.
  character(len=shr_kind_cl) :: rsf_file = 'rsf_file'
  character(len=shr_kind_cl) :: exo_coldens_file = ''
  character(len=shr_kind_cl) :: tuv_xsect_file = 'tuv_xsect_file'
  character(len=shr_kind_cl) :: o2_xsect_file = 'o2_xsect_file'
  character(len=shr_kind_cl) :: xs_coef_file = 'xs_coef_file'
  character(len=shr_kind_cl) :: xs_short_file = 'xs_short_file'
  character(len=shr_kind_cl) :: xs_long_file = 'xs_long_file'
  character(len=shr_kind_cl) :: electron_file = 'electron_file'
  character(len=shr_kind_cl) :: euvac_file = 'euvac_file'
  character(len=shr_kind_cl) :: euvacdat_file = 'euvacdat_file'

  ! solar / geomag data

  character(len=shr_kind_cl) :: solar_parms_file = ' '     ! solar variability parameters
  character(len=shr_kind_cl) :: photon_file = 'photon_file'

  ! dry dep
  
  character(len=shr_kind_cl) :: depvel_file = 'depvel_file'
  character(len=shr_kind_cl) :: depvel_lnd_file = 'depvel_lnd_file'
  character(len=shr_kind_cl) :: clim_soilw_file = 'clim_soilw_file'
  character(len=shr_kind_cl) :: season_wes_file = 'season_wes_file'

  ! emis

  character(len=shr_kind_cl) :: airpl_emis_file = '' ! airplane emissions
  character(len=shr_kind_cl) :: srf_emis_specifier(pcnst) = ''
  character(len=shr_kind_cl) :: ext_frc_specifier(pcnst) = ''

  character(len=24)  :: srf_emis_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' |  'INTERP_MISSING_MONTHS'
  integer            :: srf_emis_cycle_yr  = 0
  integer            :: srf_emis_fixed_ymd = 0
  integer            :: srf_emis_fixed_tod = 0

  character(len=24)  :: ext_frc_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' |  'INTERP_MISSING_MONTHS'
  integer            :: ext_frc_cycle_yr  = 0
  integer            :: ext_frc_fixed_ymd = 0
  integer            :: ext_frc_fixed_tod = 0

  ! fixed stratosphere
  
  character(len=shr_kind_cl) :: fstrat_file = 'fstrat_file'
  character(len=16)  :: fstrat_list(pcnst)  = ''

  ! stratospheric aerosols

  character(len=shr_kind_cl) :: sad_file = 'sad_file'

  ! trop sulf

  character(len=shr_kind_cl) :: sulf_file = 'sulf_file'

  type(time_ramp)    :: sad_timing      != time_ramp( "CYCLICAL",  19970101, 0 )

! for linoz
  character(len=shr_kind_cl) :: chlorine_loading_file = ''
  character(len=8)   :: chlorine_loading_type = 'SERIAL' ! "FIXED" or "SERIAL"
  integer            :: chlorine_loading_fixed_ymd = 0         ! YYYYMMDD for "FIXED" type
  integer            :: chlorine_loading_fixed_tod = 0         ! seconds of day for "FIXED" type

!---------------------------------------------------------------------------------
! dummy values for specific heats at constant pressure
!---------------------------------------------------------------------------------
  real(r8), parameter   :: cptmp = 666._r8

  character(len=fieldname_len) :: srcnam(gas_pcnst) ! names of source/sink tendencies

  integer :: ixcldliq, ixcldice                     ! indicies of liquid and ice cloud water
  integer :: ndx_cld
  integer :: ndx_cmfdqr
  integer :: ndx_nevapr
  integer :: ndx_prain
  integer :: ndx_cldtop
  integer :: h2o_ndx
  integer :: ixndrop             ! cloud droplet number index
  integer :: ndx_pblh

  logical :: ghg_chem = .false.      ! .true. => use ghg chem package
  logical :: chem_step = .true.
  logical :: is_active = .false.

  character(len=32) :: chem_name = 'UNSET'
  logical :: chem_rad_passive = .false.
  
  ! for MEGAN emissions
  integer, allocatable :: megan_indices_map(:) 
  real(r8),allocatable :: megan_wght_factors(:)

  ! dust constituent indices
  integer :: dust_indices(8)

!================================================================================================
contains
!================================================================================================

logical function chem_is (name)
   use phys_control,     only : cam_chempkg_is

   character(len=*), intent(in) :: name
   chem_is = cam_chempkg_is(name)

end function chem_is

!================================================================================================

  subroutine chem_register
!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents and physics buffer fields
! 
!-----------------------------------------------------------------------

    use ioFileMod,           only : getfil
    use mo_sim_dat,          only : set_sim_dat
    use chem_mods,           only : gas_pcnst, adv_mass, inv_lst, nfs, fix_mass
    use mo_tracname,         only : solsym
    use mo_chem_utls,        only : get_spc_ndx
    use progseasalts_intr,   only : progseasalts_set_idx
    use dust_intr,           only : dust_set_idx
    use short_lived_species, only : slvd_index, short_lived_map=>map, register_short_lived_species
    use cfc11star,           only : register_cfc11star
#if ( defined MODAL_AERO)
    use modal_aero_initialize_data, only : modal_aero_register
#endif
    use phys_control,        only : waccmx_is   ! WACCM-X switch query function
    use mo_photo,            only : photo_register
    use mo_aurora,           only : aurora_register
    use mo_setsoa,           only : soa_register

    implicit none

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer  :: idx                                 ! Index for adding fields to pbuf physics buffer
    integer  :: m, n                                ! tracer index
    real(r8) :: qmin                                ! min value
    logical  :: ic_from_cam2                        ! wrk variable for initial cond input
    logical  :: has_fixed_ubc                       ! wrk variable for upper bndy cond
    logical  :: has_fixed_ubflx                     ! wrk variable for upper bndy flux
    character(len=shr_kind_cl) :: locfn                     ! chemistry data filespec
    integer  :: ch4_ndx, n2o_ndx, o3_ndx 
    integer  :: co2_ndx, cfc11_ndx, cfc12_ndx, o2_1s_ndx, o2_1d_ndx, o2_ndx
    integer  :: n_ndx, no_ndx, co_ndx, h_ndx, h2_ndx, o_ndx, e_ndx, np_ndx
    integer  :: op_ndx, o1d_ndx, n2d_ndx, nop_ndx, n2p_ndx, o2p_ndx
    integer  :: cly_ndx,cl_ndx,clo_ndx,hocl_ndx,cl2_ndx,cl2o2_ndx,oclo_ndx,hcl_ndx,clono2_ndx
    integer  :: bry_ndx,br_ndx,bro_ndx,hbr_ndx,brono2_ndx,brcl_ndx,hobr_ndx

    character(len=128) :: lng_name                  ! variable long name
    logical :: cam_outfld
    character(len=128) :: mixtype
    character(len=128) :: molectype
    integer :: islvd

!-----------------------------------------------------------------------
! Set the simulation chemistry variables
!-----------------------------------------------------------------------
    call set_sim_dat

    o3_ndx    = get_spc_ndx('O3')
    ch4_ndx   = get_spc_ndx('CH4')
    n2o_ndx   = get_spc_ndx('N2O')

    co2_ndx   = get_spc_ndx('CO2')
    cfc11_ndx = get_spc_ndx('CFC11')
    cfc12_ndx = get_spc_ndx('CFC12')
    o2_1s_ndx = get_spc_ndx('O2_1S')
    o2_1d_ndx = get_spc_ndx('O2_1D')
    o2_ndx    = get_spc_ndx('O2')
    n_ndx     = get_spc_ndx('N')
    no_ndx    = get_spc_ndx('NO')
    co_ndx    = get_spc_ndx('CO')
    h_ndx     = get_spc_ndx('H')
    h2_ndx    = get_spc_ndx('H2')
    o_ndx     = get_spc_ndx('O')
    e_ndx     = get_spc_ndx('e')
    np_ndx    = get_spc_ndx('Np')
    op_ndx    = get_spc_ndx('Op')
    o1d_ndx   = get_spc_ndx('O1D_')
    n2d_ndx   = get_spc_ndx('N2D')
    n2p_ndx   = get_spc_ndx('N2p')
    nop_ndx   = get_spc_ndx('NOp')
    h2o_ndx   = get_spc_ndx('H2O')
    o2p_ndx   = get_spc_ndx('O2p')

    cly_ndx   = get_spc_ndx('CLY')
    cl_ndx    = get_spc_ndx('CL')
    clo_ndx   = get_spc_ndx('CLO')
    hocl_ndx  = get_spc_ndx('HOCL')
    cl2_ndx   = get_spc_ndx('CL2')
    cl2o2_ndx = get_spc_ndx('CL2O2')
    oclo_ndx  = get_spc_ndx('OCLO')
    hcl_ndx   = get_spc_ndx('HCL')
    clono2_ndx= get_spc_ndx('CLONO2')
    bry_ndx   = get_spc_ndx('BRY')
    br_ndx    = get_spc_ndx('BR')
    bro_ndx   = get_spc_ndx('BRO')
    hbr_ndx   = get_spc_ndx('HBR')
    brono2_ndx= get_spc_ndx('BRONO2')
    brcl_ndx  = get_spc_ndx('BRCL')
    hobr_ndx  = get_spc_ndx('HOBR')

    !-----------------------------------------------------------------------
    ! Set names of diffused variable tendencies and declare them as history variables
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    ! BAB: 2004-09-01 kludge to define a fixed ubc for water vapor
    !      required because water vapor is not declared by chemistry but
    !      has a fixed ubc only if chemistry is running.
    !-----------------------------------------------------------------------
    if (chem_is('waccm_mozart') .or. chem_is('waccm_mozart_mam3')) then

    !-----------------------------------------------------------------------------------------
    !For WACCM-X, change variable cnst_fixed_ubc(1) from .true. to .false. which eliminates 
    !the constant fixed upper boundary condition for all species
    !-----------------------------------------------------------------------------------------
      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
        cnst_fixed_ubc(1) = .false.
      else
        cnst_fixed_ubc(1) = .true.
      endif
      
    endif
    
    !-----------------------------------------------------------------------
    ! Set names of diffused variable tendencies and declare them as history variables
    !-----------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    ! For WACCM-X, change variable has_fixed_ubc from .true. to .false. which is a flag
    ! used later to check for a fixed upper boundary condition for species. 
    !----------------------------------------------------------------------------------
     do m = 1,gas_pcnst
       ic_from_cam2  = .true.
       has_fixed_ubc = .false.
       has_fixed_ubflx = .false.
       lng_name      = trim( solsym(m) )
       molectype = 'minor'

       qmin = 1.e-36_r8

       if ( m == o3_ndx ) then
          qmin          = 1.e-12_r8
       else if ( m == ch4_ndx ) then
          qmin          = 1.e-12_r8
          if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
            has_fixed_ubc = .false.   ! diffusive equilibrium at UB
          else
            has_fixed_ubc = .true.
          endif
       else if ( m == n2o_ndx .or. m == co2_ndx ) then
          qmin = 1.e-15_r8
          if( m == co2_ndx ) then
            if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
              has_fixed_ubc = .false.   ! diffusive equilibrium at UB
            else
              has_fixed_ubc = .true.
            endif
          end if
       else if( m == cfc11_ndx .or. m == cfc12_ndx ) then
          qmin = 1.e-20_r8
       else if( m == o2_1s_ndx .or. m == o2_1d_ndx ) then
          ic_from_cam2 = .false.
          if( m == o2_1d_ndx ) then
             lng_name = 'O2(1-delta)'
          else
             lng_name = 'O2(1-sigma)'
          end if
       else if( m==o2_ndx .or. m==n_ndx .or. m==no_ndx .or. m==co_ndx .or. m==h_ndx .or. m==h2_ndx .or. m==o_ndx ) then
         if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
           has_fixed_ubc = .false.   ! diffusive equilibrium at UB
           if ( m == h_ndx ) has_fixed_ubflx = .true. ! fixed flux value for H at UB
           if ( m == o2_ndx .or. m == o_ndx ) molectype = 'major'
         else
           has_fixed_ubc = .true.
         endif
       else if( m == e_ndx ) then
          lng_name = 'electron concentration'
       else if( m == np_ndx ) then
          lng_name = 'N+'
       else if( m == op_ndx ) then
          lng_name = 'O+'
       else if( m == o1d_ndx ) then
          lng_name = 'O(1D)'
       else if( m == n2d_ndx ) then
          lng_name = 'N(2D)'
       else if( m == o2p_ndx ) then
          lng_name = 'O2+'
       else if( m == n2p_ndx ) then
          lng_name = 'N2+'
       else if( m == nop_ndx ) then
          lng_name = 'NO+'
       else if( m == h2o_ndx ) then
          map2chm(1) = m
          cycle
       endif

       if ( any(solsym(m) == dust_names) ) then
          cam_outfld=.true.
       else if ( any(solsym(m) == progseasalts_names) ) then
          cam_outfld=.true.
       else
          cam_outfld=.false.
          is_active = .true.
       endif

       if ( chem_is('waccm_mozart').or.chem_is('waccm_mozart_mam3').or. &
            chem_is('waccm_ghg').or.chem_is('super_fast_llnl') ) then
          mixtype='wet'
       else
          mixtype='dry'
       endif

       islvd = slvd_index(solsym(m))

       if ( islvd > 0 ) then
          short_lived_map(islvd) = m
       else
          call cnst_add( solsym(m), adv_mass(m), cptmp, qmin, n, readiv=ic_from_cam2, cam_outfld=cam_outfld, &
               mixtype=mixtype, molectype=molectype, fixed_ubc=has_fixed_ubc, fixed_ubflx=has_fixed_ubflx, longname=trim(lng_name) )

          if ( solsym(m) == dust_names(1)) call dust_set_idx(n)
          if ( solsym(m) == progseasalts_names(1)) call progseasalts_set_idx(n)

          if( imozart == -1 ) then
             imozart = n
          end if
          map2chm(n) = m
       endif

    end do

    call register_short_lived_species()
    call soa_register()
    if ( chem_is('waccm_mozart') .or. chem_is('waccm_mozart_mam3') ) then
       call register_cfc11star()
    endif
#if ( defined MODAL_AERO)
    call modal_aero_register()
#endif
    
    if ( waccmx_is('ionosphere') ) then 
       call photo_register()
       call aurora_register()
    endif

  end subroutine chem_register

!================================================================================================
  
  subroutine chem_readnl(nlfile)

    ! Read chem namelist group.

    use abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    use mz_aerosols_intr, only: mz_aero_defaultopts, mz_aero_setopts
    use linoz_data,       only: linoz_data_defaultopts,  linoz_data_setopts
    use tracer_cnst,      only: tracer_cnst_defaultopts, tracer_cnst_setopts
    use tracer_srcs,      only: tracer_srcs_defaultopts, tracer_srcs_setopts
    use aerosol_intr,     only: aerosol_readnl
    use gas_wetdep_opts,  only: gas_wetdep_readnl

#ifdef WACCM_MOZART
    use spedata,          only: spedata_defaultopts, spedata_setopts
#endif
    use mo_sad,           only: sad_defaultopts, sad_setopts

    use upper_bc,         only: ubc_defaultopts, ubc_setopts

    use mo_drydep,       only: drydep_srf_file


   ! args

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    character(len=8)   :: sad_type = 'CYCLICAL'      ! 'CYCLICAL' | 'SERIAL' | 'FIXED'
    integer            :: sad_cycle_yr  = 0
    integer            :: sad_fixed_ymd = 0
    integer            :: sad_fixed_tod = 0

    logical            :: use_cam_sulfchem   

    real(r8)           :: sol_facti_cloud_borne

    character(len=16), dimension(pcnst) :: aer_wetdep_list = ' '

    ! linoz data
    character(len=shr_kind_cl) :: linoz_data_file               ! prescribed data file
    character(len=shr_kind_cl) :: linoz_data_filelist           ! list of prescribed data files (series of files)
    character(len=shr_kind_cl) :: linoz_data_path               ! absolute path of prescribed data files 
    character(len=24)  :: linoz_data_type               ! 'INTERP_MISSING_MONTHS' | 'CYCLICAL' | 'SERIAL' (default)
    logical            :: linoz_data_rmfile             ! remove data file from local disk (default .false.)
    integer            :: linoz_data_cycle_yr
    integer            :: linoz_data_fixed_ymd
    integer            :: linoz_data_fixed_tod

    ! trop_mozart prescribed constituent concentratons
    character(len=shr_kind_cl) :: tracer_cnst_file              ! prescribed data file
    character(len=shr_kind_cl) :: tracer_cnst_filelist          ! list of prescribed data files (series of files)
    character(len=shr_kind_cl) :: tracer_cnst_datapath          ! absolute path of prescribed data files 
    character(len=24)  :: tracer_cnst_type              ! 'INTERP_MISSING_MONTHS' | 'CYCLICAL' | 'SERIAL' (default)
    character(len=shr_kind_cl) :: tracer_cnst_specifier(MAXTRCRS) ! string array where each 
    logical            :: tracer_cnst_rmfile            ! remove data file from local disk (default .false.)
    integer            :: tracer_cnst_cycle_yr
    integer            :: tracer_cnst_fixed_ymd
    integer            :: tracer_cnst_fixed_tod

    ! trop_mozart prescribed constituent sourrces/sinks
    character(len=shr_kind_cl) :: tracer_srcs_file              ! prescribed data file
    character(len=shr_kind_cl) :: tracer_srcs_filelist          ! list of prescribed data files (series of files)
    character(len=shr_kind_cl) :: tracer_srcs_datapath          ! absolute path of prescribed data files 
    character(len=24)  :: tracer_srcs_type              ! 'INTERP_MISSING_MONTHS' | 'CYCLICAL' | 'SERIAL' (default)
    character(len=shr_kind_cl) :: tracer_srcs_specifier(MAXTRCRS) ! string array where each 
    logical            :: tracer_srcs_rmfile            ! remove data file from local disk (default .false.)
    integer            :: tracer_srcs_cycle_yr
    integer            :: tracer_srcs_fixed_ymd
    integer            :: tracer_srcs_fixed_tod

#if ( defined WACCM_MOZART || defined WACCM_GHG )
    ! Upper boundary conditions
    character(len=shr_kind_cl) :: tgcm_ubc_file
    integer            :: tgcm_ubc_cycle_yr
    integer            :: tgcm_ubc_fixed_ymd
    integer            :: tgcm_ubc_fixed_tod
    character(len=32)  :: tgcm_ubc_data_type
    character(len=shr_kind_cl) :: snoe_ubc_file
    ! Upper boundary conditions
    real(r8)           :: t_pert_ubc   ! temperature perturbation at ubc
    real(r8)           :: no_xfac_ubc  ! no multiplicative factor at ubc
#endif

#ifdef WACCM_MOZART
    ! waccm solor proton data variables
    logical            :: spe_remove_file     ! true => the offline spe file will be removed
    character(len=shr_kind_cl) :: spe_data_file       ! name of file that contains the spe data
    character(len=shr_kind_cl) :: spe_filenames_list  ! file that lists a series of spe files 
#endif
    logical            :: strat_aero_feedback ! true => radiation feed backs from strat sulfur aerosol 
#if ( defined MODAL_AERO)
    ! trop_mozart aerosol constituents that have dry deposition
    character(len=16)  :: aer_drydep_list(pcnst)
#endif

    namelist /chem_inparm/ chem_freq, airpl_emis_file, &
         solar_parms_file, euvac_file, &
         euvacdat_file, photon_file, electron_file, &
         sad_type, sad_cycle_yr, sad_fixed_ymd, sad_fixed_tod, sad_file, &
         sulf_file, depvel_file, xs_coef_file, xs_short_file, &
         exo_coldens_file, tuv_xsect_file, o2_xsect_file, &
         xs_long_file, rsf_file, &
         lght_no_prd_factor, xactive_prates, &
         depvel_lnd_file, clim_soilw_file, season_wes_file, drydep_srf_file, &
         srf_emis_type, srf_emis_cycle_yr, srf_emis_fixed_ymd, srf_emis_fixed_tod, srf_emis_specifier,  &
         fstrat_file, fstrat_list, &
         ext_frc_specifier, ext_frc_type, ext_frc_cycle_yr, ext_frc_fixed_ymd, ext_frc_fixed_tod, &
         aer_wetdep_list, use_cam_sulfchem, sol_facti_cloud_borne

    namelist /chem_inparm/ chem_rad_passive

    ! ghg chem

    namelist /chem_inparm/ bndtvg, h2orates, ghg_chem

    ! linoz inputs

    namelist /chem_inparm/ &
         linoz_data_file, linoz_data_filelist, linoz_data_path, &
         linoz_data_type, &
         linoz_data_rmfile, linoz_data_cycle_yr, linoz_data_fixed_ymd, linoz_data_fixed_tod
    namelist /chem_inparm/ &
         chlorine_loading_file, chlorine_loading_type, chlorine_loading_fixed_ymd, chlorine_loading_fixed_tod

    ! prescribed chem tracers

    namelist /chem_inparm/ &
         tracer_cnst_file, tracer_cnst_filelist, tracer_cnst_datapath, &
         tracer_cnst_type, tracer_cnst_specifier, &
         tracer_srcs_file, tracer_srcs_filelist, tracer_srcs_datapath, &
         tracer_srcs_type, tracer_srcs_specifier, &
         tracer_cnst_rmfile, tracer_cnst_cycle_yr, tracer_cnst_fixed_ymd, tracer_cnst_fixed_tod, &
         tracer_srcs_rmfile, tracer_srcs_cycle_yr, tracer_srcs_fixed_ymd, tracer_srcs_fixed_tod 

#if ( defined WACCM_MOZART || defined WACCM_GHG )
    ! upper boundary conditions
    namelist /chem_inparm/ tgcm_ubc_file, tgcm_ubc_data_type, tgcm_ubc_cycle_yr, tgcm_ubc_fixed_ymd, tgcm_ubc_fixed_tod, &
                           snoe_ubc_file, t_pert_ubc, no_xfac_ubc
#endif

#ifdef WACCM_MOZART
    ! waccm solar proton namelist variables
    namelist /chem_inparm/ spe_data_file, spe_remove_file, spe_filenames_list
#endif
    ! radiation feed back from stratospheric sulfur aerosols 
    namelist /chem_inparm/ strat_aero_feedback
#if ( defined MODAL_AERO)
    namelist /chem_inparm/ aer_drydep_list
#endif

    ! get the default settings

#if ( defined MODAL_AERO)
    call mz_aero_defaultopts( sol_facti_cloud_borne, aer_wetdep_list_out = aer_wetdep_list, &
         aer_drydep_list_out = aer_drydep_list, use_cam_sulfchem_out = use_cam_sulfchem )
#else
    call mz_aero_defaultopts( sol_facti_cloud_borne, aer_wetdep_list_out = aer_wetdep_list, &
         use_cam_sulfchem_out = use_cam_sulfchem )
#endif
    call linoz_data_defaultopts( &
         linoz_data_file_out      = linoz_data_file,      &
         linoz_data_filelist_out  = linoz_data_filelist,  &
         linoz_data_path_out      = linoz_data_path,      &
         linoz_data_type_out      = linoz_data_type,      &
         linoz_data_rmfile_out    = linoz_data_rmfile,    &
         linoz_data_cycle_yr_out  = linoz_data_cycle_yr,  &
         linoz_data_fixed_ymd_out = linoz_data_fixed_ymd, &
         linoz_data_fixed_tod_out = linoz_data_fixed_tod  ) 
    call tracer_cnst_defaultopts( &
         tracer_cnst_file_out      = tracer_cnst_file,      &
         tracer_cnst_filelist_out  = tracer_cnst_filelist,  &
         tracer_cnst_datapath_out  = tracer_cnst_datapath,  &
         tracer_cnst_type_out      = tracer_cnst_type,      &
         tracer_cnst_specifier_out = tracer_cnst_specifier, &
         tracer_cnst_rmfile_out    = tracer_cnst_rmfile,    &
         tracer_cnst_cycle_yr_out  = tracer_cnst_cycle_yr,  &
         tracer_cnst_fixed_ymd_out = tracer_cnst_fixed_ymd, &
         tracer_cnst_fixed_tod_out = tracer_cnst_fixed_tod  ) 
    call tracer_srcs_defaultopts( &
         tracer_srcs_file_out      = tracer_srcs_file,      &
         tracer_srcs_filelist_out  = tracer_srcs_filelist,  &
         tracer_srcs_datapath_out  = tracer_srcs_datapath,  &
         tracer_srcs_type_out      = tracer_srcs_type,      &
         tracer_srcs_specifier_out = tracer_srcs_specifier, &
         tracer_srcs_rmfile_out    = tracer_srcs_rmfile,    &
         tracer_srcs_cycle_yr_out  = tracer_srcs_cycle_yr,  &
         tracer_srcs_fixed_ymd_out = tracer_srcs_fixed_ymd, &
         tracer_srcs_fixed_tod_out = tracer_srcs_fixed_tod  )

#ifdef WACCM_MOZART
    ! spedata
    call spedata_defaultopts( spe_data_file_out = spe_data_file, &
         spe_remove_file_out = spe_remove_file, &
         spe_filenames_list_out = spe_filenames_list )
#endif
    call sad_defaultopts( strat_aero_feedback_out = strat_aero_feedback )

#if ( defined WACCM_MOZART || defined WACCM_GHG )
    ! Upper boundary conditions
    call ubc_defaultopts( &
         snoe_ubc_file_out =snoe_ubc_file, &
#ifdef WACCM_MOZART
         t_pert_ubc_out    =t_pert_ubc, &
         no_xfac_ubc_out   =no_xfac_ubc, &
#endif
         tgcm_ubc_file_out      = tgcm_ubc_file, &
         tgcm_ubc_data_type_out = tgcm_ubc_data_type, &
         tgcm_ubc_cycle_yr_out  = tgcm_ubc_cycle_yr, &
         tgcm_ubc_fixed_ymd_out = tgcm_ubc_fixed_ymd, &
         tgcm_ubc_fixed_tod_out = tgcm_ubc_fixed_tod )
#endif

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'chem_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, chem_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun('chem_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables

    ! control

    call mpibcast (chem_freq,         1,                               mpiint,  0, mpicom)

    call mpibcast (use_cam_sulfchem,  1,                               mpilog,  0, mpicom)
    call mpibcast (sol_facti_cloud_borne,  1,                          mpir8,   0, mpicom)
    call mpibcast (aer_wetdep_list,   len(aer_wetdep_list(1))*pcnst,   mpichar, 0, mpicom)
#if ( defined MODAL_AERO)
    call mpibcast (aer_drydep_list,   len(aer_drydep_list(1))*pcnst,   mpichar, 0, mpicom)
#endif

    call mpibcast (chem_rad_passive,  1,                               mpilog,  0, mpicom)

    ! ghg

    call mpibcast (ghg_chem,          1,                               mpilog,  0, mpicom)
    call mpibcast (bndtvg,            len(bndtvg),                     mpichar, 0, mpicom)
    call mpibcast (h2orates,          len(h2orates),                   mpichar, 0, mpicom)

    ! lightning

    call mpibcast (lght_no_prd_factor,1,                               mpir8,   0, mpicom)

    ! photolysis

    call mpibcast (rsf_file,          len(rsf_file),                   mpichar, 0, mpicom)
    call mpibcast (exo_coldens_file,  len(exo_coldens_file),           mpichar, 0, mpicom)
    call mpibcast (tuv_xsect_file,    len(tuv_xsect_file),             mpichar, 0, mpicom)
    call mpibcast (o2_xsect_file,     len(o2_xsect_file),              mpichar, 0, mpicom)
    call mpibcast (xs_coef_file,      len(xs_coef_file),               mpichar, 0, mpicom)
    call mpibcast (xs_short_file,     len(xs_short_file),              mpichar, 0, mpicom)
    call mpibcast (xs_long_file,      len(xs_long_file),               mpichar, 0, mpicom)
    call mpibcast (xactive_prates,    1,                               mpilog,  0, mpicom)
    call mpibcast (electron_file,     len(electron_file),              mpichar, 0, mpicom)
    call mpibcast (euvac_file,        len(euvac_file),                 mpichar, 0, mpicom)
    call mpibcast (euvacdat_file,     len(euvacdat_file),              mpichar, 0, mpicom)

    ! solar / geomag data

    call mpibcast (solar_parms_file,  len(solar_parms_file),           mpichar, 0, mpicom)
    call mpibcast (photon_file,       len(photon_file),                mpichar, 0, mpicom)

    ! dry dep

    call mpibcast (depvel_lnd_file,   len(depvel_lnd_file),            mpichar, 0, mpicom)
    call mpibcast (depvel_file,       len(depvel_file),                mpichar, 0, mpicom)
    call mpibcast (clim_soilw_file,   len(clim_soilw_file),            mpichar, 0, mpicom)
    call mpibcast (season_wes_file,   len(season_wes_file),            mpichar, 0, mpicom)
    call mpibcast (drydep_srf_file,   len(drydep_srf_file),            mpichar, 0, mpicom)

    ! emis

    call mpibcast (airpl_emis_file,   len(airpl_emis_file),            mpichar, 0, mpicom)
    call mpibcast (srf_emis_specifier,len(srf_emis_specifier(1))*pcnst,mpichar, 0, mpicom)
    call mpibcast (srf_emis_type,     len(srf_emis_type),              mpichar, 0, mpicom)
    call mpibcast (srf_emis_cycle_yr, 1,                               mpiint,  0, mpicom)
    call mpibcast (srf_emis_fixed_ymd,1,                               mpiint,  0, mpicom)
    call mpibcast (srf_emis_fixed_tod,1,                               mpiint,  0, mpicom)
    call mpibcast (ext_frc_specifier, len(ext_frc_specifier(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast (ext_frc_type,      len(ext_frc_type),               mpichar, 0, mpicom)
    call mpibcast (ext_frc_cycle_yr,  1,                               mpiint,  0, mpicom)
    call mpibcast (ext_frc_fixed_ymd, 1,                               mpiint,  0, mpicom)
    call mpibcast (ext_frc_fixed_tod, 1,                               mpiint,  0, mpicom)


    ! fixed stratosphere

    call mpibcast (fstrat_file,       len(fstrat_file),                mpichar, 0, mpicom)
    call mpibcast (fstrat_list,       len(fstrat_list(1))*pcnst,       mpichar, 0, mpicom)

    ! stratospheric aerosols

    call mpibcast (sad_file,          len(sad_file),                   mpichar, 0, mpicom)
    call mpibcast (sad_type,          len(sad_type),                   mpichar, 0, mpicom)
    call mpibcast (sad_cycle_yr,      1,                               mpiint,  0, mpicom)
    call mpibcast (sad_fixed_ymd,     1,                               mpiint,  0, mpicom)
    call mpibcast (sad_fixed_tod,     1,                               mpiint,  0, mpicom)

   ! trop sulf

    call mpibcast (sulf_file,         len(sulf_file),                  mpichar, 0, mpicom)

#if ( defined WACCM_MOZART || defined WACCM_GHG )
    ! upper boundary
    call mpibcast (tgcm_ubc_file,      len(tgcm_ubc_file),     mpichar, 0, mpicom)
    call mpibcast (tgcm_ubc_data_type, len(tgcm_ubc_data_type),mpichar, 0, mpicom)
    call mpibcast (tgcm_ubc_cycle_yr,  1,                      mpiint,  0, mpicom)
    call mpibcast (tgcm_ubc_fixed_ymd, 1,                      mpiint,  0, mpicom)
    call mpibcast (tgcm_ubc_fixed_tod, 1,                      mpiint,  0, mpicom)

    call mpibcast (snoe_ubc_file, len(snoe_ubc_file), mpichar, 0, mpicom)
    call mpibcast (t_pert_ubc,    1,                  mpir8,   0, mpicom)
    call mpibcast (no_xfac_ubc,   1,                  mpir8,   0, mpicom)
#endif
#ifdef WACCM_MOZART
   ! waccm solar proton variables
   call mpibcast (spe_data_file      ,len(spe_data_file)     ,mpichar,0,mpicom)
   call mpibcast (spe_remove_file    ,1                      ,mpilog, 0, mpicom )
   call mpibcast (spe_filenames_list ,len(spe_filenames_list),mpichar,0,mpicom)
#endif
   call mpibcast (strat_aero_feedback,1                      ,mpilog, 0, mpicom )

    ! linoz data

    call mpibcast (linoz_data_file,      len(linoz_data_file),                  mpichar, 0, mpicom)
    call mpibcast (linoz_data_filelist,  len(linoz_data_filelist),              mpichar, 0, mpicom)
    call mpibcast (linoz_data_path,      len(linoz_data_path),              mpichar, 0, mpicom)
    call mpibcast (linoz_data_type,      len(linoz_data_type),                  mpichar, 0, mpicom)
    call mpibcast (linoz_data_rmfile,    1,                                     mpilog, 0,  mpicom)
    call mpibcast (linoz_data_cycle_yr,       1,                                     mpiint, 0,  mpicom)
    call mpibcast (linoz_data_fixed_ymd,       1,                                     mpiint, 0,  mpicom)
    call mpibcast (linoz_data_fixed_tod,       1,                                     mpiint, 0,  mpicom)

    call mpibcast (chlorine_loading_file,len(chlorine_loading_file), mpichar, 0, mpicom)
    call mpibcast (chlorine_loading_type,len(chlorine_loading_type), mpichar, 0, mpicom)
    call mpibcast (chlorine_loading_fixed_ymd, 1,                    mpiint,  0, mpicom)
    call mpibcast (chlorine_loading_fixed_tod, 1,                    mpiint,  0, mpicom)

    ! prescribed chemical tracers

    call mpibcast (tracer_cnst_specifier, len(tracer_cnst_specifier(1))*MAXTRCRS, mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_file,      len(tracer_cnst_file),                  mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_filelist,  len(tracer_cnst_filelist),              mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_datapath,  len(tracer_cnst_datapath),              mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_type,      len(tracer_cnst_type),                  mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_rmfile,    1,                                      mpilog, 0,  mpicom)
    call mpibcast (tracer_cnst_cycle_yr,  1,                                      mpiint, 0,  mpicom)
    call mpibcast (tracer_cnst_fixed_ymd, 1,                                      mpiint, 0,  mpicom)
    call mpibcast (tracer_cnst_fixed_tod, 1,                                      mpiint, 0,  mpicom)

    call mpibcast (tracer_srcs_specifier, len(tracer_srcs_specifier(1))*MAXTRCRS, mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_file,      len(tracer_srcs_file),                  mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_filelist,  len(tracer_srcs_filelist),              mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_datapath,  len(tracer_srcs_datapath),              mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_type,      len(tracer_srcs_type),                  mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_rmfile,    1,                                      mpilog,  0, mpicom)
    call mpibcast (tracer_srcs_cycle_yr,  1,                                      mpiint,  0, mpicom)
    call mpibcast (tracer_srcs_fixed_ymd, 1,                                      mpiint,  0, mpicom)
    call mpibcast (tracer_srcs_fixed_tod, 1,                                      mpiint,  0, mpicom)

#endif

    sad_timing%type      = sad_type
    sad_timing%cycle_yr  = sad_cycle_yr
    sad_timing%fixed_ymd = sad_fixed_ymd
    sad_timing%fixed_tod = sad_fixed_tod

    ! set the options

#if ( defined MODAL_AERO)
   call mz_aero_setopts( sol_facti_cloud_borne, aer_wetdep_list_in = aer_wetdep_list, aer_drydep_list_in = aer_drydep_list, &
	use_cam_sulfchem_in = use_cam_sulfchem  )
#else
   call mz_aero_setopts( sol_facti_cloud_borne, aer_wetdep_list_in = aer_wetdep_list, use_cam_sulfchem_in = use_cam_sulfchem  )
#endif
   call linoz_data_setopts( &
        linoz_data_file_in      = linoz_data_file,      &
        linoz_data_filelist_in  = linoz_data_filelist,  &
        linoz_data_path_in      = linoz_data_path,      &
        linoz_data_type_in      = linoz_data_type,      &
        linoz_data_rmfile_in    = linoz_data_rmfile,    &
        linoz_data_cycle_yr_in  = linoz_data_cycle_yr,  &
        linoz_data_fixed_ymd_in = linoz_data_fixed_ymd, &
        linoz_data_fixed_tod_in = linoz_data_fixed_tod )
   call tracer_cnst_setopts( &
        tracer_cnst_file_in      = tracer_cnst_file,      &
        tracer_cnst_filelist_in  = tracer_cnst_filelist,  &
        tracer_cnst_datapath_in  = tracer_cnst_datapath,  &
        tracer_cnst_type_in      = tracer_cnst_type,      &
        tracer_cnst_specifier_in = tracer_cnst_specifier, &
        tracer_cnst_rmfile_in    = tracer_cnst_rmfile,    &
        tracer_cnst_cycle_yr_in  = tracer_cnst_cycle_yr,  &
        tracer_cnst_fixed_ymd_in = tracer_cnst_fixed_ymd, &
        tracer_cnst_fixed_tod_in = tracer_cnst_fixed_tod )
   call tracer_srcs_setopts( &
        tracer_srcs_file_in      = tracer_srcs_file,      &
        tracer_srcs_filelist_in  = tracer_srcs_filelist,  &
        tracer_srcs_datapath_in  = tracer_srcs_datapath,  &
        tracer_srcs_type_in      = tracer_srcs_type,      &
        tracer_srcs_specifier_in = tracer_srcs_specifier, &
        tracer_srcs_rmfile_in    = tracer_srcs_rmfile,    &
        tracer_srcs_cycle_yr_in  = tracer_srcs_cycle_yr,  &
        tracer_srcs_fixed_ymd_in = tracer_srcs_fixed_ymd, &
        tracer_srcs_fixed_tod_in = tracer_srcs_fixed_tod )

#ifdef WACCM_MOZART
   call spedata_setopts( spe_data_file_in = spe_data_file, &
        spe_remove_file_in = spe_remove_file, &
        spe_filenames_list_in = spe_filenames_list )
#endif
   call sad_setopts( strat_aero_feedback_in = strat_aero_feedback )

#if ( defined WACCM_MOZART || defined WACCM_GHG )
   ! Upper boundary conditions
   call ubc_setopts( &
        snoe_ubc_file_in =snoe_ubc_file, &
#ifdef WACCM_MOZART
        t_pert_ubc_in    =t_pert_ubc, &
        no_xfac_ubc_in   =no_xfac_ubc, &
#endif
        tgcm_ubc_file_in =tgcm_ubc_file, &
        tgcm_ubc_data_type_in = tgcm_ubc_data_type, &
        tgcm_ubc_cycle_yr_in = tgcm_ubc_cycle_yr, &
        tgcm_ubc_fixed_ymd_in = tgcm_ubc_fixed_ymd, &
        tgcm_ubc_fixed_tod_in = tgcm_ubc_fixed_tod )
#endif

   call aerosol_readnl(nlfile)     
!
   call gas_wetdep_readnl(nlfile)


  endsubroutine chem_readnl

!================================================================================================

function chem_is_active()
!----------------------------------------------------------------------- 
! Purpose: return true if this package is active
!-----------------------------------------------------------------------
   logical :: chem_is_active
!-----------------------------------------------------------------------
   chem_is_active = is_active
end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------
    use chem_mods,       only : gas_pcnst, inv_lst, nfs
    use mo_tracname,     only : solsym

!-----------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------
    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value
!-----------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------
    integer :: m
    
    chem_implements_cnst = .false.
    do m = 1,gas_pcnst
       if( trim(name) /= 'H2O' ) then
          if( trim(name) == solsym(m) ) then
             chem_implements_cnst = .true.
             exit
          end if
       end if
    end do
    do m = 1,nfs
       if( trim(name) /= 'H2O' ) then
          if( trim(name) == inv_lst(m) ) then
             chem_implements_cnst = .true.
             exit
          end if
       endif
    enddo

  end function chem_implements_cnst

  subroutine chem_init(phys_state, pbuf2d)

!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterized greenhouse gas chemistry
!          (declare history variables)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------
    use physics_buffer,      only : physics_buffer_desc, pbuf_get_index
    
    use constituents,        only : cnst_get_ind, cnst_longname
    use cam_history,         only : addfld, add_default, phys_decomp, fieldname_len
    use chem_mods,           only : gas_pcnst, nfs, inv_lst
    use mo_chemini,          only : chemini
    use mo_constants,        only : mo_constants_inti
    use mz_aerosols_intr,    only : mz_aero_initialize
    use mo_ghg_chem,         only : ghg_chem_init
    use mo_tracname,         only : solsym
    use llnl_O1D_to_2OH_adj, only : O1D_to_2OH_adj_init
    use lin_strat_chem,      only : lin_strat_chem_inti
#if ( defined MODAL_AERO )
    use modal_aero_initialize_data, only: modal_aero_initialize
#endif
    use chlorine_loading_data, only : chlorine_loading_init
    use cfc11star,             only : init_cfc11star
    use phys_control,          only : phys_getopts
    use chem_mods,             only : adv_mass
    use infnan,                only : nan, assignment(=)
    use mo_chem_utls,          only : get_spc_ndx
    use abortutils,            only : endrun

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)

    
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer :: m                                ! tracer indicies
    character(len=fieldname_len) :: spc_name
    integer :: i, n, ii

    call phys_getopts(cam_chempkg_out=chem_name)

    call mo_constants_inti
    call mz_aero_initialize()
#if ( defined MODAL_AERO )
    call modal_aero_initialize(pbuf2d)
#endif

!-----------------------------------------------------------------------
! Get liq and ice cloud water indicies
!-----------------------------------------------------------------------
    call cnst_get_ind( 'CLDLIQ', ixcldliq )
    call cnst_get_ind( 'CLDICE', ixcldice )
    call cnst_get_ind( 'NUMLIQ', ixndrop, abort=.false.  )

!-----------------------------------------------------------------------
! get pbuf indicies
!-----------------------------------------------------------------------
    ndx_cld    = pbuf_get_index('CLD')
    ndx_cmfdqr = pbuf_get_index('RPRDTOT')
    ndx_nevapr = pbuf_get_index('NEVAPR')
    ndx_prain  = pbuf_get_index('PRAIN')
    ndx_cldtop = pbuf_get_index('CLDTOP')
    ndx_pblh   = pbuf_get_index('pblh')

    call addfld( 'HEIGHT',    'm',        pverp,'A', 'geopotential height above surface at interfaces (m)', phys_decomp )
    call addfld( 'CT_H2O_GHG','kg/kg/s ', pver, 'A', 'ghg-chem h2o source/sink',                            phys_decomp )

!-----------------------------------------------------------------------
! Set names of chemistry variable tendencies and declare them as history variables
!-----------------------------------------------------------------------
    do m = 1,gas_pcnst
       spc_name = solsym(m)
       srcnam(m) = 'CT_' // spc_name ! chem tendancy (source/sink)

       call addfld( srcnam(m), 'kg/kg/s ', pver, 'A', trim(spc_name)//' source/sink', phys_decomp )
       !call add_default (srcnam(m),     1, ' ')
    end do

    if ( masterproc ) write(iulog,*) 'chem_init: addfld done'

!-----------------------------------------------------------------------
! Initialize chemistry modules
!-----------------------------------------------------------------------
    call chemini &
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

     if ( ghg_chem ) then
        call ghg_chem_init(phys_state, bndtvg, h2orates)
     endif
     
     call O1D_to_2OH_adj_init()

     call lin_strat_chem_inti(phys_state)
     call chlorine_loading_init( chlorine_loading_file, &
                                 type = chlorine_loading_type, &
                                 ymd = chlorine_loading_fixed_ymd, &
                                 tod = chlorine_loading_fixed_tod )

     if ( chem_is('waccm_mozart') .or. chem_is('waccm_mozart_mam3') ) then
        call init_cfc11star(pbuf2d)
     endif

    ! Bulk aerosols
    call cnst_get_ind('DST01',dust_indices(1), abort=.false.)
    call cnst_get_ind('DST02',dust_indices(2), abort=.false.)
    call cnst_get_ind('DST03',dust_indices(3), abort=.false.)
    call cnst_get_ind('DST04',dust_indices(4), abort=.false.)

    ! MODAL aerosols
    call cnst_get_ind('dst_a1',dust_indices(5), abort=.false.)
    call cnst_get_ind('dst_a3',dust_indices(6), abort=.false.)
    call cnst_get_ind('dst_a5',dust_indices(7), abort=.false.)
    call cnst_get_ind('dst_a7',dust_indices(8), abort=.false.)

    ! MEGAN emissions initialize
    if (shr_megan_mechcomps_n>0) then

       allocate( megan_indices_map(shr_megan_mechcomps_n) )
       allocate( megan_wght_factors(shr_megan_mechcomps_n) )
       megan_wght_factors(:) = nan

       do n=1,shr_megan_mechcomps_n
          call cnst_get_ind (shr_megan_mechcomps(n)%name,  megan_indices_map(n), abort=.false.)
          ii = get_spc_ndx(shr_megan_mechcomps(n)%name)
          if (ii>0) then
             megan_wght_factors(n) = adv_mass(ii)*1.e-3_r8 ! kg/moles (to convert moles/m2/sec to kg/m2/sec)
          else
             call endrun('gas_phase_chemdr_inti: MEGAN compound not in chemistry mechanism : '//trim(shr_megan_mechcomps(n)%name))
          endif

          ! MEGAN  history fields
          call addfld( 'MEG_'//trim(shr_megan_mechcomps(n)%name),'kg/m2/sec',1,'A',&
                       trim(shr_megan_mechcomps(n)%name)//' MEGAN emissions flux',phys_decomp)
       enddo
    endif

  end subroutine chem_init

!================================================================================
!================================================================================
  subroutine chem_reset_fluxes( fptr, cam_in )
    use camsrfexch,          only : cam_in_t     

    use constituents,        only : pcnst, cnst_get_ind
    use mo_gas_phase_chemdr, only : map2chm
    use phys_grid,           only : get_ncols_p
    use cam_cpl_indices,     only : index_x2a_Fall_flxvoc
    use cam_history,         only : outfld

    real(r8), pointer             :: fptr(:,:)        ! pointer into    array data
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)

    integer :: i,c,ig,imech, ncols  ! indices
    integer :: m,n

    ! dust_intr module modifies cflx for dust each time step
    ! -- resetting dust emissions each time step here changes answers...??

    ig=1
    do c=begchunk,endchunk
       ncols = get_ncols_p(c) 

       ! initialize chemistry constituent surface fluxes to zero
       do m = 2,pcnst
          n = map2chm(m)
          if (.not. ( any(dust_indices(:)==m) )) then
             if (n>0) cam_in(c)%cflx(:,m) = 0._r8 
          endif
       enddo
       
       if ( index_x2a_Fall_flxvoc>0 .and. shr_megan_mechcomps_n>0 ) then

          ! set MEGAN fluxes 
          do i =1,ncols
             do imech = 1,shr_megan_mechcomps_n
                cam_in(c)%cflx(i,megan_indices_map(imech)) = -fptr(index_x2a_Fall_flxvoc+imech-1, ig)*megan_wght_factors(imech)
             enddo
             ig=ig+1
          end do

          ! output MEGAN emis fluxes to history
          do n = 1,shr_megan_mechcomps_n
             call outfld('MEG_'//trim(shr_megan_mechcomps(n)%name), cam_in(c)%cflx(:ncols, megan_indices_map(n)), ncols, c)
          enddo

       endif

    end do
    
  end subroutine chem_reset_fluxes
!================================================================================

  subroutine chem_init_cnst( name, q, gcid)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Specify initial mass mixing ratios
! 
!-----------------------------------------------------------------------

    use chem_mods, only : inv_lst

    use physconst,     only : mwdry, mwch4, mwn2o, mwf11, mwf12, mwh2o, mwo3
    use chem_surfvals, only : chem_surfvals_get

    implicit none

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    character(len=*), intent(in) :: name                   !  constituent name
    real(r8), intent(inout) :: q(:,:)           !  mass mixing ratio
    integer, intent(in)     :: gcid(:)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    
    real(r8) :: rmwn2o != mwn2o/mwdry ! ratio of mol weight n2o   to dry air
    real(r8) :: rmwch4 != mwch4/mwdry ! ratio of mol weight ch4   to dry air
    real(r8) :: rmwf11 != mwf11/mwdry ! ratio of mol weight cfc11 to dry air
    real(r8) :: rmwf12 != mwf12/mwdry ! ratio of mol weight cfc12 to dry air

!-----------------------------------------------------------------------
! initialize local variables
!-----------------------------------------------------------------------

    rmwn2o = mwn2o/mwdry 
    rmwch4 = mwch4/mwdry 
    rmwf11 = mwf11/mwdry 
    rmwf12 = mwf12/mwdry 

!-----------------------------------------------------------------------
! Get initial mixing ratios
!-----------------------------------------------------------------------
    if ( any( inv_lst .eq. name ) ) then
       q(:,:) = 0.0_r8
    else
       q(:,:) = 1.e-38_r8
    endif

    if ( ghg_chem ) then
       select case (name)
       case ('N2O')
          q = rmwn2o * chem_surfvals_get('N2OVMR')
       case ('CH4')
          q = rmwch4 * chem_surfvals_get('CH4VMR')
       case ('CFC11')
          q = rmwf11 * chem_surfvals_get('F11VMR')
       case ('CFC12')
          q = rmwf12 * chem_surfvals_get('F12VMR')
       end select
    endif

  end subroutine chem_init_cnst

  subroutine chem_timestep_init(phys_state,pbuf2d)

    use time_manager,      only : get_nstep
    use time_manager,      only : get_curr_calday
    use mo_srf_emissions,  only : set_srf_emissions_time
    use mo_sulf,           only : set_sulf_time
    use mo_extfrc,         only : extfrc_timestep_init
    use mo_flbc,           only : flbc_chk
    use tracer_cnst,       only : tracer_cnst_adv
    use tracer_srcs,       only : tracer_srcs_adv
    use mo_ghg_chem,       only : ghg_chem_timestep_init

    use mo_solar_parms,    only : solar_parms_timestep_init
    use mo_jshort,         only : jshort_timestep_init
    use mo_jlong,          only : jlong_timestep_init
    use mo_aurora,         only : aurora_timestep_init
    use spedata,           only : advance_spedata
    use mo_photo,          only : photo_timestep_init
    use mo_strato_sad,     only : strato_sad_timestep_init
    use linoz_data,        only : linoz_data_adv
    use chlorine_loading_data, only : chlorine_loading_advance

    use cfc11star,         only : update_cfc11star
    use physics_buffer,    only : physics_buffer_desc

    implicit none

    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    real(r8) :: calday
    integer :: nstep

    nstep = get_nstep()
    chem_step = mod( nstep, chem_freq ) == 0

    if ( .not. chem_step ) return

    !-----------------------------------------------------------------------
    ! get current calendar day of year
    !-----------------------------------------------------------------------
    calday = get_curr_calday( )

    !-----------------------------------------------------------------------
    ! Set emissions timing factors
    !-----------------------------------------------------------------------
    call set_srf_emissions_time( pbuf2d, phys_state )

    !-----------------------------------------------------------------------
    ! Set external forcings timing factors
    !-----------------------------------------------------------------------
    call extfrc_timestep_init( pbuf2d, phys_state )

    !-----------------------------------------------------------------------
    ! Set sulf timing factors
    !-----------------------------------------------------------------------
    call set_sulf_time( pbuf2d, phys_state  )

    !-----------------------------------------------------------------------
    ! Set fixed lower boundary timing factors
    !-----------------------------------------------------------------------
    call flbc_chk

    !-----------------------------------------------------------------------
    ! Advance stratospheric aerosol densities...
    !-----------------------------------------------------------------------
    call strato_sad_timestep_init()

    !-----------------------------------------------------------------------
    ! Set fixed offline tracers
    !-----------------------------------------------------------------------
    call tracer_cnst_adv(pbuf2d, phys_state)

    !-----------------------------------------------------------------------
    ! Set fixed offline tracer sources
    !-----------------------------------------------------------------------
    call tracer_srcs_adv(pbuf2d, phys_state)

    !-----------------------------------------------------------------------
    ! Advance the linoz data
    !-----------------------------------------------------------------------
    call linoz_data_adv(pbuf2d, phys_state)
    call chlorine_loading_advance()

    if ( ghg_chem ) then
       call ghg_chem_timestep_init(phys_state)
    endif

    if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart') .or. chem_is('waccm_mozart_mam3')) then
       !-----------------------------------------------------------------------
       ! Set solar parameters
       !-----------------------------------------------------------------------
       call solar_parms_timestep_init
    endif

    !-----------------------------------------------------------------------------
    !	... setup the time interpolation for mo_photo
    !-----------------------------------------------------------------------------
    call photo_timestep_init( calday )

    !-----------------------------------------------------------------------
    ! Set jlong etf
    !-----------------------------------------------------------------------
    call jlong_timestep_init

    if (chem_is('waccm_mozart') .or. chem_is('waccm_mozart_mam3')) then
       !-----------------------------------------------------------------------
       ! Set jshort etf
       !-----------------------------------------------------------------------
       call jshort_timestep_init

       !-----------------------------------------------------------------------
       ! Set up aurora
       !-----------------------------------------------------------------------
       call aurora_timestep_init

       !-----------------------------------------------------------------------
       ! Set up solar proton data
       !-----------------------------------------------------------------------
       call advance_spedata()

       call update_cfc11star( pbuf2d, phys_state )

    endif

  end subroutine chem_timestep_init

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, &
       pbuf,  fh2o, fsds )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to parameterized greenhouse gas chemisty (source/sink).
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------
    use physics_buffer,      only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use cam_history,         only : outfld
    use time_manager,        only : get_curr_calday
    use chem_mods,           only : gas_pcnst
    use mo_gas_phase_chemdr, only : gas_phase_chemdr
    
    use spmd_utils,          only : iam
    use camsrfexch,          only : cam_in_t, cam_out_t     
    use perf_mod,            only : t_startf, t_stopf
    use tropopause,          only : tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE
    use mo_drydep,           only : drydep_update
    use mo_neu_wetdep,       only : neu_wetdep_tend, do_neu_wetdep
#if (defined MODAL_AERO)
    use modal_aero_data,     only : ntot_amode    
#endif
    use aerodep_flx,         only : aerodep_flx_prescribed
    
    implicit none

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    real(r8),            intent(in)    :: dt              ! time step
    type(physics_state), intent(in)    :: state           ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend           ! indivdual parameterization tendencies
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(inout) :: cam_out
    real(r8),            intent(out)   :: fh2o(pcols)     ! h2o flux to balance source from chemistry
    

    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8),            intent(in)    :: fsds(pcols)     ! longwave down at sfc

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer  :: i, k, m, n                         ! indicies
    integer  :: lchnk                              ! chunk identifier
    integer  :: ncol                               ! number of atmospheric columns
    real(r8) :: calday                             ! current calendar day of year
    real(r8) :: cldw(pcols,pver)                   ! cloud water (kg/kg)
    real(r8) :: chem_dt              ! time step
    real(r8) :: drydepflx(pcols,pcnst)             ! dry deposition fluxes (kg/m2/s)
    integer  :: tropLev(pcols)
    real(r8) :: ncldwtr(pcols,pver)                ! droplet number concentration (#/kg)
    real(r8), pointer :: pblh(:)
    real(r8), pointer :: prain(:,:)
    real(r8), pointer :: cldfr(:,:)
    real(r8), pointer :: cmfdqr(:,:)
    real(r8), pointer :: nevapr(:,:)
    real(r8), pointer :: cldtop(:)

    integer :: tim_ndx

    logical :: lq(pcnst)

    if ( .not. chem_step ) return

    chem_dt = chem_freq*dt

    lchnk = state%lchnk
    ncol  = state%ncol

    lq(:) = .false.
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          lq(n) = .true.
       end if
    end do
    if ( ghg_chem ) lq(1) = .true.

    call physics_ptend_init(ptend, state%psetcols, 'chemistry', lq=lq)
    
    call drydep_update( state, cam_in )

!-----------------------------------------------------------------------
! get current calendar day of year
!-----------------------------------------------------------------------
    calday = get_curr_calday()

!-----------------------------------------------------------------------
! get tropopause level
!-----------------------------------------------------------------------
    call tropopause_find(state, tropLev, primary=TROP_ALG_HYBSTOB, backup=TROP_ALG_CLIMATE)

    tim_ndx = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ndx_pblh,       pblh)
    call pbuf_get_field(pbuf, ndx_prain,      prain,  start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cld,        cldfr,  start=(/1,1,tim_ndx/), kount=(/ncol,pver,1/) )
    call pbuf_get_field(pbuf, ndx_cmfdqr,     cmfdqr, start=(/1,1/),         kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_nevapr,     nevapr, start=(/1,1/),         kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cldtop,     cldtop )

!-----------------------------------------------------------------------
! call Neu wet dep scheme
!-----------------------------------------------------------------------
    call neu_wetdep_tend(lchnk,ncol,state%q,state%pmid,state%pdel,state%zi,state%t,dt, &
         prain, nevapr, cldfr, cmfdqr, ptend%q)

!-----------------------------------------------------------------------
! compute tendencies and surface fluxes
!-----------------------------------------------------------------------
    call t_startf( 'chemdr' )
    do k = 1,pver
       cldw(:ncol,k) = state%q(:ncol,k,ixcldliq) + state%q(:ncol,k,ixcldice)
       if (ixndrop>0) &
            ncldwtr(:ncol,k) = state%q(:ncol,k,ixndrop)
    end do

    call gas_phase_chemdr(lchnk, ncol, imozart, state%q, ptend%q, &
                          cam_in%cflx, state%phis, state%zm, state%zi, calday, &
                          state%t, state%pmid, state%pdel, state%pint, &
                          cldw, tropLev, ncldwtr, state%u, state%v, &
                          chem_dt, state%ps, xactive_prates, &
                          fsds, cam_in%ts, cam_in%asdir, cam_in%ocnfrac, cam_in%icefrac, &
                          cam_out%precc, cam_out%precl, cam_in%snowhland, ghg_chem, state%latmapback, &
                          chem_name, drydepflx, pbuf)

    call t_stopf( 'chemdr' )

!-----------------------------------------------------------------------
! set flags for tracer tendencies (water and gas phase constituents)
! record tendencies on history files
!-----------------------------------------------------------------------
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          call outfld( srcnam(m), ptend%q(:,:,n), pcols, lchnk )
       end if

       ! if the user has specified prescribed aerosol dep fluxes then 
       ! do not set cam_out dep fluxes according to the prognostic aerosols
       if (.not.aerodep_flx_prescribed()) then
          ! set deposition fluxes in the export state
          select case (trim(cnst_name(n)))
          case('CB1')
             do i = 1, ncol
                cam_out%bcphodry(i) = max(drydepflx(i,n), 0._r8)
             end do
          case('CB2')
             do i = 1, ncol
                cam_out%bcphidry(i) = max(drydepflx(i,n), 0._r8)
             end do
          case('OC1')
             do i = 1, ncol
                cam_out%ocphodry(i) = max(drydepflx(i,n), 0._r8)
             end do
          case('OC2')
             do i = 1, ncol
                cam_out%ocphidry(i) = max(drydepflx(i,n), 0._r8)
             end do
          end select
       endif
    end do
    if ( ghg_chem ) then
       ptend%lq(1) = .true.
       call outfld( 'CT_H2O_GHG', ptend%q(:,:,1), pcols, lchnk )
    endif

    call outfld( 'HEIGHT', state%zi(:ncol,:),  ncol, lchnk )

!-----------------------------------------------------------------------
!  turn off water vapor tendency if radiatively passive
!-----------------------------------------------------------------------
    if (chem_rad_passive) then
       ptend%lq(1) = .false.
       ptend%q(:ncol,:,1) = 0._r8
    endif

!-----------------------------------------------------------------------
! Compute water vapor flux required to make conservation check
!-----------------------------------------------------------------------
    fh2o(:ncol) = 0._r8
    do k = 1,pver
       fh2o(:ncol) = fh2o(:ncol) + ptend%q(:ncol,k,1)*state%pdel(:ncol,k)/gravit
    end do
  end subroutine chem_timestep_tend

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine chem_final
  end subroutine chem_final

!-------------------------------------------------------------------
!-------------------------------------------------------------------

  subroutine chem_init_restart( File )
    use pio, only : file_desc_t
    use tracer_cnst,      only: init_tracer_cnst_restart
    use tracer_srcs,      only: init_tracer_srcs_restart
    use linoz_data, only : init_linoz_data_restart
    implicit none
    type(file_desc_t),intent(inout) :: File     ! pio File pointer

    !
    ! data for offline tracers
    !
    call init_tracer_cnst_restart(File)
    call init_tracer_srcs_restart(File)
    call init_linoz_data_restart(File)
  end subroutine chem_init_restart
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine chem_write_restart( File )
    use tracer_cnst, only: write_tracer_cnst_restart
    use tracer_srcs, only: write_tracer_srcs_restart
    use linoz_data,  only: write_linoz_data_restart
    use pio, only : file_desc_t
    implicit none
    type(file_desc_t) :: File

    !
    ! data for offline tracers
    !
    call write_tracer_cnst_restart(File)
    call write_tracer_srcs_restart(File)
    call write_linoz_data_restart(File)
  end subroutine chem_write_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine chem_read_restart( File )
    use tracer_cnst, only: read_tracer_cnst_restart
    use tracer_srcs, only: read_tracer_srcs_restart
    use linoz_data,  only: read_linoz_data_restart

    use pio, only : file_desc_t
    implicit none
    type(file_desc_t) :: File

    !
    ! data for offline tracers
    !
    call read_tracer_cnst_restart(File)
    call read_tracer_srcs_restart(File)
    call read_linoz_data_restart(File)
  end subroutine chem_read_restart

end module chemistry
