#define CAM_VERSION_IS_ACME


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_amicphys --- does modal aerosol gas-aerosol exchange
!
! !INTERFACE:
   module modal_aero_amicphys

! !USES:
  use shr_kind_mod,    only:  r8 => shr_kind_r8
  use cam_abortutils,  only:  endrun
  use cam_logfile,     only:  iulog
  use chem_mods,       only:  gas_pcnst
  use physconst,       only:  pi
  use ppgrid,          only:  pcols, pver
  use modal_aero_data, only:  ntot_aspectype, ntot_amode, nsoa, npoa, nbc
! use ref_pres,        only:  top_lev => clim_modal_aero_top_lev  ! this is for gg02a
  use ref_pres,        only:  top_lev => trop_cloud_top_lev       ! this is for ee02c

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public modal_aero_amicphys_intr, modal_aero_amicphys_init

! !PUBLIC DATA MEMBERS:
  type :: misc_vars_aa_type
! using this derived type reduces the number of changes needed to add more mosaic diagnostics to history
     real(r8) :: ncluster_tend_nnuc_1grid
#if ( defined( MOSAIC_SPECIES ) )
     real(r8) :: cnvrg_fail_1grid
     real(r8) :: max_kelvin_iter_1grid
     real(r8), dimension(5,4) :: xnerr_astem_negative_1grid
#endif
  end type misc_vars_aa_type

  logical, public :: mosaic = .true. !BSINGH -  Added logical for mosaic model

  integer, parameter :: pcnstxx = gas_pcnst

  real(r8), public, protected :: n_so4_monolayers_pcage = huge(1.0_r8)
! number of so4(+nh4) monolayers needed to "age" a carbon particle

  real(r8), public, protected :: dr_so4_monolayers_pcage = huge(1.0_r8)
! thickness of the so4 monolayers (m)
! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3,
!    --> 1 mol so4(+nh4)  = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
! aging criterion is approximate so do not try to distinguish
!    sulfuric acid, bisulfate, ammonium sulfate

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
  integer, public :: cldy_rh_sameas_clear = 0
! this is only used for some specific box model tests
#endif

  integer, public :: mdo_gaexch_cldy_subarea = 0
! controls if gas condensation is done in cloudy subarea
!    1 = yes ; 0 = no

  integer, public :: gaexch_h2so4_uptake_optaa = 2
! controls treatment of h2so4 condensation in mam_gasaerexch_1subarea
!    1 = sequential   calc. of gas-chem prod then condensation loss
!    2 = simultaneous calc. of gas-chem prod and  condensation loss

  integer, public :: newnuc_h2so4_conc_optaa = 2
! controls treatment of h2so4 concentrationin mam_newnuc_1subarea
!    1 = use average value calculated in standard cam5.2.10 and earlier
!    2 = use average value calculated in mam_gasaerexch_1subarea
!   11 = use average of initial and final values from mam_gasaerexch_1subarea
!   12 = use final value from mam_gasaerexch_1subarea

  integer, public :: rename_method_optaa = 40
! controls renaming parameterization

  integer, public :: update_qaerwat = 0
  integer, public :: update_dgncur_a = 0
  integer, public :: update_dgncur_awet = 0
! controls updating of qaerwat
! controls updating of dgncur_a
! controls updating of dgncur_awet and wetdens_host

  real (r8) :: newnuc_adjust_factor_dnaitdt = 1.0_r8
  real (r8) :: newnuc_adjust_factor_pbl     = 1.0_r8


#if ( defined MODAL_AERO_3MODE )
  integer, parameter :: max_gas = nsoa + 1
  ! the +3 in max_aer are dst, ncl, so4
  integer, parameter :: max_aer = nsoa + npoa + nbc + 3
#elif ( defined MODAL_AERO_4MODE )
  integer, parameter :: max_gas = nsoa + 1
  ! the +3 in max_aer are dst, ncl, so4
  integer, parameter :: max_aer = nsoa + npoa + nbc + 3
#elif ( defined MODAL_AERO_4MODE_MOM )
  integer, parameter :: max_gas = nsoa + 1
  ! the +4 in max_aer are dst, ncl, so4, mom
  integer, parameter :: max_aer = nsoa + npoa + nbc + 4
#elif ( ( defined MODAL_AERO_7MODE ) && ( defined MOSAIC_SPECIES ) )
  integer, parameter :: max_gas = nsoa + 4
  ! the +8 in max_aer are dst, ncl(=na), so4, no3, cl, nh4, ca, co3 
  integer, parameter :: max_aer = nsoa + npoa + nbc + 8
#elif ( defined MODAL_AERO_7MODE )
  integer, parameter :: max_gas = nsoa + 2
  ! the +4 in max_aer are dst, ncl, so4, nh4
  integer, parameter :: max_aer = nsoa + npoa + nbc + 4
#elif ( defined MODAL_AERO_8MODE )
  integer, parameter :: max_gas = nsoa + 2
  ! the +4 in max_aer are dst, ncl, so4, mom ???
  integer, parameter :: max_aer = nsoa + npoa + nbc + 4
#elif ( defined MODAL_AERO_9MODE )
  integer, parameter :: max_gas = nsoa + 2
  ! the +4+5 in max_aer are dst, ncl, so4, nh4 and 5 marine organics
  integer, parameter :: max_aer = nsoa + npoa + nbc + 4 + 5
#endif

#if (( defined MODAL_AERO_8MODE ) || ( defined MODAL_AERO_4MODE ) || ( defined MODAL_AERO_4MODE_MOM ))
  integer, parameter :: ntot_amode_extd = ntot_amode
#else
  integer, parameter :: ntot_amode_extd = ntot_amode + 1
! integer, parameter :: ntot_amode_extd = ntot_amode
#endif

  integer, parameter :: max_mode_fresh = 1

  integer, parameter :: max_mode = ntot_amode_extd + max_mode_fresh
  public max_mode !BSINGH - used in module_mosaic_cam_init.F90

  integer, parameter :: max_coagpair = 100

#if ( defined MODAL_AERO_9MODE )
  integer, parameter :: max_agepair = 3
#else
  integer, parameter :: max_agepair = 1
#endif

  integer, parameter :: maxsubarea = 2

  integer, parameter :: nqtendaa = 4
  integer, parameter :: iqtend_cond = 1
  integer, parameter :: iqtend_rnam = 2
  integer, parameter :: iqtend_nnuc = 3
  integer, parameter :: iqtend_coag = 4
  integer, parameter :: nqqcwtendaa = 1
  integer, parameter :: iqqcwtend_rnam = 1

  integer, parameter :: iqqcwtend_match_iqtend(nqtendaa) = (/ 0, iqqcwtend_rnam, 0, 0 /)

  logical, parameter :: aging_include_seasalt = .false.
                      ! when .true., aging (by coagulation) includes contribution of seasalt
                      ! early versions of mam neglected the seasalt contribution

  ! species indices for various qgas_--- arrays
  integer :: igas_soa, igas_h2so4, igas_nh3, igas_hno3, igas_hcl
  ! species indices for various qaer_--- arrays
  !    when nsoa > 1, igas_soa and iaer_soa are indices of the first soa species
  !    when nbc  > 1, iaer_bc  is index of the first bc  species
  !    when npom > 1, iaer_pom is index of the first pom species
  integer :: iaer_bc, iaer_dst, iaer_ncl, iaer_nh4, iaer_pom, iaer_soa, iaer_so4, &
             iaer_mpoly, iaer_mprot, iaer_mlip, iaer_mhum, iaer_mproc, iaer_mom, &
             iaer_no3, iaer_cl, iaer_ca, iaer_co3
  integer :: i_agepair_pca, i_agepair_macc, i_agepair_mait
  integer :: lmap_gas(max_gas)
  integer :: lmap_aer(max_aer,max_mode), lmapbb_aer(max_aer,max_mode), &
             lmap_aercw(max_aer,max_mode)
  integer :: lmap_num(max_mode), lmap_numcw(max_mode)
  integer :: lmapcc_all(gas_pcnst)
  integer, parameter :: lmapcc_val_gas = 1, lmapcc_val_aer = 2, lmapcc_val_num = 3
  integer :: ngas, naer
  integer :: nacc, nait, npca, nufi, nmacc, nmait

  integer :: n_agepair, n_coagpair
  integer :: modefrm_agepair(max_agepair), modetoo_agepair(max_agepair)
  integer :: mode_aging_optaa(max_mode)
  integer :: modefrm_coagpair(max_coagpair), modetoo_coagpair(max_coagpair), &
             modeend_coagpair(max_coagpair)

  integer :: lun82,   lun97,   lun98,   lun13n,   lun15n
  logical :: ldiag82, ldiag97, ldiag98, ldiag13n, ldiag15n
  logical :: ldiagd1

  real(r8) :: accom_coef_gas(max_gas)
  real(r8) :: alnsg_aer(max_mode)
  real(r8) :: dgnum_aer(max_mode), dgnumhi_aer(max_mode), dgnumlo_aer(max_mode)
  real(r8) :: dens_aer(max_aer)
  real(r8) :: dens_so4a_host
  real(r8) :: fac_m2v_aer(max_aer)        ! converts (mol-aero/mol-air) to (m3-aero/mol-air)
  real(r8) :: fac_eqvso4hyg_aer(max_aer)  ! converts a species volume to a volume of so4
                                          !    (or nh4hso4) having same hygroscopicity
  real(r8) :: fac_m2v_eqvhyg_aer(max_aer) ! = fac_m2v_aer * fac_eqvso4hyg_aer

  real(r8) :: fcvt_gas(max_gas), fcvt_aer(max_aer), fcvt_num, fcvt_wtr
  real(r8) :: fcvt_dgnum_dvolmean(max_mode)
  real(r8) :: hygro_aer(max_aer)
  real(r8) :: mw_gas(max_gas), mw_aer(max_aer)
  real(r8) :: mwhost_gas(max_gas), mwhost_aer(max_aer), mwhost_num
  real(r8) :: mw_nh4a_host, mw_so4a_host
  real(r8) :: mwuse_soa(nsoa), mwuse_poa(npoa)
  real(r8) :: sigmag_aer(max_mode)
  real(r8) :: vol_molar_gas(max_gas)

! following were used in aging calcs but are no longer needed
!    fac_m2v_so4, fac_m2v_nh4, fac_m2v_soa(:)
!    fac_m2v_pcarbon(:)
!    soa_equivso4_factor(:)


  character(len=16) :: name_gas(max_gas), name_aerpfx(max_aer), &
     name_aer(max_aer,max_mode), name_aercw(max_aer,max_mode), &
     name_num(max_mode), name_numcw(max_mode)

  character(len=8) :: suffix_q_coltendaa(nqtendaa) = &
     (/ '_sfgaex1', '_sfgaex2', '_sfnnuc1', '_sfcoag1' /)
  character(len=8) :: suffix_qqcw_coltendaa(nqqcwtendaa) = &
                    '_sfgaex2'

  logical :: do_q_coltendaa(gas_pcnst,nqtendaa) = .false.
  logical :: do_qqcw_coltendaa(gas_pcnst,nqqcwtendaa) = .false.

! *** following 3 variables should eventually be in modal_aero_data
  real(r8) :: specmw2_amode(ntot_aspectype,ntot_amode)
  real(r8) :: specdens2_amode(ntot_aspectype,ntot_amode)
  real(r8) :: spechygro2(ntot_aspectype,ntot_amode)


! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------


  contains


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine modal_aero_amicphys_intr(                             &
                        mdo_gasaerexch,     mdo_rename,          &
                        mdo_newnuc,         mdo_coag,            &
                        lchnk,    ncol,     nstep,               &
                        loffset,  deltat,                        &
                        latndx,   lonndx,                        &
                        t,        pmid,     pdel,                &
                        zm,       pblh,                          &
                        qv,       cld,                           &
                        q,                  qqcw,                &
                        q_pregaschem,                            &
                        q_precldchem,       qqcw_precldchem,     &
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
                        nqtendbb,           nqqcwtendbb,         &
                        q_tendbb,           qqcw_tendbb,         &
#endif
                        dgncur_a,           dgncur_awet,         &
                        wetdens_host,                            &
                        qaerwat                                  )


! !USES:
use cam_history,       only:  outfld, fieldname_len
use chem_mods,         only:  adv_mass
use constituents,      only:  cnst_name
use physconst,         only:  gravit, mwdry, r_universal
use wv_saturation,     only:  qsat
use phys_control,      only:  phys_getopts

use modal_aero_data,   only:  &
    cnst_name_cw, &
    lmassptr_amode, lmassptrcw_amode, lptr2_soa_g_amode, &
    nspec_amode, &
    numptr_amode, numptrcw_amode
use modal_aero_newnuc, only:  adjust_factor_pbl_ratenucl


implicit none

! !PARAMETERS:
   integer,  intent(in)    :: mdo_gasaerexch, mdo_rename, mdo_newnuc, mdo_coag
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric columns in the chunk
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   integer,  intent(in)    :: latndx(pcols), lonndx(pcols)
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
   integer,  intent(in)    :: nqtendbb             ! dimension for q_tendbb
   integer,  intent(in)    :: nqqcwtendbb          ! dimension for qqcw_tendbb
#endif

   real(r8), intent(in)    :: deltat               ! time step (s)

   real(r8), intent(inout) :: q(ncol,pver,pcnstxx) ! current tracer mixing ratios (TMRs)
                                                   ! these values are updated (so out /= in)
                                                   ! *** MUST BE  #/kmol-air for number
                                                   ! *** MUST BE mol/mol-air for mass
                                                   ! *** NOTE ncol dimension
   real(r8), intent(inout) :: qqcw(ncol,pver,pcnstxx) 
                                                   ! like q but for cloud-borner tracers
                                                   ! these values are updated
   real(r8), intent(in)    :: q_pregaschem(ncol,pver,pcnstxx)    ! q TMRs    before gas-phase chemistry
   real(r8), intent(in)    :: q_precldchem(ncol,pver,pcnstxx)    ! q TMRs    before cloud chemistry
   real(r8), intent(in)    :: qqcw_precldchem(ncol,pver,pcnstxx) ! qqcw TMRs before cloud chemistry
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
   real(r8), intent(inout) :: q_tendbb(ncol,pver,pcnstxx,nqtendbb)    ! TMR tendencies for box-model diagnostic output
   real(r8), intent(inout) :: qqcw_tendbb(ncol,pver,pcnstxx,nqqcwtendbb) 
#endif

   real(r8), intent(in)    :: t(pcols,pver)        ! temperature at model levels (K)
   real(r8), intent(in)    :: pmid(pcols,pver)     ! pressure at model level centers (Pa)
   real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: zm(pcols,pver)       ! altitude (above ground) at level centers (m)
   real(r8), intent(in)    :: pblh(pcols)          ! planetary boundary layer depth (m)
   real(r8), intent(in)    :: qv(pcols,pver)       ! specific humidity (kg/kg)
   real(r8), intent(in)    :: cld(ncol,pver)       ! cloud fraction (-) *** NOTE ncol dimension
   real(r8), intent(inout) :: dgncur_a(pcols,pver,ntot_amode)
   real(r8), intent(inout) :: dgncur_awet(pcols,pver,ntot_amode)
                                 ! dry & wet geo. mean dia. (m) of number distrib.
   real(r8), intent(inout) :: wetdens_host(pcols,pver,ntot_amode)
                                 ! interstitial aerosol wet density (kg/m3)
   real(r8), intent(inout), optional :: &
                              qaerwat(pcols,pver,ntot_amode)
                                 ! aerosol water mixing ratio (kg/kg, NOT mol/mol)

! !DESCRIPTION: 
! calculates changes to gas and aerosol TMRs (tracer mixing ratios) from
!    gas-aerosol exchange (condensation/evaporation)
!    growth from smaller to larger modes (renaming) due to both 
!       condensation and cloud chemistry
!    new particle nucleation
!    coagulation
!    transfer of particles from hydrophobic modes to hydrophilic modes (aging)
!       due to condensation and coagulation
!
! the incoming mixing ratios (q and qqcw) are updated before output
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from earlier version of CAM5 modal aerosol routines
!                  for these processes
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
   integer, parameter :: method_soa = 2
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr modal_aero_soaexch

   integer :: i, icol_diag, ipass, iq
   integer :: itmpa, itmpb, itmpc, itmpd
   integer :: iqtend, iqqcwtend
   integer :: iaer, igas
   integer :: j, jac, jsoa, jsub
   integer :: jclea, jcldy
   integer :: k
   integer :: l, l2, l3, la, lb, lc, lmz, lsfrm, lstoo
   integer :: lun, lund
   integer :: m
   integer :: n, niter, niter_max, ntot_soamode
   integer :: nsubarea, ncldy_subarea

   logical :: do_cond, do_rename, do_newnuc, do_coag
   logical :: iscldy_subarea(maxsubarea)

   character(len=fieldname_len+3) :: fieldname
   character(len=6)   :: tmpch6a, tmpch6c
   character(len=200) :: tmp_str

   real (r8) :: pdel_fac

!----------------------------------------------------------------------
   logical   :: history_aerocom    ! Output the aerocom history
!-----------------------------------------------------------------------
 

      real(8), parameter :: fcld_locutoff = 1.0e-5_r8
! cloud chemistry is only on when cld(i,k) >= 1.0e-5_r8
! it may be that the macrophysics has a higher threshold that this
      real(8), parameter :: fcld_hicutoff = 0.999_r8

      real(r8) :: afracsub(maxsubarea)
      real(r8) :: dgn_a(max_mode), dgn_awet(max_mode)
      real(r8) :: ev_sat(pcols,pver)
      real(r8) :: fclea, fcldy, fcldybb
      real(r8) :: nufine_3dtend_nnuc(pcols,pver)
      real(r8) :: ncluster_3dtend_nnuc(pcols,pver)
      real(r8) :: qv_sat(pcols,pver)
      real(r8) :: relhumgcm, relhumsub(maxsubarea)
      real(r8) :: soag_3dtend_cond(pcols,pver,nsoa)
      real(r8) :: tmpa, tmpb, tmpc
      real(r8) :: tmp_qa_clea, tmp_qa_cldy, tmp_qa_gcav
      real(r8) :: tmp_qc_cldy, tmp_qc_gcav
      real(r8) :: tmp_aa, tmp_aa_clea, tmp_aa_cldy
      real(r8) :: tmp_kxt, tmp_kxt2, tmp_pxt, tmp_pok
      real(r8) :: tmp_q1, tmp_q2, tmp_q3, tmp_q4, tmp_q5, tmp_qdot4
      real(r8) :: wetdens(max_mode)


! qgcmN and qqcwgcmN (N=1:4) are grid-cell mean tracer mixing ratios (TMRs, mol/mol or #/kmol)
!    N=1 - before gas-phase chemistry
!    N=2 - before cloud chemistry
!    N=3 - incoming values (before gas-aerosol exchange, newnuc, coag)
!    N=4 - outgoing values (after  gas-aerosol exchange, newnuc, coag)
      real(r8), dimension( 1:gas_pcnst ) :: &
         qgcm1, qgcm2, qgcm3, qgcm4, &
         qqcwgcm1, qqcwgcm2, qqcwgcm3, qqcwgcm4
      real(r8), dimension( 1:gas_pcnst, 1:nqtendaa ) :: &
         qgcm_tendaa
      real(r8), dimension( 1:gas_pcnst, 1:nqqcwtendaa ) :: &
         qqcwgcm_tendaa
      real(r8), dimension( 1:ntot_amode_extd ) :: &
         qaerwatgcm3, qaerwatgcm4   ! aerosol water mixing ratios (mol/mol)

! qsubN and qqcwsubN (N=1:4) are TMRs in sub-areas
!    currently there are just clear and cloudy sub-areas
!    the N=1:4 have same meanings as for qgcmN
      real(r8), dimension( 1:gas_pcnst, 1:maxsubarea ) :: &
         qsub1, qsub2, qsub3, qsub4, &
         qqcwsub1, qqcwsub2, qqcwsub3, qqcwsub4
      real(r8), dimension( 1:gas_pcnst, 1:nqtendaa, 1:maxsubarea ) :: &
         qsub_tendaa
      real(r8), dimension( 1:gas_pcnst, 1:nqqcwtendaa, 1:maxsubarea ) :: &
         qqcwsub_tendaa
      real(r8), dimension( 1:ntot_amode_extd, 1:maxsubarea ) :: &
         qaerwatsub3, qaerwatsub4   ! aerosol water mixing ratios (mol/mol)

! q_coltendaa and qqcw_coltendaa are column-integrated tendencies
!    for different processes, which are output to history
! the processes are condensation/evaporation (and associated aging), 
!    renaming, coagulation, and nucleation
      real(r8), dimension( 1:pcols, 1:gas_pcnst, 1:nqtendaa ) :: &
         q_coltendaa
      real(r8), dimension( 1:pcols, 1:gas_pcnst, 1:nqqcwtendaa ) :: &
         qqcw_coltendaa

#if ( defined( MOSAIC_SPECIES ) )
      real(r8) :: cnvrg_fail(pcols,pver) !BSINGH -  For tracking MOSAIC convergence failures
      real(r8) :: max_kelvin_iter(pcols,pver)  !BSINGH -  For tracking when max is hit for kelvin iterations
      real(r8) :: xnerr_astem_negative(pcols,pver,5,4)
#endif

      type ( misc_vars_aa_type ) :: misc_vars_aa




      adjust_factor_pbl_ratenucl = newnuc_adjust_factor_pbl

#if ( defined CAM_VERSION_IS_ACME ) 
      history_aerocom = .false.
#else
      call phys_getopts( history_aerocom_out        = history_aerocom )
#endif


      icol_diag = -1
      if (ldiag1 > 0) then
      if (nstep < 3) then
         do i = 1, ncol
!           if ((latndx(i) == 23) .and. (lonndx(i) == 37)) icol_diag = i
            if ((latndx(i) == 47) .and. (lonndx(i) ==121)) icol_diag = i  ! amazon
         end do
      end if
      end if

      do_cond   = ( mdo_gasaerexch > 0 )
      do_rename = ( mdo_rename > 0 )
      do_newnuc = ( mdo_newnuc > 0 )
      do_coag   = ( mdo_coag > 0 )

      q_coltendaa = 0.0_r8 ; qqcw_coltendaa = 0.0_r8
      nufine_3dtend_nnuc = 0.0_r8
      ncluster_3dtend_nnuc = 0.0_r8
      soag_3dtend_cond = 0.0_r8

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
! these variables otherwise undefined
      q_tendbb = 0.0_r8 ; qqcw_tendbb = 0.0_r8
#endif

#if ( defined( MOSAIC_SPECIES ) )
      cnvrg_fail(1:pcols,1:pver) = 0.0_r8
      max_kelvin_iter(1:pcols,1:pver)  = 0.0_r8
      xnerr_astem_negative(1:pcols,1:pver,1:5,1:4) = 0.0_r8
#endif

! turn off history selectively for comparison with dd06f
      if ( (.not. do_cond) .and. (.not. do_rename) ) then
         do_q_coltendaa(:,iqtend_cond) = .false.
         do_q_coltendaa(:,iqtend_rnam) = .false.
         do_qqcw_coltendaa(:,iqqcwtend_rnam) = .false.
      end if
      if ( .not. do_newnuc ) then
         do_q_coltendaa(:,iqtend_nnuc) = .false.
      end if
      if ( .not. do_coag ) then
         do_q_coltendaa(:,iqtend_coag) = .false.
      end if

! get saturation mixing ratio
      call qsat( t(1:ncol,1:pver), pmid(1:ncol,1:pver), &
                 ev_sat(1:ncol,1:pver), qv_sat(1:ncol,1:pver) )

main_k_loop: &
      do k = top_lev, pver
main_i_loop: &
      do i = 1, ncol

      if ( ldiag13n ) lun13n = 129 + i


!
! determine the number of sub-areas, their fractional areas, and relative humidities
!
! if cloud fraction ~= 0, the grid-cell has a single clear  sub-area      (nsubarea = 1)
! if cloud fraction ~= 1, the grid-cell has a single cloudy sub-area      (nsubarea = 1)
! otherwise,              the grid-cell has a clear and a cloudy sub-area (nsubarea = 2)
!
      if (cld(i,k) < fcld_locutoff) then
! note that cloud chemistry is only on when cld(i,k) >= 1.0e-5_r8
! it may be that the macrophysics has a higher threshold that this
         fcldy = 0.0_r8
         nsubarea = 1 ; ncldy_subarea = 0
         jclea = 1 ; jcldy = 0
      else if (cld(i,k) > fcld_hicutoff) then
         fcldy = 1.0_r8
         nsubarea = 1 ; ncldy_subarea = 1
         jclea = 0 ; jcldy = 1
      else
         fcldy = cld(i,k)
         nsubarea = 2 ; ncldy_subarea = 1
         jclea = 1 ; jcldy = 2
      end if
      fclea = 1.0_r8 - fcldy
      fcldybb = max( cld(i,k), 1.0e-6_r8 )

      iscldy_subarea(:) = .false.
      if (jcldy > 0) iscldy_subarea(jcldy) = .true.

      afracsub(:) = 0.0_r8
      if (jclea > 0) afracsub(jclea) = fclea
      if (jcldy > 0) afracsub(jcldy) = fcldy

      relhumgcm = max( 0.0_r8, min( 1.0_r8, qv(i,k)/qv_sat(i,k) ) )
      if (ncldy_subarea <= 0) then
         relhumsub(:) = relhumgcm
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      else if (cldy_rh_sameas_clear > 0) then
         relhumsub(:) = relhumgcm
#endif
      else
         relhumsub(jcldy) = 1.0_r8
         if (jclea > 0) then
            tmpa = (relhumgcm - afracsub(jcldy))/afracsub(jclea)
            relhumsub(jclea) = max( 0.0_r8, min( 1.0_r8, tmpa ) )
         end if
      end if

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( ldiag13n ) then
      write(lun13n,'(/a,3i5)') 'modal_aero_amicphys_intr mapping at nstep, i, k', nstep, i, k
      write(lun13n,'(a,1p,5i12  )') 'jclea, jcldy, ncldy, nsubc', &
         jclea, jcldy, ncldy_subarea, nsubarea
      write(lun13n,'(a,1p,5e12.4)') 'cld, fcldy, fcldybb, fclea', &
         cld(i,k), fcldy, fcldybb, fclea
      write(lun13n,'(a,1p,5e12.4)') 'relhumav, relhumsub(1:2)  ', &
         relhumgcm, relhumsub(1:2)
      end if
#endif


      do lmz = 1, gas_pcnst
         qgcm1(lmz)    = max( 0.0_r8, q_pregaschem(i,k,lmz) )
         qgcm2(lmz)    = max( 0.0_r8, q_precldchem(i,k,lmz) )
         qqcwgcm2(lmz) = max( 0.0_r8, qqcw_precldchem(i,k,lmz) )
         qgcm3(lmz)    = max( 0.0_r8, q(i,k,lmz) )
         qqcwgcm3(lmz) = max( 0.0_r8, qqcw(i,k,lmz) )
      end do
      qaerwatgcm3(:) = 0.0_r8
      if ( present( qaerwat ) ) then
         qaerwatgcm3(1:ntot_amode) = max( 0.0_r8, qaerwat(i,k,1:ntot_amode) )
      end if

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      n = min( maxsubarea, nsubarea+1 )
#else
      n = nsubarea
#endif
      qsub1(:,1:n) = 0.0_r8
      qsub2(:,1:n) = 0.0_r8
      qsub3(:,1:n) = 0.0_r8
      qsub4(:,1:n) = 0.0_r8
      qqcwsub1(:,1:n) = 0.0_r8
      qqcwsub2(:,1:n) = 0.0_r8
      qqcwsub3(:,1:n) = 0.0_r8
      qqcwsub4(:,1:n) = 0.0_r8
      qaerwatsub3(:,1:n) = 0.0_r8
      qaerwatsub4(:,1:n) = 0.0_r8

!
! calculate initial (i.e., before cond/rnam/nnuc/coag) tracer mixing ratios within the sub-areas
!    for all-clear or all-cloudy cases, the sub-area TMRs are equal to the grid-cell means
!    for partly cloudy case, they are different.  This is primarily because the
!       interstitial aerosol mixing ratios are assumed lower in the cloudy sub-area than in
!       the clear sub-area, because much of the aerosol is activated in the cloudy sub-area.
!
      if ( (jclea > 0) .and. (jcldy > 0) .and. &
           (jclea+jcldy == 3) .and. (nsubarea == 2) ) then
! partly cloudy case

! set gas mixing ratios in sub-areas (for the condensing gases only!!)
         do lmz = 1, gas_pcnst
            if (lmapcc_all(lmz) /= lmapcc_val_gas) cycle

            ! assume gas in both sub-areas before gas-chem and cloud-chem equal grid-cell mean
            qsub1(lmz,1:nsubarea) = qgcm1(lmz)
            qsub2(lmz,1:nsubarea) = qgcm2(lmz)

            ! assume gas in clear sub-area after cloud-chem equals before cloud-chem value
            qsub3(lmz,jclea) = qsub2(lmz,jclea)
            ! gas in cloud sub-area then determined by grid-cell mean and clear values
            qsub3(lmz,jcldy) = (qgcm3(lmz) - fclea*qsub3(lmz,jclea))/fcldy
            ! check that this does not produce a negative value
            if (qsub3(lmz,jcldy) < 0.0_r8) then
               qsub3(lmz,jcldy) = 0.0_r8
               qsub3(lmz,jclea) = qgcm3(lmz)/fclea
            end if
         end do

! set aerosol mixing ratios in sub-areas
         do n = 1, ntot_amode

         do l2 = 0, nspec_amode(n)

            if (l2 <= 1) then
            ! calculcate partitioning factors
               if (l2 == 0) then
                  la = numptr_amode(n) - loffset
                  lc = numptrcw_amode(n) - loffset
                  tmp_qa_gcav = qgcm2(la)
                  tmp_qc_gcav = qqcwgcm2(lc)
               else
                  tmp_qa_gcav = 0.0_r8
                  tmp_qc_gcav = 0.0_r8
                  do l3 = 1, nspec_amode(n)
                     la = lmassptr_amode(l3,n) - loffset
                     tmp_qa_gcav = tmp_qa_gcav + qgcm2(la)
                     lc = lmassptrcw_amode(l3,n) - loffset
                     tmp_qc_gcav = tmp_qc_gcav + qqcwgcm2(lc)
                  end do
               end if

               tmp_qc_cldy = tmp_qc_gcav/fcldy
               tmp_qa_cldy = max( 0.0_r8, ((tmp_qa_gcav+tmp_qc_gcav) - tmp_qc_cldy) )
               tmp_qa_clea = (tmp_qa_gcav - fcldy*tmp_qa_cldy)/fclea

               ! *** question ***
               !    use same tmp_aa_clea/cldy for everything ?
               !    use one for number and one for all masses (based on total mass) ?
               !    use separate ones for everything ?
               ! maybe one for number and one for all masses is best,
               !    because number and mass have different activation fractions
               ! *** question ***
               tmp_aa = max( 1.e-35_r8, tmp_qa_clea*fclea ) / max( 1.e-35_r8, tmp_qa_gcav )
               tmp_aa = max( 0.0_r8, min( 1.0_r8, tmp_aa ) )
               tmp_aa_clea = tmp_aa/fclea
               tmp_aa_cldy = (1.0_r8-tmp_aa)/fcldy

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
               if ( n <= 2 .and. ldiag13n ) then
                  if (n==1 .and. l2==0) write(lun13n,'(a)')
                  write(lun13n,'(a,2i3, 1p,6e12.4)') 'n, l2, tmp_aa, tmp_aa_clea, tmp_aa_cldy', &
                     n, l2, tmp_aa, tmp_aa_clea, tmp_aa_cldy
                  tmpa = 1.0e-6_r8*mwhost_num/mwdry
                  if (l2 > 0) tmpa = 1.0e9_r8
                  write(lun13n,'(a, 6x, 1p,6e12.4)') 'qct, qcy, qat, qay, qax, qtt           ', &
                     tmpa*tmp_qc_gcav, tmpa*tmp_qc_cldy, tmpa*tmp_qa_gcav, &
                     tmpa*tmp_qa_cldy, tmpa*tmp_qa_clea, tmpa*(tmp_qc_gcav+tmp_qa_gcav)
               end if
#endif
            end if ! (l2 <= 1)

            if (l2 == 0) then
               la = numptr_amode(n) - loffset
               lc = numptrcw_amode(n) - loffset
            else
               la = lmassptr_amode(l2,n) - loffset
               lc = lmassptrcw_amode(l2,n) - loffset
            end if

            qsub2(la,jclea) = qgcm2(la)*tmp_aa_clea
            qsub2(la,jcldy) = qgcm2(la)*tmp_aa_cldy
            qqcwsub2(lc,jclea) = 0.0_r8
            qqcwsub2(lc,jcldy) = qqcwgcm2(lc)/fcldy

            qsub3(la,jclea) = qgcm3(la)*tmp_aa_clea
            qsub3(la,jcldy) = qgcm3(la)*tmp_aa_cldy
            qqcwsub3(lc,jclea) = 0.0_r8
            qqcwsub3(lc,jcldy) = qqcwgcm3(lc)/fcldy

            end do ! l2
         end do ! n

      else if ((jclea == 1) .and. (jcldy == 0) .and. (nsubarea == 1)) then
! all clear, or cld < 1e-5
! in this case, fclea=1 and fcldy=0
!
! put all the gases and interstitial aerosols in the clear sub-area
!    and set mix-ratios = 0 in cloudy sub-area
! for cloud-borne aerosol, do nothing
!    because the grid-cell-mean cloud-borne aerosol will be left unchanged
!    (i.e., this routine only changes qqcw when cld >= 1e-5)
!
         do lmz = 1, gas_pcnst
            if (lmapcc_all(lmz) <= 0) cycle
            qsub1(lmz,jclea) = qgcm1(lmz)
            qsub2(lmz,jclea) = qgcm2(lmz)
            qsub3(lmz,jclea) = qgcm3(lmz)
            qqcwsub2(lmz,jclea) = qqcwgcm2(lmz)
            qqcwsub3(lmz,jclea) = qqcwgcm3(lmz)
         end do

      else if ((jclea == 0) .and. (jcldy == 1) .and. (nsubarea == 1)) then
! all cloudy, or cld > 0.999
! in this case, fcldy= and fclea=0
!
! put all the gases and interstitial aerosols in the cloudy sub-area
!    and set mix-ratios = 0 in clear sub-area
!
         do lmz = 1, gas_pcnst
            if (lmapcc_all(lmz) <= 0) cycle
            qsub1(lmz,jcldy) = qgcm1(lmz)
            qsub2(lmz,jcldy) = qgcm2(lmz)
            qsub3(lmz,jcldy) = qgcm3(lmz)
            qqcwsub2(lmz,jcldy) = qqcwgcm2(lmz)
            qqcwsub3(lmz,jcldy) = qqcwgcm3(lmz)
         end do

      else
! this should not happen
         write(tmp_str,'(a,3(1x,i10))') &
            '*** modal_aero_amicphys - bad jclea, jcldy, nsubarea', &
            jclea, jcldy, nsubarea
         call endrun( tmp_str )
      end if

! aerosol water -- how to treat this in sub-areas needs more work/thinking
! currently modal_aero_water_uptake calculates qaerwat using
!    the grid-cell mean interstital-aerosol mix-rats and the clear-area rh
      do jsub = 1, nsubarea
         qaerwatsub3(:,jsub) = qaerwatgcm3(:)
      end do

      if (nsubarea == 1) then
! the j=2 subarea is used for some diagnostics
! but is not used in actual calculations
         j = 2
         qsub1(:,j) = 0.0_r8
         qsub2(:,j) = 0.0_r8
         qsub3(:,j) = 0.0_r8
         qqcwsub2(:,j) = 0.0_r8
         qqcwsub3(:,j) = 0.0_r8
      end if


! diagnostics after forming sub-areas
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( ldiag13n ) then
      do l2 = 1, 4
         if (l2 == 1) then
            igas = igas_h2so4
         else if (l2 == 3) then
            igas = igas_nh3
            if (igas <= 0) cycle
         else if (l2 == 4) then
            igas = -3
         else
            igas = 1
         end if
         if (igas > 0) then
            l = lmap_gas(igas)
            tmpch6a = name_gas(igas)
         else
            l = -igas
            tmpch6a = cnst_name(l+loffset)
         end if
         tmpa = 1.0e9
         write(lun13n,'(a)')
         write(lun13n,'(2a,1p,4e12.4)') tmpch6a, ' host  1-3', &
            q_pregaschem(i,k,l)*tmpa, q_precldchem(i,k,l)*tmpa, q(i,k,l)*tmpa
         write(lun13n,'(2a,1p,4e12.4)') tmpch6a, ' gm    1-3', &
            qgcm1(l)*tmpa, qgcm2(l)*tmpa, qgcm3(l)*tmpa
         j = jclea ; if (j <= 0) j = nsubarea+1
         write(lun13n,'(2a,1p,4e12.4)') tmpch6a, ' clear 1-3', &
            qsub1(l,j)*tmpa, qsub2(l,j)*tmpa, qsub3(l,j)*tmpa
         j = jcldy ; if (j <= 0) j = nsubarea+1
         write(lun13n,'(2a,1p,4e12.4)') tmpch6a, ' cloud 1-3', &
            qsub1(l,j)*tmpa, qsub2(l,j)*tmpa, qsub3(l,j)*tmpa
      end do ! l2

      n = 1
      do l2 = 1, 3
         if (l2 == 1) then
            tmpa = 1.0e-6_r8/28.966_r8
            la = lmap_num(n)
            lc = lmap_numcw(n)
            tmpch6a = name_num(n)
            tmpch6c = name_numcw(n)
         else 
            if (l2 == 2) then
               iaer = iaer_so4
            else
               iaer = iaer_soa
            end if
            tmpa = 1.0e9_r8
            la = lmap_aer(iaer,n)
            lc = lmap_aercw(iaer,n)
            tmpch6a = name_aer(iaer,n)
            tmpch6c = name_aercw(iaer,n)
         end if
         write(lun13n,'(a)')
         write(lun13n,'(4a,1p,2(2x,2e12.4))') tmpch6a, ' host  2-3; ', tmpch6c, ' ...', &
            q_precldchem(i,k,la)*tmpa, q(i,k,la)*tmpa, &
            qqcw_precldchem(i,k,lc)*tmpa, qqcw(i,k,lc)*tmpa
         write(lun13n,'(4a,1p,2(2x,2e12.4))') tmpch6a, ' gm    2-3; ', tmpch6c, ' ...', &
            qgcm2(la)*tmpa, qgcm3(la)*tmpa, &
            qqcwgcm2(lc)*tmpa, qqcwgcm3(lc)*tmpa
         j = jclea ; if (j <= 0) j = nsubarea+1
         write(lun13n,'(4a,1p,2(2x,2e12.4))') tmpch6a, ' clear 2-3; ', tmpch6c, ' ...', &
            qsub2(la,j)*tmpa, qsub3(la,j)*tmpa, &
            qqcwsub2(lc,j)*tmpa, qqcwsub3(lc,j)*tmpa
         j = jcldy ; if (j <= 0) j = nsubarea+1
         write(lun13n,'(4a,1p,2(2x,2e12.4))') tmpch6a, ' cloud 2-3; ', tmpch6c, ' ...', &
            qsub2(la,j)*tmpa, qsub3(la,j)*tmpa, &
            qqcwsub2(lc,j)*tmpa, qqcwsub3(lc,j)*tmpa
      end do ! l2
      end if ! ( ldiag13n )
#endif

!
! start integration
!
      do n = 1, max_mode
         if (n <= ntot_amode) then
            dgn_a(n)           = dgncur_a(i,k,n)
            dgn_awet(n)        = dgncur_awet(i,k,n)
            wetdens(n)         = max( 1000.0_r8, wetdens_host(i,k,n) )
         else
            dgn_a(n) = 0.0_r8
            dgn_awet(n) = 0.0_r8
            wetdens(n) = 1000.0_r8
         end if
      end do

      misc_vars_aa%ncluster_tend_nnuc_1grid = ncluster_3dtend_nnuc(i,k)
#if ( defined ( MOSAIC_SPECIES ) )
      misc_vars_aa%cnvrg_fail_1grid = cnvrg_fail(i,k)
      misc_vars_aa%max_kelvin_iter_1grid = max_kelvin_iter(i,k)
      misc_vars_aa%xnerr_astem_negative_1grid(1:5,1:4) = xnerr_astem_negative(pcols,pver,1:5,1:4)
#endif


      lund = iulog  ! for cambox, iulog=93 at this point

!      ubroutine mam_amicphys_1gridcell(          &
!        do_cond,            do_rename,           &
!        do_newnuc,          do_coag,             &
!        nstep,    lchnk,    i,         k,        &
!        latndx,   lonndx,   lund,                &
!        loffset,  deltat,                        &
!        nsubarea,  ncldy_subarea,                &
!        iscldy_subarea,     afracsub,            &
!        temp,     pmid,     pdel,                &
!        zmid,     pblh,     relhumsub,           &
!        dgn_a,    dgn_awet, wetdens,             &
!        qsub1,                                   &
!        qsub2, qqcwsub2,                         &
!        qsub3, qqcwsub3,                         &
!        qsub4, qqcwsub4,                         &
!        qsub_tendaa, qqcwsub_tendaa              )

      call mam_amicphys_1gridcell(                &
         do_cond,             do_rename,          &
         do_newnuc,           do_coag,            &
         nstep,    lchnk,     i,         k,       &
         latndx(i),           lonndx(i), lund,    &
         loffset,  deltat,                        &
         nsubarea,  ncldy_subarea,                &
         iscldy_subarea,      afracsub,           &
         t(i,k),   pmid(i,k), pdel(i,k),          &
         zm(i,k),  pblh(i),   relhumsub,          &
         dgn_a,    dgn_awet,  wetdens,            &
         qsub1,                                   &
         qsub2, qqcwsub2,                         &
         qsub3, qqcwsub3, qaerwatsub3,            &
         qsub4, qqcwsub4, qaerwatsub4,            &
         qsub_tendaa, qqcwsub_tendaa,             &
         misc_vars_aa                             )


!
! form new grid-mean mix-ratios
!
      if (nsubarea == 1) then
         qgcm4(:) = qsub4(:,1)
         qgcm_tendaa(:,:) = qsub_tendaa(:,:,1)
         qaerwatgcm4(1:ntot_amode) = qaerwatsub4(1:ntot_amode,1)
      else
         qgcm4(:) = 0.0_r8
         qgcm_tendaa(:,:) = 0.0_r8
         do j = 1, nsubarea
            qgcm4(:) = qgcm4(:) + qsub4(:,j)*afracsub(j)
            qgcm_tendaa(:,:) = qgcm_tendaa(:,:) + qsub_tendaa(:,:,j)*afracsub(j)
         end do
         ! for aerosol water use the clear sub-area value
         qaerwatgcm4(1:ntot_amode) = qaerwatsub4(1:ntot_amode,jclea)
      end if

      if (ncldy_subarea <= 0) then
         qqcwgcm4(:) = qqcwgcm3(:)
         qqcwgcm_tendaa(:,:) = 0.0_r8
      else if (nsubarea == 1) then
         qqcwgcm4(:) = qqcwsub4(:,1)
         qqcwgcm_tendaa(:,:) = qqcwsub_tendaa(:,:,1)
      else
         qqcwgcm4(:) = 0.0_r8
         qqcwgcm_tendaa(:,:) = 0.0_r8
         do j = 1, nsubarea
            if ( .not. iscldy_subarea(j) ) cycle
            qqcwgcm4(:) = qqcwgcm4(:) + qqcwsub4(:,j)*afracsub(j)
            qqcwgcm_tendaa(:,:) = qqcwgcm_tendaa(:,:) + qqcwsub_tendaa(:,:,j)*afracsub(j)
         end do
      end if

      do lmz = 1, gas_pcnst
         if (lmapcc_all(lmz) > 0) then
            q(i,k,lmz) = max(qgcm4(lmz),0.0_r8)  ! HW, to ensure non-negative
            if (lmapcc_all(lmz) >= lmapcc_val_aer) then
               qqcw(i,k,lmz) = max(qqcwgcm4(lmz),0.0_r8)  !HW, to ensure non-negative
            end if
         end if
      end do

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if (iqtend_cond <= nqtendbb) q_tendbb(i,k,:,iqtend_cond) = qgcm_tendaa(:,iqtend_cond)
      if (iqtend_rnam <= nqtendbb) q_tendbb(i,k,:,iqtend_rnam) = qgcm_tendaa(:,iqtend_rnam)
      if (iqtend_nnuc <= nqtendbb) q_tendbb(i,k,:,iqtend_nnuc) = qgcm_tendaa(:,iqtend_nnuc)
      if (iqtend_coag <= nqtendbb) q_tendbb(i,k,:,iqtend_coag) = qgcm_tendaa(:,iqtend_coag)
      if (iqqcwtend_rnam <= nqqcwtendbb) qqcw_tendbb(i,k,:,iqqcwtend_rnam) = qqcwgcm_tendaa(:,iqqcwtend_rnam)
#endif
      if ( update_qaerwat > 0 .and. present( qaerwat ) ) then
         qaerwat(i,k,1:ntot_amode) = qaerwatgcm4(1:ntot_amode)
      end if


! diagnostics after forming sub-areas
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( ldiag13n ) then
      do l2 = 1, 4
         if (l2 == 1) then
            igas = igas_h2so4
         else if (l2 == 3) then
            igas = igas_nh3
            if (igas <= 0) cycle
         else if (l2 == 4) then
            igas = -3
         else
            igas = 1
         end if
         if (igas > 0) then
            l = lmap_gas(igas)
            tmpch6a = name_gas(igas)
         else
            l = -igas
            tmpch6a = cnst_name(l+loffset)
         end if
         tmpa = 1.0e9
         write(lun13n,'(a)')
         write(lun13n,'(2a,1p,4e12.4)') tmpch6a, ' host  1-4', &
            q_pregaschem(i,k,l)*tmpa, q_precldchem(i,k,l)*tmpa, 0.0, q(i,k,l)*tmpa
         write(lun13n,'(2a,1p,4e12.4)') tmpch6a, ' gm    1-4', &
            qgcm1(l)*tmpa, qgcm2(l)*tmpa, qgcm3(l)*tmpa, qgcm4(l)*tmpa
         j = jclea ; if (j <= 0) j = nsubarea+1
         write(lun13n,'(2a,1p,4e12.4)') tmpch6a, ' clear 1-4', &
            qsub1(l,j)*tmpa, qsub2(l,j)*tmpa, qsub3(l,j)*tmpa, qsub4(l,j)*tmpa
         j = jcldy ; if (j <= 0) j = nsubarea+1
         write(lun13n,'(2a,1p,4e12.4)') tmpch6a, ' cloud 1-4', &
            qsub1(l,j)*tmpa, qsub2(l,j)*tmpa, qsub3(l,j)*tmpa, qsub4(l,j)*tmpa
      end do ! l2
      end if ! ( ldiag13n )
#endif


! increment column tendencies
      pdel_fac = pdel(i,k)/gravit
      do iqtend = 1, nqtendaa
      do l = 1, gas_pcnst
         if ( do_q_coltendaa(l,iqtend) ) then
            q_coltendaa(i,l,iqtend) = q_coltendaa(i,l,iqtend) + qgcm_tendaa(l,iqtend)*pdel_fac
         end if
         if (iqtend <= nqqcwtendaa) then
         if ( do_qqcw_coltendaa(l,iqtend) ) then
            qqcw_coltendaa(i,l,iqtend) = qqcw_coltendaa(i,l,iqtend) + qqcwgcm_tendaa(l,iqtend)*pdel_fac
         end if
         end if
      end do ! l
      end do ! iqtend

      if ( history_aerocom ) then
         ! 3d soa tendency for aerocom
         ! note that flux units (kg/m2/s) are used here instead of tendency units (kg/kg/s or kg/m3/s)
         do jsoa = 1, nsoa
            l = lptr2_soa_g_amode(jsoa) - loffset
            soag_3dtend_cond(i,k,jsoa) = qgcm_tendaa(l,iqtend_cond)*(adv_mass(l)/mwdry)*(pdel(i,k)/gravit)
         end do
         ! 3d number nucleation tendency for aerocom - units are (#/m3/s)
         ! so multiply qgcm_tendaa (#/kmol/s) by air molar density (kmol/m3)
         l = numptr_amode(nait) - loffset
         nufine_3dtend_nnuc(i,k) = qgcm_tendaa(l,iqtend_nnuc) * (pmid(i,k)/(r_universal*t(i,k)))
      end if


      ncluster_3dtend_nnuc(i,k) = misc_vars_aa%ncluster_tend_nnuc_1grid
#if ( defined ( MOSAIC_SPECIES ) )
      cnvrg_fail(i,k) = misc_vars_aa%cnvrg_fail_1grid 
      max_kelvin_iter(i,k) = misc_vars_aa%max_kelvin_iter_1grid 
      xnerr_astem_negative(pcols,pver,1:5,1:4) = misc_vars_aa%xnerr_astem_negative_1grid(1:5,1:4)
#endif

      end do main_i_loop

      end do main_k_loop


! output column tendencies to history
! the ordering here is to allow comparison of fort.90 files from box model testing
!    but is not important for regular cam simulations
      do ipass = 1, 3

         if (ipass == 1) then
            itmpa = iqtend_cond ; itmpb = iqtend_rnam
            itmpc = iqqcwtend_rnam ; itmpd = iqqcwtend_rnam
         else if (ipass == 2) then
            itmpa = iqtend_nnuc ; itmpb = iqtend_nnuc
            itmpc = 0 ; itmpd = 0
         else
            itmpa = iqtend_coag ; itmpb = iqtend_coag
            itmpc = 0 ; itmpd = 0
         end if

         do l = 1, gas_pcnst
            do iqtend = itmpa, itmpb
               if (iqtend <= 0) cycle
               if ( do_q_coltendaa(l,iqtend) ) then
                  q_coltendaa(1:ncol,l,iqtend) = q_coltendaa(1:ncol,l,iqtend)*(adv_mass(l)/mwdry)
                  fieldname = trim(cnst_name(l+loffset)) // suffix_q_coltendaa(iqtend)
                  call outfld( fieldname, q_coltendaa(1:ncol,l,iqtend), ncol, lchnk )
               end if
            end do ! iqtend
            do iqqcwtend = itmpc, itmpd
               if (iqqcwtend <= 0) cycle
               if ( do_qqcw_coltendaa(l,iqqcwtend) ) then
                  qqcw_coltendaa(1:ncol,l,iqqcwtend) = qqcw_coltendaa(1:ncol,l,iqqcwtend)* (adv_mass(l)/mwdry)
                  fieldname = trim(cnst_name_cw(l+loffset)) // suffix_qqcw_coltendaa(iqqcwtend)
                  call outfld( fieldname, qqcw_coltendaa(1:ncol,l,iqqcwtend), ncol, lchnk )
               end if
            end do ! iqqcwtend
         end do ! l

         if ( ipass==1 .and. history_aerocom ) then
            do jsoa = 1, nsoa
               l = lptr2_soa_g_amode(jsoa)
               fieldname = trim(cnst_name(l)) // '_sfgaex3d'
               call outfld( fieldname, soag_3dtend_cond(1:ncol,:,jsoa), ncol, lchnk )
            end do
            l = numptr_amode(nait)
            fieldname = trim(cnst_name(l)) // '_nuc1'
            call outfld( fieldname, nufine_3dtend_nnuc(1:ncol,:), ncol, lchnk )
            fieldname = trim(cnst_name(l)) // '_nuc2'
            call outfld( fieldname, ncluster_3dtend_nnuc(1:ncol,:), ncol, lchnk )
         end if

      end do ! ipass

#if ( defined( MOSAIC_SPECIES ) )
      if ( mosaic ) then
         !BSINGH - output MOSAIC convergence fail tracking:
         call outfld( 'convergence_fail', cnvrg_fail(1:ncol,:), ncol, lchnk )
         call outfld( 'max_kelvin_iter' , max_kelvin_iter(1:ncol,:),  ncol, lchnk )

         do n = 1, 4
         do m = 1, 5
            fieldname = ' '
            write( fieldname(1:16), '(a,i1,a,i1)') 'astem_negval_', m, '_', n
            call outfld( fieldname, xnerr_astem_negative(1:ncol,1:pver,m,n), ncol, lchnk )
         end do
         end do
      end if
#endif

      return
!EOC
      end subroutine modal_aero_amicphys_intr


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine mam_amicphys_1gridcell(          &
         do_cond,            do_rename,           &
         do_newnuc,          do_coag,             &
         nstep,    lchnk,    i,        k,         &
         latndx,   lonndx,   lund,                &
         loffset,  deltat,                        &
         nsubarea,  ncldy_subarea,                &
         iscldy_subarea,     afracsub,            &
         temp,     pmid,     pdel,                &
         zmid,     pblh,     relhumsub,           &
         dgn_a,    dgn_awet, wetdens,             &
         qsub1,                                   &
         qsub2, qqcwsub2,                         &
         qsub3, qqcwsub3, qaerwatsub3,            &
         qsub4, qqcwsub4, qaerwatsub4,            &
         qsub_tendaa, qqcwsub_tendaa,             &
         misc_vars_aa                             )
!
! calculates changes to gas and aerosol sub-area TMRs (tracer mixing ratios)
!    for the current grid cell (with indices = lchnk,i,k)
! qsub3 and qqcwsub3 are the incoming current TMRs
! qsub4 and qqcwsub4 are the outgoing updated TMRs
!
      logical,  intent(in)    :: do_cond, do_rename, do_newnuc, do_coag
      logical,  intent(in)    :: iscldy_subarea(maxsubarea)

      integer,  intent(in)    :: nstep                 ! model time-step number
      integer,  intent(in)    :: lchnk                 ! chunk identifier
      integer,  intent(in)    :: i, k                  ! column and level indices
      integer,  intent(in)    :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in)    :: lund                  ! logical unit for diagnostic output
      integer,  intent(in)    :: loffset
      integer,  intent(in)    :: nsubarea, ncldy_subarea

      real(r8), intent(in)    :: deltat                ! time step (s)
      real(r8), intent(in)    :: afracsub(maxsubarea)   ! sub-area fractional area (0-1)

      real(r8), intent(in)    :: temp                  ! temperature at model levels (K)
      real(r8), intent(in)    :: pmid                  ! pressure at layer center (Pa)
      real(r8), intent(in)    :: pdel                  ! pressure thickness of layer (Pa)
      real(r8), intent(in)    :: zmid                  ! altitude (above ground) at layer center (m)
      real(r8), intent(in)    :: pblh                  ! planetary boundary layer depth (m)
      real(r8), intent(in)    :: relhumsub(maxsubarea)  ! sub-area relative humidity (0-1)
      real(r8), intent(inout) :: dgn_a(max_mode)
      real(r8), intent(inout) :: dgn_awet(max_mode)
                                    ! dry & wet geo. mean dia. (m) of number distrib.
      real(r8), intent(inout) :: wetdens(max_mode)
                                    ! interstitial aerosol wet density (kg/m3)
                                    ! dry & wet geo. mean dia. (m) of number distrib.

! qsubN and qqcwsubN (N=1:4) are tracer mixing ratios (TMRs, mol/mol or #/kmol) in sub-areas
!    currently there are just clear and cloudy sub-areas
!    the N=1:4 have same meanings as for qgcmN
!    N=1 - before gas-phase chemistry
!    N=2 - before cloud chemistry
!    N=3 - incoming values (before gas-aerosol exchange, newnuc, coag)
!    N=4 - outgoing values (after  gas-aerosol exchange, newnuc, coag)
      real(r8), intent(in   ), dimension( 1:gas_pcnst, 1:maxsubarea ) :: &
         qsub1, qsub2, qsub3, qqcwsub2, qqcwsub3 
      real(r8), intent(inout), dimension( 1:gas_pcnst, 1:maxsubarea ) :: &
         qsub4, qqcwsub4
      real(r8), intent(inout), dimension( 1:ntot_amode_extd, 1:maxsubarea ) :: &
         qaerwatsub3, qaerwatsub4   ! aerosol water mixing ratios (mol/mol)
! qsub_tendaa and qqcwsub_tendaa are TMR tendencies
!    for different processes, which are used to produce history output
! the processes are condensation/evaporation (and associated aging), 
!    renaming, coagulation, and nucleation
      real(r8), intent(inout), dimension( 1:gas_pcnst, 1:nqtendaa, 1:maxsubarea ) :: &
         qsub_tendaa
      real(r8), intent(inout), dimension( 1:gas_pcnst, 1:nqqcwtendaa, 1:maxsubarea ) :: &
         qqcwsub_tendaa
      type ( misc_vars_aa_type ), intent(inout) :: misc_vars_aa

! local
      integer :: iaer, igas
      integer :: jsub
      integer :: l
      integer :: n
      logical :: do_cond_sub, do_rename_sub, do_newnuc_sub, do_coag_sub
      logical :: do_map_gas_sub

      real(r8), dimension( 1:max_gas ) :: &
         qgas1, qgas2, qgas3, qgas4
      real(r8), dimension( 1:max_mode ) :: &
         qnum2, qnum3, qnum4, &
         qnumcw2, qnumcw3, qnumcw4
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer2, qaer3, qaer4, &
         qaercw2, qaercw3, qaercw4
      real(r8), dimension( 1:max_mode ) :: &
         qwtr3, qwtr4

      real(r8), dimension( 1:max_gas, 1:nqtendaa ) :: &
         qgas_delaa
      real(r8), dimension( 1:max_mode, 1:nqtendaa ) :: &
         qnum_delaa
      real(r8), dimension( 1:max_mode, 1:nqqcwtendaa ) :: &
         qnumcw_delaa
      real(r8), dimension( 1:max_aer, 1:max_mode, 1:nqtendaa ) :: &
         qaer_delaa
      real(r8), dimension( 1:max_aer, 1:max_mode, 1:nqqcwtendaa ) :: &
         qaercw_delaa

      real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf, tmpn

      type ( misc_vars_aa_type ), dimension(nsubarea) :: misc_vars_aa_sub


! the q--4 values will be equal to q--3 values unless they get changed
      qsub4(:,1:nsubarea) = qsub3(:,1:nsubarea)
      qqcwsub4(:,1:nsubarea) = qqcwsub3(:,1:nsubarea)
      qaerwatsub4(:,1:nsubarea) = qaerwatsub3(:,1:nsubarea)

      qsub_tendaa(:,:,1:nsubarea) = 0.0_r8
      qqcwsub_tendaa(:,:,1:nsubarea) = 0.0_r8

      do jsub = 1, nsubarea
          misc_vars_aa_sub(jsub) = misc_vars_aa
      end do


main_jsub_loop: &
      do jsub = 1, nsubarea

      if ( iscldy_subarea(jsub) .eqv. .true. ) then
         do_cond_sub   = do_cond
         do_rename_sub = do_rename
         do_newnuc_sub = .false.
         do_coag_sub   = .false.
         if (mdo_gaexch_cldy_subarea <= 0) do_cond_sub = .false.
      else
         do_cond_sub   = do_cond
         do_rename_sub = do_rename
         do_newnuc_sub = do_newnuc
         do_coag_sub   = do_coag
      end if
      do_map_gas_sub = do_cond_sub .or. do_newnuc_sub
      

! map incoming sub-area mix-ratios to gas/aer/num arrays

      qgas1(:) = 0.0_r8
      qgas2(:) = 0.0_r8
      qgas3(:) = 0.0_r8
      qgas4(:) = 0.0_r8
      if ( do_map_gas_sub .eqv. .true. ) then
! for cldy subarea, only do gases if doing gaexch
      do igas = 1, ngas
         l = lmap_gas(igas)
         qgas1(igas) = qsub1(l,jsub)*fcvt_gas(igas)
         qgas2(igas) = qsub2(l,jsub)*fcvt_gas(igas)
         qgas3(igas) = qsub3(l,jsub)*fcvt_gas(igas)
         qgas4(igas) = qgas3(igas)
      end do
      end if

      qaer2(:,:) = 0.0_r8
      qnum2(:)   = 0.0_r8
      qaer3(:,:) = 0.0_r8
      qnum3(:)   = 0.0_r8
      qaer4(:,:) = 0.0_r8
      qnum4(:)   = 0.0_r8
      qwtr3(:)   = 0.0_r8
      qwtr4(:)   = 0.0_r8
      do n = 1, ntot_amode
         l = lmap_num(n)
         qnum2(n) = qsub2(l,jsub)*fcvt_num
         qnum3(n) = qsub3(l,jsub)*fcvt_num
         qnum4(n) = qnum3(n)
         do iaer = 1, naer
            l = lmap_aer(iaer,n)
            if (l > 0) then
               qaer2(iaer,n) = qsub2(l,jsub)*fcvt_aer(iaer)
               qaer3(iaer,n) = qsub3(l,jsub)*fcvt_aer(iaer)
               qaer4(iaer,n) = qaer3(iaer,n)
            end if
         end do
         qwtr3(n) = qaerwatsub3(n,jsub)*fcvt_wtr
         qwtr4(n) = qwtr3(n)
      end do ! n

      if ( iscldy_subarea(jsub) .eqv. .true. ) then
! only do cloud-borne for cloudy
      qaercw2(:,:) = 0.0_r8
      qnumcw2(:)   = 0.0_r8
      qaercw3(:,:) = 0.0_r8
      qnumcw3(:)   = 0.0_r8
      qaercw4(:,:) = 0.0_r8
      qnumcw4(:)   = 0.0_r8
      do n = 1, ntot_amode
         l = lmap_numcw(n)
         qnumcw2(n) = qqcwsub2(l,jsub)*fcvt_num
         qnumcw3(n) = qqcwsub3(l,jsub)*fcvt_num
         qnumcw4(n) = qnumcw3(n)
         do iaer = 1, naer
            l = lmap_aercw(iaer,n)
            if (l > 0) then
               qaercw2(iaer,n) = qqcwsub2(l,jsub)*fcvt_aer(iaer)
               qaercw3(iaer,n) = qqcwsub3(l,jsub)*fcvt_aer(iaer)
               qaercw4(iaer,n) = qaercw3(iaer,n)
            end if
         end do
      end do ! n
      end if


      if ( iscldy_subarea(jsub) .eqv. .true. ) then

      call mam_amicphys_1subarea_cloudy(             &
         do_cond_sub,            do_rename_sub,      &
         do_newnuc_sub,          do_coag_sub,        &
         nstep,      lchnk,      i,        k,        &
         latndx,     lonndx,     lund,               &
         loffset,    deltat,                         &
         jsub,                   nsubarea,           &
         iscldy_subarea(jsub),   afracsub(jsub),     &
         temp,       pmid,       pdel,               &
         zmid,       pblh,       relhumsub(jsub),    &
         dgn_a,      dgn_awet,   wetdens,            &
         qgas1,      qgas3,      qgas4,              &
         qgas_delaa,                                 &
         qnum3,      qnum4,                          &
         qnum_delaa,                                 &
         qaer2,      qaer3,      qaer4,              &
         qaer_delaa,                                 &
         qwtr3,      qwtr4,                          &
         qnumcw3,    qnumcw4,                        &
         qnumcw_delaa,                               &
         qaercw2,    qaercw3,    qaercw4,            &
         qaercw_delaa,                               &
         misc_vars_aa_sub(jsub)                      )

      else

      call mam_amicphys_1subarea_clear(              &
         do_cond_sub,            do_rename_sub,      &
         do_newnuc_sub,          do_coag_sub,        &
         nstep,      lchnk,      i,        k,        &
         latndx,     lonndx,     lund,               &
         loffset,    deltat,                         &
         jsub,                   nsubarea,           &
         iscldy_subarea(jsub),   afracsub(jsub),     &
         temp,       pmid,       pdel,               &
         zmid,       pblh,       relhumsub(jsub),    &
         dgn_a,      dgn_awet,   wetdens,            &
         qgas1,      qgas3,      qgas4,              &
         qgas_delaa,                                 &
         qnum3,      qnum4,      qnum_delaa,         &
         qaer3,      qaer4,      qaer_delaa,         &
         qwtr3,      qwtr4,                          &
         misc_vars_aa_sub(jsub)                      )

      end if

      if ((nsubarea == 1) .or. (iscldy_subarea(jsub) .eqv. .false.)) then
         misc_vars_aa%ncluster_tend_nnuc_1grid = misc_vars_aa%ncluster_tend_nnuc_1grid &
                                               + misc_vars_aa_sub(jsub)%ncluster_tend_nnuc_1grid*afracsub(jsub)
#if ( defined ( MOSAIC_SPECIES ) )
         misc_vars_aa%cnvrg_fail_1grid      = misc_vars_aa_sub(jsub)%cnvrg_fail_1grid 
         misc_vars_aa%max_kelvin_iter_1grid = misc_vars_aa_sub(jsub)%max_kelvin_iter_1grid 
         misc_vars_aa%xnerr_astem_negative_1grid(1:5,1:4) = misc_vars_aa_sub(jsub)%xnerr_astem_negative_1grid(1:5,1:4)
#endif
      end if



! map gas/aer/num arrays (mix-ratio and del=change) back to sub-area arrays

      if ( do_map_gas_sub .eqv. .true. ) then
      do igas = 1, ngas
         l = lmap_gas(igas)
         qsub4(l,jsub) = qgas4(igas)/fcvt_gas(igas)
         qsub_tendaa(l,:,jsub) = qgas_delaa(igas,:)/(fcvt_gas(igas)*deltat)
      end do
      end if

      do n = 1, ntot_amode
         l = lmap_num(n)
         qsub4(l,jsub) = qnum4(n)/fcvt_num
         qsub_tendaa(l,:,jsub) = qnum_delaa(n,:)/(fcvt_num*deltat)
         do iaer = 1, naer
            l = lmap_aer(iaer,n)
            if (l > 0) then
               qsub4(l,jsub) = qaer4(iaer,n)/fcvt_aer(iaer)
               qsub_tendaa(l,:,jsub) = qaer_delaa(iaer,n,:)/(fcvt_aer(iaer)*deltat)
            end if
         end do
         qaerwatsub4(n,jsub) = qwtr4(n)/fcvt_wtr

         if ( iscldy_subarea(jsub) ) then
         l = lmap_numcw(n)
         qqcwsub4(l,jsub) = qnumcw4(n)/fcvt_num
         qqcwsub_tendaa(l,:,jsub) = qnumcw_delaa(n,:)/(fcvt_num*deltat)
         do iaer = 1, naer
            l = lmap_aercw(iaer,n)
            if (l > 0) then
               qqcwsub4(l,jsub) = qaercw4(iaer,n)/fcvt_aer(iaer)
               qqcwsub_tendaa(l,:,jsub) = qaercw_delaa(iaer,n,:)/(fcvt_aer(iaer)*deltat)
            end if
         end do
         end if
      end do ! n


      end do main_jsub_loop



      return
      end subroutine mam_amicphys_1gridcell


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine mam_amicphys_1subarea_cloudy(       &
         do_cond,                do_rename,          &
         do_newnuc,              do_coag,            &
         nstep,      lchnk,      i,        k,        &
         latndx,     lonndx,     lund,               &
         loffset,    deltat,                         &
         jsub,                   nsubarea,           &
         iscldy_subarea,         afracsub,           &
         temp,       pmid,       pdel,               &
         zmid,       pblh,       relhum,             &
         dgn_a,      dgn_awet,   wetdens,            &
         qgas1,      qgas3,      qgas4,              &
         qgas_delaa,                                 &
         qnum3,      qnum4,                          &
         qnum_delaa,                                 &
         qaer2,      qaer3,      qaer4,              &
         qaer_delaa,                                 &
         qwtr3,      qwtr4,                          &
         qnumcw3,    qnumcw4,                        &
         qnumcw_delaa,                               &
         qaercw2,    qaercw3,    qaercw4,            &
         qaercw_delaa,                               &
         misc_vars_aa_sub                            )
!
! calculates changes to gas and aerosol sub-area TMRs (tracer mixing ratios)
!    for a single cloudy sub-area (with indices = lchnk,i,k,jsub)
! qgas3, qaer3, qaercw3, qnum3, qnumcw3 are the current incoming TMRs
! qgas4, qaer4, qaercw4, qnum4, qnumcw4 are the updated outgoing TMRs
!
! when do_cond = false, this routine only calculates changes involving
!    growth from smaller to larger modes (renaming) following cloud chemistry
!    so gas TMRs are not changed
! when do_cond = true, this routine also calculates changes involving
!    gas-aerosol exchange (condensation/evaporation)
!    transfer of particles from hydrophobic modes to hydrophilic modes (aging)
!       due to condensation
! currently this routine does not do
!    new particle nucleation - because h2so4 gas conc. should be very low in cloudy air
!    coagulation - because cloud-borne aerosol would need to be included
!
      use physconst, only:  r_universal

      logical,  intent(in)    :: do_cond, do_rename, do_newnuc, do_coag
      logical,  intent(in)    :: iscldy_subarea        ! true if sub-area is cloudy
      integer,  intent(in)    :: lchnk                 ! chunk identifier
      integer,  intent(in)    :: nstep                 ! model time-step number
      integer,  intent(in)    :: i, k                  ! column and level indices
      integer,  intent(in)    :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in)    :: lund                  ! logical unit for diagnostic output
      integer,  intent(in)    :: loffset
      integer,  intent(in)    :: jsub, nsubarea        ! sub-area index, number of sub-areas

      real(r8), intent(in)    :: afracsub              ! fractional area of sub-area (0-1)
      real(r8), intent(in)    :: deltat                ! time step (s)

      real(r8), intent(in)    :: temp                  ! temperature at model levels (K)
      real(r8), intent(in)    :: pmid                  ! pressure at layer center (Pa)
      real(r8), intent(in)    :: pdel                  ! pressure thickness of layer (Pa)
      real(r8), intent(in)    :: zmid                  ! altitude (above ground) at layer center (m)
      real(r8), intent(in)    :: pblh                  ! planetary boundary layer depth (m)
      real(r8), intent(in)    :: relhum                ! relative humidity (0-1)

      real(r8), intent(inout) :: dgn_a(max_mode)
      real(r8), intent(inout) :: dgn_awet(max_mode)
                                    ! dry & wet geo. mean dia. (m) of number distrib.
      real(r8), intent(inout) :: wetdens(max_mode)
                                    ! interstitial aerosol wet density (kg/m3)
                                    ! dry & wet geo. mean dia. (m) of number distrib.

! qXXXN (X=gas,aer,wat,num; N=1:4) are sub-area mixing ratios
!    XXX=gas - gas species
!    XXX=aer - aerosol mass  species (excluding water)
!    XXX=wat - aerosol water
!    XXX=num - aerosol number
!    N=1 - before gas-phase chemistry
!    N=2 - before cloud chemistry
!    N=3 - current incoming values (before gas-aerosol exchange, newnuc, coag)
!    N=4 - updated outgoing values (after  gas-aerosol exchange, newnuc, coag)
!
! qXXX_delaa are TMR changes (not tendencies)
!    for different processes, which are used to produce history output
! for a clear sub-area, the processes are condensation/evaporation (and associated aging),
!    renaming, coagulation, and nucleation
      real(r8), intent(in   ), dimension( 1:max_gas ) :: &
         qgas1, qgas3
      real(r8), intent(inout), dimension( 1:max_gas ) :: &
         qgas4
      real(r8), intent(inout), dimension( 1:max_gas, 1:nqtendaa ) :: &
         qgas_delaa

      real(r8), intent(in   ), dimension( 1:max_mode ) :: &
         qnum3
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum4
      real(r8), intent(inout), dimension( 1:max_mode, 1:nqtendaa ) :: &
         qnum_delaa

      real(r8), intent(in   ), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer2, qaer3
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer4
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode, 1:nqtendaa ) :: &
         qaer_delaa

      real(r8), intent(in   ), dimension( 1:max_mode ) :: &
         qwtr3
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr4

      real(r8), intent(in   ), dimension( 1:max_mode ) :: &
         qnumcw3
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnumcw4
      real(r8), intent(inout), dimension( 1:max_mode, 1:nqqcwtendaa ) :: &
         qnumcw_delaa

      real(r8), intent(in   ), dimension( 1:max_aer, 1:max_mode ) :: &
         qaercw2, qaercw3
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaercw4
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode, 1:nqqcwtendaa ) :: &
         qaercw_delaa

      type ( misc_vars_aa_type ), intent(inout) :: misc_vars_aa_sub

! local
      integer, parameter :: ntot_poaspec = npoa
      integer, parameter :: ntot_soaspec = nsoa

      integer :: iaer, igas, ip
      integer :: jtsubstep
      integer :: ll
      integer :: modefrm, modetoo
! if mtoo_renamexf(n) >  0, then mode n gets renamed into mode mtoo_renamexf(n)
! if mtoo_renamexf(n) <= 0, then mode n does not have renaming
      integer :: mtoo_renamexf(max_mode)
      integer :: n, ntsubstep
      integer :: n_mode
      integer :: ntot_soamode

      logical, parameter :: flag_pcarbon_opoa_frac_zero   = .true.
      logical, parameter :: flag_nh4_lt_2so4_each_step  = .false.

      logical :: skip_soamode(max_mode)   ! true if this mode does not have soa

      real(r8), dimension( 1:max_gas ) :: &
         qgas_cur, qgas_sv1, qgas_avg
      real(r8), dimension( 1:max_gas ) :: &
         qgas_del_cond, qgas_del_nnuc, qgas_netprod_otrproc
                ! qgas_netprod_otrproc = gas net production rate from other processes
                !    such as gas-phase chemistry and emissions (mol/mol/s)
                ! this allows the condensation (gasaerexch) routine to apply production and condensation loss
                !    together, which is more accurate numerically
                ! NOTE - must be >= zero, as numerical method can fail when it is negative
                ! NOTE - currently only the values for h2so4 and nh3 should be non-zero

! qxxx_del_yyyy    are mix-ratio changes over full time step (deltat)
! qxxx_delsub_yyyy are mix-ratio changes over time sub-step (dtsubstep)
      real(r8), dimension( 1:max_mode ) :: &
         qnum_cur, qnum_sv1
      real(r8), dimension( 1:max_mode ) :: &
         qnum_del_cond, qnum_del_rnam, qnum_del_nnuc, qnum_del_coag, &
         qnum_delsub_cond, qnum_delsub_coag

      real(r8), dimension( 1:max_mode ) :: &
         qnumcw_cur, qnumcw_sv1
      real(r8), dimension( 1:max_mode ) :: &
         qnumcw_del_rnam

      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur, qaer_sv1
      real(r8), dimension( 1:max_aer, 1:max_agepair ) :: &
         qaer_delsub_coag_in
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_del_cond, qaer_del_rnam, qaer_del_nnuc, qaer_del_coag, &
         qaer_delsub_grow4rnam, &
         qaer_delsub_cond, qaer_delsub_coag

      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaercw_cur, qaercw_sv1
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaercw_del_rnam, &
         qaercw_delsub_grow4rnam

      real(r8), dimension( 1:max_mode ) :: &
         qwtr_cur

      real(r8) :: aircon                           ! air molar density (kmol/m3)
      real(r8) :: del_h2so4_gasprod
      real(r8) :: del_h2so4_aeruptk
      real(r8) :: dnclusterdt
      real(r8) :: dtsubstep                        ! time sub-step
      real(r8) :: gas_diffus(max_gas)              ! gas diffusivity at current temp and pres (m2/s)
      real(r8) :: gas_freepath(max_gas)            ! gas mean free path at current temp and pres (m)

      real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf
      real(r8) :: tmp_relhum
      real(r8) :: uptkaer(max_gas,max_mode)
      real(r8) :: uptkrate_h2so4



! air molar density (kmol/m3)
      aircon = pmid/(r_universal*temp)

      n_mode = ntot_amode

      qgas_cur = qgas3
      qaer_cur = qaer3
      qnum_cur = qnum3
      qwtr_cur = qwtr3
      qnumcw_cur = qnumcw3
      qaercw_cur = qaercw3


      qgas_netprod_otrproc(:) = 0.0_r8
      if ( ( do_cond                         ) .and. &
           ( gaexch_h2so4_uptake_optaa == 2 ) ) then
         do igas = 1, ngas
            if ((igas == igas_h2so4) .or. (igas == igas_nh3)) then
! if gaexch_h2so4_uptake_optaa == 2, then
!    if qgas increases from pre-gaschem to post-cldchem,
!       start from the pre-gaschem mix-ratio and add in the production
!       during the integration
!    if it decreases, 
!       start from post-cldchem mix-ratio
! *** currently just do this for h2so4 and nh3
               qgas_netprod_otrproc(igas) = (qgas3(igas) - qgas1(igas))/deltat
               if ( qgas_netprod_otrproc(igas) >= 0.0_r8 ) then
                  qgas_cur(igas) = qgas1(igas)
               else
                  qgas_netprod_otrproc(igas) = 0.0_r8
               end if
            end if
         end do ! igas
      end if


      qgas_del_cond = 0.0_r8
      qgas_del_nnuc = 0.0_r8

      qaer_del_cond = 0.0_r8
      qaer_del_rnam = 0.0_r8
      qaer_del_nnuc = 0.0_r8
      qaer_del_coag = 0.0_r8
      qaer_delsub_cond = 0.0_r8

      qaercw_del_rnam = 0.0_r8

      qnum_del_cond = 0.0_r8
      qnum_del_rnam = 0.0_r8
      qnum_del_nnuc = 0.0_r8
      qnum_del_coag = 0.0_r8
      qnum_delsub_cond = 0.0_r8

      qnumcw_del_rnam = 0.0_r8

      dnclusterdt = 0.0_r8


      ntsubstep = 1
      dtsubstep = deltat
      if (ntsubstep > 1) dtsubstep = deltat/ntsubstep

      del_h2so4_gasprod = max( qgas3(igas_h2so4)-qgas1(igas_h2so4), 0.0_r8 )/ntsubstep

!
!
! loop over multiple time sub-steps
!
!
jtsubstep_loop: &
      do jtsubstep = 1, ntsubstep


!
!
! gas-aerosol exchange
!
!
      uptkrate_h2so4 = 0.0_r8
do_cond_if_block10: &
      if ( do_cond ) then

      qgas_sv1 = qgas_cur
      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

#if ( defined( MOSAIC_SPECIES ) )
      if ( mosaic ) then
         tmp_relhum = min( relhum, 0.98_r8 )
         call mosaic_gasaerexch_1subarea_intr(     nstep,                &!Intent(ins)
              lchnk,             i,                k,           jsub,    &
              temp,              tmp_relhum,       pmid,                 &
              aircon,            dtsubstep,        n_mode,               &
              dgn_a,             dgn_awet,         qaer_cur,             &!Intent(inouts)
              qgas_cur,          qnum_cur,         qwtr_cur,             &
              qgas_avg,          qgas_netprod_otrproc,                   &
              uptkrate_h2so4,    misc_vars_aa_sub                        )
      else
#endif
         call mam_gasaerexch_1subarea(                                &
           nstep,             lchnk,                                  &
           i,                 k,                jsub,                 &
           jtsubstep,         ntsubstep,                              &
           latndx,            lonndx,           lund,                 &
           dtsubstep,                                                 &
           temp,              pmid,             aircon,               &
           n_mode,                                                    &
           qgas_cur,          qgas_avg,                               &
           qgas_netprod_otrproc,                                      &
           qaer_cur,                                                  &
           qnum_cur,                                                  &
           qwtr_cur,                                                  &
           dgn_a,             dgn_awet,         wetdens,              &
           uptkaer,           uptkrate_h2so4                          )
#if ( defined( MOSAIC_SPECIES ) )
      end if
#endif

      if (newnuc_h2so4_conc_optaa == 11) then
         qgas_avg(igas_h2so4) = 0.5_r8*(qgas_sv1(igas_h2so4) + qgas_cur(igas_h2so4))
      else if (newnuc_h2so4_conc_optaa == 12) then
         qgas_avg(igas_h2so4) = qgas_cur(igas_h2so4)
      end if

      qgas_del_cond = qgas_del_cond + (qgas_cur - (qgas_sv1 + qgas_netprod_otrproc*dtsubstep))
      qnum_delsub_cond = qnum_cur - qnum_sv1
      qaer_delsub_cond = qaer_cur - qaer_sv1
! qaer_delsub_grow4rnam = change in qaer_del_cond during latest condensation calculations
      qaer_delsub_grow4rnam = qaer_cur - qaer_sv1

      del_h2so4_aeruptk = qgas_cur(igas_h2so4) &
                       - (qgas_sv1(igas_h2so4) + qgas_netprod_otrproc(igas_h2so4)*dtsubstep)

      else ! do_cond_if_block10

      qgas_avg(1:ngas) = qgas_cur(1:ngas)
      qaer_delsub_grow4rnam(:,:) = 0.0_r8

      del_h2so4_aeruptk = 0.0_r8

      end if do_cond_if_block10


!
!
! renaming after "continuous growth"
!
!
do_rename_if_block30: &
      if ( do_rename ) then

      mtoo_renamexf(:) = 0
      mtoo_renamexf(nait) = nacc

! qaer_delsub_grow4rnam   = change in qaer from cloud chemistry and gas condensation
! qaercw_delsub_grow4rnam = change in qaercw from cloud chemistry
      qaer_delsub_grow4rnam   = (qaer3 - qaer2)/ntsubstep  + qaer_delsub_grow4rnam
      qaercw_delsub_grow4rnam = (qaercw3 - qaercw2)/ntsubstep

      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur
      qnumcw_sv1 = qnumcw_cur
      qaercw_sv1 = qaercw_cur

      call mam_rename_1subarea(                                      &
         nstep,             lchnk,                                   &
         i,                 k,                jsub,                  &
         latndx,            lonndx,           lund,                  &
         iscldy_subarea,                                             &
         mtoo_renamexf,                                              &
         n_mode,                                                     &
         qnum_cur,                                                   &
         qaer_cur,          qaer_delsub_grow4rnam,                   &
         qwtr_cur,                                                   &
         qnumcw_cur,                                                 &
         qaercw_cur,        qaercw_delsub_grow4rnam                  )

      qnum_del_rnam = qnum_del_rnam + (qnum_cur - qnum_sv1)
      qaer_del_rnam = qaer_del_rnam + (qaer_cur - qaer_sv1)
      qnumcw_del_rnam = qnumcw_del_rnam + (qnumcw_cur - qnumcw_sv1)
      qaercw_del_rnam = qaercw_del_rnam + (qaercw_cur - qaercw_sv1)

      end if do_rename_if_block30


!
!
! primary carbon aging
!
!
      if ( ( n_agepair > 0        ) .and. &
           ( do_cond .eqv. .true. ) ) then

      qaer_delsub_coag_in = 0.0_r8
      qaer_delsub_coag = 0.0_r8
      qnum_delsub_coag = 0.0_r8

      call mam_pcarbon_aging_1subarea(                              &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         dtsubstep,         dgn_a,            do_cond,              &
         n_mode,                                                    &
         qnum_cur,          qnum_delsub_cond, qnum_delsub_coag,     &
         qaer_cur,          qaer_delsub_cond, qaer_delsub_coag,     &
         qaer_delsub_coag_in,                                       &
         qwtr_cur                                                   )

      end if


! accumulate sub-step q-dels
      if ( do_cond ) then
         qnum_del_cond = qnum_del_cond + qnum_delsub_cond
         qaer_del_cond = qaer_del_cond + qaer_delsub_cond
      end if

      end do jtsubstep_loop


!
!
! final mix ratios
!
!
      qgas4 = qgas_cur
      qaer4 = qaer_cur
      qnum4 = qnum_cur
      qwtr4 = qwtr_cur
      qnumcw4 = qnumcw_cur
      qaercw4 = qaercw_cur

! final mix ratio changes

      qgas_delaa(:,iqtend_cond) = qgas_del_cond(:)
      qgas_delaa(:,iqtend_rnam) = 0.0_r8
      qgas_delaa(:,iqtend_nnuc) = 0.0_r8
      qgas_delaa(:,iqtend_coag) = 0.0_r8

      qnum_delaa(:,iqtend_cond) = qnum_del_cond(:)
      qnum_delaa(:,iqtend_rnam) = qnum_del_rnam(:)
      qnum_delaa(:,iqtend_nnuc) = 0.0_r8
      qnum_delaa(:,iqtend_coag) = 0.0_r8

      qaer_delaa(:,:,iqtend_cond) = qaer_del_cond(:,:)
      qaer_delaa(:,:,iqtend_rnam) = qaer_del_rnam(:,:)
      qaer_delaa(:,:,iqtend_nnuc) = 0.0_r8
      qaer_delaa(:,:,iqtend_coag) = 0.0_r8

      qnumcw_delaa(:,iqqcwtend_rnam) = qnumcw_del_rnam(:)

      qaercw_delaa(:,:,iqqcwtend_rnam) = qaercw_del_rnam(:,:)

      misc_vars_aa_sub%ncluster_tend_nnuc_1grid = dnclusterdt

      return
      end subroutine mam_amicphys_1subarea_cloudy


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine mam_amicphys_1subarea_clear(        &
         do_cond,                do_rename,          &
         do_newnuc,              do_coag,            &
         nstep,      lchnk,      i,        k,        &
         latndx,     lonndx,     lund,               &
         loffset,    deltat,                         &
         jsub,                   nsubarea,           &
         iscldy_subarea,         afracsub,           &
         temp,       pmid,       pdel,               &
         zmid,       pblh,       relhum,             &
         dgn_a,      dgn_awet,   wetdens,            &
         qgas1,      qgas3,      qgas4,              &
         qgas_delaa,                                 &
         qnum3,      qnum4,      qnum_delaa,         &
         qaer3,      qaer4,      qaer_delaa,         &
         qwtr3,      qwtr4,                          &
         misc_vars_aa_sub                            )
!
! calculates changes to gas and aerosol sub-area TMRs (tracer mixing ratios)
!    for a single clear sub-area (with indices = lchnk,i,k,jsub)
! qgas3, qaer3, qnum3 are the current incoming TMRs
! qgas4, qaer4, qnum4 are the updated outgoing TMRs
!
! this routine calculates changes involving 
!    gas-aerosol exchange (condensation/evaporation)
!    growth from smaller to larger modes (renaming) due to condensation
!    new particle nucleation
!    coagulation
!    transfer of particles from hydrophobic modes to hydrophilic modes (aging)
!       due to condensation and coagulation
!
      use physconst, only:  r_universal

      logical,  intent(in)    :: do_cond, do_rename, do_newnuc, do_coag
      logical,  intent(in)    :: iscldy_subarea        ! true if sub-area is cloudy
      integer,  intent(in)    :: lchnk                 ! chunk identifier
      integer,  intent(in)    :: nstep                 ! model time-step number
      integer,  intent(in)    :: i, k                  ! column and level indices
      integer,  intent(in)    :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in)    :: lund                  ! logical unit for diagnostic output
      integer,  intent(in)    :: loffset
      integer,  intent(in)    :: jsub, nsubarea        ! sub-area index, number of sub-areas

      real(r8), intent(in)    :: afracsub              ! fractional area of sub-area (0-1)
      real(r8), intent(in)    :: deltat                ! time step (s)

      real(r8), intent(in)    :: temp                  ! temperature at model levels (K)
      real(r8), intent(in)    :: pmid                  ! pressure at layer center (Pa)
      real(r8), intent(in)    :: pdel                  ! pressure thickness of layer (Pa)
      real(r8), intent(in)    :: zmid                  ! altitude (above ground) at layer center (m)
      real(r8), intent(in)    :: pblh                  ! planetary boundary layer depth (m)
      real(r8), intent(in)    :: relhum                ! relative humidity (0-1)

      real(r8), intent(inout) :: dgn_a(max_mode)
      real(r8), intent(inout) :: dgn_awet(max_mode)
                                    ! dry & wet geo. mean dia. (m) of number distrib.
      real(r8), intent(inout) :: wetdens(max_mode)
                                    ! interstitial aerosol wet density (kg/m3)
                                    ! dry & wet geo. mean dia. (m) of number distrib.

! qXXXN (X=gas,aer,wat,num; N=1:4) are sub-area mixing ratios
!    XXX=gas - gas species
!    XXX=aer - aerosol mass  species (excluding water)
!    XXX=wat - aerosol water
!    XXX=num - aerosol number
!    N=1 - before gas-phase chemistry
!    N=2 - before cloud chemistry
!    N=3 - current incoming values (before gas-aerosol exchange, newnuc, coag)
!    N=4 - updated outgoing values (after  gas-aerosol exchange, newnuc, coag)
!
! qXXX_delaa are TMR changes (not tendencies)
!    for different processes, which are used to produce history output
! for a clear sub-area, the processes are condensation/evaporation (and associated aging), 
!    renaming, coagulation, and nucleation
      real(r8), intent(in   ), dimension( 1:max_gas ) :: &
         qgas1, qgas3
      real(r8), intent(inout), dimension( 1:max_gas ) :: &
         qgas4
      real(r8), intent(inout), dimension( 1:max_gas, 1:nqtendaa ) :: &
         qgas_delaa

      real(r8), intent(in   ), dimension( 1:max_mode ) :: &
         qnum3
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum4
      real(r8), intent(inout), dimension( 1:max_mode, 1:nqtendaa ) :: &
         qnum_delaa

      real(r8), intent(in   ), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer3
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer4
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode, 1:nqtendaa ) :: &
         qaer_delaa

      real(r8), intent(in   ), dimension( 1:max_mode ) :: &
         qwtr3
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr4

      type ( misc_vars_aa_type ), intent(inout) :: misc_vars_aa_sub

! local
      integer, parameter :: ntot_poaspec = npoa
      integer, parameter :: ntot_soaspec = nsoa

      integer :: iaer, igas, ip
      integer :: jtsubstep
      integer :: ll
      integer :: modefrm, modetoo
! if mtoo_renamexf(n) >  0, then mode n gets renamed into mode mtoo_renamexf(n)
! if mtoo_renamexf(n) <= 0, then mode n does not have renaming
      integer :: mtoo_renamexf(max_mode)
      integer :: n, ntsubstep
      integer :: n_mode
      integer :: ntot_soamode

      logical, parameter :: flag_pcarbon_opoa_frac_zero   = .true.
      logical, parameter :: flag_nh4_lt_2so4_each_step  = .false.

      logical :: skip_soamode(max_mode)   ! true if this mode does not have soa

      real(r8), dimension( 1:max_gas ) :: &
         qgas_cur, qgas_sv1, qgas_avg
      real(r8), dimension( 1:max_gas ) :: &
         qgas_del_cond, qgas_del_nnuc, qgas_netprod_otrproc
                ! qgas_netprod_otrproc = gas net production rate from other processes
                !    such as gas-phase chemistry and emissions (mol/mol/s)
                ! this allows the condensation (gasaerexch) routine to apply production and condensation loss 
                !    together, which is more accurate numerically
                ! NOTE - must be >= zero, as numerical method can fail when it is negative
                ! NOTE - currently only the values for h2so4 and nh3 should be non-zero

! qxxx_del_yyyy    are mix-ratio changes over full time step (deltat)
! qxxx_delsub_yyyy are mix-ratio changes over time sub-step (dtsubstep)
      real(r8), dimension( 1:max_mode ) :: &
         qnum_cur, qnum_sv1
      real(r8), dimension( 1:max_mode ) :: &
         qnum_del_cond, qnum_del_rnam, qnum_del_nnuc, qnum_del_coag, &
         qnum_delsub_cond, qnum_delsub_coag

      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur, qaer_sv1
      real(r8), dimension( 1:max_aer, 1:max_agepair ) :: &
         qaer_delsub_coag_in
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_del_cond, qaer_del_rnam, qaer_del_nnuc, qaer_del_coag, &
         qaer_delsub_grow4rnam, &
         qaer_delsub_cond, qaer_delsub_coag

      real(r8), dimension( 1:max_mode ) :: &
         qwtr_cur

      real(r8) :: aircon                           ! air molar density (kmol/m3)
      real(r8) :: del_h2so4_gasprod
      real(r8) :: del_h2so4_aeruptk
      real(r8) :: dnclusterdt, dnclusterdt_substep
      real(r8) :: dtsubstep                        ! time sub-step
      real(r8) :: gas_diffus(max_gas)              ! gas diffusivity at current temp and pres (m2/s)
      real(r8) :: gas_freepath(max_gas)            ! gas mean free path at current temp and pres (m)

      real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf
      real(r8) :: uptkaer(max_gas,max_mode)
      real(r8) :: uptkrate_h2so4



! air molar density (kmol/m3)
      aircon = pmid/(r_universal*temp)

      n_mode = ntot_amode

      qgas_cur = qgas3
      qaer_cur = qaer3
      qnum_cur = qnum3
      qwtr_cur = qwtr3

      qgas_netprod_otrproc(:) = 0.0_r8
      if ( ( do_cond                         ) .and. &
           ( gaexch_h2so4_uptake_optaa == 2 ) ) then
         do igas = 1, ngas
            if ((igas == igas_h2so4) .or. (igas == igas_nh3)) then
! if gaexch_h2so4_uptake_optaa == 2, then
!    if qgas increases from pre-gaschem to post-cldchem,
!       start from the pre-gaschem mix-ratio and add in the production
!       during the integration
!    if it decreases, 
!       start from post-cldchem mix-ratio
! *** currently just do this for h2so4 and nh3
               qgas_netprod_otrproc(igas) = (qgas3(igas) - qgas1(igas))/deltat
               if ( qgas_netprod_otrproc(igas) >= 0.0_r8 ) then
                  qgas_cur(igas) = qgas1(igas)
               else
                  qgas_netprod_otrproc(igas) = 0.0_r8
               end if
            end if
         end do ! igas
      end if

      qgas_del_cond = 0.0_r8
      qgas_del_nnuc = 0.0_r8

      qaer_del_cond = 0.0_r8
      qaer_del_rnam = 0.0_r8
      qaer_del_nnuc = 0.0_r8
      qaer_del_coag = 0.0_r8
      qaer_delsub_coag_in = 0.0_r8
      qaer_delsub_cond = 0.0_r8
      qaer_delsub_coag = 0.0_r8

      qnum_del_cond = 0.0_r8
      qnum_del_rnam = 0.0_r8
      qnum_del_nnuc = 0.0_r8
      qnum_del_coag = 0.0_r8
      qnum_delsub_cond = 0.0_r8
      qnum_delsub_coag = 0.0_r8

      dnclusterdt = 0.0_r8


      ntsubstep = 1
      dtsubstep = deltat
      if (ntsubstep > 1) dtsubstep = deltat/ntsubstep

      del_h2so4_gasprod = max( qgas3(igas_h2so4)-qgas1(igas_h2so4), 0.0_r8 )/ntsubstep

!
!
! loop over multiple time sub-steps
!
!
jtsubstep_loop: &
      do jtsubstep = 1, ntsubstep


!
!
! gas-aerosol exchange
!
!
      uptkrate_h2so4 = 0.0_r8
do_cond_if_block10: &
      if ( do_cond ) then

      qgas_sv1 = qgas_cur
      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

#if ( defined( MOSAIC_SPECIES ) )
      if ( mosaic ) then
         call mosaic_gasaerexch_1subarea_intr(     nstep,                &!Intent(ins)
              lchnk,             i,                k,           jsub,    &
              temp,              relhum,           pmid,                 &
              aircon,            dtsubstep,        n_mode,               &
              dgn_a,             dgn_awet,         qaer_cur,             &!Intent(inouts)
              qgas_cur,          qnum_cur,         qwtr_cur,             &
              qgas_avg,          qgas_netprod_otrproc,                   &
              uptkrate_h2so4,    misc_vars_aa_sub                        )
      else
#endif
         call mam_gasaerexch_1subarea(                                &
           nstep,             lchnk,                                  &
           i,                 k,                jsub,                 &
           jtsubstep,         ntsubstep,                              &
           latndx,            lonndx,           lund,                 &
           dtsubstep,                                                 &
           temp,              pmid,             aircon,               &
           n_mode,                                                    &
           qgas_cur,          qgas_avg,                               &
           qgas_netprod_otrproc,                                      &
           qaer_cur,                                                  &
           qnum_cur,                                                  &
           qwtr_cur,                                                  &
           dgn_a,             dgn_awet,         wetdens,              &
           uptkaer,           uptkrate_h2so4                          )
#if ( defined( MOSAIC_SPECIES ) )
      end if
#endif
         
      if (newnuc_h2so4_conc_optaa == 11) then
         qgas_avg(igas_h2so4) = 0.5_r8*(qgas_sv1(igas_h2so4) + qgas_cur(igas_h2so4))
      else if (newnuc_h2so4_conc_optaa == 12) then
         qgas_avg(igas_h2so4) = qgas_cur(igas_h2so4)
      end if

      qgas_del_cond = qgas_del_cond + (qgas_cur - (qgas_sv1 + qgas_netprod_otrproc*dtsubstep))
      qnum_delsub_cond = qnum_cur - qnum_sv1
      qaer_delsub_cond = qaer_cur - qaer_sv1
! qaer_del_grow4rnam = change in qaer_del_cond during latest condensation calculations
      qaer_delsub_grow4rnam = qaer_cur - qaer_sv1

      del_h2so4_aeruptk = qgas_cur(igas_h2so4) &
                       - (qgas_sv1(igas_h2so4) + qgas_netprod_otrproc(igas_h2so4)*dtsubstep)

      else ! do_cond_if_block10

      qgas_avg(1:ngas) = qgas_cur(1:ngas)
      qaer_delsub_grow4rnam(:,:) = 0.0_r8

      del_h2so4_aeruptk = 0.0_r8

      end if do_cond_if_block10


!
!
! renaming after "continuous growth"
!
!
do_rename_if_block30: &
      if ( do_rename ) then

      mtoo_renamexf(:) = 0
      mtoo_renamexf(nait) = nacc

      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

      call mam_rename_1subarea(                                    &
         nstep,             lchnk,                                 &
         i,                 k,                jsub,                &
         latndx,            lonndx,           lund,                &
         iscldy_subarea,                                           &
         mtoo_renamexf,                                            &
         n_mode,                                                   &
         qnum_cur,                                                 &
         qaer_cur,          qaer_delsub_grow4rnam,                 &
         qwtr_cur                                                  )

      qnum_del_rnam = qnum_del_rnam + (qnum_cur - qnum_sv1)
      qaer_del_rnam = qaer_del_rnam + (qaer_cur - qaer_sv1)

      end if do_rename_if_block30


!
!
! new particle formation (nucleation)
!
!
do_newnuc_if_block50: &
      if ( do_newnuc ) then

      qgas_sv1 = qgas_cur
      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

      call mam_newnuc_1subarea(                                     &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         dtsubstep,                                                 &
         temp,              pmid,             aircon,               &
         zmid,              pblh,             relhum,               &
         uptkrate_h2so4,   del_h2so4_gasprod, del_h2so4_aeruptk,    &
         n_mode,                                                    &
         qgas_cur,          qgas_avg,                               &
         qnum_cur,                                                  &
         qaer_cur,                                                  &
         qwtr_cur,                                                  &
         dnclusterdt_substep                                        )

      qgas_del_nnuc = qgas_del_nnuc + (qgas_cur - qgas_sv1)
      qnum_del_nnuc = qnum_del_nnuc + (qnum_cur - qnum_sv1)
      qaer_del_nnuc = qaer_del_nnuc + (qaer_cur - qaer_sv1)
      dnclusterdt = dnclusterdt + dnclusterdt_substep*(dtsubstep/deltat)

      end if do_newnuc_if_block50


!
!
! coagulation part
!
!
      if ( do_coag ) then

      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

      call mam_coag_1subarea(                                       &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         dtsubstep,                                                 &
         temp,              pmid,             aircon,               &
         dgn_a,             dgn_awet,         wetdens,              &
         n_mode,                                                    &
         qnum_cur,                                                  &
         qaer_cur,          qaer_delsub_coag_in,                    &
         qwtr_cur                                                   )

      qnum_delsub_coag = qnum_cur - qnum_sv1
      qaer_delsub_coag = qaer_cur - qaer_sv1

      end if


!
!
! primary carbon aging
!
!
      if ( n_agepair > 0 ) then

      call mam_pcarbon_aging_1subarea(                              &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         dtsubstep,         dgn_a,            do_cond,              &
         n_mode,                                                    &
         qnum_cur,          qnum_delsub_cond, qnum_delsub_coag,     &
         qaer_cur,          qaer_delsub_cond, qaer_delsub_coag,     &
         qaer_delsub_coag_in,                                       &
         qwtr_cur                                                   )

      end if


! accumulate sub-step q-dels
      if ( do_coag ) then
         qnum_del_coag = qnum_del_coag + qnum_delsub_coag
         qaer_del_coag = qaer_del_coag + qaer_delsub_coag
      end if
      if ( do_cond ) then
         qnum_del_cond = qnum_del_cond + qnum_delsub_cond
         qaer_del_cond = qaer_del_cond + qaer_delsub_cond
      end if

      end do jtsubstep_loop


!
!
! final mix ratios
!
!
      qgas4 = qgas_cur
      qaer4 = qaer_cur
      qnum4 = qnum_cur
      qwtr4 = qwtr_cur

! final mix ratio changes

      qgas_delaa(:,iqtend_cond) = qgas_del_cond(:)
      qgas_delaa(:,iqtend_rnam) = 0.0_r8
      qgas_delaa(:,iqtend_nnuc) = qgas_del_nnuc(:)
      qgas_delaa(:,iqtend_coag) = 0.0_r8

      qnum_delaa(:,iqtend_cond) = qnum_del_cond(:)
      qnum_delaa(:,iqtend_rnam) = qnum_del_rnam(:)
      qnum_delaa(:,iqtend_nnuc) = qnum_del_nnuc(:)
      qnum_delaa(:,iqtend_coag) = qnum_del_coag(:)

      qaer_delaa(:,:,iqtend_cond) = qaer_del_cond(:,:)
      qaer_delaa(:,:,iqtend_rnam) = qaer_del_rnam(:,:)
      qaer_delaa(:,:,iqtend_nnuc) = qaer_del_nnuc(:,:)
      qaer_delaa(:,:,iqtend_coag) = qaer_del_coag(:,:)

      misc_vars_aa_sub%ncluster_tend_nnuc_1grid = dnclusterdt

      return
      end subroutine mam_amicphys_1subarea_clear


!---------------------------------------------------------------------
!---------------------------------------------------------------------
#if ( defined( MOSAIC_SPECIES ) )
      subroutine mosaic_gasaerexch_1subarea_intr(  nstep,                &!Intent(ins)
              lchnk,             i_in,             k_in,        jsub_in, &
              temp,              relhum,           pmid,                 &
              aircon,            dtsubstep,        n_mode,               &
              dgn_a,             dgn_awet,         qaer_cur,             &!Intent(inouts)
              qgas_cur,          qnum_cur,         qwtr_cur,             &
              qgas_avg,          qgas_netprod_otrproc,                   &
              uptkrate_h2so4,    misc_vars_aa_sub                        )
        !------------------------------------------------------------------------------!
        !Purpose: This routine acts as an interface between Mosaic and CAM
        !Future work:
        !===========
        !1. Clean Mosaic code and get rid of the arguments which stays constant
        !   for the entire simulation
        !3. Please handle the Mosaic counters, either use pbuf or make them internal to 
        !   Mosaic
        !4. Use get_nstep() for it_mosaic or pull out the it_mosaic .eq. 1 computation 
        !   to the init routines
        !5. SOA from CAM is stored in LIM2 of Mosaic. Rest of the 7 SOA species in 
        !   Mosaic are populated with zeros
        !6. Some variables in Mosaic had to be initialized to zero. Please revisit and 
        !  fix whatever is necessary 
        !7. jhyst_leg is constant for now and is equal to jhyst_up
        !
        !Author: Balwinder Singh (PNNL)
        !------------------------------------------------------------------------------!
        !Use statements
        use module_mosaic_box_aerchem, only: mosaic_box_aerchemistry
        use infnan,                    only: nan, bigint
        use physconst,                 only: mwh2o
        use module_data_mosaic_aero,   only: naer_mosaic => naer, &
             inh4_a, ilim2_a, iso4_a, ina_a, icl_a, ibc_a, ioin_a, ioc_a, &
             ino3_a, icl_a,   ica_a,  ico3_a, &
             ilim2_g, ih2so4_g, inh3_g, ihno3_g, ihcl_g, &
             jhyst_up, jtotal, &
             nbin_a, nbin_a_max, ngas_volatile, nmax_astem, nmax_mesa, nsalt, &
             mosaic_vars_aa_type
#ifdef SPMD
        use spmd_dyn,                  only: mpicom_xy, iam
        use units,                     only: getunit, freeunit
#endif

        implicit none
        
        !Args: intent(in)
        integer,  intent(in) :: lchnk                 ! chunk identifier
        integer,  intent(in) :: nstep                 ! model time-step number
        integer,  intent(in) :: i_in, k_in            ! column and level indices
        integer,  intent(in) :: jsub_in               ! subarea index
      
        real(r8), intent(in) :: temp             !Temperature at model levels (K)
        real(r8), intent(in) :: relhum           !Relative humidity (0-1)
        real(r8), intent(in) :: pmid             !Pressure at layer center (Pa)
        real(r8), intent(in) :: aircon           !Air molar density (kmol/m3)
        real(r8), intent(in) :: dtsubstep        !Time sub-step (s)
        integer,  intent(in) :: n_mode           !current number of active modes

        !Args: intent(inout)
        real(r8), intent(inout) :: dgn_a(max_mode)            !Dry geo. mean dia. (m) of number distrib.
        real(r8), intent(inout) :: dgn_awet(max_mode)         !Wet geo. mean dia. (m) of number distrib.
        real(r8), intent(inout) :: qaer_cur(max_aer,max_mode) !Current aerosol mass mix ratios (mol/mol)
        real(r8), intent(inout) :: qgas_cur(max_gas)          !Current gas mix ratios (mol/mol)
        real(r8), intent(inout) :: qnum_cur(max_mode)         !Current aerosol number mix ratios (#/kmol)
        real(r8), intent(inout) :: qwtr_cur(max_mode)         !Current aerosol water mix ratios (mol/mol)
        real(r8), intent(inout) :: qgas_avg(max_gas)          !average gas conc. over dtchem time step (mol/mol)
        real(r8), intent(inout) :: qgas_netprod_otrproc(max_gas)
                  ! qgas_netprod_otrproc = gas net production rate from other processes
                  !    such as gas-phase chemistry and emissions (mol/mol/s)
                  ! this allows the condensation (gasaerexch) routine to apply production and condensation loss 
                  !    together, which is more accurate numerically
                  ! NOTE - must be >= zero, as numerical method can fail when it is negative
                  ! NOTE - currently only the values for h2so4 and nh3 should be non-zero
        real(r8), intent(inout) :: uptkrate_h2so4  ! rate of h2so4 uptake by aerosols (1/s)
        type ( misc_vars_aa_type ), intent(inout) :: misc_vars_aa_sub

        !Local Variables - [To be sent as args to Mosaic code]
        integer  :: ierr
!       integer  :: it_mosaic     !Time step counter for Mosaic
!       integer  :: jASTEM_fail   !Counter to indicate if the ASTEM convergence failed in Mosaic
        real(r8) :: dtchem        !Timestep in seconds 
        real(r8) :: T_K           !Temperature in K

        
        integer :: mcall_load_mosaic_parameters !Flag to decide whether to call 'load_mosaic_parameters' or not(*BALLI not used anymore)
        integer :: mcall_print_aer_in           !Flag to decide whether to call 'print_aer' or not

        integer, dimension(nbin_a_max) :: jaerosolstate     !Aerosol state (solid, liquid, gas)
!       integer, dimension(nbin_a_max) :: iter_mesa         !MESA iterations counters
        integer, dimension(nbin_a_max) :: jaerosolstate_bgn !Aerosol state at the begining (solid, liquid, gas)
        integer, dimension(nbin_a_max) :: jhyst_leg


        real(r8) :: aH2O              !Relative humidity in fraction(variaes between 0 and 1)
        real(r8) :: P_atm             !Pressure in atm units
        real(r8) :: RH_pc             !Relative humidity in %age(variaes between 0 and 100)
        real(r8) :: cair_mol_m3       !Air molar density (mol/m3)

        real(r8), dimension(nbin_a_max) :: water_a    !Current aerosol water mix ratios (kg/m3)
        real(r8), dimension(nbin_a_max) :: sigmag_a   !Geometric standard deviation for aerosol mode 
        real(r8), dimension(nbin_a_max) :: Dp_dry_a   !Dry geo. mean dia. (cm) of number distrib. 

        real(r8), dimension(nbin_a_max) :: num_a            !Current aerosol number mix ratios (#/cm3)
        real(r8), dimension(nbin_a_max) :: dp_wet_a         !Diameter of aerosol in (cm) 
        real(r8), dimension(nbin_a_max) :: mass_dry_a_bgn   !g/cc(air) **BALLI*** comment missing
        real(r8), dimension(nbin_a_max) :: mass_dry_a       !g/cc(air) **BALLI*** comment missing
        real(r8), dimension(nbin_a_max) :: dens_dry_a_bgn   !g/cc      **BALLI*** comment missing
        real(r8), dimension(nbin_a_max) :: dens_dry_a       !g/cc      **BALLI*** comment missing
        real(r8), dimension(nbin_a_max) :: water_a_hyst     !kg(water)/m^3(air) hysteresis (at 60% RH) **BALLI*** comment missing
        real(r8), dimension(nbin_a_max) :: aH2O_a           !Relative humidity in fraction(variaes between 0 and 1)
        real(r8), dimension(nbin_a_max) :: gam_ratio
        
        real(r8), dimension(ngas_volatile) :: gas            !Current gas mix ratios (nano mol/m3)
        real(r8), dimension(ngas_volatile) :: gas_avg          ! average gas conc. over dtchem time step (nmol/m3)
        real(r8), dimension(ngas_volatile) :: gas_netprod_otrproc
                  ! gas_netprod_otrproc = gas net production rate from other processes
                  !    such as gas-phase chemistry and emissions (mol/mol/s)
                  ! NOTE - must be >= zero, as numerical method can fail when it is negative
                  ! NOTE - currently for mosaic, only the value for h2so4 can be non-zero
        real(r8), dimension(naer_mosaic,3,nbin_a_max) :: aer !Current aerosol mass mix ratios (nano mol/m3)

        !Local Variables - [other local variables]
        character(len=500) :: tmp_str, nlfile, sballi
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
        character(len=500) :: infile
        logical, parameter :: debug_mosaic = .false.
        integer, parameter :: iam = 0
#endif
!       logical  :: zero_water_flag, flag_itr_kel
        integer  :: imode, ibin, iaer, igas, istate, isalt, ibin_in, iaer_in, istate_in
        integer  :: unitn

        !BALLI -  Following should be in the modules as parameter
        real(r8), parameter ::  oneatminv = 1.0_r8/1.01325e5_r8  
        !BALLI -  Following should be in the modules as parameter -  ENDS

        !BSINGH - For converting CAM units to Mosaic units
        real(r8) :: nano_mult_cair, aer_tmp
        real(r8) :: num_cam_to_mos_units, wtr_cam_to_mos_units

        !BSINGH - For converting Mosaic units to CAM units
        real(r8) :: nano_mult_cair_inv
        real(r8) :: num_mos_to_cam_units, wtr_mos_to_cam_units

        !BSINGH -  For debugging mosiac code
!       integer, dimension(6)          :: hostgridinfo
        integer, dimension(nbin_a_max) :: jaerosolstate_in
        integer, dimension(nbin_a_max) :: jhyst_leg_in

        real(r8), dimension(nbin_a_max) :: num_a_in 
        real(r8), dimension(nbin_a_max) :: dp_wet_a_in
        real(r8), dimension(nbin_a_max) :: water_a_in   
        real(r8), dimension(nbin_a_max) :: sigmag_a_in
        real(r8), dimension(nbin_a_max) :: Dp_dry_a_in  

        real(r8), dimension(ngas_volatile) :: gas_in, gas_netprod_otrproc_in, gas_avg_in

        real(r8), dimension(naer_mosaic,3,nbin_a_max) :: aer_in
        
        real(r8), dimension(naer_mosaic) :: kappa_nonelectro

        type (mosaic_vars_aa_type) :: mosaic_vars_aa

        !BSINGH - For storing points having trouble converging
        logical,parameter ::  convergence_pt_trk = .true. !For tracking points where convergence failed, let the run proceed
!       logical :: f_neg_vol_tmp


        ! allocate the allocatable parts of mosaic_vars_aa
        allocate( mosaic_vars_aa%iter_mesa(nbin_a_max), stat=ierr )
        if (ierr /= 0) then
           print *, '*** subr mosaic_gasaerexch_1subarea_intr - allocate error for mosaic_vars_aa%iter_mesa'
           stop
        end if


        !if(nstep>17 .and. i_in == 14)write(202,*)'AMICPHYS I K:',i_in,k_in,nstep,lchnk
        !------------------------------------------------------------!
        !------------------------------------------------------------!
        !Populate MOSAIC variables
        !------------------------------------------------------------!
        !------------------------------------------------------------!
        
        !Counters:
        !BSINGH - This counter is internal to Mosaic model.
        !         It indicates if ASTEM convergence failed in Mosaic
!       jASTEM_fail = 0
        mosaic_vars_aa%jastem_fail = 0
        
        !BSINGH - This is time step number in Mosaic
!       it_mosaic   = nstep
        mosaic_vars_aa%it_mosaic = nstep
        
        !Inputs for Mosaic model (Should be intent-ins for Mosaic model)
        aH2O         = relhum               !Relative humidity [fraction between 0 and 1]
        T_K          = temp                 !Temperature in K
        P_atm        = pmid   * oneatminv   !Pressure (atm)
        RH_pc        = aH2O   * 100.0_r8    !Relative humidity [%age between 0 and 100]
        cair_mol_m3  = aircon * 1000.0_r8   !Air molar density (mol/m3){units conversion: aircon[kmol/m3] * 1.0e3[mol/kmol]}
        dtchem       = dtsubstep            !timestep (s)
        
        jhyst_leg(1:nbin_a_max)      = jhyst_up

        !Flags to control Mosaic model
        mcall_load_mosaic_parameters = 1    !**BALLI.. This flag is not used anymore
        mcall_print_aer_in           = 0    !**BALLI...insert a dummy call to print_aer
        
        
        !Populate aersols
        nbin_a = n_mode   ! current number of modes

        aer(:,:,:) = 0.0_r8 !initialized to zero

        !Populate aerosol numbers and water species
        num_a(:)   = 0.0_r8 !Initialized to zero
        water_a(:) = 0.0_r8 !initialized to zero
        
        !BSINGH - units of qnum_cur in CAM are #/kmol of air. In Mosaic, units are #/cm3
        !Units conversion: qnum_cur[#/kmol] * 1.0e-3[kmol/mol] * cair_mol_m3[mol/m3] * 1.0e-6[m3/cm3]

        num_cam_to_mos_units = 1.0e-3_r8 * cair_mol_m3 * 1.0e-6_r8
        
        !BSINGH - units for water in CAM are mol/mol. In Mosaic, units are kg/m3 
        !Units conversion: qwtr_cur[mol/mol] * mwh2o[g/mol] * cair_mol_m3[mol/m3] * 1.0e-3[kg/g]
        
        wtr_cam_to_mos_units = mwh2o * cair_mol_m3 * 1.0e-3_r8
        
        
        nano_mult_cair  =  cair_mol_m3 * 1.0e9_r8
        
        do imode = 1, n_mode
           !Notes:
           !1. NCL(sea salt) of CAM is mapped in NA and CL of MOSAIC
           !2. SOA of CAM is lumped into LIM2 species of MOSAIC !BALLI *ASK RAHUL and Dick
           !3. Species NO3, MSA, CO3, Ca do not exist in CAM therefore not mapped here
           !4. Species ARO1, ARO2, ALK1, OLE1, API1, API2, LIM1 are SOA species in MOSAIC
           !   which are not used in CAM-MOSAIC framework as of now
           !5. CAM units are (mol/mol of air) which are converted to Mosaic units (nano mol/m3).
           
           !Units conversion:qaer_cur[mol/mol] * cair_mol_m3[mol/m3] * 1.0e9[nmol/mol] 
           aer(inh4_a,  jtotal, imode)  = qaer_cur(iaer_nh4, imode) * nano_mult_cair 
           aer(ilim2_a, jtotal, imode)  = qaer_cur(iaer_soa, imode) * nano_mult_cair
           aer(iso4_a,  jtotal, imode)  = qaer_cur(iaer_so4, imode) * nano_mult_cair 
           aer(ina_a,   jtotal, imode)  = qaer_cur(iaer_ncl, imode) * nano_mult_cair 
           if (iaer_cl  > 0) then
              aer(icl_a,   jtotal, imode)  = qaer_cur(iaer_cl,  imode) * nano_mult_cair 
           else
              aer(icl_a,   jtotal, imode)  = qaer_cur(iaer_ncl, imode) * nano_mult_cair 
           end if
           if (iaer_no3 > 0) &
                aer(ino3_a,  jtotal, imode)  = qaer_cur(iaer_no3, imode) * nano_mult_cair 
           if (iaer_ca  > 0) &
                aer(ica_a,   jtotal, imode)  = qaer_cur(iaer_ca,  imode) * nano_mult_cair 
           if (iaer_co3 > 0) &
                aer(ico3_a,  jtotal, imode)  = qaer_cur(iaer_co3, imode) * nano_mult_cair 
           
           !Units of BC, OC and DST in CAM are (mol/mol of air) and nano-g/m3 in MOSAIC
           !Units conversion:qaer_cur[mol/mol] * mw_aer[g/mol] * cair_mol_m3[mol/m3] * 1.0e9[nano-g/g] 
           aer(ibc_a,  jtotal, imode)  = qaer_cur(iaer_bc,  imode) * mw_aer(iaer_bc)  * nano_mult_cair 
           aer(ioin_a, jtotal, imode)  = qaer_cur(iaer_dst, imode) * mw_aer(iaer_dst) * nano_mult_cair !BSINGH - "Other inorganic(oin)" in Mosaic is DST in CAM
           aer(ioc_a,  jtotal, imode)  = qaer_cur(iaer_pom, imode) * mw_aer(iaer_pom) * nano_mult_cair
           
           !Populate aerosol number and water species
           num_a(imode)   = qnum_cur(imode) * num_cam_to_mos_units
           water_a(imode) = qwtr_cur(imode) * wtr_cam_to_mos_units
        end do
        
        !Populate gases
        gas(:) = 0.0_r8  !Initialized to zero
        !BSINGH - only 3 gases are avialble in CAM (SOAG, H2SO4, NH3). 
        !SOAG is stored in LIM2 gas species as of now
        !CAM units are (mol/mol of air) which are converted to Mosaic units (nano mol/m3).
        gas_avg(:) = 0.0_r8
        
        !Units conversion:qgas_cur[mol/mol] * cair_mol_m3[mol/m3] * 10.0e9[nmol/mol] 
        gas(ilim2_g)  = qgas_cur(igas_soa)   * nano_mult_cair
        gas(ih2so4_g) = qgas_cur(igas_h2so4) * nano_mult_cair 
        gas(inh3_g)   = qgas_cur(igas_nh3)   * nano_mult_cair 
        if (igas_hno3 > 0) &
             gas(ihno3_g)   = qgas_cur(igas_hno3)   * nano_mult_cair 
        if (igas_hcl > 0) &
             gas(ihcl_g)   = qgas_cur(igas_hcl)   * nano_mult_cair 
        
        !Populate gas_netprod_otrproc
        gas_netprod_otrproc(:) = 0.0_r8
        gas_netprod_otrproc(ih2so4_g) = qgas_netprod_otrproc(igas_h2so4) * nano_mult_cair
        ! nh3 gas-phase chem production should be zero (unless we include emissions)
        !    and doing simultaneous production and condensation in mosaic is more complicated
        !    for nh3 that for h2so4
        ! so for now, just add in the production here
        gas(inh3_g)   = gas(inh3_g) + max( qgas_netprod_otrproc(igas_nh3)*dtchem, 0.0_r8 ) * nano_mult_cair
        
        
        !BSINGH - Initialize the following variables as 'nan' and then assign values to a subset of their dimesions
        Dp_dry_a(:)            = nan
        sigmag_a(:)            = nan
        dp_wet_a(:)            = nan
        
        sigmag_a(1:n_mode) = sigmag_aer(1:n_mode)       !Geometric standard deviation for aerosol mode 
        !Dry geo. mean dia.(cm) of number distrib [convert from m to cm]!**BALLI: check if it meant to be in,inout or out only and units also
        Dp_dry_a(1:n_mode) = dgn_a(1:n_mode)    * 100.0_r8 * fcvt_dgnum_dvolmean(1:n_mode)
        !Wet geo. mean dia.(cm) of number distrib [convert from m to cm]!**BALLI: check if it meant to be in,inout or out only and untis also
        dp_wet_a(1:n_mode) = dgn_awet(1:n_mode) * 100.0_r8 * fcvt_dgnum_dvolmean(1:n_mode)
        
        
        !BSINGH - These are output variables from Mosaic. 
        !Declared as nan to make sure that they are not inadvertently 'used' before assignment.
!       iter_mesa(:)         = bigint
        mosaic_vars_aa%iter_mesa(1:nbin_a_max) = 0
        jaerosolstate_bgn(:) = bigint
        jaerosolstate(:)     = bigint

        mass_dry_a_bgn(:)    = nan
        mass_dry_a(:)        = nan
        dens_dry_a_bgn(:)    = nan
        dens_dry_a(:)        = nan
        water_a_hyst(:)      = nan
        aH2O_a(:)            = nan
        gam_ratio(:)         = nan
        
        
        !------------------------------------------------------------!
        !------------------------------------------------------------!
        !END [Populate MOSAIC variables]
        !------------------------------------------------------------!
        !------------------------------------------------------------!


#if ( defined( CAMBOX_ACTIVATE_THIS ) )
        !BSINGH - This section is required ONLY for the MAM box model
        !         to see if it can reproduce errors encountered by the 
        !         CAM model during runtime(e.g. convergence errors).
        !         This block repopulate all the information which is 
        !         going into the mosaic box (intent-ins and intent-inouts). 
        !         It is a binary read to preserve the accuracy.
  
        if(debug_mosaic) then
           !Read a binary file which has all the inputs to the mosaic box
           !and stop the model
           
           unitn = 101
           infile = 'mosaic_error_7.bin'
           open( unitn, file=trim(infile), status='old', form='unformatted', CONVERT = 'BIG_ENDIAN' )
           
           read(unitn)aH2O
           read(unitn)T_K
           read(unitn)P_atm
           read(unitn)RH_pc
           read(unitn)dtchem
           
           do ibin = 1, ntot_amode !nbin_a_max
              read(unitn)num_a(ibin),water_a(ibin),Dp_dry_a(ibin),        &
                   sigmag_a(ibin),dp_wet_a(ibin),jhyst_leg(ibin),          &
                   jaerosolstate(ibin)
           end do
           
           
           do igas = 1, ngas_volatile
              read(unitn) gas(igas), gas_avg(igas), gas_netprod_otrproc(igas)
           enddo
           
           do ibin = 1, ntot_amode !nbin_a_max
              do istate = 1, 3
                 do iaer = 1 , naer
                    read(unitn)iaer_in,istate_in,ibin_in, aer_tmp
                    aer(iaer_in,istate_in,ibin_in) = aer_tmp                    
                 end do
              end do
           end do
           close(unitn)
        endif
        !BSINGH -----xxx ENDS reading file for debugging mosaic xxxx----
#endif
        
        
        !Store the variables which are intent(inout) to Mosaic box model
        !for debuging purposes
        aer_in(:,:,:)             = aer(:,:,:)
        num_a_in(:)               = num_a(:)
        water_a_in(:)             = water_a(:)
        gas_in(:)                 = gas(:)
        Dp_dry_a_in(:)            = Dp_dry_a(:)
        sigmag_a_in(:)            = sigmag_a(:)
        dp_wet_a_in(:)            = dp_wet_a(:)
        jhyst_leg_in(:)           = jhyst_leg(:)
        jaerosolstate_in(:)       = jaerosolstate(:)
        gas_netprod_otrproc_in(:) = gas_netprod_otrproc(:)
        gas_avg_in(:)             = gas_avg(:)
        

        !BSINGH - zero_water_flag becomes .true. if water is zero in liquid phase
!       zero_water_flag = .false. 
        mosaic_vars_aa%zero_water_flag = .false.
        !BSINGH - flag_itr_kel becomes true when kelvin iteration in mdofule_mosaic_ext.F90 are greater then 100
!       flag_itr_kel    = .false.        
        mosaic_vars_aa%flag_itr_kel = .false.


        !Store grid info
!       hostgridinfo(1)   = i_in
!       hostgridinfo(2)   = k_in
!       hostgridinfo(3)   = lchnk
!       hostgridinfo(4:6) = bigint
        mosaic_vars_aa%hostgridinfo(1)   = i_in
        mosaic_vars_aa%hostgridinfo(2)   = k_in
        mosaic_vars_aa%hostgridinfo(3)   = lchnk
        mosaic_vars_aa%hostgridinfo(4:6) = bigint
        mosaic_vars_aa%it_host = 0


        ! *** maybe these should be bigint or nan ???
        mosaic_vars_aa%f_mos_fail = -1
        mosaic_vars_aa%isteps_astem = 0
        mosaic_vars_aa%isteps_astem_max = 0
        mosaic_vars_aa%jastem_call = 0
        mosaic_vars_aa%jmesa_call = 0
        mosaic_vars_aa%jmesa_fail = 0
        mosaic_vars_aa%niter_mesa_max = 0
        mosaic_vars_aa%nmax_astem = nmax_astem
        mosaic_vars_aa%nmax_mesa = nmax_mesa
        mosaic_vars_aa%cumul_steps_astem = 0.0_r8
        mosaic_vars_aa%niter_mesa = 0.0_r8
        mosaic_vars_aa%xnerr_astem_negative(:,:) = 0.0_r8


        ! set kappa values for non-electrolyte species
        ! reason for doing this here is that if cam eventually has multiple varieties of dust and/or pom, 
        !    then the dust hygroscopicity may vary spatially and temporally,
        !    and the kappa values cannot be constants
        kappa_nonelectro(:) = 0.0_r8
        kappa_nonelectro(ibc_a  ) = 0.0001  ! previously kappa_poa = 0.0001
        kappa_nonelectro(ioc_a  ) = 0.0001  ! previously kappa_bc  = 0.0001
        kappa_nonelectro(ilim2_a) = 0.1     ! previously kappa_soa = 0.1
        kappa_nonelectro(ioin_a ) = 0.06    ! previously kappa_oin = 0.06


        !Call MOSAIC parameterization
        !BSINGH - jASTEM_fail is in arg list to know if the mosiac model converged or not
        !BSINGH - Following variables are not required by CAM but they still exist in the
        !         calling arguments as intent-outs as Mosaic model needs them to be in the 
        !         arg list:
        !         gam_ratio, iter_mesa, aH2O_a,jaerosolstate, mass_dry_a_bgn, mass_dry_a, 
        !         dens_dry_a_bgn, dens_dry_a, water_a_hyst, jaerosolstate_bgn

! *** ff03h version ***
!       call mosaic_box_aerchemistry(                                                   &
!            hostgridinfo,            it_mosaic,    aH2O,               T_K,            &!Intent-ins
!            P_atm,                   RH_pc,        dtchem,                             &
!            mcall_load_mosaic_parameters,          mcall_print_aer_in, sigmag_a,       &
!            jaerosolstate,           aer,                                              &!Intent-inouts
!            num_a,                   water_a,      gas,                                &
!            gas_avg,                 gas_netprod_otrproc,              Dp_dry_a,       &
!            dp_wet_a,                jhyst_leg,    zero_water_flag,    flag_itr_kel,   &
!            mass_dry_a_bgn,          mass_dry_a,                                       &!Intent-outs
!            dens_dry_a_bgn,          dens_dry_a,   water_a_hyst,       aH2O_a,         &
!            gam_ratio,               jaerosolstate_bgn,                jASTEM_fail,    &
!            iter_MESA,               f_neg_vol_tmp                                     )

! *** ff04a version ***
        call mosaic_box_aerchemistry(               aH2O,               T_K,            &!Intent-ins
             P_atm,                   RH_pc,        dtchem,                             &
             mcall_load_mosaic_parameters,          mcall_print_aer_in, sigmag_a,       &
             kappa_nonelectro,                                                          &
             jaerosolstate,           aer,                                              &!Intent-inouts
             num_a,                   water_a,      gas,                                &
             gas_avg,                 gas_netprod_otrproc,              Dp_dry_a,       &
             dp_wet_a,                jhyst_leg,                                        &
             mosaic_vars_aa,                                                            &
             mass_dry_a_bgn,          mass_dry_a,                                       &!Intent-outs
             dens_dry_a_bgn,          dens_dry_a,   water_a_hyst,       aH2O_a,         &
             uptkrate_h2so4,          gam_ratio,    jaerosolstate_bgn                   )

! *** ff04a version ***
!  subr       mosaic_box_aerchemistry(        aH2O,               T_K,            &!Intent-ins
!       P_atm,                  RH_pc,        dtchem,                             &
!       mcall_load_mosaic_parameters,         mcall_print_aer_in, sigmag_a,       &
!       kappa_nonelectro,                                                         &
!       jaerosolstate,          aer,                                              &!Intent-inouts
!       num_a,                  water_a,      gas,                                &
!       gas_avg,                gas_netprod_otrproc,              Dp_dry_a,       &
!       dp_wet_a,               jhyst_leg,                                        &
!       mosaic_vars_aa,                                                           &
!       mass_dry_a_bgn,         mass_dry_a,                                       &!Intent-outs
!       dens_dry_a_bgn,         dens_dry_a,   water_a_hyst,       aH2O_a,         &
!       uptkrate_h2so4,         gam_ratio,    jaerosolstate_bgn                   )

        if (mosaic_vars_aa%flag_itr_kel) then
           misc_vars_aa_sub%max_kelvin_iter_1grid = misc_vars_aa_sub%max_kelvin_iter_1grid + 1.0_r8
        endif
        
        if (mosaic_vars_aa%jASTEM_fail > 0 .or. mosaic_vars_aa%zero_water_flag .or. mosaic_vars_aa%f_mos_fail > 0 ) then !solver in ASTEM didn't converge
           
           !Let the run proceed and track the points(i,k) where the run fails convergence 
           if(convergence_pt_trk .and. mosaic_vars_aa%jASTEM_fail > 0 ) then
              misc_vars_aa_sub%cnvrg_fail_1grid = misc_vars_aa_sub%cnvrg_fail_1grid + 1.0_r8
           else
              !Printout a binary file which has all the inputs to the mosaic box
              !and stop the model
              
              !Generate a unit number and form file name based on process number
#ifdef SPMD
              unitn  = getunit()
              write(tmp_str,*)iam
              write(nlfile,*)'mosaic_error_',trim(adjustl(tmp_str)),'.bin'               
#else
              unitn = 101
              nlfile = 'mosiac_error.txt'
#endif
              !Open a binary file, remember it is written out as BIG ENDIAN
              open( unitn, file=trim(nlfile), status='unknown', form = 'unformatted' )
              
              write(unitn)aH2O    !Write relative humidity
              write(unitn)T_K     !Write relative temp 
              write(unitn)P_atm   !Write relative pressure
              write(unitn)RH_pc
              write(unitn)dtchem
              !Write variables with 'nbin_a_max' dimension
              do ibin = 1, ntot_amode!nbin_a_max
                 write(unitn)num_a_in(ibin),water_a_in(ibin),Dp_dry_a_in(ibin),        &
                      sigmag_a_in(ibin),dp_wet_a_in(ibin),jhyst_leg_in(ibin),          &
                      jaerosolstate_in(ibin)
              end do
              
              !Write gas array
              do igas = 1, ngas_volatile
                 write(unitn) gas_in(igas), gas_avg_in(igas), gas_netprod_otrproc_in(igas)
              enddo
              
              !Write aerosols
              do ibin = 1, ntot_amode !nbin_a_max
                 do istate = 1, 3
                    do iaer = 1 , naer_mosaic
                       write(unitn)iaer,istate,ibin,aer_in(iaer,istate,ibin)
                    end do
                 end do
              end do
              !Close the file 
              close(unitn)
#ifdef SPMD
              !free unit number
              call freeunit(unitn)
#endif
              !Write error message and stop the model.
              write(tmp_str,*) 'Error in Mosaic, jASTEM_fail= ', mosaic_vars_aa%jASTEM_fail, &
                 ' zero_water_flag: ', mosaic_vars_aa%zero_water_flag, &
                 ' f_mos_fail:', mosaic_vars_aa%f_mos_fail
              call endrun (tmp_str)
           endif
        endif
        

        ! copy other diagnostic outputs (that are written to history) from mosaic_vars_aa to misc_vars_aa_sub
        misc_vars_aa_sub%xnerr_astem_negative_1grid(:,:) = mosaic_vars_aa%xnerr_astem_negative(:,:)

        ! deallocate the allocatable parts of mosaic_vars_aa
        deallocate( mosaic_vars_aa%iter_mesa, stat=ierr )
        if (ierr /= 0) then
           print *, '*** subr mosaic_gasaerexch_1subarea_intr - deallocate error for mosaic_vars_aa%iter_mesa'
           stop
        end if


        !------------------------------------------------------------!
        !------------------------------------------------------------!
        !Process MOSAIC output and store it in CAM data structures
        !------------------------------------------------------------!
        !------------------------------------------------------------! 
        !BSINGH - units of qnum_cur in CAM are #/kmol of air. In Mosaic, units are #/cm3
        num_mos_to_cam_units = 1.0_r8/num_cam_to_mos_units !Take inverse of cam_to_mos units
        num_cam_to_mos_units = nan                         !To avoid inadvertent use

        !BSINGH - units for water in CAM are mol/mol. In Mosaic, units are kg/m3 
        wtr_mos_to_cam_units = 1.0_r8/wtr_cam_to_mos_units !Take inverse of cam_to_mos units
        wtr_cam_to_mos_units = nan                         !To avoid inadvertent use

        nano_mult_cair_inv  =  1.0_r8/nano_mult_cair       !Take inverse of cam to mosaic units
        nano_mult_cair      = nan                          !To avoid inadvertent use
        
        do imode = 1, n_mode
           !Notes:
           !1. NCL(sea salt) of CAM is mapped in NA and CL of MOSAIC
           !2. SOA of CAM is lumped into LIM2 species of MOSAIC !BALLI *ASK RAHUL and Dick
           !3. Species NO3, MSA, CO3, Ca do not exist in CAM therefore not mapped here
           !4. Species ARO1, ARO2, ALK1, OLE1, API1, API2, LIM1 are SOA species in MOSAIC
           !   which are not used in CAM-MOSAIC framework as of now
           !5. CAM units are (mol/mol of air) and  Mosaic units are (nano mol/m3).

           qaer_cur(iaer_nh4, imode) = aer(inh4_a,  jtotal , imode) * nano_mult_cair_inv
           qaer_cur(iaer_soa, imode) = aer(ilim2_a, jtotal , imode) * nano_mult_cair_inv
           qaer_cur(iaer_so4, imode) = aer(iso4_a,  jtotal , imode) * nano_mult_cair_inv
           qaer_cur(iaer_ncl, imode) = aer(ina_a,   jtotal , imode) * nano_mult_cair_inv
           if (iaer_cl  > 0) &
           qaer_cur(iaer_cl,  imode) = aer(icl_a,   jtotal , imode) * nano_mult_cair_inv
           if (iaer_no3 > 0) &
           qaer_cur(iaer_no3, imode) = aer(ino3_a,  jtotal , imode) * nano_mult_cair_inv
           if (iaer_ca  > 0) &
           qaer_cur(iaer_ca,  imode) = aer(ica_a,   jtotal , imode) * nano_mult_cair_inv
           if (iaer_co3 > 0) &
           qaer_cur(iaer_co3, imode) = aer(ico3_a,  jtotal , imode) * nano_mult_cair_inv

           !Units of BC, OC and DST in CAM are (mol/mol of air) and nano-g/m3 in MOSAIC
           qaer_cur(iaer_bc,  imode) = (aer(ibc_a,  jtotal , imode)/mw_aer(iaer_bc))  * nano_mult_cair_inv
           qaer_cur(iaer_dst, imode) = (aer(ioin_a, jtotal , imode)/mw_aer(iaer_dst)) * nano_mult_cair_inv !BSINGH - "Other inorganic" in Mosaic is DST in CAM
           qaer_cur(iaer_pom, imode) = (aer(ioc_a,  jtotal , imode)/mw_aer(iaer_pom)) * nano_mult_cair_inv

           !Populate aerosol number and water species
           qnum_cur(imode) = num_a(imode)   * num_mos_to_cam_units
           qwtr_cur(imode) = water_a(imode) * wtr_mos_to_cam_units
        end do

        !BSINGH - only 3 gases are avialble in CAM (SOAG, H2SO4, NH3). 
        !SOAG is stored in LIM2 gas species as of now

        qgas_cur(igas_soa)   = gas(ilim2_g)  * nano_mult_cair_inv
        qgas_cur(igas_h2so4) = gas(ih2so4_g) * nano_mult_cair_inv 
        qgas_cur(igas_nh3)   = gas(inh3_g)   * nano_mult_cair_inv 

        qgas_avg(igas_soa)   = gas_avg(ilim2_g)  * nano_mult_cair_inv
        qgas_avg(igas_h2so4) = gas_avg(ih2so4_g) * nano_mult_cair_inv
        qgas_avg(igas_nh3)   = gas_avg(inh3_g)   * nano_mult_cair_inv

        if (igas_hno3 > 0) then
           qgas_cur(igas_hno3) = gas(    ihno3_g) * nano_mult_cair_inv
           qgas_avg(igas_hno3) = gas_avg(ihno3_g) * nano_mult_cair_inv
        end if
        if (igas_hcl > 0) then
           qgas_cur(igas_hcl) = gas(    ihcl_g) * nano_mult_cair_inv
           qgas_avg(igas_hcl) = gas_avg(ihcl_g) * nano_mult_cair_inv
        end if

        !update mode  diameters 
        !Dry geo. mean dia.(m) of number distrib [convert from cm to m]
        dgn_a(1:n_mode)          = Dp_dry_a(1:n_mode) * 0.01_r8 / fcvt_dgnum_dvolmean(1:n_mode)
        !Wet geo. mean dia.(m) of number distrib [convert from cm to m]
        dgn_awet(1:n_mode)       = dp_wet_a(1:n_mode) * 0.01_r8 / fcvt_dgnum_dvolmean(1:n_mode)

        !------------------------------------------------------------!
        !------------------------------------------------------------!
        !END [Process MOSAIC output ....]
        !------------------------------------------------------------!
        !------------------------------------------------------------!


      end subroutine mosaic_gasaerexch_1subarea_intr
#endif


!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine mam_gasaerexch_1subarea(                           &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         jtsubstep,         ntsubstep,                              &
         latndx,            lonndx,           lund,                 &
         dtsubstep,                                                 &
         temp,              pmid,             aircon,               &
         n_mode,                                                    &
         qgas_cur,          qgas_avg,                               &
         qgas_netprod_otrproc,                                      &
         qaer_cur,                                                  &
         qnum_cur,                                                  &
         qwtr_cur,                                                  &
         dgn_a,             dgn_awet,         wetdens,              &
         uptkaer,           uptkrate_h2so4                          )

! uses

      implicit none

! arguments
      integer,  intent(in) :: nstep                 ! model time-step number
      integer,  intent(in) :: lchnk                 ! chunk identifier
      integer,  intent(in) :: i, k                  ! column and level indices
      integer,  intent(in) :: jsub                  ! sub-area index
      integer,  intent(in) :: jtsubstep, ntsubstep  ! time substep info from calling routine
      integer,  intent(in) :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in) :: lund                  ! logical unit for diagnostic output
      integer,  intent(in) :: n_mode                ! current number of modes (including temporary)

      real(r8), intent(in) :: dtsubstep        ! integration timestep (s)
      real(r8), intent(in) :: temp             ! air temperature (K)
      real(r8), intent(in) :: pmid             ! air pressure at model levels (Pa)
      real(r8), intent(in) :: aircon           ! air molar concentration (kmol/m3)

      real(r8), intent(inout), dimension( 1:max_gas ) :: &
         qgas_cur, &  ! current gas mix ratios (mol/mol)
         qgas_avg     ! average gas mix ratios over the dtsubstep integration
      real(r8), intent(in   ), dimension( 1:max_gas ) :: &
         qgas_netprod_otrproc
                ! qgas_netprod_otrproc = gas net production rate from other processes
                !    such as gas-phase chemistry and emissions (mol/mol/s)
                ! this allows the condensation (gasaerexch) routine to apply production and condensation loss 
                !    together, which is more accurate numerically
                ! NOTE - must be >= zero, as numerical method can fail when it is negative
                ! NOTE - currently only the values for h2so4 and nh3 should be non-zero
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur     ! current aerosol mass mix ratios (mol/mol)
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur     ! current aerosol number mix ratios (#/kmol)
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur     ! current aerosol water mix ratios (mol/mol)
! qgas/aer/num/wtr_cur values are updated during the dtsubstep integration

      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         dgn_a,    &  ! dry geo. mean dia. (m) of number distrib.
         dgn_awet, &  ! wet geo. mean dia. (m) of number distrib.
         wetdens      ! interstitial aerosol wet density (kg/m3)
      real(r8), intent(inout), dimension( 1:max_gas, 1:max_mode ) :: &
         uptkaer      ! gas to aerosol mass transfer rate (1/s)
      real(r8), intent(inout) :: uptkrate_h2so4
                      ! h2so4(g) to aerosol mass transfer rate, summed over all modes (1/s)
                      ! this is needed by the nucleation routine (mam_newnuc_1subarea)

! local
      integer :: iaer, igas, ip
      integer :: ll
      integer :: n

      logical, parameter :: flag_nh4_lt_2so4_each_step  = .false.

      real(r8), dimension( 1:max_gas ) :: &
         gas_diffus,   &  ! gas diffusivity at current temp and pres (m2/s)
         gas_freepath     ! gas mean free path at current temp and pres (m)

      real(r8), dimension( 1:max_gas ) :: &
         qgas_prv

      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_prv

      real(r8) :: tmpa, tmpb, tmpc
      real(r8) :: tmp_kxt, tmp_kxt2, tmp_pxt, tmp_pok
      real(r8) :: tmp_q1, tmp_q2, tmp_q3, tmp_q4, tmp_q5
      real(r8) :: tmp_qdel_cond
      real(r8) :: uptkrate(max_mode)


      qgas_avg(1:ngas) = 0.0_r8


! calc gas uptake (mass transfer) rates
      if (jtsubstep == 1) then

      tmpa = pmid/1.013e5_r8
      do igas = 1, ngas
         gas_diffus(igas) = gas_diffusivity( &
                               temp, tmpa, mw_gas(igas), vol_molar_gas(igas) )

         tmpb = mean_molecular_speed( temp, mw_gas(igas) )

         gas_freepath(igas) = 3.0_r8 * gas_diffus(igas) / tmpb

!        subr gas_aer_uptkrates_1box1gas( &
!           accom, gasdiffus, gasfreepath, &
!           beta, nmode, dgncur_awet, lnsg, uptkrate )
         call gas_aer_uptkrates_1box1gas( &
            accom_coef_gas(igas), gas_diffus(igas), gas_freepath(igas), &
            0.0_r8, ntot_amode, dgn_awet, alnsg_aer, uptkrate )

         iaer = igas
         do n = 1, ntot_amode
            if ( lmap_aer(iaer,n) > 0 .or. & 
                 mode_aging_optaa(n) > 0 ) then
               ! uptkrate is for number = 1 #/m3, so mult. by number conc. (#/m3)
               uptkaer(igas,n) = uptkrate(n) * (qnum_cur(n) * aircon)
            else
               ! mode does not contain this species
               uptkaer(igas,n) = 0.0_r8
            end if
         end do
      end do ! igas

      do igas = 1, ngas
         ! use cam5.1.00 uptake rates
         if (igas <= nsoa    ) uptkaer(igas,1:ntot_amode) = uptkaer(igas_h2so4,1:ntot_amode)*0.81
         if (igas == igas_nh3) uptkaer(igas,1:ntot_amode) = uptkaer(igas_h2so4,1:ntot_amode)*2.08
      end do ! igas
      uptkrate_h2so4 = sum( uptkaer(igas_h2so4,1:ntot_amode) )

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( k == pver .and. ldiagd1 ) write(lund,'(a,2i4,1p,10e11.3)') 'i,k,h2so4_uprt', i, k, uptkaer(igas_h2so4,1:ntot_amode)
!     if (i==1 .and. k==4) then
!        write(*,*) 'uptake rates at i=1, k=4, igas down, nmode across'
!        do igas = 1, ngas
!           write(*,'(1p,10e10.2)') uptkaer(igas,1:ntot_amode)
!        end do
!        write(*,*) 'dgn_awet then sigmag then qnum'
!        write(*,'(1p,10e10.2)') dgn_awet(1:ntot_amode)
!        write(*,'(1p,10e10.2)') sigmag_aer(1:ntot_amode)
!        write(*,'(1p,10e10.2)') qnum_cur(1:ntot_amode)
!     end if
#endif

      end if ! (jtsubstep == 1)


! do soa
      call mam_soaexch_1subarea(                                    &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         dtsubstep,                                                 &
         temp,              pmid,             aircon,               &
         n_mode,                                                    &
         qgas_cur,          qgas_avg,                               &
         qaer_cur,                                                  &
         qnum_cur,                                                  &
         qwtr_cur,                                                  &
         uptkaer                                                    )


! do other gases (that are assumed non-volatile) with no time sub-stepping
      do igas = nsoa+1, ngas
         iaer = igas
         qgas_prv(igas)          = qgas_cur(igas)
         qaer_prv(iaer,1:n_mode) = qaer_cur(iaer,1:n_mode)
      end do

      do igas = nsoa+1, ngas
         iaer = igas
         if ( (igas == igas_hno3) .or. &
              (igas == igas_hcl ) ) cycle

         tmpa = sum( uptkaer(igas,1:n_mode) )
         tmp_kxt = tmpa*dtsubstep
         tmp_pxt = qgas_netprod_otrproc(igas)*dtsubstep
         tmp_q1 = qgas_prv(igas)
         ! tmp_q1 = mix-rat at t=tcur
         ! tmp_q3 = mix-rat at t=tcur+dtsubstep
         ! tmp_q4 = avg mix-rat between t=tcur and t=tcur+dtsubstep
         if (tmp_kxt >= 1.0e-20_r8) then
            if (tmp_kxt > 0.001_r8) then
               tmp_pok = tmp_pxt/tmp_kxt
               tmp_q3 = (tmp_q1 - tmp_pok)*exp(-tmp_kxt) + tmp_pok
               tmp_q4 = (tmp_q1 - tmp_pok)*(1.0_r8 - exp(-tmp_kxt))/tmp_kxt + tmp_pok
            else
               tmp_kxt2 = tmp_kxt*tmp_kxt
               tmp_q3 = tmp_q1 *(1.0_r8 - tmp_kxt        + tmp_kxt2*0.5_r8) &
                      + tmp_pxt*(1.0_r8 - tmp_kxt*0.5_r8 + tmp_kxt2/6.0_r8)
               tmp_q4 = tmp_q1 *(1.0_r8 - tmp_kxt*0.5_r8 + tmp_kxt2/6.0_r8) &
                      + tmp_pxt*(0.5_r8 - tmp_kxt/6.0_r8 + tmp_kxt2/24.0_r8)
            end if
            qgas_cur(igas) = tmp_q3
            tmp_qdel_cond = (tmp_q1 + tmp_pxt) - tmp_q3
            qgas_avg(igas) = tmp_q4
            do n = 1, n_mode
               if (uptkaer(igas,n) <= 0.0_r8) cycle
               tmpc = tmp_qdel_cond*(uptkaer(igas,n)/tmpa)
               qaer_cur(iaer,n) = qaer_prv(iaer,n) + tmpc
            end do

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
            if ( ldiag82 ) then
            if (i==1 .and. k==pver .and. igas==igas_h2so4) then
               tmp_q2 = tmp_q1 + tmp_pxt

               write(lun82,'(/a,2i5,1p,8e17.9)') 'gasaer - i, k, sum_uprt_so4, qav', &
                  i, k, tmpa, -1.0
               tmp_q2 = max( 1.0e-30_r8, tmp_q2 )
               write(lun82,'(/a,2i5,1p,8e17.9)') 'gasaer - i, k, q1, q2, q3, q4   ', &
                  i, k, tmp_q1, tmp_q2, tmp_q3, tmp_q4
               write(lun82,'(/a,2i5,1p,8e17.9)') 'gasaer - i, k, k*t, p*t, p/k, t ', &
                  i, k, tmp_kxt, tmp_pxt, tmp_pok, dtsubstep, tmp_qdel_cond
            end if
            end if
#endif

         else
         ! tmp_kxt < 1.0e-20_r8 so uptake to aerosols ~= 0.0
         ! in this case, do not bother to update qaer_cur
            tmp_q3 = tmp_q1 + tmp_pxt
            tmp_q4 = tmp_q1 + tmp_pxt*0.5_r8
            qgas_cur(igas) = tmp_q3
            qgas_avg(igas) = tmp_q4
         end if
      end do ! igas

      if ( igas_nh3 > 0 ) then
! do not allow nh4 to exceed 2*so4 (molar basis)
         iaer = iaer_nh4 ; igas = igas_nh3
         do n = 1, n_mode
            if (uptkaer(igas,n) <= 0.0_r8) cycle
            tmpa = qaer_cur(iaer,n) - 2.0_r8*qaer_cur(iaer_so4,n)
            if (tmpa > 0.0_r8) then
               qaer_cur(iaer,n) = qaer_cur(iaer,n) - tmpa
               qgas_cur(igas)   = qgas_cur(igas)   + tmpa
               qgas_avg(igas) = qgas_avg(igas) + tmpa*0.5_r8
            end if
         end do
      end if


      return
      end subroutine mam_gasaerexch_1subarea


!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine mam_soaexch_1subarea(                              &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         dtsubstep,                                                 &
         temp,              pmid,             aircon,               &
         n_mode,                                                    &
         qgas_cur,          qgas_avg,                               &
         qaer_cur,                                                  &
         qnum_cur,                                                  &
         qwtr_cur,                                                  &
         uptkaer                                                    )
!
! calculate soa condensation/evaporation for i,k,jsub over time dtsubstep
!

! uses
      use modal_aero_data, only: lptr2_soa_a_amode


      implicit none

! arguments
      integer,  intent(in) :: nstep                 ! model time-step number
      integer,  intent(in) :: lchnk                 ! chunk identifier
      integer,  intent(in) :: i, k                  ! column and level indices
      integer,  intent(in) :: jsub                  ! sub-area index
      integer,  intent(in) :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in) :: lund                  ! logical unit for diagnostic output
      integer,  intent(in) :: n_mode                ! current number of modes (including temporary)

      real(r8), intent(in) :: dtsubstep        ! current integration timestep (s)
      real(r8), intent(in) :: temp             ! temperature (K)
      real(r8), intent(in) :: pmid             ! pressure at model levels (Pa)
      real(r8), intent(in) :: aircon           ! air molar concentration (kmol/m3)

      real(r8), intent(inout), dimension( 1:max_gas ) :: &
         qgas_cur, qgas_avg
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur
      real(r8), intent(in   ), dimension( 1:max_gas, 1:max_mode ) :: &
         uptkaer

! local
      integer, parameter :: ntot_poaspec = npoa
      integer, parameter :: ntot_soaspec = nsoa

      integer :: iaer, igas, ip
      integer :: ll
      integer :: n, niter, niter_max
      integer :: ntot_soamode

      logical, parameter :: flag_pcarbon_opoa_frac_zero   = .true.

      logical :: skip_soamode(max_mode)   ! true if this mode does not have soa

      real(r8), dimension( 1:max_gas ) :: &
         qgas_prv

      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_prv

      real(r8) :: uptkaer_soag_tmp(nsoa,max_mode)

      real(r8), parameter :: a_min1 = 1.0e-20
      real(r8), parameter :: g_min1 = 1.0e-20
      real(r8), parameter :: alpha_astem = 0.05_r8 ! parameter used in calc of time step
      real(r8), parameter :: dtsub_fixed = -1.0    ! fixed sub-step for time integration (s)
!     real(r8), parameter :: dtsub_fixed = 10.0    ! fixed sub-step for time integration (s)
      real(r8), parameter :: rgas = 8.3144_r8      ! gas constant in J/K/mol

      real(r8) :: a_ooa_sum_tmp(max_mode)          ! total ooa (=soa+opoa) in a mode
      real(r8) :: a_opoa(max_mode)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa(ntot_soaspec,max_mode)     ! soa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa_tmp(ntot_soaspec,max_mode) ! temporary soa aerosol mixrat (mol/mol)
      real(r8) :: beta(ntot_soaspec,max_mode)      ! dtcur*xferrate
      real(r8) :: delh_vap_soa(ntot_soaspec)       ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
      real(r8) :: del_g_soa_tmp(ntot_soaspec)
      real(r8) :: dtcur                            ! current time step (s)
      real(r8) :: dtfull                           ! full time step (s)
      real(r8) :: dtmax                            ! = (dtfull-tcur)
      real(r8) :: dtsum_qgas_avg
      real(r8) :: g0_soa(ntot_soaspec)             ! ambient soa gas equilib mixrat (mol/mol at actual mw)
      real(r8) :: g_soa(ntot_soaspec)              ! soa gas mixrat (mol/mol at actual mw)
      real(r8) :: g_star(ntot_soaspec,max_mode)    ! soa gas mixrat that is in equilib
                                                   ! with each aerosol mode (mol/mol)
      real(r8) :: mw_poa(ntot_poaspec)             ! actual molec wght of poa
      real(r8) :: mw_soa(ntot_soaspec)             ! actual molec wght of soa
      real(r8) :: opoa_frac(ntot_poaspec,max_mode) ! fraction of poa that is opoa
      real(r8) :: phi(ntot_soaspec,max_mode)       ! "relative driving force"
      real(r8) :: p0_soa(ntot_soaspec)             ! soa gas equilib vapor presssure (atm)
      real(r8) :: p0_soa_298(ntot_soaspec)         ! p0_soa_298 = soa gas equilib vapor presssure (atm) at 298 k
      real(r8) :: sat(ntot_soaspec,max_mode)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                   !    used by the numerical integration scheme -- it is not a saturation rato!
      real(r8) :: tcur                             ! current integration time (from 0 s)

      real(r8) :: tmpa, tmpb, tmpc

      real(r8) :: tot_soa(ntot_soaspec)            ! g_soa + sum( a_soa(:) )


! calc ntot_soamode = "last" mode on which soa is allowed to condense
      ntot_soamode = 0
      do n = 1, ntot_amode
         if (n == nufi) cycle
         if (mode_aging_optaa(n) > 0) ntot_soamode = n
         if (lptr2_soa_a_amode(n,1) > 0) ntot_soamode = n
      end do
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( i*k == top_lev .and. ldiagd1 ) write(lund,'(/a,5i5)') &
         'ntot_amode, ntot_amode_extd, n_mode, ntot_soamode', &
          ntot_amode, ntot_amode_extd, n_mode, ntot_soamode
#endif

      opoa_frac = 0.1_r8
! for primary carbon mode, set opoa_frac=0 for consistency with older code
! (this could be changed)
      if ( flag_pcarbon_opoa_frac_zero ) then
         if (npca > 0) opoa_frac(:,npca) = 0.0_r8
      end if

      delh_vap_soa = 156.0e3
!     delh_vap_soa =  30.0e3  ! 11-jun-2012
      p0_soa_298 = 1.0e-10

! calc ambient equilibrium soa gas
      do ll = 1, ntot_soaspec
         p0_soa(ll) = p0_soa_298(ll) * &
                  exp( -(delh_vap_soa(ll)/rgas)*((1.0/temp)-(1.0/298.0)) )
         g0_soa(ll) = 1.01325e5*p0_soa(ll)/pmid
      end do

      niter_max = 1000
      niter = 0
      dtfull = dtsubstep
      tcur = 0.0
      dtcur = 0.0
      phi(:,:) = 0.0
      g_star(:,:) = 0.0
      g_soa(:) = 0.0
      a_opoa(:) = 0.0
      a_soa(:,:) = 0.0

!
! main integration loop -- does multiple substeps to reach dtfull
!
      qgas_avg(1:nsoa) = 0.0_r8
      dtsum_qgas_avg = 0.0_r8

time_loop: &
      do while (tcur < dtfull-1.0e-3_r8 )

      niter = niter + 1
      if (niter > niter_max) exit


! set qxxx_prv to be current value
      qgas_prv(1:nsoa) = qgas_cur(1:nsoa)
      qaer_prv = qaer_cur
!     qaer_num = qnum_cur


! determine which modes have non-zero transfer rates
!    and are involved in the soa gas-aerosol transfer
! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
      do n = 1, ntot_soamode
         skip_soamode(n) = .true.
         do ll = 1, ntot_soaspec
            if (uptkaer(ll,n) > 1.0e-15_r8) then
               uptkaer_soag_tmp(ll,n) = uptkaer(ll,n)
               skip_soamode(n) = .false.
            else
               uptkaer_soag_tmp(ll,n) = 0.0_r8
            end if
         end do
      end do

! load incoming soag and soaa into temporary arrays
! force things to be non-negative
! calc tot_soa(ll)
! calc a_opoa (always slightly >0)
!
! *** questions ***
! > why not use qgas and qaer instead of g_soa and a_soa
! > why not calc the following on every substep because 
!      nuc and coag may change things:
!      skip)soamode, uptkaer_soag_tmp, tot_soa, a_opoa
! > include gasprod for soa ??
! > create qxxx_bgn = qxxx_cur at the very beginning (is it needed)
!
      do ll = 1, ntot_soaspec
         g_soa(ll) = max( qgas_prv(ll), 0.0_r8 )
         tot_soa(ll) = g_soa(ll)
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            a_soa(ll,n) = max( qaer_prv(ll,n), 0.0_r8 )
            tot_soa(ll) = tot_soa(ll) + a_soa(ll,n)
         end do
      end do

      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         a_opoa(n) = 0.0_r8
         do ll = 1, ntot_poaspec
            a_opoa(n) = a_opoa(n) + opoa_frac(ll,n) * max( qaer_prv(iaer_pom+ll-1,n), 0.0_r8 )
         end do
      end do


! determine time step
      tmpa = 0.0  ! time integration parameter for all soa species
      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         a_ooa_sum_tmp(n) = a_opoa(n) + sum( a_soa(1:ntot_soaspec,n) )
      end do
      do ll = 1, ntot_soaspec
         tmpb = 0.0  ! time integration parameter for a single soa species
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            sat(ll,n) = g0_soa(ll)/max( a_ooa_sum_tmp(n), a_min1 )
            g_star(ll,n) = sat(ll,n)*a_soa(ll,n)
            phi(ll,n) = (g_soa(ll) - g_star(ll,n))/max( g_soa(ll), g_star(ll,n), g_min1 )
            tmpb = tmpb + uptkaer_soag_tmp(ll,n)*abs(phi(ll,n))
         end do
         tmpa = max( tmpa, tmpb )
      end do

      if (dtsub_fixed > 0.0_r8) then
         dtcur = dtsub_fixed
         tcur = tcur + dtcur
      else
         dtmax = dtfull-tcur
         if (dtmax*tmpa <= alpha_astem) then
! here alpha_astem/tmpa >= dtmax, so this is final substep
            dtcur = dtmax
            tcur = dtfull
         else
            dtcur = alpha_astem/tmpa
            tcur = tcur + dtcur
         end if
      end if


! step 1 - for modes where soa is condensing, estimate "new" a_soa(ll,n)
!    using an explicit calculation with "old" g_soa
!    and g_star(ll,n) calculated using "old" a_soa(ll,n)
! do this to get better estimate of "new" a_soa(ll,n) and sat(ll,n)
      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         do ll = 1, ntot_soaspec
            ! first ll loop calcs a_soa_tmp(ll,n) & a_ooa_sum_tmp
            a_soa_tmp(ll,n) = a_soa(ll,n)
            beta(ll,n) = dtcur*uptkaer_soag_tmp(ll,n)
            del_g_soa_tmp(ll) = g_soa(ll) - g_star(ll,n)
            if (del_g_soa_tmp(ll) > 0.0_r8) then
               a_soa_tmp(ll,n) = a_soa(ll,n) + beta(ll,n)*del_g_soa_tmp(ll)
            end if
         end do
         a_ooa_sum_tmp(n) = a_opoa(n) + sum( a_soa_tmp(1:ntot_soaspec,n) )
         do ll = 1, ntot_soaspec
            ! second ll loop calcs sat & g_star
            if (del_g_soa_tmp(ll) > 0.0_r8) then
               sat(ll,n) = g0_soa(ll)/max( a_ooa_sum_tmp(n), a_min1 )
               g_star(ll,n) = sat(ll,n)*a_soa_tmp(ll,n)   ! this just needed for diagnostics
            end if
         end do
      end do


! step 2 - implicit in g_soa and semi-implicit in a_soa,
!    with g_star(ll,n) calculated semi-implicitly
      do ll = 1, ntot_soaspec
         tmpa = 0.0
         tmpb = 0.0
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            tmpa = tmpa + a_soa(ll,n)/(1.0_r8 + beta(ll,n)*sat(ll,n))
            tmpb = tmpb + beta(ll,n)/(1.0_r8 + beta(ll,n)*sat(ll,n))
         end do

         g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_r8 + tmpb)
         g_soa(ll) = max( 0.0_r8, g_soa(ll) )
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            a_soa(ll,n) = (a_soa(ll,n) + beta(ll,n)*g_soa(ll))/   &
                       (1.0_r8 + beta(ll,n)*sat(ll,n))
         end do
      end do


! update mix ratios for soa species
      do igas = 1, nsoa
         iaer = igas
         qgas_cur(igas) = g_soa(igas)
         tmpc = qgas_cur(igas) - qgas_prv(igas)
         qgas_avg(igas) = qgas_avg(igas) + dtcur*(qgas_prv(igas) + 0.5_r8*tmpc)
         do n = 1, ntot_soamode
            qaer_cur(iaer,n) = a_soa(iaer,n)
            tmpc = qaer_cur(iaer,n) - qaer_prv(iaer,n)
         end do
      end do


      dtsum_qgas_avg = dtsum_qgas_avg + dtcur

      end do time_loop

! convert qgas_avg from sum_over[ qgas*dt_cut ] to an average
      do igas = 1, nsoa
         qgas_avg(igas) = max( 0.0_r8, qgas_avg(igas)/dtsum_qgas_avg )
      end do


      return
      end subroutine mam_soaexch_1subarea


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine mam_rename_1subarea(                               &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         iscldy_subarea,                                            &
         mtoo_renamexf,                                             &
         n_mode,                                                    &
         qnum_cur,                                                  &
         qaer_cur,          qaer_del_grow4rnam,                     &
         qwtr_cur,                                                  &
         qnumcw_cur,                                                &
         qaercw_cur,        qaercw_del_grow4rnam                    )

#if ( defined CAM_VERSION_IS_ACME ) 
      use shr_spfn_mod, only: erfc => shr_spfn_erfc  ! acme version of cam
#else
      use error_function,  only: erfc                ! mozart-mosaic version of cam
#endif

      logical,  intent(in)    :: iscldy_subarea        ! true if sub-area is cloudy
      integer,  intent(in)    :: nstep                 ! model time-step number
      integer,  intent(in)    :: lchnk                 ! chunk identifier
      integer,  intent(in)    :: i, k                  ! column and level indices
      integer,  intent(in)    :: jsub                  ! sub-area index
      integer,  intent(in)    :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in)    :: lund                  ! logical unit for diagnostic output
      integer,  intent(in)    :: mtoo_renamexf(max_mode)
      integer,  intent(in)    :: n_mode                ! current number of modes (including temporary)

      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur
      real(r8), intent(in   ), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_del_grow4rnam
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur

      real(r8), intent(inout), optional, dimension( 1:max_mode ) :: &
         qnumcw_cur
      real(r8), intent(inout), optional, dimension( 1:max_aer, 1:max_mode ) :: &
         qaercw_cur
      real(r8), intent(in   ), optional, dimension( 1:max_aer, 1:max_mode ) :: &
         qaercw_del_grow4rnam


! !DESCRIPTION: 
! computes TMR (tracer mixing ratio) tendencies for "mode renaming"
!    during a continuous growth process
! currently this transfers number and mass (and surface) from the aitken
!    to accumulation mode after gas condensation or stratiform-cloud
!    aqueous chemistry
! (convective cloud aqueous chemistry not yet implemented)
!
! !REVISION HISTORY:
!

! local variables
      integer :: iaer
      integer :: mfrm, mtoo
      integer :: n, npair

      integer, parameter :: ldiag1 = 0

      real(r8), parameter :: frelax = 27.0_r8
      real(r8), parameter :: onethird = 1.0_r8/3.0_r8

      real(r8) :: deldryvol_a(ntot_amode)
      real(r8) :: deldryvol_c(ntot_amode)
      real(r8) :: dp_belowcut(max_mode)
      real(r8) :: dp_cut(max_mode)
      real(r8) :: dgn_aftr, dgn_xfer
      real(r8) :: dgn_t_new, dgn_t_old, dgn_t_oldaa
      real(r8) :: dryvol_t_del, dryvol_t_new
      real(r8) :: dryvol_t_old, dryvol_t_oldaa, dryvol_t_oldbnd
      real(r8) :: dryvol_a(ntot_amode)
      real(r8) :: dryvol_c(ntot_amode)
      real(r8) :: dryvol_smallest(ntot_amode)
      real(r8) :: factoraa(ntot_amode)
      real(r8) :: factoryy(ntot_amode)
      real(r8) :: lndp_cut(max_mode)
      real(r8) :: lndgn_new, lndgn_old
      real(r8) :: lndgv_new, lndgv_old
      real(r8) :: num_t_old, num_t_oldbnd
      real(r8) :: tailfr_volnew, tailfr_volold
      real(r8) :: tailfr_numnew, tailfr_numold
      real(r8) :: tmpa, tmpb, tmpd
      real(r8) :: tmp_alnsg2(max_mode)
      real(r8) :: v2nhirlx(ntot_amode), v2nlorlx(ntot_amode)
      real(r8) :: xfercoef, xfertend
      real(r8) :: xferfrac_vol, xferfrac_num, xferfrac_max
      real(r8) :: yn_tail, yv_tail


      xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

! calculate variable used in the renamingm mode" of each renaming pair
! also compute dry-volume change during the continuous growth process
      npair = 0
      do n = 1, ntot_amode
         mtoo = mtoo_renamexf(n)
         if (mtoo <= 0) cycle

         npair = npair + 1
         mfrm = n
         factoraa(mfrm) = (pi/6.)*exp(4.5*(alnsg_aer(mfrm)**2))
         factoraa(mtoo) = (pi/6.)*exp(4.5*(alnsg_aer(mtoo)**2))
         factoryy(mfrm) = sqrt( 0.5 )/alnsg_aer(mfrm)
! dryvol_smallest is a very small volume mixing ratio (m3-AP/kmol-air)
! used for avoiding overflow.  it corresponds to dp = 1 nm
! and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air
         dryvol_smallest(mfrm) = 1.0e-25
!        v2nlorlx(mfrm) = voltonumblo_amode(mfrm)*frelax
!        v2nhirlx(mfrm) = voltonumbhi_amode(mfrm)/frelax
         v2nlorlx(mfrm) = ( 1._r8 / ( (pi/6._r8)* &
            (dgnumlo_aer(mfrm)**3._r8)*exp(4.5_r8*alnsg_aer(mfrm)**2._r8) ) ) * frelax
         v2nhirlx(mfrm) = ( 1._r8 / ( (pi/6._r8)* &
            (dgnumhi_aer(mfrm)**3._r8)*exp(4.5_r8*alnsg_aer(mfrm)**2._r8) ) ) / frelax

         tmp_alnsg2(mfrm) = 3.0 * (alnsg_aer(mfrm)**2)
         dp_cut(mfrm) = sqrt(   &
            dgnum_aer(mfrm)*exp(1.5*(alnsg_aer(mfrm)**2)) *   &
            dgnum_aer(mtoo)*exp(1.5*(alnsg_aer(mtoo)**2)) )
         lndp_cut(mfrm) = log( dp_cut(mfrm) )
         dp_belowcut(mfrm) = 0.99*dp_cut(mfrm)
      end do
      if (npair <= 0) return

! compute aerosol dry-volume for the "from mode" of each renaming pair
! also compute dry-volume change during the continuous growth process
      do n = 1, ntot_amode
         mtoo = mtoo_renamexf(n)
         if (mtoo <= 0) cycle

         tmpa = 0.0_r8 ; tmpb = 0.0_r8
         do iaer = 1, naer
!   fac_m2v_aer converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
            tmpa = tmpa + qaer_cur(iaer,n)*fac_m2v_aer(iaer)
            tmpb = tmpb + qaer_del_grow4rnam(iaer,n)*fac_m2v_aer(iaer)
         end do
         dryvol_a(n) = tmpa-tmpb ! dry volume before growth
         deldryvol_a(n) = tmpb   ! change to dry volume due to growth

         if ( iscldy_subarea ) then
         tmpa = 0.0_r8 ; tmpb = 0.0_r8
         do iaer = 1, naer
!   fac_m2v_aer converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
            tmpa = tmpa + qaercw_cur(iaer,n)*fac_m2v_aer(iaer)
            tmpb = tmpb + qaercw_del_grow4rnam(iaer,n)*fac_m2v_aer(iaer)
         end do
         dryvol_c(n) = tmpa-tmpb 
         deldryvol_c(n) = tmpb
         end if ! ( iscldy_subarea ) then

      end do


!
!   loop over renaming pairs
!
mainloop1_ipair:  do n = 1, ntot_amode

      mfrm = n
      mtoo = mtoo_renamexf(n)
      if (mtoo <= 0) cycle mainloop1_ipair

!   dryvol_t_old is the old total (a+c) dry-volume for the "from" mode 
!      in m^3-AP/kmol-air
!   dryvol_t_new is the new total dry-volume
!      (old/new = before/after the continuous growth)
!   num_t_old is total number in particles/kmol-air
      if ( iscldy_subarea ) then
         dryvol_t_old = dryvol_a(mfrm) + dryvol_c(mfrm)
         dryvol_t_del = deldryvol_a(mfrm) + deldryvol_c(mfrm)
         num_t_old = (qnum_cur(mfrm) + qnumcw_cur(mfrm))
      else
         dryvol_t_old = dryvol_a(mfrm)
         dryvol_t_del = deldryvol_a(mfrm)
         num_t_old = qnum_cur(mfrm)
      end if
      dryvol_t_new = dryvol_t_old + dryvol_t_del

!   no renaming if dryvol_t_new ~ 0 or dryvol_t_del ~ 0
      if (dryvol_t_new .le. dryvol_smallest(mfrm)) cycle mainloop1_ipair
      dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
      if (rename_method_optaa .ne. 40) then
         if (dryvol_t_del .le. 1.0e-6*dryvol_t_oldbnd) cycle mainloop1_ipair
      end if

      num_t_old = max( 0.0_r8, num_t_old )
      dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
      num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
      num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )

!   no renaming if dgnum < "base" dgnum, 
      dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
      if (dgn_t_new .le. dgnum_aer(mfrm)) cycle mainloop1_ipair

!   compute new fraction of number and mass in the tail (dp > dp_cut)
      lndgn_new = log( dgn_t_new )
      lndgv_new = lndgn_new + tmp_alnsg2(mfrm)
      yn_tail = (lndp_cut(mfrm) - lndgn_new)*factoryy(mfrm)
      yv_tail = (lndp_cut(mfrm) - lndgv_new)*factoryy(mfrm)
      tailfr_numnew = 0.5_r8*erfc( yn_tail )
      tailfr_volnew = 0.5_r8*erfc( yv_tail )

!   compute old fraction of number and mass in the tail (dp > dp_cut)
      dgn_t_old =   &
            (dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird
      dgn_t_oldaa = dgn_t_old
      dryvol_t_oldaa = dryvol_t_old
      
      if (rename_method_optaa .eq. 40) then
         if (dgn_t_old .gt. dp_belowcut(mfrm)) then
            ! this revised volume corresponds to dgn_t_old == dp_belowcut, and same number conc
            dryvol_t_old = dryvol_t_old * (dp_belowcut(mfrm)/dgn_t_old)**3
            dgn_t_old = dp_belowcut(mfrm)
         end if
         if ((dryvol_t_new-dryvol_t_old) .le. 1.0e-6_r8*dryvol_t_oldbnd) cycle mainloop1_ipair
      else if (dgn_t_new .ge. dp_cut(mfrm)) then
!         if dgn_t_new exceeds dp_cut, use the minimum of dgn_t_old and 
!         dp_belowcut to guarantee some transfer
          dgn_t_old = min( dgn_t_old, dp_belowcut(mfrm) )
      end if
      lndgn_old = log( dgn_t_old )
      lndgv_old = lndgn_old + tmp_alnsg2(mfrm)
      yn_tail = (lndp_cut(mfrm) - lndgn_old)*factoryy(mfrm)
      yv_tail = (lndp_cut(mfrm) - lndgv_old)*factoryy(mfrm)
      tailfr_numold = 0.5_r8*erfc( yn_tail )
      tailfr_volold = 0.5_r8*erfc( yv_tail )

!   transfer fraction is difference between new and old tail-fractions
!   transfer fraction for number cannot exceed that of mass
      tmpa = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_old
      if (tmpa .le. 0.0_r8) cycle mainloop1_ipair

      xferfrac_vol = min( tmpa, dryvol_t_new )/dryvol_t_new
      xferfrac_vol = min( xferfrac_vol, xferfrac_max ) 
      xferfrac_num = tailfr_numnew - tailfr_numold
      xferfrac_num = max( 0.0_r8, min( xferfrac_num, xferfrac_vol ) )
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( ldiag98 ) write(lun98,'(/a,2i3,1p,10e11.3)') &
         'rename i,k, xf n/v', i, k, xferfrac_num, xferfrac_vol
#endif

#if ( defined( CAMBOX_NEVER_ACTIVATE_THIS ) )
!   diagnostic output start ----------------------------------------
       if (ldiag1 > 0) then
       icol_diag = -1
       if ((lonndx(i) == 37) .and. (latndx(i) == 23)) icol_diag = i
       if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
 !      write(lund,97010) fromwhere, nstep, lchnk, i, k, ipair
       write(lund,97010) fromwhere, nstep, latndx(i), lonndx(i), k, ipair
       write(lund,97020) 'drv olda/oldbnd/old/new/del',   &
             dryvol_t_oldaa, dryvol_t_oldbnd, dryvol_t_old, dryvol_t_new, dryvol_t_del
       write(lund,97020) 'num old/oldbnd, dgnold/new ',   &
             num_t_old, num_t_oldbnd, dgn_t_old, dgn_t_new
       write(lund,97020) 'tailfr v_old/new, n_old/new',   &
             tailfr_volold, tailfr_volnew, tailfr_numold, tailfr_numnew
       tmpa = max(1.0d-10,xferfrac_vol) / max(1.0d-10,xferfrac_num)
       dgn_xfer = dgn_t_new * tmpa**onethird
       tmpa = max(1.0d-10,(1.0d0-xferfrac_vol)) /   &
               max(1.0d-10,(1.0d0-xferfrac_num))
       dgn_aftr = dgn_t_new * tmpa**onethird
       write(lund,97020) 'xferfrac_v/n; dgn_xfer/aftr',   &
             xferfrac_vol, xferfrac_num, dgn_xfer, dgn_aftr
 !97010      format( / 'RENAME ', a, '  nx,lc,i,k,ip', i8, 4i4 )
 97010      format( / 'RENAME ', a, '  nx,lat,lon,k,ip', i8, 4i4 )
 97020      format( a, 6(1pe15.7) )
       end if
       end if ! (ldiag1 > 0)
!   diagnostic output end   ------------------------------------------
#endif


!
!   compute changes to number and species masses
!
      tmpd = qnum_cur(mfrm)*xferfrac_num
      qnum_cur(mfrm) = qnum_cur(mfrm) - tmpd
      qnum_cur(mtoo) = qnum_cur(mtoo) + tmpd
      do iaer = 1, naer
         tmpd = qaer_cur(iaer,mfrm)*xferfrac_vol
         qaer_cur(iaer,mfrm) = qaer_cur(iaer,mfrm) - tmpd
         qaer_cur(iaer,mtoo) = qaer_cur(iaer,mtoo) + tmpd
      end do ! iaer

      if ( iscldy_subarea ) then
      tmpd = qnumcw_cur(mfrm)*xferfrac_num
      qnumcw_cur(mfrm) = qnumcw_cur(mfrm) - tmpd
      qnumcw_cur(mtoo) = qnumcw_cur(mtoo) + tmpd
      do iaer = 1, naer
         tmpd = qaercw_cur(iaer,mfrm)*xferfrac_vol
         qaercw_cur(iaer,mfrm) = qaercw_cur(iaer,mfrm) - tmpd
         qaercw_cur(iaer,mtoo) = qaercw_cur(iaer,mtoo) + tmpd
      end do ! iaer
      end if ! ( iscldy_subarea ) then


#if ( defined( CAMBOX_NEVER_ACTIVATE_THIS ) )
!   diagnostic output start ----------------------------------------
                if (ldiag1 > 0) then
                if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
                  if (lstooa .gt. 0) then
                    write(lund,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend),   &
                        deltat*dqdt(i,k,lstooa), deltat*(dqdt(i,k,lstooa) + xfertend)
                  else
                    write(lund,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend)
                  end if
                end if
                end if
!   diagnostic output end   ------------------------------------------
#endif


      end do mainloop1_ipair


      return
      end subroutine mam_rename_1subarea


!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine mam_newnuc_1subarea(                               &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         deltat,                                                    &
         temp,              pmid,             aircon,               &
         zmid,              pblh,             relhum,               &
         uptkrate_h2so4,   del_h2so4_gasprod, del_h2so4_aeruptk,    &
         n_mode,                                                    &
         qgas_cur,          qgas_avg,                               &
         qnum_cur,                                                  &
         qaer_cur,                                                  &
         qwtr_cur,                                                  &
         dnclusterdt                                                )

! uses
      use chem_mods,     only: adv_mass

      use modal_aero_newnuc, only: &
         mer07_veh02_nuc_mosaic_1box, qh2so4_cutoff

      implicit none

! arguments
      integer,  intent(in) :: nstep                 ! model time-step number
      integer,  intent(in) :: lchnk                 ! chunk identifier
      integer,  intent(in) :: i, k                  ! column and level indices
      integer,  intent(in) :: jsub                  ! sub-area index
      integer,  intent(in) :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in) :: lund                  ! logical unit for diagnostic output
      integer,  intent(in) :: n_mode                ! current number of modes (including temporary)

      real(r8), intent(in) :: deltat           ! model timestep (s)
      real(r8), intent(in) :: temp             ! temperature (K)
      real(r8), intent(in) :: pmid             ! pressure at model levels (Pa)
      real(r8), intent(in) :: aircon           ! air molar concentration (kmol/m3)
      real(r8), intent(in) :: zmid             ! midpoint height above surface (m)
      real(r8), intent(in) :: pblh             ! pbl height (m)
      real(r8), intent(in) :: relhum           ! relative humidity (0-1)
      real(r8), intent(in) :: uptkrate_h2so4
      real(r8), intent(in) :: del_h2so4_gasprod
      real(r8), intent(in) :: del_h2so4_aeruptk

      real(r8), intent(inout) :: dnclusterdt   ! cluster nucleation rate (#/m3/s)

      real(r8), intent(inout), dimension( 1:max_gas ) :: &
         qgas_cur
      real(r8), intent(in   ), dimension( 1:max_gas ) :: &
         qgas_avg
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur

! DESCRIPTION: 
!   computes changes due to aerosol nucleation (new particle formation)
!       treats both nucleation and subsequent growth of new particles
!          to aitken mode size
!   uses the following parameterizations
!       vehkamaki et al. (2002) parameterization for binary
!           homogeneous nucleation (h2so4-h2o) plus
!       kerminen and kulmala (2002) parameterization for
!           new particle loss during growth to aitken size
!
! REVISION HISTORY:
!   R.Easter 2007.09.14:  Adapted from MIRAGE2 code and CMAQ V4.6 code
!

! local variables
      integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
      integer, parameter :: newnuc_method_flagaa = 11
!     integer, parameter :: newnuc_method_flagaa = 12
        !  1=merikanto et al (2007) ternary   2=vehkamaki et al (2002) binary
        ! 11=merikanto ternary + first-order boundary layer
        ! 12=merikanto ternary + second-order boundary layer

      integer :: itmp
      integer :: l
      integer :: ldiagveh02
      integer :: m

      real(r8) :: dens_nh4so4a
      real(r8) :: dmdt_ait, dmdt_aitsv1, dmdt_aitsv2, dmdt_aitsv3
      real(r8) :: dndt_ait, dndt_aitsv1, dndt_aitsv2, dndt_aitsv3
      real(r8) :: dnh4dt_ait, dso4dt_ait
      real(r8) :: dpnuc
      real(r8) :: dplom_mode(1), dphim_mode(1)
      real(r8) :: mass1p
      real(r8) :: mass1p_aithi, mass1p_aitlo 
      real(r8) :: qh2so4_cur, qh2so4_avg, qh2so4_del
      real(r8) :: qnh3_cur, qnh3_del, qnh4a_del
      real(r8) :: qnuma_del
      real(r8) :: qso4a_del
      real(r8) :: relhumnn
      real(r8) :: tmpa, tmpb, tmpc
      real(r8) :: tmp_q2, tmp_q3
      real(r8) :: tmp_q_del
      real(r8) :: tmp_frso4, tmp_uptkrate

      character(len=1) :: tmpch1, tmpch2, tmpch3


! begin
      dnclusterdt = 0.0_r8

! qh2so4_cur = current qh2so4, after aeruptk
! qh2so4_avg = average qh2so4 over time-step
      qh2so4_cur = qgas_cur(igas_h2so4)

      if ( (gaexch_h2so4_uptake_optaa == 1) .and. &
           (newnuc_h2so4_conc_optaa   == 1) ) then
! estimate qh2so4_avg using the method in standard cam5.2 modal_aero_newnuc

         ! skip if h2so4 vapor < qh2so4_cutoff
         if (qh2so4_cur <= qh2so4_cutoff) goto 80000

         tmpa = max( 0.0_r8, del_h2so4_gasprod )
         tmp_q3 = qh2so4_cur
         ! tmp_q2 = qh2so4 before aeruptk
         ! (note tmp_q3, tmp_q2 both >= 0.0)
         tmp_q2 = tmp_q3 + max( 0.0_r8, -del_h2so4_aeruptk )

         ! tmpb = log( tmp_q2/tmp_q3 ) BUT with some checks added
         if (tmp_q2 <= tmp_q3) then
            tmpb = 0.0_r8
         else
            tmpc = tmp_q2 * exp( -20.0_r8 )
            if (tmp_q3 <= tmpc) then
               tmp_q3 = tmpc
               tmpb = 20.0_r8
            else
               tmpb = log( tmp_q2/tmp_q3 )
            end if
         end if
         ! d[ln(qh2so4)]/dt (1/s) from uptake (condensation) to aerosol
         tmp_uptkrate = tmpb/deltat

!   qh2so4_avg = estimated average qh2so4
!   when production & loss are done simultaneously
         if (tmpb <= 0.1_r8) then
            qh2so4_avg = tmp_q3*(1.0_r8 + 0.5_r8*tmpb) - 0.5_r8*tmpa
         else
            tmpc = tmpa/tmpb
            qh2so4_avg = (tmp_q3 - tmpc)*((exp(tmpb)-1.0_r8)/tmpb) + tmpc
         end if
      else
! use qh2so4_avg and first-order loss rate calculated in mam_gasaerexch_1subarea
         qh2so4_avg = qgas_avg(igas_h2so4)
         tmp_uptkrate = uptkrate_h2so4
      end if

      if (qh2so4_avg <= qh2so4_cutoff) goto 80000

      if (igas_nh3 > 0) then
          qnh3_cur = max( 0.0_r8, qgas_cur(igas_nh3) )
      else
          qnh3_cur = 0.0_r8
      end if

!   dry-diameter limits for "grown" new particles
      dplom_mode(1) = exp( 0.67_r8*log(dgnumlo_aer(nait))   &
                         + 0.33_r8*log(dgnum_aer(nait)) )
      dphim_mode(1) = dgnumhi_aer(nait)

!   mass1p_... = mass (kg) of so4 & nh4 in a single particle of diameter ...
!                (assuming same dry density for so4 & nh4)
!      mass1p_aitlo - dp = dplom_mode(1)
!      mass1p_aithi - dp = dphim_mode(1)
      tmpa = dens_so4a_host*pi/6.0_r8
      mass1p_aitlo = tmpa*(dplom_mode(1)**3)
      mass1p_aithi = tmpa*(dphim_mode(1)**3)

!   limit RH to between 0.1% and 99%
      relhumnn = max( 0.01_r8, min( 0.99_r8, relhum ) )


!   call ... routine to get nucleation rates
      ldiagveh02 = -1
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if (ldiag2 > 0) then
      if ((lonndx == 37) .and. (latndx == 23)) then
      if ((k >= 24) .or. (mod(k,4) == 0)) then
         ldiagveh02 = +1
         write(lund,'(/a,i8,3i4,f8.2,1p,4e10.2)')   &
            'veh02 call - nstep,lat,lon,k; tk,rh,p,cair',   &
            nstep, latndx, lonndx, k,   &
            temp, relhumnn, pmid, aircon*1.0e3_r8
         ! output aircon at (mol/m3)
      end if
      end if
      end if ! (ldiag2 > 0)
#endif

      call mer07_veh02_nuc_mosaic_1box(   &
           newnuc_method_flagaa,   &
           deltat, temp, relhumnn, pmid,   &
           zmid, pblh,   &
           qh2so4_cur, qh2so4_avg, qnh3_cur, tmp_uptkrate,   &
           mw_so4a_host,   &
           1, 1, dplom_mode, dphim_mode,   &
           itmp, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a,   &
           ldiagveh02, dnclusterdt )
!----------------------------------------------------------------------
!       subr mer07_veh02_nuc_mosaic_1box(   &
!          newnuc_method_flagaa,   &
!          dtnuc, temp_in, rh_in, press_in,   &
!          qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
!          nsize, maxd_asize, dplom_sect, dphim_sect,   &
!          isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
!          qh2so4_del, qnh3_del, dens_nh4so4a )
!
!! subr arguments (in)
!        real(r8), intent(in) :: dtnuc             ! nucleation time step (s)
!        real(r8), intent(in) :: temp_in           ! temperature, in k
!        real(r8), intent(in) :: rh_in             ! relative humidity, as fraction
!        real(r8), intent(in) :: press_in          ! air pressure (pa)
!
!        real(r8), intent(in) :: qh2so4_cur, qh2so4_avg
!                                                  ! gas h2so4 mixing ratios (mol/mol-air)
!        real(r8), intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
!             ! qxxx_cur = current value (after gas chem and condensation)
!             ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
!        real(r8), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)

!
!        integer, intent(in) :: nsize                    ! number of aerosol size bins
!        integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
!        real(r8), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
!        real(r8), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)
!
!! subr arguments (out)
!        integer, intent(out) :: isize_nuc         ! size bin into which new particles go
!        real(r8), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
!        real(r8), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
!        real(r8), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
!        real(r8), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
!        real(r8), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
!                                                  ! aerosol changes are > 0; gas changes are < 0
!        real(r8), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)
!----------------------------------------------------------------------


!   convert qnuma_del from (#/mol-air) to (#/kmol-air)
      qnuma_del = qnuma_del*1.0e3_r8
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
        if ( ldiag97 ) then
        write(lun97,'(/a,2i5,1p,3e11.3,2x,3e11.3)') &
           'newnuc - i, k, qavg s/n, uprt, rh 1/2', &
           i, k, qh2so4_avg, qh2so4_cur, qnh3_cur, &
           tmp_uptkrate, relhum, relhumnn
        write(lun97,'( a,10x,1p,3e11.3,2x,3e11.3)') &
           '               del qn, qso4a, qnh4a  ', &
           qnuma_del, qso4a_del, qnh4a_del
        end if
#endif

!   number nuc rate (#/kmol-air/s) from number nuc amt
      dndt_ait = qnuma_del/deltat

!   fraction of mass nuc going to so4
      tmpa = qso4a_del*mw_so4a_host
      if (igas_nh3 > 0) then
         tmpb = tmpa + qnh4a_del*mw_nh4a_host
         tmp_frso4 = max( tmpa, 1.0e-35_r8 )/max( tmpb, 1.0e-35_r8 )
      else
         tmpb = tmpa
         tmp_frso4 = 1.0_r8
      end if

!   mass nuc rate (kg/kmol-air/s) from mass nuc amts
      dmdt_ait = max( 0.0_r8, (tmpb/deltat) ) 

      dndt_aitsv1 = dndt_ait
      dmdt_aitsv1 = dmdt_ait
      dndt_aitsv2 = 0.0
      dmdt_aitsv2 = 0.0
      dndt_aitsv3 = 0.0
      dmdt_aitsv3 = 0.0
      tmpch1 = ' '
      tmpch2 = ' '

      if (dndt_ait < 1.0e2) then
!   ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
         dndt_ait = 0.0
         dmdt_ait = 0.0
         tmpch1 = 'A'

      else
         dndt_aitsv2 = dndt_ait
         dmdt_aitsv2 = dmdt_ait
         tmpch1 = 'B'

!   mirage2 code checked for complete h2so4 depletion here,
!   but this is now done in mer07_veh02_nuc_mosaic_1box
         mass1p = dmdt_ait/dndt_ait
         dndt_aitsv3 = dndt_ait
         dmdt_aitsv3 = dmdt_ait

!   apply particle size constraints
         if (mass1p < mass1p_aitlo) then
!   reduce dndt to increase new particle size
            dndt_ait = dmdt_ait/mass1p_aitlo
            tmpch1 = 'C'
         else if (mass1p > mass1p_aithi) then
!   reduce dmdt to decrease new particle size
            dmdt_ait = dndt_ait*mass1p_aithi
            tmpch1 = 'E'
         end if
      end if

! *** apply adjustment factor to avoid unrealistically high
!     aitken number concentrations in mid and upper troposphere
      dndt_ait = dndt_ait * newnuc_adjust_factor_dnaitdt
      dmdt_ait = dmdt_ait * newnuc_adjust_factor_dnaitdt

      tmp_q_del = dndt_ait*deltat
      qnum_cur(     nait) = qnum_cur(     nait) + tmp_q_del

!   dso4dt_ait, dnh4dt_ait are (kmol/kmol-air/s)
      dso4dt_ait = dmdt_ait*tmp_frso4/mw_so4a_host
      dnh4dt_ait = dmdt_ait*(1.0_r8 - tmp_frso4)/mw_nh4a_host

      if (dso4dt_ait > 0.0_r8) then
         tmp_q_del = dso4dt_ait*deltat
         qaer_cur(     iaer_so4,nait) = qaer_cur(     iaer_so4,nait) + tmp_q_del

         tmp_q_del = min( tmp_q_del, qgas_cur(igas_h2so4) )
         qgas_cur(     igas_h2so4) = qgas_cur(     igas_h2so4) - tmp_q_del
      end if

      if ((igas_nh3 > 0) .and. (dnh4dt_ait > 0.0_r8)) then
         tmp_q_del = dnh4dt_ait*deltat
         qaer_cur(     iaer_nh4,nait) = qaer_cur(     iaer_nh4,nait) + tmp_q_del

         tmp_q_del = min( tmp_q_del, qgas_cur(igas_nh3) )
         qgas_cur(     igas_nh3) = qgas_cur(     igas_nh3) - tmp_q_del
      end if

!!   temporary diagnostic
!        if (ldiag3 > 0) then
!        if ((dndt_ait /= 0.0_r8) .or. (dmdt_ait /= 0.0_r8)) then
!           write(lund,'(3a,1x,i7,3i5,1p,5e12.4)')   &
!              'newnucxx', tmpch1, tmpch2, nstep, lchnk, i, k,   &
!              dndt_ait, dmdt_ait, cldx
!!          call endrun( 'modal_aero_newnuc_sub' )
!        end if
!        end if


#if ( defined( CAMBOX_NEVER_ACTIVATE_THIS ) )
!   diagnostic output start ----------------------------------------
       if (ldiag4 > 0) then
       if ((lonndx == 37) .and. (latndx == 23)) then
       if ((k >= 24) .or. (mod(k,4) == 0)) then
        write(lund,97010) nstep, latndx, lonndx, k, temp, aircon*1.0e3_r8
        write(lund,97020) 'pmid                         ',   &
                pmid
        write(lund,97030) 'qv,qvsw,      rh_av, rh_clr  ',   &
                qv(i,k), qvswtr,       relhumav, relhum
        write(lund,97020) 'h2so4_cur,       _av, nh3_cur',   &
             qh2so4_cur,         qh2so4_avg, qnh3_cur
        write(lund,97020) 'del_h2so4_gasprod, _aeruptk  ',   &
             del_h2so4_gasprod(i,k), del_h2so4_aeruptk(i,k),   &
             tmp_uptkrate*3600.0
        write(lund,97020) ' '
        write(lund,97050) 'tmpch1, tmpch2               ', tmpch1, tmpch2
        write(lund,97020) 'dndt_, dmdt_aitsv1           ',   &
                          dndt_aitsv1, dmdt_aitsv1
        write(lund,97020) 'dndt_, dmdt_aitsv2           ',   &
                          dndt_aitsv2, dmdt_aitsv2
        write(lund,97020) 'dndt_, dmdt_aitsv3           ',   &
                          dndt_aitsv3, dmdt_aitsv3
        write(lund,97020) 'dndt_, dmdt_ait              ',   &
                          dndt_ait, dmdt_ait
        write(lund,97020) 'dso4dt_, dnh4dt_ait          ',   &
                          dso4dt_ait, dnh4dt_ait
        write(lund,97020) 'qso4a_del, qh2so4_del        ',   &
                          qso4a_del, qh2so4_del
        write(lund,97020) 'qnh4a_del, qnh3_del          ',   &
                          qnh4a_del, qnh3_del
        write(lund,97020) 'dqdt(h2so4), (nh3)           ',   &
              dqdt(i,k,l_h2so4), dqdt(i,k,l_nh3) 
        write(lund,97020) 'dqdt(so4a), (nh4a), (numa)   ',   &
              dqdt(i,k,lso4ait), dqdt(i,k,lnh4ait), dqdt(i,k,lnumait)
 
       dpnuc = 0.0
       if (dndt_aitsv1 > 1.0e-5) dpnuc = (6.0*dmdt_aitsv1/   &
                   (pi*dens_so4a_host*dndt_aitsv1))**0.3333333
       if (dpnuc > 0.0) then
       write(lund,97020) 'dpnuc,      dp_aitlo, _aithi ',   &
                    dpnuc, dplom_mode(1), dphim_mode(1)
       write(lund,97020) 'mass1p, mass1p_aitlo, _aithi ',   &
                    mass1p, mass1p_aitlo, mass1p_aithi
       end if
 
97010  format( / 'NEWNUC nstep,lat,lon,k,tk,cair', i8, 3i4, f8.2, 1pe12.4 )
97020  format( a, 1p, 6e12.4 )
97030  format( a, 1p, 2e12.4, 0p, 5f10.6 )
97040  format( 29x, 1p, 6e12.4 )
97050  format( a, 2(3x,a) )
       end if ! ((k >= 24) .or. (mod(k,4) == 0))
       end if ! ((lonndx == 37) .and. (latndx == 23))
       end if ! (ldiag4 > 0)
!   diagnostic output end   ------------------------------------------
#endif


80000 continue


      return
      end subroutine mam_newnuc_1subarea


!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine mam_coag_1subarea(                                 &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         deltat,                                                    &
         temp,              pmid,             aircon,               &
         dgn_a,             dgn_awet,         wetdens,              &
         n_mode,                                                    &
         qnum_cur,                                                  &
         qaer_cur,          qaer_del_coag_in,                       &
         qwtr_cur                                                   )

! coag between aitken, pcarbon, and accum modes
! inter-modal coag of ultrafine mode

! uses
      use modal_aero_coag, only: getcoags_wrapper_f

      implicit none

! arguments
      integer,  intent(in) :: nstep                 ! model time-step number
      integer,  intent(in) :: lchnk                 ! chunk identifier
      integer,  intent(in) :: i, k                  ! column and level indices
      integer,  intent(in) :: jsub                  ! sub-area index
      integer,  intent(in) :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in) :: lund                  ! logical unit for diagnostic output
      integer,  intent(in) :: n_mode                ! current number of modes (including temporary)

      real(r8), intent(in) :: deltat                ! model timestep (s)
      real(r8), intent(in) :: temp                  ! temperature at model levels (K)
      real(r8), intent(in) :: pmid                  ! pressure at layer center (Pa)
      real(r8), intent(in) :: aircon                ! air molar concentration (kmol/m3)
      real(r8), intent(in) :: dgn_a(max_mode)
      real(r8), intent(in) :: dgn_awet(max_mode)
                                    ! dry & wet geo. mean dia. (m) of number distrib.
      real(r8), intent(in) :: wetdens(max_mode)
                                    ! interstitial aerosol wet density (kg/m3)
                                    ! dry & wet geo. mean dia. (m) of number distrib.

      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur
      real(r8), intent(out), dimension( 1:max_aer, 1:max_agepair ) :: &
         qaer_del_coag_in
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur

! local variables
      integer :: iaer, ip
      integer :: modefrm, modetoo
      integer :: n

      real(r8), parameter :: epsilonx1 = epsilon( 1.0_r8 )
      real(r8), parameter :: epsilonx2 = epsilonx1*2.0_r8

      real(r8) :: tmp1, tmp2, tmp3, tmp4
      real(r8) :: tmpa, tmpb, tmpc, tmpn
      real(r8) :: tmp_dq, tmp_xf
      real(r8) :: xbetaij2i, xbetaij2j, xbetaii2, xbetajj2
      real(r8) :: ybetaij0(max_coagpair), ybetaij3(max_coagpair), &
                  ybetaii0(max_coagpair), ybetajj0(max_coagpair)

      real(r8), dimension( 1:max_mode ) :: &
         qnum_tmpa, qnum_tmpb, qnum_tmpc
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_tmpa, qaer_tmpb, qaer_tmpc

! DESCRIPTION: 

      qnum_tmpa = max( 0.0_r8, qnum_cur )
      qaer_tmpa = max( 0.0_r8, qaer_cur )
      qnum_tmpb = qnum_tmpa
      qaer_tmpb = qaer_tmpa
      qaer_del_coag_in = 0.0_r8

!
! compute coagulation rates using cmaq "fast" method
!    (based on E. Whitby's approximation approach)
! here subr. arguments are all in mks unit
!
      lun15n = 149 + i
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( ldiag15n ) write(lun15n,'(//a,3i5,1p,e11.3)') 'coag - nstep,i,k', nstep, i, k, aircon
#endif

      do ip = 1, n_coagpair
         modefrm = modefrm_coagpair(ip)
         modetoo = modetoo_coagpair(ip)

!        call getcoags_wrapper_f(              &
!           airtemp, airprs,                        &
!           dgatk, dgacc,                           &
!           sgatk, sgacc,                           &
!           xxlsgat, xxlsgac,                       &
!           pdensat, pdensac,                       &
!           betaij0, betaij2i, betaij2j, betaij3,   &
!           betaii0, betaii2, betajj0, betajj2      )
         call getcoags_wrapper_f(                                       &
            temp,                        pmid,                          &
            dgn_awet(modefrm),           dgn_awet(modetoo),             &
            sigmag_aer(modefrm),         sigmag_aer(modetoo),           &
            alnsg_aer(modefrm),          alnsg_aer(modetoo),            &
            wetdens(modefrm),            wetdens(modetoo),              &
            ybetaij0(ip), xbetaij2i,     xbetaij2j,     ybetaij3(ip),   &
            ybetaii0(ip), xbetaii2,      ybetajj0(ip),  xbetajj2        )

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
!   short diagnostics for coag coefficients
         if ( ldiag15n ) write(lun15n,'(a,i5,1p,10e11.3)') 'ip, ybeta       ', ip, &
            ybetaij0(ip), ybetaij3(ip), ybetaii0(ip), ybetajj0(ip)

!   long diagnostics for coag coefficients
!        if ( ldiag15n ) then
!        if ( ip == 1 ) then
!        write(lun15n,'(a,1p,2e12.4)') 'temp, pmid', temp, pmid
!        write(lun15n,'(a/a/a/a/a)') &
!           'modefrm, modetoo, ip,                                  ', &
!           'ybetaij0(ip), ybetaij3(ip), ybetaii0(ip), ybetajj0(ip), [m3/s]', &
!           'dgn_awet(modefrm)*1.0e6,     dgn_awet(modetoo)*1.0e6,  ', &
!           'sigmag_aer(modefrm),         sigmag_aer(modetoo),      ', &
!           'wetdens(modefrm)*1.0e-3,     wetdens(modetoo)*1.0e-3   ' 
!        end if
!        write(lun15n,'(a,2i3,i5,1p,4e11.3, 0p,2x,2f7.4,2(2x,2f6.3))') &
!           'ip, ybeta ', modefrm, modetoo, ip, &
!           ybetaij0(ip), ybetaij3(ip), ybetaii0(ip), ybetajj0(ip), &
!           dgn_awet(modefrm)*1.0e6,     dgn_awet(modetoo)*1.0e6,   &
!           sigmag_aer(modefrm),         sigmag_aer(modetoo),       &
!           wetdens(modefrm)*1.0e-3,     wetdens(modetoo)*1.0e-3
!        end if
#endif

         ! convert coag coefficients from (m3/s) to (kmol-air/s)
         ybetaij0(ip) = ybetaij0(ip)*aircon
         ybetaij3(ip) = ybetaij3(ip)*aircon
         ybetaii0(ip) = ybetaii0(ip)*aircon
         ybetajj0(ip) = ybetajj0(ip)*aircon
      end do ! ip


! first calculate changes to number
! use the following order because
!    accum   number loss depends on accum number
!    pcarbon number loss depends on pcarbon and accum number
!    maccum  number loss depends on maccum, pcarbon, and accum number
!    aitken  number loss depends on aitken, maccum, pcarbon, and accum number
!    maitken number loss depends on maitken, aitken, maccum, pcarbon, and accum number
! the average number concencentrations (over current time step)
!    of other modes can thus be used to calculate the number loss of a mode

! accum mode number loss - analytical solution
      tmpa = max( 0.0_r8, deltat*ybetajj0(1) )
      qnum_tmpb(nacc) = qnum_tmpa(nacc) / &
         ( 1.0_r8 + ybetajj0(1)*deltat*qnum_tmpa(nacc) ) 
      qnum_tmpc(nacc) = (qnum_tmpa(nacc) + qnum_tmpb(nacc))*0.5_r8

! pcarbon mode number loss - approximate analytical solution
!    using average number conc. for accum mode
      if (npca > 0) then
         tmpa = max( 0.0_r8, deltat*ybetaij0(2)*qnum_tmpc(nacc) )
         tmpb = max( 0.0_r8, deltat*ybetaii0(2) )
         tmpn = qnum_tmpa(npca)
         if (tmpa < 1.0e-5_r8) then
            qnum_tmpb(npca) = tmpn / &
               ( 1.0_r8 + (tmpa+tmpb*tmpn)*(1.0_r8 + 0.5_r8*tmpa) )
         else
            tmpc = exp(-tmpa)
            qnum_tmpb(npca) = tmpn*tmpc / &
               ( 1.0_r8 + (tmpb*tmpn/tmpa)*(1.0_r8-tmpc) )
         end if
         qnum_tmpc(npca) = (qnum_tmpa(npca) + qnum_tmpb(npca))*0.5_r8
      end if

! marine-organics accum mode number loss - approximate analytical solution
!    using average number conc. for accum and pcarbon modes
      if (nmacc > 0) then
         tmpa = ybetaij0( 9)*qnum_tmpc(nacc) & 
              + ybetaij0(10)*qnum_tmpc(npca)
         tmpa = max( 0.0_r8, deltat*tmpa )
         tmpb = max( 0.0_r8, deltat*ybetaii0(9) )
         tmpn = qnum_tmpa(nmacc)
         if (tmpa < 1.0e-5_r8) then
            qnum_tmpb(nmacc) = tmpn / &
               ( 1.0_r8 + (tmpa+tmpb*tmpn)*(1.0_r8 + 0.5_r8*tmpa) )
         else
            tmpc = exp(-tmpa)
            qnum_tmpb(nmacc) = tmpn*tmpc / &
               ( 1.0_r8 + (tmpb*tmpn/tmpa)*(1.0_r8-tmpc) )
         end if
         qnum_tmpc(nmacc) = (qnum_tmpa(nmacc) + qnum_tmpb(nmacc))*0.5_r8
      end if

! aitken mode number loss - approximate analytical solution
!    using average number conc. for accum, pcarbon, and marine-org accum modes
      tmpa = ybetaij0(1)*qnum_tmpc(nacc)
      if (npca  > 0) tmpa = tmpa + ybetaij0(3)*qnum_tmpc(npca)
      if (nmacc > 0) tmpa = tmpa + ybetaij0(4)*qnum_tmpc(nmacc)
      tmpa = max( 0.0_r8, deltat*tmpa )
      tmpb = max( 0.0_r8, deltat*ybetaii0(1) )
      tmpn = qnum_tmpa(nait)
      if (tmpa < 1.0e-5_r8) then
         qnum_tmpb(nait) = tmpn / &
            ( 1.0_r8 + (tmpa+tmpb*tmpn)*(1.0_r8 + 0.5_r8*tmpa) )
      else
         tmpc = exp(-tmpa)
         qnum_tmpb(nait) = tmpn*tmpc / &
            ( 1.0_r8 + (tmpb*tmpn/tmpa)*(1.0_r8-tmpc) )
      end if
      qnum_tmpc(nait) = (qnum_tmpa(nait) + qnum_tmpb(nait))*0.5_r8

! marine-organics aitken mode number loss - approximate analytical solution
!    using average number conc. for accum, pcarbon, aitken, and marine-org accum modes
      if (nmait > 0) then
         tmpa = ybetaij0(5)*qnum_tmpc(nacc) + ybetaij0(7)*qnum_tmpc(nait) &
              + ybetaij0(6)*qnum_tmpc(npca) + ybetaij0(8)*qnum_tmpc(nmacc)
         tmpa = max( 0.0_r8, deltat*tmpa )
         tmpb = max( 0.0_r8, deltat*ybetaii0(5) )
         tmpn = qnum_tmpa(nmait)
         if (tmpa < 1.0e-5_r8) then
            qnum_tmpb(nmait) = tmpn / &
               ( 1.0_r8 + (tmpa+tmpb*tmpn)*(1.0_r8 + 0.5_r8*tmpa) )
         else
            tmpc = exp(-tmpa)
            qnum_tmpb(nmait) = tmpn*tmpc / &
               ( 1.0_r8 + (tmpb*tmpn/tmpa)*(1.0_r8-tmpc) )
         end if
         qnum_tmpc(nmait) = (qnum_tmpa(nmait) + qnum_tmpb(nmait))*0.5_r8
      end if


! now calculate mass transfers between modes
! the transfer amounts are calculated using as an exponential decay of
!     the initial mass concentrations, 
!     where the decay rate is calculated using the average (over time step) 
!     number concentrations for each mode
! the mass transfer calculations are first-order accurate in time, 
!     because the mass transferred out of a mode does not
!     include any mass transferred in during the time step
! with this approach, the ordering is not important, but the mass transfer
!     calculations are done in the reverse order of the number loss calculations

! mass transfer out of marine-organics aitken mode
!    uses average number conc. for accum, aitken, pcarbon, and marine-org accum modes
      if (nmait > 0) then
         tmp1 = max( 0.0_r8, ybetaij3(5)*qnum_tmpc(nacc) )
         tmp2 = max( 0.0_r8, ybetaij3(6)*qnum_tmpc(npca) )
         tmp3 = max( 0.0_r8, ybetaij3(7)*qnum_tmpc(nait) )
         tmp4 = max( 0.0_r8, ybetaij3(8)*qnum_tmpc(nmacc) )
         tmpa = tmp1 + tmp2 + tmp3 + tmp4
         tmpc = deltat*tmpa
         if (tmpc > epsilonx2) then
            ! calc coag change only when it is not ~= zero
            tmp_xf = 1.0_r8 - exp(-tmpc)
            tmp2 = tmp2/tmpa
            tmp3 = tmp3/tmpa
            tmp4 = tmp4/tmpa
            tmp1 = 1.0_r8 - (tmp2 + tmp3 + tmp4)
            do iaer = 1, naer
               tmp_dq = tmp_xf*qaer_tmpa(iaer,nmait)
               qaer_tmpb(iaer,nmait) = qaer_tmpb(iaer,nmait) - tmp_dq
               qaer_tmpb(iaer,nacc ) = qaer_tmpb(iaer,nacc ) + tmp_dq*tmp1
               qaer_tmpb(iaer,npca ) = qaer_tmpb(iaer,npca ) + tmp_dq*tmp2
               qaer_tmpb(iaer,nait ) = qaer_tmpb(iaer,nait ) + tmp_dq*tmp3
               qaer_tmpb(iaer,nmacc) = qaer_tmpb(iaer,nmacc) + tmp_dq*tmp4
               qaer_del_coag_in(iaer,i_agepair_pca ) &
                                    = qaer_del_coag_in(iaer,i_agepair_pca ) + tmp_dq*tmp2
               qaer_del_coag_in(iaer,i_agepair_macc) &
                                    = qaer_del_coag_in(iaer,i_agepair_macc) + tmp_dq*tmp4
            end do
         end if
      end if

! (ip == 2)  modefrm = npca  ; modetoo = nacc
! (ip == 1)  modefrm = nait  ; modetoo = nacc
! (ip == 3)  modefrm = nait  ; modetoo = npca
! (ip == 4)  modefrm = nait  ; modetoo = nmacc
! (ip == 5)  modefrm = nmait ; modetoo = nacc
! (ip == 6)  modefrm = nmait ; modetoo = npca
! (ip == 7)  modefrm = nmait ; modetoo = nait
! (ip == 8)  modefrm = nmait ; modetoo = nmacc
! (ip == 9)  modefrm = nmacc ; modetoo = nacc
! (ip ==10)  modefrm = nmacc ; modetoo = npca

! mass transfer out of aitken mode
!    uses average number conc. for accum, pcarbon, and marine-org accum modes
      tmp1 = max( 0.0_r8, ybetaij3(1)*qnum_tmpc(nacc) )
      if (nmacc > 0) then
         tmp2 = max( 0.0_r8, ybetaij3(3)*qnum_tmpc(npca) )
         tmp3 = max( 0.0_r8, ybetaij3(4)*qnum_tmpc(nmacc) )
      else if (npca > 0) then
         tmp2 = max( 0.0_r8, ybetaij3(3)*qnum_tmpc(npca) )
         tmp3 = 0.0_r8
      else
         tmp2 = 0.0_r8
         tmp3 = 0.0_r8
      end if
      tmpa = tmp1 + tmp2 + tmp3
      tmpc = deltat*tmpa
      if (tmpc > epsilonx2) then
         ! calc coag change only when it is not ~= zero
         tmp_xf = 1.0_r8 - exp(-tmpc)
         if (nmacc > 0) then
            tmp2 = tmp2/tmpa
            tmp3 = tmp3/tmpa
            tmp1 = 1.0_r8 - (tmp2 + tmp3)
            do iaer = 1, naer
               tmp_dq = tmp_xf*qaer_tmpa(iaer,nait)
               qaer_tmpb(iaer,nait)  = qaer_tmpb(iaer,nait)  - tmp_dq
               qaer_tmpb(iaer,nacc)  = qaer_tmpb(iaer,nacc)  + tmp_dq*tmp1
               qaer_tmpb(iaer,npca)  = qaer_tmpb(iaer,npca)  + tmp_dq*tmp2
               qaer_tmpb(iaer,nmacc) = qaer_tmpb(iaer,nmacc) + tmp_dq*tmp3
               qaer_del_coag_in(iaer,i_agepair_pca) &
                                    = qaer_del_coag_in(iaer,i_agepair_pca) + tmp_dq*tmp2
               qaer_del_coag_in(iaer,i_agepair_macc) &
                                    = qaer_del_coag_in(iaer,i_agepair_macc) + tmp_dq*tmp3
            end do
         else if (npca > 0) then
            tmp2 = tmp2/tmpa
            tmp1 = 1.0_r8 - tmp2
            do iaer = 1, naer
               tmp_dq = tmp_xf*qaer_tmpa(iaer,nait)
               qaer_tmpb(iaer,nait) = qaer_tmpb(iaer,nait) - tmp_dq
               qaer_tmpb(iaer,nacc) = qaer_tmpb(iaer,nacc) + tmp_dq*tmp1
               qaer_tmpb(iaer,npca) = qaer_tmpb(iaer,npca) + tmp_dq*tmp2
               qaer_del_coag_in(iaer,i_agepair_pca) &
                                    = qaer_del_coag_in(iaer,i_agepair_pca) + tmp_dq*tmp2
            end do
         else
            do iaer = 1, naer
               tmp_dq = tmp_xf*qaer_tmpa(iaer,nait)
               qaer_tmpb(iaer,nait) = qaer_tmpb(iaer,nait) - tmp_dq
               qaer_tmpb(iaer,nacc) = qaer_tmpb(iaer,nacc) + tmp_dq
            end do
         end if
      end if

!! old version for 3 and 7 mode only
!! mass transfer out of aitken mode mass
!!    uses average number conc. for accum and pcarbon modes
!      tmpa = max( 0.0_r8, ybetaij3(1)*qnum_tmpc(nacc) )
!      if (npca > 0) then
!         tmpb = max( 0.0_r8, ybetaij3(3)*qnum_tmpc(npca) )
!         tmpc = tmpa + tmpb
!      else
!         tmpc = tmpa
!      end if
!      tmpc = deltat*tmpc
!      if (tmpc > epsilonx2) then
!         ! calc coag change only when it is not ~= zero
!         tmp_xf = 1.0_r8 - exp(-tmpc)
!         if (npca > 0) then
!            tmp2 = tmpb/(tmpa + tmpb + epsilonx1)
!            tmp1 = 1.0_r8 - tmp2
!            do iaer = 1, naer
!               tmp_dq = tmp_xf*qaer_tmpa(iaer,nait)
!               qaer_tmpb(iaer,nait) = qaer_tmpb(iaer,nait) - tmp_dq
!               qaer_tmpb(iaer,nacc) = qaer_tmpb(iaer,nacc) + tmp_dq*tmp1
!               qaer_tmpb(iaer,npca) = qaer_tmpb(iaer,npca) + tmp_dq*tmp2
!               qaer_del_coag_in(iaer,i_agepair_pca) &
!                                    = qaer_del_coag_in(iaer,i_agepair_pca) + tmp_dq*tmp2
!            end do
!         else
!            do iaer = 1, naer
!               tmp_dq = tmp_xf*qaer_tmpa(iaer,nait)
!               qaer_tmpb(iaer,nait) = qaer_tmpb(iaer,nait) - tmp_dq
!               qaer_tmpb(iaer,nacc) = qaer_tmpb(iaer,nacc) + tmp_dq
!            end do
!         end if
!      end if

! mass transfer out of marine-organics accum mode
!    uses average number conc. for accum and pcarbon modes
      if (nmacc > 0) then
         tmp1 = max( 0.0_r8, ybetaij3( 9)*qnum_tmpc(nacc) )
         tmp2 = max( 0.0_r8, ybetaij3(10)*qnum_tmpc(npca) )
         tmpa = tmp1 + tmp2
         tmpc = deltat*tmpa
         if (tmpc > epsilonx2) then
            ! calc coag change only when it is not ~= zero
            tmp_xf = 1.0_r8 - exp(-tmpc)
            tmp2 = tmp2/tmpa
            tmp1 = 1.0_r8 - tmp2
            do iaer = 1, naer
               tmp_dq = tmp_xf*qaer_tmpa(iaer,nmacc)
               qaer_tmpb(iaer,nmacc) = qaer_tmpb(iaer,nmacc) - tmp_dq
               qaer_tmpb(iaer,nacc ) = qaer_tmpb(iaer,nacc ) + tmp_dq*tmp1
               qaer_tmpb(iaer,npca ) = qaer_tmpb(iaer,npca ) + tmp_dq*tmp2
               qaer_del_coag_in(iaer,i_agepair_pca ) &
                                    = qaer_del_coag_in(iaer,i_agepair_pca ) + tmp_dq*tmp2
            end do
         end if
      end if

! mass transfer out of pcarbon mode 
!    uses average number conc. for accum mode
      if (npca > 0) then
         tmpc = max( 0.0_r8, ybetaij3(2)*qnum_tmpc(nacc) )
         tmpc = deltat*tmpc
         if (tmpc > epsilonx2) then
            tmp_xf = 1.0_r8 - exp(-tmpc)
            do iaer = 1, naer
               tmp_dq = tmp_xf*qaer_tmpa(iaer,npca)
               qaer_tmpb(iaer,npca) = qaer_tmpb(iaer,npca) - tmp_dq
               qaer_tmpb(iaer,nacc) = qaer_tmpb(iaer,nacc) + tmp_dq
            end do
         end if
      end if

! mass transfer out of accum mode - there is no transfer out of this mode

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( ldiag15n ) then
! diagnostics
      do n = 1, 3
         write(lun15n,'(a,i5,1p,10e11.3)') 'n, qnum_tmpa/b/c', n, &
         qnum_tmpa(n), qnum_tmpc(n), qnum_tmpb(n), &
         qnum_tmpb(n)-qnum_tmpa(n), &
         1.0_r8-qnum_tmpb(n)/max(1.0e-5_r8,qnum_tmpa(n))
      end do
      do n = 1, 3
         write(lun15n,'(a,i5,1p,10e11.3)') 'n, dgnd/w, densw', n, &
         dgn_a(n), dgn_awet(n), wetdens(n)
      end do
      end if
#endif

      qnum_cur = qnum_tmpb
      qaer_cur = qaer_tmpb

      return
      end subroutine mam_coag_1subarea


!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine mam_pcarbon_aging_1subarea(                        &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         deltat,            dgn_a,            do_cond,              &
         n_mode,                                                    &
         qnum_cur,          qnum_del_cond,    qnum_del_coag,        &
         qaer_cur,          qaer_del_cond,    qaer_del_coag,        &
         qaer_del_coag_in,                                          &
         qwtr_cur                                                   )

! uses

      implicit none

! arguments
      integer,  intent(in) :: nstep                 ! model time-step number
      integer,  intent(in) :: lchnk                 ! chunk identifier
      integer,  intent(in) :: i, k                  ! column and level indices
      integer,  intent(in) :: jsub                  ! sub-area index
      integer,  intent(in) :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in) :: lund                  ! logical unit for diagnostic output
      integer,  intent(in) :: n_mode                ! current number of modes (including temporary)

      logical,  intent(in) :: do_cond               ! true if condensation (gas-aerosol exch) is on

      real(r8), intent(in) :: deltat                ! model timestep (s)
      real(r8), intent(in), dimension( 1:max_mode ) :: &
                              dgn_a                 ! dgnum_dry of mode

      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur, qnum_del_cond, qnum_del_coag
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur, qaer_del_cond, qaer_del_coag
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_agepair ) :: &
         qaer_del_coag_in
! *** need to add qaer_del_coag_nmoacc_in
      real(r8), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur

! local variables
      integer :: iaer, ipair, itmpa
      integer :: nfrm, ntoo

      real(r8) :: fac_volsfc
      real(r8) :: tmpa, tmp1, tmp2, tmp3, tmp4
      real(r8) :: vol_core, vol_shell
      real(r8) :: xferfrac_max, xferfrac_pcage


! 
agepair_loop1: &
      do ipair = 1, n_agepair

      nfrm = modefrm_agepair(ipair)
      ntoo = modetoo_agepair(ipair)

      vol_shell = qaer_cur(iaer_so4,nfrm)*fac_m2v_aer(iaer_so4)
      tmp3 = qaer_del_cond(iaer_so4,nfrm)  *fac_m2v_aer(iaer_so4)
      tmp4 = qaer_del_coag_in(iaer_so4,ipair)*fac_m2v_aer(iaer_so4)

      do iaer = 1, naer
!     species that contribute to aging are 
!        so4 (but it is already done)
!        soa, nh4 and no3
!        ncl and cl (when aging_include_seasalt == .true.)
         if ( (iaer <= nsoa    ) .or. &
              (iaer == iaer_nh4) .or. &
              (iaer == iaer_no3) .or. &
              (iaer == iaer_cl ) ) then
            continue
         else if (iaer == iaer_ncl) then
            if (aging_include_seasalt .eqv. .false.) cycle
         else
            cycle
         end if

         if ( (iaer == iaer_cl                    ) .and. &
              (aging_include_seasalt .eqv. .false.) ) then
            ! special case - only include the cl from condensation
            tmp1 = max( qaer_del_cond(iaer,nfrm), 0.0_r8 )
            tmp2 = max( qaer_del_coag_in(iaer,ipair), 0.0_r8 ) + tmp1
            if (tmp2 >= 1.0e-35_r8) then
               vol_shell = vol_shell + qaer_cur(iaer,nfrm)*fac_m2v_eqvhyg_aer(iaer)*(tmp1/tmp2)
               tmp3 = tmp3 + qaer_del_cond(iaer,nfrm)*fac_m2v_eqvhyg_aer(iaer)
            end if
         else
            vol_shell = vol_shell + qaer_cur(iaer,nfrm)*fac_m2v_eqvhyg_aer(iaer)
            tmp3 = tmp3 + qaer_del_cond(iaer,nfrm)*fac_m2v_eqvhyg_aer(iaer)
            tmp4 = tmp4 + qaer_del_coag_in(iaer,ipair)*fac_m2v_eqvhyg_aer(iaer)
         end if
      end do

      if ( do_cond ) then
         tmp3 = max( tmp3, 1.0e-35_r8 )
         tmp3 = tmp3/(tmp3 + max( tmp4, 0.0_r8 ))
      else
         tmp3 = 0.0_r8
      end if
      tmp4 = 1.0_r8 - tmp3

      vol_core = 0.0
      do iaer = 1, naer
         ! for core volume, only include the mapped species 
         !    which are primary and low hygroscopicity
         if (lmap_aer(iaer,nfrm) > 0) &
         vol_core = vol_core + qaer_cur(iaer,nfrm)*fac_m2v_aer(iaer)
      end do
!   ratio1 = vol_shell/vol_core = 
!      actual hygroscopic-shell-volume/carbon-core-volume after gas uptake
!   ratio2 = 6.0_r8*dr_so4_monolayers_pcage/(dgncur_a*fac_volsfc)
!      = (shell-volume corresponding to n_so4_monolayers_pcage)/core-volume 
!      The 6.0/(dgncur_a*fac_volsfc) = (mode-surface-area/mode-volume)
!   Note that vol_shell includes both so4+nh4 AND soa as "equivalent so4",
!      The soa_equivso4_factor accounts for the lower hygroscopicity of soa.
!
!   Define xferfrac_pcage = min( 1.0, ratio1/ratio2)
!   But ratio1/ratio2 == tmp1/tmp2, and coding below avoids possible overflow 
!
      fac_volsfc = exp( 2.5*(alnsg_aer(nfrm)**2) )
      xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

      tmp1 = vol_shell*dgn_a(nfrm)*fac_volsfc
      tmp2 = max( 6.0_r8*dr_so4_monolayers_pcage*vol_core, 0.0_r8 )
      if (tmp1 >= tmp2) then
         xferfrac_pcage = xferfrac_max
      else
         xferfrac_pcage = min( tmp1/tmp2, xferfrac_max )
      end if

      do iaer = 1, naer
         if (lmap_aer(iaer,nfrm) > 0) then
            ! species is pom or bc
            ! transfer the aged fraction to accum mode
            ! include this transfer change in the cond and/or coag change (for mass budget)
            tmpa = qaer_cur(iaer,nfrm)*xferfrac_pcage
            qaer_cur(iaer,nfrm)      = qaer_cur(iaer,nfrm) - tmpa
            qaer_cur(iaer,ntoo)      = qaer_cur(iaer,ntoo) + tmpa
            qaer_del_cond(iaer,nfrm) = qaer_del_cond(iaer,nfrm) - tmpa*tmp3
            qaer_del_cond(iaer,ntoo) = qaer_del_cond(iaer,ntoo) + tmpa*tmp3
            qaer_del_coag(iaer,nfrm) = qaer_del_coag(iaer,nfrm) - tmpa*tmp4
            qaer_del_coag(iaer,ntoo) = qaer_del_coag(iaer,ntoo) + tmpa*tmp4
         else
            ! species is soa, so4, or nh4 produced by condensation or coagulation
            ! transfer all of it to accum mode
            ! also transfer the condensation and coagulation changes
            !    to accum mode (for mass budget)
            qaer_cur(iaer,ntoo)      = qaer_cur(iaer,ntoo) &
                                     + qaer_cur(iaer,nfrm)
            qaer_del_cond(iaer,ntoo) = qaer_del_cond(iaer,ntoo) &
                                     + qaer_del_cond(iaer,nfrm)
            qaer_del_coag(iaer,ntoo) = qaer_del_coag(iaer,ntoo) &
                                     + qaer_del_coag(iaer,nfrm)
            qaer_cur(iaer,nfrm)      = 0.0_r8
            qaer_del_cond(iaer,nfrm) = 0.0_r8
            qaer_del_coag(iaer,nfrm) = 0.0_r8
         end if
      end do
      ! number - transfer the aged fraction to accum mode
      ! include this transfer change in the cond and/or coag change (for mass budget)
      tmpa = qnum_cur(nfrm)*xferfrac_pcage
      qnum_cur(nfrm)      = qnum_cur(nfrm) - tmpa
      qnum_cur(ntoo)      = qnum_cur(ntoo) + tmpa
      qnum_del_cond(nfrm) = qnum_del_cond(nfrm) - tmpa*tmp3
      qnum_del_cond(ntoo) = qnum_del_cond(ntoo) + tmpa*tmp3
      qnum_del_coag(nfrm) = qnum_del_coag(nfrm) - tmpa*tmp4
      qnum_del_coag(ntoo) = qnum_del_coag(ntoo) + tmpa*tmp4

      end do agepair_loop1

      return
      end subroutine mam_pcarbon_aging_1subarea


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      function mean_molecular_speed( temp, rmw )
      implicit none
      real(8) :: mean_molecular_speed  ! (m/s)
      real(8) :: temp                  ! temperature (K)
      real(8) :: rmw                   ! molec. weight (g/mol)
      mean_molecular_speed = 145.5_8 * sqrt(temp/rmw)
      return
      end function mean_molecular_speed


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      function gas_diffusivity( t_k, p_atm, rmw, vm )
      implicit none
      real(8) :: gas_diffusivity       ! (m2/s)
      real(8) :: t_k                   ! temperature (K)
      real(8) :: p_atm                 ! pressure (atmospheres)
      real(8) :: rmw                   ! molec. weight (g/mol)
      real(8) :: vm                    ! molar volume (units = ??)

      real(8) :: dgas

      dgas = (1.0e-3_8 * t_k**1.75_8 * sqrt(1./rmw + 0.035_8))/   &
             (p_atm * (vm**0.3333333333333333_8 + 2.7189_8)**2)
      gas_diffusivity = dgas*1.0e-4_8
      return
      end function gas_diffusivity


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine gas_aer_uptkrates_1box1gas( &
         accom, gasdiffus, gasfreepath, &
         beta_inp, n_mode, dgncur_awet, lnsg, uptkrate )
!
!                         /
!   computes   uptkrate = | dx  dN/dx  gas_conden_rate(Dp(x))
!                         /
!   using Gauss-Hermite quadrature of order nghq=2
!
!       Dp = particle diameter (cm)
!       x = ln(Dp)
!       dN/dx = log-normal particle number density distribution
!       gas_conden_rate(Dp) = 2 * pi * gasdiffus * Dp * F(Kn,ac)
!           F(Kn,ac) = Fuchs-Sutugin correction factor
!           Kn = Knudsen number
!           ac = accomodation coefficient
!
      implicit none

      integer, parameter :: r8 = 8

      integer,  intent(in)  :: n_mode                ! number of modes

      real(r8), intent(in)  :: accom                ! accomodation coefficient (--)
      real(r8), intent(in)  :: gasdiffus            ! gas diffusivity (m2/s)
      real(r8), intent(in)  :: gasfreepath          ! gas mean free path (m)
      real(r8), intent(in)  :: beta_inp             ! quadrature parameter (--)
      real(r8), intent(in)  :: dgncur_awet(n_mode)
                               ! mode-median wet diameter of number distribution (m)
      real(r8), intent(in)  :: lnsg(n_mode)
                               ! ln( sigmag )  (--)
      real(r8), intent(out) :: uptkrate(n_mode)  
                               ! gas-to-aerosol mass transfer rates (1/s)
                               ! for number concentration = 1 #/m3


! local 
      integer, parameter :: nghq = 2
      integer :: i, iq, k, l1, l2, la, n

      real(r8), parameter :: tworootpi = 3.5449077018110320_r8
      real(r8), parameter :: root2 =     1.4142135623730950_r8
      real(r8), parameter :: one = 1.0_r8
      real(r8), parameter :: two = 2.0_r8

      real(r8) :: accomxp283, accomxp75
      real(r8) :: beta
      real(r8) :: const
      real(r8) :: dp, dum_m2v
      real(r8) :: fuchs_sutugin
      real(r8) :: knudsen
      real(r8) :: lndp, lndpgn
      real(r8) :: sumghq
      real(r8) :: tmpa
      real(r8), save :: xghq(nghq), wghq(nghq) ! quadrature abscissae and weights

      data xghq / 0.70710678, -0.70710678 /
      data wghq / 0.88622693,  0.88622693 /


      accomxp283 = accom * 0.283_r8
      accomxp75  = accom * 0.75_r8

! outermost loop over all modes
      do n = 1, n_mode

         lndpgn = log( dgncur_awet(n) )   ! (m)

! beta = dln(uptake_rate)/dln(dp)
!      = 2.0 in free molecular regime, 1.0 in continuum regime
! if uptake_rate ~= a * (dp**beta), then the 2 point quadrature is very accurate
         if (abs(beta_inp-1.5_r8) > 0.5_r8) then
!           dp = dgncur_awet(n) * exp( 1.5_r8*(lnsg(n)**2) )
            dp = dgncur_awet(n)
            knudsen = two*gasfreepath/dp
            ! tmpa = dln(fuchs_sutugin)/d(knudsen)
            tmpa = one/(one+knudsen) - (two*knudsen + one + accomxp283) / &
                      ( knudsen*( knudsen + one + accomxp283 ) + accomxp75 )
            beta = one - knudsen*tmpa
            beta = max( one, min( two, beta ) )
         else
            beta = beta_inp
         end if

         const  = tworootpi * exp( beta*lndpgn + 0.5_r8*(beta*lnsg(n))**2 )
         
!   sum over gauss-hermite quadrature points
         sumghq = 0.0
         do iq = 1, nghq
            lndp = lndpgn + beta*lnsg(n)**2 + root2*lnsg(n)*xghq(iq)
            dp = exp(lndp)

            knudsen = two*gasfreepath/dp

!           fkn = ( 0.75*accomcoef*(1. + xkn) ) / &
!                 ( xkn**2 + xkn + 0.283*xkn*accomcoef + 0.75*accomcoef )
            fuchs_sutugin = &
                  ( accomxp75*(one + knudsen) ) / &
                  ( knudsen*( knudsen + one + accomxp283 ) + accomxp75 )

            sumghq = sumghq + wghq(iq)*dp*fuchs_sutugin/(dp**beta)
         end do
         uptkrate(n) = const * gasdiffus * sumghq    

      end do   ! "do n = 1, ntot_soamode"


      return
      end subroutine gas_aer_uptkrates_1box1gas



!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine modal_aero_amicphys_init( imozart, species_class,n_so4_monolayers_pcage_in)

!-----------------------------------------------------------------------
!
! Purpose:
!    set do_adjust and do_aitken flags
!    create history fields for column tendencies associated with
!       modal_aero_calcsize
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

use cam_history, only  :  fieldname_len
use cam_logfile, only  :  iulog
use chem_mods, only    :  adv_mass
use constituents, only :  pcnst, cnst_get_ind, cnst_name
use mo_chem_utls, only :  get_spc_ndx
use mo_tracname,  only :  solsym
use physconst, only    :  mwdry, mwh2o
use spmd_utils, only   :  masterproc
use phys_control,only  :  phys_getopts

use modal_aero_data, only : &
    cnst_name_cw, &
    dgnum_amode, dgnumlo_amode, dgnumhi_amode, &
    lmassptr_amode, lmassptrcw_amode, &
    modeptr_accum, modeptr_aitken, modeptr_pcarbon, modeptr_ufine, &
    modeptr_maccum, modeptr_maitken, &
    nspec_amode, &
    numptr_amode, numptrcw_amode, sigmag_amode

implicit none

!-----------------------------------------------------------------------
! arguments
   integer, intent(in)  :: imozart
   integer, intent(in)  :: species_class(:)
   real(r8), intent(in) :: n_so4_monolayers_pcage_in

!-----------------------------------------------------------------------
! local
   integer, parameter :: big_neg_int = -999888777
   integer, parameter :: big_pos_int =  999888777
   integer  :: iaer, igas, ip, ipair, itmpa
   integer  :: j, jac, jsoa
   integer  :: l, l1, l2, lac
   integer  :: lmz, lmz2, loffset
   integer  :: l_so4g, l_nh4g, l_msag
   integer  :: m
   integer  :: n, na, nb, nc
   integer  :: nspec

   real(r8) :: tmp1, tmp2

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(128)                 :: msg, fmtaa
   character(2)                   :: tmpch2
   !-----------------------------------------------------------------------
 
!namelist variables
n_so4_monolayers_pcage  = n_so4_monolayers_pcage_in
dr_so4_monolayers_pcage = n_so4_monolayers_pcage * 4.76e-10


#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      ldiag82  = .true.  ; lun82  = 82
      ldiag97  = .true.  ; lun97  = 97
      ldiag98  = .true.  ; lun98  = 98
      ldiag13n = .true.  ; lun13n = 130
      ldiag15n = .true.  ; lun15n = 150
      ldiagd1  = .true. 
#else
      ldiag82  = .false. ; lun82  = iulog
      ldiag97  = .false. ; lun97  = iulog
      ldiag98  = .false. ; lun98  = iulog
      ldiag13n = .false. ; lun13n = iulog
      ldiag15n = .false. ; lun15n = iulog
      ldiagd1  = .false.
#endif


      call mam_set_lptr2_and_specxxx2


      mwuse_soa(:) = 150.0_r8
      mwuse_poa(:) = 150.0_r8

! set ngas, name_gas, and igas_xxx
! set naer, name_aerpfx, and iaer_xxx
      name_gas    = "???"
      name_aerpfx = "???"
      name_aer    = "???"
      name_aercw  = "???"
      name_num    = "???"
      name_numcw  = "???"

      igas_h2so4 = 0 ; igas_nh3 = 0
      iaer_bc  = 0 ; iaer_dst = 0 
      iaer_ncl = 0 ; iaer_nh4 = 0 
      iaer_pom = 0 ; iaer_soa = 0 
      iaer_so4 = 0
      iaer_no3 = 0 ; iaer_cl  = 0 
      iaer_ca  = 0 ; iaer_co3 = 0 
      iaer_mpoly = 0 ; iaer_mprot = 0 
      iaer_mlip  = 0 ; iaer_mhum = 0 
      iaer_mproc = 0 ; iaer_mom = 0

      if (nsoa == 1) then
         name_gas(1) = 'SOAG'
         name_aerpfx(1) = 'soa'
      else if (nsoa == 2) then
         jsoa =      1 ; name_gas(jsoa) = 'SOAGa'  ; name_aerpfx(jsoa) = 'soaa'  ! jsoa=1
         jsoa = jsoa+1 ; name_gas(jsoa) = 'SOAGb'  ; name_aerpfx(jsoa) = 'soab'  ! jsoa=2
      else if (nsoa == 6) then
         jsoa =      1 ; name_gas(jsoa) = 'SOAGa1' ; name_aerpfx(jsoa) = 'soaa1' ! jsoa=1
         jsoa = jsoa+1 ; name_gas(jsoa) = 'SOAGa2' ; name_aerpfx(jsoa) = 'soaa2' ! jsoa=2
         jsoa = jsoa+1 ; name_gas(jsoa) = 'SOAGa3' ; name_aerpfx(jsoa) = 'soaa3' ! jsoa=3
         jsoa = jsoa+1 ; name_gas(jsoa) = 'SOAGb1' ; name_aerpfx(jsoa) = 'soab1' ! jsoa=4
         jsoa = jsoa+1 ; name_gas(jsoa) = 'SOAGb2' ; name_aerpfx(jsoa) = 'soab2' ! jsoa=5
         jsoa = jsoa+1 ; name_gas(jsoa) = 'SOAGb3' ; name_aerpfx(jsoa) = 'soab3' ! jsoa=6
      else
         call endrun( 'modal_aero_amicphys_init ERROR - bad nsoa' )
      end if
      ngas = nsoa
      naer = nsoa
      igas_soa = 1
      iaer_soa = 1

      ngas = ngas + 1
      name_gas(ngas) = 'H2SO4'
      naer = naer + 1
      name_aerpfx(naer) = 'so4'
      igas_h2so4 = ngas
      iaer_so4 = naer

      if ( (ntot_amode==7) .or. &
           (ntot_amode==8) .or. &
           (ntot_amode==9) ) then
         ngas = ngas + 1
         name_gas(ngas) = 'NH3'
         naer = naer + 1
         name_aerpfx(naer) = 'nh4'
         igas_nh3 = ngas
         iaer_nh4 = naer
      end if

#if ( ( defined MODAL_AERO_7MODE ) && ( defined MOSAIC_SPECIES ) )
      ngas = ngas + 1
      name_gas(ngas) = 'HNO3'
      naer = naer + 1
      name_aerpfx(naer) = 'no3'
      igas_hno3 = ngas
      iaer_no3 = naer

      ngas = ngas + 1
      name_gas(ngas) = 'HCL'
      naer = naer + 1
      name_aerpfx(naer) = 'cl'
      igas_hcl = ngas
      iaer_cl = naer
#endif

      iaer_pom = naer + 1
      if (npoa == 1) then
         naer = naer + 1
         name_aerpfx(naer) = 'pom'
      else if (npoa == 2) then
         naer = naer + 1
         name_aerpfx(naer) = 'poma'
         naer = naer + 1
         name_aerpfx(naer) = 'pomb'
      else
         call endrun( 'modal_aero_amicphys_init ERROR - bad npoa' )
      end if

      iaer_bc = naer + 1
      if (nbc == 1) then
         naer = naer + 1
         name_aerpfx(naer) = 'bc'
      else if (nbc == 2) then
         naer = naer + 1
         name_aerpfx(naer) = 'bca'
         naer = naer + 1
         name_aerpfx(naer) = 'bcb'
      else
         call endrun( 'modal_aero_amicphys_init ERROR - bad nbc' )
      end if

      naer = naer + 1
      name_aerpfx(naer) = 'ncl'
      iaer_ncl = naer
      naer = naer + 1
      name_aerpfx(naer) = 'dst'
      iaer_dst = naer

#if ( ( defined MODAL_AERO_7MODE ) && ( defined MOSAIC_SPECIES ) )
      naer = naer + 1
      name_aerpfx(naer) = 'ca'
      iaer_ca = naer
      naer = naer + 1
      name_aerpfx(naer) = 'co3'
      iaer_co3 = naer
#endif

#if ( defined MODAL_AERO_4MODE_MOM )
      naer = naer + 1
      name_aerpfx(naer) = 'mom'
      iaer_mom = naer
#endif

      if (ntot_amode==9) then
         naer = naer + 1
         name_aerpfx(naer) = 'mpoly'
         iaer_mpoly = naer
         naer = naer + 1
         name_aerpfx(naer) = 'mprot'
         iaer_mprot = naer
         naer = naer + 1
         name_aerpfx(naer) = 'mlip'
         iaer_mlip = naer
         naer = naer + 1
         name_aerpfx(naer) = 'mhum'
         iaer_mhum = naer
         naer = naer + 1
         name_aerpfx(naer) = 'mproc'
         iaer_mproc = naer
      end if

      if ((ngas /= max_gas) .or. (naer /= max_aer)) then
         write(iulog,'(a,4i10)') 'ngas, max_gas, naer, max_aer', &
            ngas, max_gas, naer, max_aer
         call endrun( 'modal_aero_amicphys_init ERROR - bad ngas or naer' )
      end if

      lmapcc_all(:) = 0

! set gas mapping
      loffset = imozart - 1
      lmap_gas(:) = 0
      mwhost_gas(:) = 1.0_r8
      mw_gas(:) = 1.0_r8
      fcvt_gas(:) = 1.0_r8

      vol_molar_gas = 42.88_r8  ! value for h2so4
      accom_coef_gas = 0.65_r8  ! value for h2so4

      do igas = 1, ngas
         call cnst_get_ind( name_gas(igas), l, .false. )
         if (l < 1 .or. l > pcnst) then
            msg = 'modal_aero_amicphys_init ERROR - lmap_gas for ' // name_gas(igas)
            call endrun( msg )
         end if
         lmz = l - loffset 
         lmz2 = get_spc_ndx( name_gas(igas) )
         if (lmz /= lmz2 .or. lmz <= 0) then
            msg = 'modal_aero_amicphys_init ERROR - lmz /= lmz2 for ' // name_gas(igas)
            call endrun( msg )
         end if
         lmapcc_all(lmz) = lmapcc_val_gas
         lmap_gas(igas) = lmz

         mwhost_gas(igas) = adv_mass(lmz)
         mw_gas(igas) = mwhost_gas(igas)
         if (igas <= nsoa) mw_gas(igas) = mwuse_soa(igas)
         fcvt_gas(igas) = mwhost_gas(igas)/mw_gas(igas)

         if (igas <= nsoa) then
            vol_molar_gas(igas) = vol_molar_gas(igas_h2so4) * (mw_gas(igas)/98.0_r8)
         else if (igas == igas_nh3) then
            vol_molar_gas(igas) = 14.90_r8
         else if (igas == igas_hno3) then
            vol_molar_gas(igas) = 24.11_r8
         else if (igas == igas_hcl) then
            vol_molar_gas(igas) = 21.48_r8
         end if
! values from mosaic code
!           v_molar(iv)= 42.88_r8  ! h2so4
!           v_molar(iv)= 24.11_r8  ! hno3
!           v_molar(iv)= 21.48_r8  ! hcl
!           v_molar(iv)= 14.90_r8  ! nh3
!           v_molar(iv)= 65.0_r8   ! soa

      end do ! igas

! set aerosol mass and number mapping
      lmap_aer(:,:) = 0
      lmap_aercw(:,:) = 0
      lmapbb_aer(:,:) = 0
      dens_aer(:) = 1.0_r8
      hygro_aer(:) = 1.0_r8
      mw_aer(:) = 1.0_r8
      fcvt_aer(:) = 1.0_r8

      lmap_num(:) = 0
      lmap_numcw(:) = 0
      fcvt_num = 1.0_r8  ! leave number mix-ratios unchanged (#/kmol-air)

      fcvt_wtr = mwdry/mwh2o  ! convert aerosol water mix-ratios from (kg/kg) to (mol/mol)

      do n = 1, ntot_amode
         do lac = 1, 2
         do l1 = 1, nspec_amode(n)
            if (lac == 1) then
               l = lmassptr_amode(l1,n)
               lmz = l - loffset 
               tmpnamea = cnst_name(l)
               tmpch2 = '_a'
            else
               l = lmassptrcw_amode(l1,n)
               lmz = 0
               tmpnamea = cnst_name_cw(l)
               tmpch2 = '_c'
            end if
            iaer = 0
            do j = 1, naer
               if (n <= 9) then
                  write(tmpnameb,'(2a,i1)') trim(name_aerpfx(j)), tmpch2, n
               else
                  write(tmpnameb,'(2a,i2)') trim(name_aerpfx(j)), tmpch2, n
               end if
               if (tmpnamea == tmpnameb) then
                  iaer = j
                  exit
               end if
            end do
            if (iaer <= 0) then
               msg = 'modal_aero_amicphys_init ERROR - lmap_aer for ' // tmpnamea
               call endrun( msg )
            end if
            if (lac == 1) then
               name_aer(iaer,n) = tmpnamea
               lmz2 = get_spc_ndx( tmpnamea )
               if (lmz /= lmz2 .or. lmz <= 0) then
                  msg = 'modal_aero_amicphys_init ERROR - lmz /= lmz2 for ' // tmpnamea
                  call endrun( msg )
               end if
               lmapcc_all(lmz) = lmapcc_val_aer
               lmap_aer(iaer,n) = l - loffset
               lmapbb_aer(iaer,n) = l1

               dens_aer(iaer) = specdens2_amode(l1,n)
               hygro_aer(iaer) = spechygro2(l1,n)
               mwhost_aer(iaer) = specmw2_amode(l1,n)
               mw_aer(iaer) = mwhost_aer(iaer)

               itmpa = iaer - iaer_pom + 1
               if (iaer <= nsoa) then
                  mw_aer(iaer) = mwuse_soa(iaer)
               else if ((1 <= itmpa) .and. (itmpa <= npoa)) then
                  mw_aer(iaer) = mwuse_poa(itmpa)
               end if
               fcvt_aer(iaer) = mwhost_aer(iaer)/mw_aer(iaer)
               fac_m2v_aer(iaer) = mw_aer(iaer)/dens_aer(iaer)
            else
               name_aercw(iaer,n) = tmpnamea
               lmap_aercw(iaer,n) = l - loffset
            end if

         end do ! l1
         end do ! lac

         lmap_num(n) = numptr_amode(n) - loffset
         name_num(n) = cnst_name(numptr_amode(n))
         lmz = lmap_num(n)
         lmz2 = get_spc_ndx( name_num(n) )
         if (lmz /= lmz2 .or. lmz <= 0) then
            msg = 'modal_aero_amicphys_init ERROR - lmz /= lmz2 for ' // name_num(n)
            call endrun( msg )
         end if
         lmapcc_all(lmz) = lmapcc_val_num

         lmap_numcw(n) = numptrcw_amode(n) - loffset
         name_numcw(n) = cnst_name_cw(numptrcw_amode(n))
         mwhost_num = 1.0_r8

      end do ! n

      do iaer = 1, naer
         fac_eqvso4hyg_aer(iaer) = hygro_aer(iaer)/hygro_aer(iaer_so4)
         fac_m2v_eqvhyg_aer(iaer) = fac_m2v_aer(iaer) * fac_eqvso4hyg_aer(iaer)
      end do ! naer

      sigmag_aer(:) = 1.8_r8
      sigmag_aer(1:ntot_amode) = sigmag_amode(1:ntot_amode)
      alnsg_aer(1:max_mode) = log(sigmag_aer(1:max_mode))

      dgnum_aer(:)   =  3.0e-9_r8
      dgnumlo_aer(:) =  1.0e-9_r8
      dgnumhi_aer(:) = 10.0e-9_r8
      dgnum_aer(1:ntot_amode)   = dgnum_amode(1:ntot_amode)
      dgnumhi_aer(1:ntot_amode) = dgnumhi_amode(1:ntot_amode)
      dgnumlo_aer(1:ntot_amode) = dgnumlo_amode(1:ntot_amode)

      ! converts number geometric_mean diameter to volume-mean diameter
      fcvt_dgnum_dvolmean(1:max_mode) = exp( 1.5_r8*(alnsg_aer(1:max_mode)**2) )

      dens_so4a_host = dens_aer(iaer_so4)
      mw_so4a_host = mwhost_aer(iaer_so4)
      if (iaer_nh4 > 0) then
         mw_nh4a_host = mwhost_aer(iaer_nh4)
      else
         mw_nh4a_host = mw_so4a_host
      end if


      nacc = modeptr_accum
      nait = modeptr_aitken
      npca = modeptr_pcarbon
      nufi = modeptr_ufine
#if ( defined MODAL_AERO_9MODE )
      nmacc = modeptr_maccum
      nmait = modeptr_maitken
#else
      nmacc = big_neg_int
      nmait = big_neg_int
#endif
      if ( nufi <= 0 .and. &
           ntot_amode_extd > ntot_amode ) nufi = ntot_amode_extd

! aging pairs
      ipair = 0
      modefrm_agepair(:) = big_neg_int
      modetoo_agepair(:) = big_neg_int
      mode_aging_optaa(:) = 0
      i_agepair_pca = big_neg_int ; i_agepair_macc = big_neg_int ; i_agepair_mait = big_neg_int ;
      if (npca > 0 .and. nacc > 0) then
         ipair = ipair + 1
         modefrm_agepair(ipair) = npca
         modetoo_agepair(ipair) = nacc
         i_agepair_pca = ipair
         mode_aging_optaa(npca) = 1
      end if
      if (nmacc > 0 .and. nacc > 0) then
         ipair = ipair + 1
         modefrm_agepair(ipair) = nmacc
         modetoo_agepair(ipair) = nacc
         i_agepair_macc = ipair
         mode_aging_optaa(nmacc) = 1
      end if
      if (nmait > 0 .and. nait > 0) then
         ipair = ipair + 1
         modefrm_agepair(ipair) = nmait
         modetoo_agepair(ipair) = nait
         i_agepair_mait = ipair
         mode_aging_optaa(nmait) = 1
      end if
      n_agepair = ipair

! coagulation pairs
!
! mam version	modes involved in coagulation          # of coag pairs
! -----------   -----------------------------          ---------------
! 3 mode	accum, aitken                              1
! 4,7 mode	accum, aitken, pcarbon                     3
! 9 mode	accum, aitken, pcarbon, maccum, maitken   10
! (pcarbon = primary carbon)
! (maccum  = primary marine-organics accum)
! (maitken = primary marine-organics aitken)
!
! 9 mode -- 5 participating modes and 10 possible coagulation pairs
!    6 possible coagulation pairs involve a smaller and a larger sized mode
!       the resulting particle is placed in the larger-sized mode
!          aitken  + [ accum, pcarbon, maccum ]
!          maitken + [ accum, pcarbon, maccum ]
!    4 possible coagulation pairs involve similar sized modes
!       the resulting particle is placed in the mode that is aged 
!       or contains the largest number of species
!          pcarbon + accum   --> accum   (aged)
!          maitken + aitken  --> aitken  (aged)
!          maccum  + accum   --> accum   (aged)
!          maccum  + pcarbon --> pcarbon (largest number of species)
!    note that 2 of the coagulation pairs results in aging, so
!          aitken  + pcarbon --> pcarbon (temporary) --> accum
!          aitken  + maccum  --> maccum  (temporary) --> accum
!    each mode also has self-coagulation which only affects number
      ipair = 0
      modefrm_coagpair(:) = big_neg_int
      modetoo_coagpair(:) = big_neg_int
      modeend_coagpair(:) = big_neg_int
      do ip = 1, 11
         na = big_neg_int ; nb = big_neg_int
         nc = big_pos_int
         if (ip == 1) then
            na = nait ; nb = nacc
         else if (ip == 2) then
            na = npca ; nb = nacc
         else if (ip == 3) then
            na = nait ; nb = npca
            nc = nacc
         else if (ip == 4) then
            na = nait ; nb = nmacc
            nc = nacc
         else if (ip == 5) then
            na = nmait ; nb = nacc
         else if (ip == 6) then
            na = nmait ; nb = npca
         else if (ip == 7) then
            na = nmait ; nb = nait
         else if (ip == 8) then
            na = nmait ; nb = nmacc
         else if (ip == 9) then
            na = nmacc ; nb = nacc
         else if (ip == 10) then
            na = nmacc ; nb = npca
         end if
         if (nc == big_pos_int) nc = nb

         if (na < 1 .or. nb < 1 .or. nc < 1) cycle
         ipair = ipair + 1
         modefrm_coagpair(ipair) = na
         modetoo_coagpair(ipair) = nb
         modeend_coagpair(ipair) = nc
      end do
      n_coagpair = ipair


! diagnostics
      if ( masterproc ) then
         write(iulog,'(/a)') 'modal_aero_amicphys_init start'

         write(iulog,'(/a,i12)') &
            'mdo_gaexch_cldy_subarea    ', mdo_gaexch_cldy_subarea
         write(iulog,'( a,i12)') &
            'gaexch_h2so4_uptake_optaa  ', gaexch_h2so4_uptake_optaa
         write(iulog,'( a,i12)') &
            'newnuc_h2so4_conc_optaa    ', newnuc_h2so4_conc_optaa
         write(iulog,'( a,i12)') &
            'rename_method_optaa        ', rename_method_optaa
         write(iulog,'( a,1p,e12.4)') &
            'newnuc_adjust_factor_pbl   ', newnuc_adjust_factor_pbl

         write(iulog,'(/a56,10i5)') &
           'ngas, max_gas, naer, max_aer', &
            ngas, max_gas, naer, max_aer
         write(iulog,'(/a56,10i5)') &
           'nsoa, npoa, nbc', &
            nsoa, npoa, nbc
         write(iulog,'(/a56,10i5)') &
           'igas_soa, igas_h2so4, igas_nh3, igas_hno3, igas_hcl', &
            igas_soa, igas_h2so4, igas_nh3, igas_hno3, igas_hcl
         write(iulog,'(/a56,10i5)') &
           'iaer_soa, iaer_so4, iaer_nh4, iaer_no3, iaer_cl', &
            iaer_soa, iaer_so4, iaer_nh4, iaer_no3, iaer_cl
         write(iulog,'(/a56,10i5)') &
           'iaer_pom, iaer_bc, iaer_ncl, iaer_dst, iaer_ca, iaer_co3', &
            iaer_pom, iaer_bc, iaer_ncl, iaer_dst, iaer_ca, iaer_co3
         write(iulog,'(/a56,10i5)') &
           'iaer_mom, ...mpoly, ...mprot, ...mlip, ...mhum, ...mproc', &
            iaer_mom, iaer_mpoly, iaer_mprot, iaer_mlip, iaer_mhum, iaer_mproc
         write(iulog,'(/a)') &
           'fac_eqvso4hyg_aer(1:naer)'
         write(iulog,'(4(a,1pe10.3,3x))') &
           ( name_aerpfx(iaer)(1:6), fac_eqvso4hyg_aer(iaer), iaer=1,naer )

         write(iulog,'(/a)') 'igas, lmap, name, mwhost, mw, fcvt, accom, vmol'
         do igas = 1, ngas
            write(iulog,'(2i4,2x,a,2f10.4,1p,3e12.4)') &
               igas, lmap_gas(igas), name_gas(igas), &
               mwhost_gas(igas), mw_gas(igas), fcvt_gas(igas), &
               accom_coef_gas(igas), vol_molar_gas(igas)
         end do

         do n = 1, ntot_amode
            write(iulog,'(/a,i5)') &
               'iaer, lmap, name, mwhost, mw, fcvt, dens, fac_m2v, hygro for mode', n
            write(iulog,'(2i4,2x,a,20x,1p,e12.4)') &
               0, lmap_num(n), name_num(n), fcvt_num
            write(iulog,'(2i4,2x,a,20x,1p,e12.4)') &
               0, lmap_numcw(n), name_numcw(n)
            do iaer = 1, naer
               if (lmap_aer(iaer,n) > 0) then
                  if (max(mwhost_aer(iaer),mw_aer(iaer)) <= 9999.9999_r8) then
                     fmtaa = '(2i4,2x,a,2f10.4,1p,5e12.4)'
                  else
                     fmtaa = '(2i4,2x,a,2f10.2,1p,5e12.4)'
                  end if
                  write(iulog,fmtaa) &
                     iaer, lmap_aer(iaer,n), name_aer(iaer,n), &
                     mwhost_aer(iaer), mw_aer(iaer), fcvt_aer(iaer), &
                     dens_aer(iaer), fac_m2v_aer(iaer), hygro_aer(iaer)
                  write(iulog,'(2i4,2x,a,2f10.4,1p,4e12.4)') &
                     iaer, lmap_aercw(iaer,n), name_aercw(iaer,n)
               end if
            end do
         end do ! n

         write(iulog,'(/a)') 'l, lmz, lmapcc_all, species_class, name'
         do lmz = 1, gas_pcnst
            l = lmz + loffset
            j = -99
            if (l <= pcnst) j = species_class(l)
            write(iulog,'(4i5,2x,a)') &
               lmz+loffset, lmz, lmapcc_all(lmz), j, solsym(lmz)
         end do

      end if ! ( masterproc )


      call m_a_amicphys_init_history( loffset )


      if ( masterproc ) write(iulog,'(/a)') 'modal_aero_amicphys_init end'

      return
      end subroutine modal_aero_amicphys_init

 
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine mam_set_lptr2_and_specxxx2
!
! initializes the following
!    lptr2_soa_a_amode, lptr2_soa_g_amode, &
!    specdens2_amode, spechygro2, specmw2_amode
! when the multiple nbc/npoa/nsoa flavors is implemented,
!    this can be done in modal_aero_initialize_data
!
      use cam_logfile, only  :  iulog
      use constituents, only :  pcnst, cnst_get_ind, cnst_name

      use modal_aero_data, only : &
          lspectype_amode, &
          lptr_soa_a_amode, lptr2_soa_a_amode, lptr2_soa_g_amode, &
          nspec_amode, ntot_amode, &
          specdens_amode, spechygro, specmw_amode

      implicit none

      integer :: jsoa
      integer :: l1, l2
      integer :: n


      if (nsoa == 1) then
         jsoa = 1
         call cnst_get_ind( 'SOAG', l1, .false. )
         if (l1 < 1 .or. l1 > pcnst) &
            call endrun( 'mam_set_lptr2_and_specxxx2 ERROR - no SOAG' )
         lptr2_soa_g_amode(jsoa) = l1
         do n = 1, ntot_amode
            lptr2_soa_a_amode(n,jsoa) = lptr_soa_a_amode(n)
         end do
      else
         call endrun( 'mam_set_lptr2_and_specxxx2 ERROR - expecting nsoa = 1' )
      end if

      do n = 1, ntot_amode
         do l1 = 1, nspec_amode(n)
            l2 = lspectype_amode(l1,n)
            specmw2_amode(l1,n) = specmw_amode(l2)
            specdens2_amode(l1,n) = specdens_amode(l2)
            spechygro2(l1,n) = spechygro(l2)
         end do
      end do


      return
      end subroutine mam_set_lptr2_and_specxxx2


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine m_a_amicphys_init_history( loffset )

!-----------------------------------------------------------------------
!
! Purpose:
!    set do_adjust and do_aitken flags
!    create history fields for column tendencies associated with
!       modal_aero_calcsize
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

use cam_history, only  :  addfld, horiz_only, add_default, fieldname_len
use cam_logfile, only  :  iulog
use constituents, only :  pcnst, cnst_get_ind, cnst_name
use spmd_utils, only   :  masterproc
use phys_control,only  :  phys_getopts

use modal_aero_data, only : &
    cnst_name_cw, &
    modeptr_accum, modeptr_aitken, modeptr_pcarbon, modeptr_ufine
!use modal_aero_rename

implicit none

!-----------------------------------------------------------------------
! arguments
   integer, intent(in)  :: loffset

!-----------------------------------------------------------------------
! local
   integer  :: iaer, igas, ipair, iok
   integer  :: lmz, lmza, lmzb, lmzc
   integer  :: m
   integer  :: n, na, nb, nc

   real(r8) :: tmp1, tmp2

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(128)                 :: msg
   character(8)                   :: unit
   character(2)                   :: tmpch2

   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
   logical                        :: history_aerocom    ! Output the aerocom history
   !-----------------------------------------------------------------------
 

      call phys_getopts( history_aerosol_out        = history_aerosol   )
#if ( defined CAM_VERSION_IS_ACME ) 
      history_aerocom = .false.
#else
      call phys_getopts( history_aerocom_out        = history_aerocom )
#endif


!
! set the do_q_coltendaa
!
      do_q_coltendaa(:,:) = .false.


! gas-->aer condensation and resulting aging
      do igas = 1, ngas
         lmz = lmap_gas(igas)
         if (lmz <= 0) cycle
         do_q_coltendaa(lmz,iqtend_cond) = .true.
         iaer = igas
         do n = 1, ntot_amode
            lmz = lmap_aer(iaer,n)
            if (lmz <= 0) cycle
            do_q_coltendaa(lmz,iqtend_cond) = .true.
         end do ! n
      end do ! igas

      do ipair = 1, n_agepair
         na = modefrm_agepair(ipair)
         nb = modetoo_agepair(ipair)
         if (na < 1 .or. nb < 1) cycle

         lmza = lmap_num(na)
         lmzb = lmap_num(nb)
         do_q_coltendaa(lmza,iqtend_cond) = .true.
         do_q_coltendaa(lmzb,iqtend_cond) = .true.
         do iaer = 1, naer
            lmza = lmap_aer(iaer,na)
            lmzb = lmap_aer(iaer,nb)
            if (lmza > 0) then
               do_q_coltendaa(lmza,iqtend_cond) = .true.
               if (lmzb > 0) do_q_coltendaa(lmzb,iqtend_cond) = .true.
            end if
         end do ! iaer
      end do ! ipair

!  define history fields for gas-->aer condensation and resulting aging
      do lmz = 1, gas_pcnst
         if ( do_q_coltendaa(lmz,iqtend_cond)) then
            tmpnamea = cnst_name(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfgaex1'
            long_name = trim(tmpnamea) // ' gas-aerosol-exchange primary column tendency'
            unit = 'kg/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
         end if
      end do ! lmz

!  define history fields for 3d soa production for aerocom
      do igas = 1, nsoa
         lmz = lmap_gas(igas)
         if (lmz <= 0) cycle
         if ( .not. do_q_coltendaa(lmz,iqtend_cond)) cycle
         if ( .not. history_aerocom ) cycle

         tmpnamea = cnst_name(lmz+loffset)
         fieldname = trim(tmpnamea) // '_sfgaex3d'
         long_name = trim(tmpnamea) // ' gas-aerosol-exchange primary 3d tendency'
         unit = 'kg/m2/s'
         call addfld( fieldname, (/ 'lev' /), 'A', unit, long_name )
         call add_default( fieldname, 1, ' ' )
         if ( masterproc ) write(iulog,'(3(a,3x),2i5)') &
            'gasaerexch addfld', fieldname, unit, igas, lmz+loffset
      end do


! renaming during gas-->aer condensation or cloud chemistry
      na = modeptr_aitken
      nb = modeptr_accum
      if (na > 0 .and. nb > 0) then
         lmza = lmap_num(na)
         lmzb = lmap_num(nb)
         do_q_coltendaa(lmza,iqtend_rnam) = .true.
         do_q_coltendaa(lmzb,iqtend_rnam) = .true.
         lmza = lmap_numcw(na)
         lmzb = lmap_numcw(nb)
         do_qqcw_coltendaa(lmza,iqqcwtend_rnam) = .true.
         do_qqcw_coltendaa(lmzb,iqqcwtend_rnam) = .true.
         do iaer = 1, naer
            lmza = lmap_aer(iaer,na)
            lmzb = lmap_aer(iaer,nb)
            if (lmza > 0) then
               do_q_coltendaa(lmza,iqtend_rnam) = .true.
               if (lmzb > 0) do_q_coltendaa(lmzb,iqtend_rnam) = .true.
            end if
            lmza = lmap_aercw(iaer,na)
            lmzb = lmap_aercw(iaer,nb)
            if (lmza > 0) then
               do_qqcw_coltendaa(lmza,iqqcwtend_rnam) = .true.
               if (lmzb > 0) do_qqcw_coltendaa(lmzb,iqqcwtend_rnam) = .true.
            end if
         end do ! iaer
      end if ! (na > 0 .and. nb > 0)

!  define history fields for renaming during gas-->aer condensation or cloud chemistry
      do lmz = 1, gas_pcnst
         if ( do_q_coltendaa(lmz,iqtend_rnam)) then
            tmpnamea = cnst_name(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfgaex2'
            long_name = trim(tmpnamea) // ' gas-aerosol-exchange renaming column tendency'
            unit = 'kg/m2/s'
            if (tmpnamea(1:4) == 'num_' .or. tmpnamea(1:4) == 'NUM_') unit = '#/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
         end if
         if ( do_qqcw_coltendaa(lmz,iqqcwtend_rnam)) then
            tmpnamea = cnst_name_cw(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfgaex2'
            long_name = trim(tmpnamea) // ' gas-aerosol-exchange renaming column tendency'
            unit = 'kg/m2/s'
            if (tmpnamea(1:4) == 'num_' .or. tmpnamea(1:4) == 'NUM_') unit = '#/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
         end if
      end do ! lmz



! coagulation
      do ipair = 1, n_coagpair
         na = modefrm_coagpair(ipair)
         nb = modetoo_coagpair(ipair)
         nc = modeend_coagpair(ipair)
         if (na < 1 .or. nb < 1 .or. nc < 1) cycle

         lmza = lmap_num(na)
         lmzb = lmap_num(nb)
         lmzc = lmap_num(nc)
         do_q_coltendaa(lmza,iqtend_coag) = .true.
         do_q_coltendaa(lmzb,iqtend_coag) = .true.
         do_q_coltendaa(lmzc,iqtend_coag) = .true.
         do iaer = 1, naer
            lmza = lmap_aer(iaer,na)
            lmzb = lmap_aer(iaer,nb)
            lmzc = lmap_aer(iaer,nc)
            if (lmza > 0) then
               do_q_coltendaa(lmza,iqtend_coag) = .true.
               if (lmzc > 0) do_q_coltendaa(lmzc,iqtend_coag) = .true.
            end if
            if (nb == nc) cycle
            if (lmzb > 0) then
               do_q_coltendaa(lmzb,iqtend_coag) = .true.
               if (lmzc > 0) do_q_coltendaa(lmzc,iqtend_coag) = .true.
            end if
         end do ! iaer
      end do ! ipair

!  define history fields for coagulation
      do lmz = 1, gas_pcnst
         if ( do_q_coltendaa(lmz,iqtend_coag)) then
            tmpnamea = cnst_name(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfcoag1'
            long_name = trim(tmpnamea) // ' modal_aero coagulation column tendency'
            unit = 'kg/m2/s'
            if (tmpnamea(1:4) == 'num_' .or. tmpnamea(1:4) == 'NUM_') unit = '#/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'modal_aero_coag_init addfld', fieldname, unit
         end if
      end do ! lmz


! nucleation
      n = modeptr_aitken
      do igas = 1, ngas
         iok = 0
         if (igas == igas_h2so4) iok = 1
         if (igas == igas_nh3  ) iok = 1
         if (iok <= 0) cycle
         lmz = lmap_gas(igas)
         if (lmz > 0) then
            do_q_coltendaa(lmz,iqtend_nnuc) = .true.
            iaer = igas
            lmz = lmap_aer(iaer,n)
            if (lmz > 0) do_q_coltendaa(lmz,iqtend_nnuc) = .true.
          end if
      end do ! igas
      lmzc = lmap_num(n)
      do_q_coltendaa(lmzc,iqtend_nnuc) = .true.

!  define history fields for nucleation
      do lmz = 1, gas_pcnst
         if ( do_q_coltendaa(lmz,iqtend_nnuc)) then
            tmpnamea = cnst_name(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfnnuc1'
            long_name = trim(tmpnamea) // ' modal_aero new particle nucleation column tendency'
            unit = 'kg/m2/s'
            if (tmpnamea(1:4) == 'num_' .or. tmpnamea(1:4) == 'NUM_') unit = '#/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'modal_aero_newnuc_init addfld', fieldname, unit
         end if
      end do ! lmz

      if ( history_aerocom ) then
            tmpnamea = cnst_name(lmzc+loffset)
            fieldname = trim(tmpnamea) // '_nuc1'
            long_name = trim(tmpnamea) // ' modal_aero new particle nucleation tendency'
            unit = '#/m3/s'
            call addfld( fieldname, (/ 'lev' /), 'A', unit, long_name )
            call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,2x))') &
                'modal_aero_newnuc_init addfld', fieldname, unit

            fieldname = trim(tmpnamea) // '_nuc2'
            long_name = trim(tmpnamea) // ' modal_aero cluster nucleation rate'
            unit = '#/m3/s'
            call addfld( fieldname, (/ 'lev' /), 'A', unit, long_name )
            call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,2x))') &
                'modal_aero_newnuc_init addfld', fieldname, unit
      endif


#if ( defined( MOSAIC_SPECIES ) )
      if ( mosaic ) then
         !BSINGH - Adding addfld and add_default call for tracking convergence failures
         call addfld('convergence_fail', (/ 'lev' /), 'A', 'no units', 'For tracking MOSAIC convergence failure' )
         call addfld('max_kelvin_iter', (/ 'lev' /), 'A', 'no units', 'For tracking when MOSAIC kelvin iterations hit max ' )
         call add_default( 'convergence_fail', 1, ' ' )
         call add_default( 'max_kelvin_iter', 1, ' ' )

         do n = 1, 4
         do m = 1, 5
            fieldname = ' '
            write( fieldname(1:16), '(a,i1,a,i1)') 'astem_negval_', m, '_', n
            call addfld( fieldname, (/ 'lev' /), 'A', 'no units', 'For tracking ASTEM negative values' )
            call add_default( fieldname, 1, ' ' )
         end do
         end do
      end if
#endif

      return
      end subroutine m_a_amicphys_init_history

 
!----------------------------------------------------------------------

end module modal_aero_amicphys


