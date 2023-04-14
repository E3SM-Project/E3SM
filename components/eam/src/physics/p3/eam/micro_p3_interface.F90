module micro_p3_interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Interface between E3SM and P3 microphysics
  !!
  !! Author: Peter Caldwell
  !!
  !! Last updated: 2018-09-12
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use shr_kind_mod,   only: rtype=>shr_kind_r8
  use ppgrid,         only: pcols,pver,pverp

!comment: I think Kai added handle_errmsg. It would be better to 
!use standard E3SM libraries if possible.
  use error_messages, only: handle_errmsg

  use physics_types,  only: physics_state, &
                            physics_ptend, &
                            physics_ptend_init
  use physconst,      only: mwdry, cpair, mwh2o, gravit, rair, cpliq, pi, &
                            rh2o, latvap, latice, tmelt, rhoh2o, rairv 
  use constituents,   only: cnst_add, pcnst, sflxnam, apcnst, bpcnst, pcnst,&
                            cnst_name, cnst_get_ind,cnst_longname
  use physics_buffer, only: physics_buffer_desc, dtype_r8, &
                            pbuf_get_field, pbuf_add_field,dyn_time_lvls,dtype_i4, &
                            pbuf_set_field, pbuf_get_index, &
                            pbuf_old_tim_idx
  use ref_pres,       only: top_lev=>trop_cloud_top_lev
  use phys_control,   only: phys_getopts,use_hetfrz_classnuc
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use time_manager,   only: is_first_step, get_curr_date
  use perf_mod,       only: t_startf, t_stopf
  use micro_p3_utils, only: do_Cooper_inP3
  use pio,            only: file_desc_t, pio_nowrite
  use cam_pio_utils,    only: cam_pio_openfile,cam_pio_closefile
  use cam_grid_support, only: cam_grid_check, cam_grid_id, cam_grid_get_dim_names
  use ncdio_atm,       only: infld
  use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp,psubcols
  use cam_history_support, only: add_hist_coord
       
  implicit none
  save

  public :: micro_p3_init, micro_p3_register, micro_p3_tend, &
            micro_p3_init_cnst, micro_p3_implements_cnst,    &
            micro_p3_readnl                               
            

  character(len=16), parameter :: unset_str = 'UNSET'

  private

  !Define indices for state%q constituents at module level so
  !defining them in micro_p3_register makes them permanently 
  !available.
  CHARACTER(len=16) :: precip_frac_method = 'max_overlap'  ! AaronDonahue, Hard-coded for now, should be fixed in the future

  integer, public ::    &
       ixcldliq = -1,   & ! cloud liquid amount index
       ixcldice = -1,      & ! ice index
       ixnumliq = -1,   & ! cloud liquid number index
       ixnumice = -1,   & ! cloud ice number index
       ixrain   = -1,   & ! rain index
       ixnumrain= -1,   & ! rain number index
       ixcldrim = -1,      & ! rime index ??
       ixrimvol  = -1,  & ! rime volume index ??
       ixqm  = -1      ! ?? index ??

!! pbuf 
   integer :: &
      cldo_idx,           &
      qme_idx,            &
      precip_total_tend_idx,          &
      nevapr_idx,         &
      dei_idx,            &
      rate1_cw2pr_st_idx, &
      mu_idx,             &
      lambdac_idx,        &
      rei_idx,            &
      rel_idx,            &
      ls_flxprc_idx,      &
      ls_flxsnw_idx,      &
      ls_reffrain_idx,    &
      ls_reffsnow_idx,    &
      cv_reffliq_idx,     &
      cv_reffice_idx,     &
      qr_evap_tend_idx,      &
      cmeliq_idx,         &
      relvar_idx,         &
      qv_prev_idx,        &
      t_prev_idx,         &
      accre_enhan_idx,    &
      mon_ccn_1_idx,      &
      mon_ccn_2_idx,      &
      current_month      !Needed for prescribed CCN option         

! Physics buffer indices for fields registered by other modules
   integer :: &
      ast_idx = -1            

   integer :: &
      ni_activated_idx = -1,           &
      npccn_idx = -1,          &
      prec_str_idx = -1,       &
      prec_pcw_idx = -1,       &
      prec_sed_idx = -1,       &
      snow_str_idx = -1,       &
      snow_pcw_idx = -1,       &
      snow_sed_idx = -1
   
   integer :: &
      frzimm_idx = -1, &
      frzcnt_idx = -1, &
      frzdep_idx = -1

   real(rtype) :: &
      micro_mg_accre_enhan_fac = huge(1.0_rtype), & !Accretion enhancement factor from namelist
      prc_coef1_in             = huge(1.0_rtype), &
      prc_exp_in               = huge(1.0_rtype), &
      prc_exp1_in              = huge(1.0_rtype), &
      p3_autocon_coeff         = huge(1.0_rtype), &
      p3_accret_coeff          = huge(1.0_rtype), &
      p3_qc_autocon_expon      = huge(1.0_rtype), &
      p3_nc_autocon_expon      = huge(1.0_rtype), &
      p3_qc_accret_expon       = huge(1.0_rtype), &
      p3_wbf_coeff             = huge(1.0_rtype), &
      p3_max_mean_rain_size    = huge(1.0_rtype), &
      p3_embryonic_rain_size   = huge(1.0_rtype)
   

   integer :: ncnst

   character(len=8), parameter :: &      ! Constituent names
      cnst_names(8) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE', &
                      'RAINQM', 'CLDRIM','NUMRAI','BVRIM '/)

   character(len=128) :: micro_p3_lookup_dir     = unset_str ! location of p3 input files
   character(len=16)  :: micro_p3_tableversion   = unset_str ! P3 table version
   logical            :: micro_aerosolactivation = .false.   ! Use aerosol activation
   logical            :: micro_subgrid_cloud     = .false.   ! Use subgrid cloudiness
   logical            :: micro_tend_output       = .false.   ! Default microphysics tendencies to output file
   logical            :: do_prescribed_CCN       = .false.   ! Use prescribed CCN

   contains
!===============================================================================
subroutine micro_p3_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'micro_p3_cam_readnl'


end subroutine micro_p3_readnl
  !================================================================================================

  subroutine micro_p3_register()

  logical :: prog_modal_aero ! prognostic aerosols

  end subroutine micro_p3_register

  !================================================================================================
  function micro_p3_implements_cnst(name)

    ! Return true if specified constituent is implemented by the
    ! microphysics package

    character(len=*), intent(in) :: name        ! constituent name
    logical :: micro_p3_implements_cnst    ! return value


  end function micro_p3_implements_cnst


  !================================================================================================

  subroutine micro_p3_init_cnst(name, q)

    ! Initialize the microphysics constituents, if they are
    ! not read from the initial file.

    character(len=*), intent(in) :: name     ! constituent name
    real(rtype), intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)


  end subroutine micro_p3_init_cnst

  !================================================================================================

  subroutine get_prescribed_CCN_from_file(micro_p3_lookup_dir,month_int, ccn_values)

      character*(*), intent(in)  :: micro_p3_lookup_dir !directory of the lookup tables
      real(rtype), intent(inout) :: ccn_values(pcols,pver,begchunk:endchunk)    
      integer, intent(in) :: month_int

      !internal variables
      character(len=100) :: base_file_name
      character(len=500) :: filename
      character(len=20) :: dim1name, dim2name
      character(len=20) :: mon_str
      type(file_desc_t) :: nccn_ncid
      integer :: year, month, day, tod, next_month, grid_id
      logical :: found = .false.

   
  end subroutine get_prescribed_CCN_from_file

  !================================================================================================

  subroutine micro_p3_init(pbuf2d)
    use micro_p3,       only: p3_init
    use cam_history,    only: addfld, add_default, horiz_only 
    use micro_p3_utils, only: micro_p3_utils_init

    type(physics_buffer_desc),  pointer :: pbuf2d(:,:)
    integer        :: m, mm
    integer        :: ierr
    logical :: history_amwg         ! output the variables used by the AMWG diag package
    logical :: history_verbose      ! produce verbose history output
    logical :: history_budget       ! Output tendencies and state variables for CAM4
    integer :: budget_histfile      ! output history file number for budget fields
                                   ! temperature, water vapor, cloud ice and cloud

    !needed for prescribed CCN option:
    character(len=20) :: base_file_name
    character(len=500) :: filename, filename_next_month
    character(len=20) :: dim1name, dim2name
    type(file_desc_t) :: nccn_ncid
    integer :: year, month, day, tod, next_month, grid_id
    logical :: found = .false.
    real(rtype), pointer :: ccn_values(:,:,:)


  end subroutine micro_p3_init

  !================================================================================================
    subroutine get_cloud_fraction(its,ite,kts,kte,ast,qc,qr,qi,method, &
                  cld_frac_i,cld_frac_l,cld_frac_r)
      
       use micro_p3_utils, only: mincld, qsmall

       integer,intent(in)                                 :: its,ite,kts,kte
       real(rtype),dimension(its:ite,kts:kte),intent(in)  :: ast, qc, qr, qi
       character(len=16),intent(in)                       :: method
       real(rtype),dimension(its:ite,kts:kte),intent(out) :: cld_frac_i, cld_frac_l, cld_frac_r
       real(rtype),dimension(its:ite,kts:kte)             :: cldm

       integer  :: i,k
       integer  :: ktop, kbot, kdir

    end subroutine get_cloud_fraction

  !================================================================================================

    subroutine get_prescribed_CCN(nccn_prescribed,micro_p3_lookup_dir,its,ite,kts,kte,pbuf,lchnk)

      !INOUT/OUTPUT VARIABLES
      integer,intent(in) :: its,ite,kts,kte,lchnk
      real(rtype),dimension(its:ite,kts:kte),intent(inout)  :: nccn_prescribed
      character*(*), intent(in)    :: micro_p3_lookup_dir !directory of the lookup tables
      type(physics_buffer_desc),   pointer       :: pbuf(:)



    end subroutine get_prescribed_CCN

  !================================================================================================
  subroutine micro_p3_tend(state, ptend, dtime, pbuf)

    use phys_grid,      only: get_rlat_all_p, get_rlon_all_p, get_gcol_all_p
    use time_manager,   only: get_nstep
    use cam_history,    only: outfld
    use time_manager,   only: get_nstep
    use micro_p3,       only: p3_main
    use micro_p3_utils, only: avg_diameter, &
                              rho_h2o, &
                              rho_h2os, &
                              qsmall, &
                              mincld, & 
                              inv_cp 

    !INPUT/OUTPUT VARIABLES
    type(physics_state),         intent(in)    :: state
    type(physics_ptend),         intent(out)   :: ptend
    real(rtype),                 intent(in)    :: dtime
    type(physics_buffer_desc),   pointer       :: pbuf(:)
    logical :: lq(pcnst)   !list of what constituents to update

    !INTERNAL VARIABLES
    real(rtype) :: dz(pcols,pver)        !geometric layer thickness              m
    real(rtype) :: cldliq(pcols,pver)     !cloud liquid water mixing ratio        kg/kg
    real(rtype) :: numliq(pcols,pver)     !cloud liquid water drop concentraiton  #/kg
    real(rtype) :: rain(pcols,pver)       !rain water mixing ratio                kg/kg
    real(rtype) :: numrain(pcols,pver)    !rain water number concentration        #/kg
    real(rtype) :: qv(pcols,pver)         !water vapor mixing ratio               kg/kg
    real(rtype) :: ice(pcols,pver)        !total ice water mixing ratio           kg/kg
    real(rtype) :: qm(pcols,pver)      !rime ice mixing ratio                  kg/kg
    real(rtype) :: numice(pcols,pver)     !total ice crystal number concentration #/kg
    real(rtype) :: rimvol(pcols,pver)     !rime volume mixing ratio               m3/kg
    real(rtype) :: temp(pcols,pver)       !temperature copy needed for tendency   K
    real(rtype) :: th(pcols,pver)         !potential temperature                  K
    real(rtype) :: precip_liq_surf(pcols)         !precipitation rate, liquid             m s-1
    real(rtype) :: precip_ice_surf(pcols)         !precipitation rate, solid              m s-1

    real(rtype) :: rho_qi(pcols,pver)  !bulk density of ice                    kg m-1
    real(rtype) :: pres(pcols,pver)       !pressure at midlevel                   hPa
    real(rtype) :: qv2qi_depos_tend(pcols,pver)
    real(rtype) :: precip_liq_flux(pcols,pver+1)     !grid-box average rain flux (kg m^-2s^-1) pverp
    real(rtype) :: precip_ice_flux(pcols,pver+1)     !grid-box average ice/snow flux (kg m^-2s^-1) pverp
    real(rtype) :: rflx(pcols,pver+1)     !grid-box average rain flux (kg m^-2s^-1) pverp
    real(rtype) :: sflx(pcols,pver+1)     !grid-box average ice/snow flux (kg m^-2s^-1) pverp
    real(rtype) :: cflx(pcols,pver+1)     !grid-box average cloud flux (kg m^-2s^-1) pverp
    real(rtype) :: exner(pcols,pver)      !exner formula for converting between potential and normal temp
    real(rtype) :: cld_frac_r(pcols,pver)      !rain cloud fraction
    real(rtype) :: cld_frac_l(pcols,pver)      !liquid cloud fraction
    real(rtype) :: cld_frac_i(pcols,pver)      !ice cloud fraction
    real(rtype) :: tend_out(pcols,pver,49) !microphysical tendencies
    real(rtype), dimension(pcols,pver) :: liq_ice_exchange ! sum of liq-ice phase change tendenices
    real(rtype), dimension(pcols,pver) :: vap_liq_exchange ! sum of vap-liq phase change tendenices
    real(rtype), dimension(pcols,pver) :: vap_ice_exchange ! sum of vap-ice phase change tendenices
    real(rtype), dimension(pcols,pver) :: diag_equiv_reflectivity,diag_ze_rain,diag_ze_ice ! equivalent reflectivity [dBz]

    !Prescribed CCN concentration
    real(rtype), dimension(pcols,pver) :: nccn_prescribed

    ! PBUF Variables
    real(rtype), pointer :: ast(:,:)      ! Relative humidity cloud fraction
    real(rtype), pointer :: ni_activated(:,:)     ! ice nucleation number
    real(rtype), pointer :: npccn(:,:)    ! liquid activation number tendency
    real(rtype), pointer :: cmeliq(:,:)
    !!
    real(rtype), pointer :: prec_str(:)    ! [Total] Sfc flux of precip from stratiform [ m/s ]
    real(rtype), pointer :: prec_sed(:)    ! Surface flux of total cloud water from sedimentation
    real(rtype), pointer :: prec_pcw(:)    ! Sfc flux of precip from microphysics [ m/s ]
    real(rtype), pointer :: snow_str(:)    ! [Total] Sfc flux of snow from stratiform   [ m/s ]
    real(rtype), pointer :: snow_pcw(:)    ! Sfc flux of snow from microphysics [ m/s ]
    real(rtype), pointer :: snow_sed(:)    ! Surface flux of cloud ice from sedimentation
    real(rtype), pointer :: relvar(:,:)    ! cloud liquid relative variance [-]
    real(rtype), pointer :: cldo(:,:)      ! Old cloud fraction
    real(rtype), pointer :: qr_evap_tend(:,:) ! precipitation evaporation rate
    real(rtype), pointer :: qv_prev(:,:)   ! qv from previous p3_main call
    real(rtype), pointer :: t_prev(:,:)    ! t from previous p3_main call
    !! wetdep 
    real(rtype), pointer :: qme(:,:)
    real(rtype), pointer :: precip_total_tend(:,:)        ! Total precipitation (rain + snow)
    real(rtype), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
    !! COSP simulator
    real(rtype), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
    real(rtype), pointer :: rei(:,:)          ! Ice effective drop size (microns)
    real(rtype), pointer :: flxprc(:,:)     ! P3 grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
    real(rtype), pointer :: flxsnw(:,:)     ! P3 grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
    real(rtype), pointer :: reffrain(:,:)   ! P3 diagnostic rain effective radius (um)
    real(rtype), pointer :: reffsnow(:,:)   ! P3 diagnostic snow effective radius (um)
    real(rtype), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
    real(rtype), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)
    !! radiation 
    real(rtype), pointer :: dei(:,:)          ! Ice effective diameter (um)
    real(rtype), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
    real(rtype), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
    ! DONE PBUF
    ! For recording inputs/outputs to p3_main
    real(rtype) :: p3_main_inputs(pcols,pver+1,17) ! Record of inputs for p3_main
    real(rtype) :: p3_main_outputs(pcols,pver+1,31) ! Record of outputs for p3_main

    ! Derived Variables
    real(rtype) :: icimrst(pcols,pver) ! stratus ice mixing ratio - on grid
    real(rtype) :: icwmrst(pcols,pver) ! stratus water mixing ratio - on grid
    real(rtype) :: rho(pcols,pver)
    real(rtype) :: drout2(pcols,pver)
    real(rtype) :: reff_rain(pcols,pver)
    real(rtype) :: col_location(pcols,3),tmp_loc(pcols)  ! Array of column lon (index 1) and lat (index 2)
    integer     :: tmpi_loc(pcols) ! Global column index temp array

    ! Variables used for microphysics output
    real(rtype) :: aqrain(pcols,pver)
    real(rtype) :: anrain(pcols,pver)
    real(rtype) :: nfice(pcols,pver)
    real(rtype) :: efcout(pcols,pver)      
    real(rtype) :: efiout(pcols,pver)      
    real(rtype) :: ncout(pcols,pver)      
    real(rtype) :: niout(pcols,pver)      
    real(rtype) :: freqr(pcols,pver)      
    real(rtype) :: freql(pcols,pver)      
    real(rtype) :: freqi(pcols,pver)      
    real(rtype) :: cdnumc(pcols)      
    real(rtype) :: icinc(pcols,pver) 
    real(rtype) :: icwnc(pcols,pver) 

    ! variables for the CNT primary / heterogeneous freezing
    real(rtype), pointer :: frzimm(:,:)
    real(rtype), pointer :: frzcnt(:,:)
    real(rtype), pointer :: frzdep(:,:)
    real(rtype) :: frzimm_in(pcols,pver)
    real(rtype) :: frzcnt_in(pcols,pver)
    real(rtype) :: frzdep_in(pcols,pver)
 
    integer :: it                      !timestep counter                       -
    integer :: its, ite                !horizontal bounds (column start,finish)
    integer :: kts                     !closest level to TOM                   -
    integer :: kte                     !near surface level                     -

    logical :: do_predict_nc           !prognostic droplet concentration or not?
    logical :: do_subgrid_clouds       !use subgrid cloudiness in tendency calculations?
    integer :: icol, ncol, k
    integer :: psetcols, lchnk
    integer :: itim_old

    ! For rrtmg optics. specified distribution.
    real(rtype), parameter :: dcon   = 25.e-6_rtype         ! Convective size distribution effective radius (um)
    real(rtype), parameter :: mucon  = 5.3_rtype            ! Convective size distribution shape parameter
    real(rtype), parameter :: deicon = 50._rtype            ! Convective ice effective diameter (um)

    end subroutine micro_p3_tend

end module micro_p3_interface
