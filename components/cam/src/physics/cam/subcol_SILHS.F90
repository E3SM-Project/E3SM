module subcol_SILHS
   !---------------------------------------------------------------------------
   ! Purpose:
   !
   ! Implement a subcolumn scheme based on the Subgrid Importance Latin Hypercube 
   ! Sampling (SILHS) functionality of the CLUBB moist turbulence parameterization.
   !
   !---------------------------------------------------------------------------

   use shr_kind_mod,     only: r8=>shr_kind_r8, r4=>shr_kind_r4, i4=>shr_kind_i4
   use physics_types,    only: physics_state, physics_tend, physics_ptend
   use ppgrid,           only: pcols, psubcols, pver, pverp
   use constituents,     only: pcnst, cnst_get_ind
!   use abortutils,       only: endrun
   use shr_sys_mod,      only: endrun => shr_sys_abort
   use cam_logfile,      only: iulog
   use cam_history,      only: addfld, add_default, outfld, horiz_only
#ifdef CLUBB_SGS
#ifdef SILHS
   use clubb_api_module, only: hmp2_ip_on_hmm2_ip_ratios_type
#endif
#endif
   use physconst,     only: cpair, gravit, latvap, latice, rair

   implicit none

      private
      save

   public :: subcol_register_SILHS  ! 
   public :: subcol_init_SILHS      ! Initialize 
   public :: subcol_gen_SILHS       ! Generate subcolumn fields by calling SILHS 
   public :: subcol_readnl_SILHS    ! SILHS namelist reader
   public :: subcol_ptend_avg_SILHS
#ifdef SILHS
   public :: subcol_SILHS_var_covar_driver, subcol_SILHS_massless_droplet_destroyer
   private :: Abs_Temp_profile
   private :: StaticEng_profile
   ! Calc subcol mean ! Calc subcol variance
   private :: meansc
   private :: stdsc
#endif

   !-----
   ! Private module vars
   !-----
   ! Pbuf indicies
   integer :: thlm_idx, num_subcols_idx, pdf_params_idx, rcm_idx, rtm_idx, ice_supersat_idx, &
              alst_idx, cld_idx, wp2_idx, qrain_idx, qsnow_idx, nrain_idx, nsnow_idx, ztodt_idx

   logical :: subcol_SILHS_weight  ! if set, sets up weights for averaging subcolumns for SILHS
   integer :: subcol_SILHS_numsubcol ! number of subcolumns for this run
   logical :: docldfracscaling = .false. ! Weight tendencies by cloud fraction

   character(len=256) :: subcol_SILHS_corr_file_path
   character(len=16)  :: subcol_SILHS_corr_file_name

   logical :: subcol_SILHS_q_to_micro, &
              subcol_SILHS_n_to_micro, &
              subcol_SILHS_use_clear_col, &
              subcol_SILHS_meanice, &
              subcol_SILHS_constrainmn

   logical :: subcol_SILHS_var_covar_src
   logical :: subcol_SILHS_destroy_massless_droplets

   real(r8) :: subcol_SILHS_ncnp2_on_ncnm2

   ! There may or may not be a better place to put this.
   real(r8), parameter :: p0_clubb = 100000._r8


!   real(r8) :: subcol_SILHS_c6rt, subcol_SILHS_c7, subcol_SILHS_c8, subcol_SILHS_c11, &
!               subcol_SILHS_c11b, subcol_SILHS_gamma_coef, &
!               subcol_SILHS_mult_coef, subcol_SILHS_mu

   real(r8) :: ztodt  ! model timestep
#ifdef CLUBB_SGS
#ifdef SILHS
   type(hmp2_ip_on_hmm2_ip_ratios_type) :: hmp2_ip_on_hmm2_ip_ratios
#endif
#endif

contains

   subroutine subcol_register_SILHS()

      !--------------------------------
      ! Register fields needed by SILHS in the physics buffer
      ! Currently, most fields needed by SILHS but calculated by CLUBB are registered
      ! by clubb in clubb_intr.F90.
      ! 
      !--------------------------------
#ifdef CLUBB_SGS
#ifdef SILHS
      use physics_buffer,  only: pbuf_add_field, dtype_r8
      use ppgrid,          only: pver, pverp, pcols, psubcols


     ! Diagnostic fields needed for subcol_SILHS, need to be grid-only
      call pbuf_add_field('QRAIN',      'global',dtype_r8,(/pcols,pver/), qrain_idx)
      call pbuf_add_field('QSNOW',      'global',dtype_r8,(/pcols,pver/), qsnow_idx)
      call pbuf_add_field('NRAIN',      'global',dtype_r8,(/pcols,pver/), nrain_idx)
      call pbuf_add_field('NSNOW',      'global',dtype_r8,(/pcols,pver/), nsnow_idx)

#endif
#endif
   end subroutine subcol_register_SILHS


   subroutine subcol_readnl_SILHS(nlfile)
#ifdef CLUBB_SGS
#ifdef SILHS
      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use spmd_utils,      only: masterproc, masterprocid, mpicom
      use spmd_utils,      only: mpi_integer, mpi_logical, mpi_character, mpir8
      use clubb_api_module,only: core_rknd
#endif
#endif
      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
#ifdef CLUBB_SGS
#ifdef SILHS
      namelist /subcol_SILHS_nl/ subcol_SILHS_weight, &
                                 subcol_SILHS_numsubcol, &
                                 subcol_SILHS_corr_file_path, &
                                 subcol_SILHS_corr_file_name, &
                                 subcol_SILHS_q_to_micro, &
                                 subcol_SILHS_n_to_micro, &
                                 subcol_SILHS_ncnp2_on_ncnm2, &
                                 hmp2_ip_on_hmm2_ip_ratios, &
                                 subcol_SILHS_meanice, &
                                 subcol_SILHS_use_clear_col, &
                                 subcol_SILHS_constrainmn
!                                 subcol_SILHS_c6rt, subcol_SILHS_c7, &
!                                 subcol_SILHS_c8, subcol_SILHS_c11, subcol_SILHS_c11b, &
!                                 subcol_SILHS_gamma_coef, subcol_SILHS_mult_coef, subcol_SILHS_mu

      real( kind = core_rknd ), parameter :: &
      init_value = -999._core_rknd ! Initial value for the parameters, used to detect missing values


      !-----------------------------------------------------------------------------
      ! Set defaults

      subcol_SILHS_var_covar_src = .true.  ! TODO: put this in namelist
      subcol_SILHS_destroy_massless_droplets = .true.  ! TODO: put this in namelist

      ! Eric Raut changed a default.
      ! hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip = 0.5_core_rknd
      ! initialize the namelist varriables in hmp2_ip_on_hmm2_ip_ratios
      ! If not specified via namelist, all set to be equal to rrp2_ip_on_rrm2_ip
      ! after processing the namelist file.

      hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip = 1.0_core_rknd
      hmp2_ip_on_hmm2_ip_ratios%Nrp2_ip_on_Nrm2_ip = init_value
      hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip = init_value
      hmp2_ip_on_hmm2_ip_ratios%Nsp2_ip_on_Nsm2_ip = init_value
      hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip = init_value
      hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip = init_value

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'subcol_SILHS_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, subcol_SILHS_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun('subcol_readnl_SILHS: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
         if (hmp2_ip_on_hmm2_ip_ratios%Nrp2_ip_on_Nrm2_ip == init_value) &
             hmp2_ip_on_hmm2_ip_ratios%Nrp2_ip_on_Nrm2_ip =  &
             hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip

         if (hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip == init_value) &
             hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip = &
             hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip

         if (hmp2_ip_on_hmm2_ip_ratios%Nsp2_ip_on_Nsm2_ip == init_value) &
             hmp2_ip_on_hmm2_ip_ratios%Nsp2_ip_on_Nsm2_ip = &
             hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip

         if (hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip == init_value) &
             hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip = &
             hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip

         if (hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip == init_value) &
             hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip = &
             hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip

      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpi_bcast(subcol_SILHS_weight,    1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_numsubcol, 1, mpi_integer, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_corr_file_path, len(subcol_SILHS_corr_file_path), &
                     mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_corr_file_name, len(subcol_SILHS_corr_file_name), &
                     mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_use_clear_col, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_constrainmn, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_meanice, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_q_to_micro, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_n_to_micro, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_var_covar_src,1,mpi_logical,masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_destroy_massless_droplets, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_ncnp2_on_ncnm2, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_ratios%nrp2_ip_on_nrm2_ip, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_ratios%nip2_ip_on_nim2_ip, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_ratios%nsp2_ip_on_nsm2_ip, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c6rt, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c7, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c8, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c11, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c11b, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_gamma_coef, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_mult_coef, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_mu, 1, mpir8, masterprocid, mpicom, ierr)

! SPMD
#endif
! SILHS
#endif
! CLUBB_SGS
#endif
   end subroutine subcol_readnl_SILHS


   subroutine subcol_init_SILHS(pbuf2d)

      !--------------------------------
      ! Read in parameters and initialize SILHS PDF fields.
      ! Set up indexes into Pbuf fields.
      ! Register history outputs.
      !--------------------------------

      use physics_buffer,          only: physics_buffer_desc, pbuf_get_field, &
                                         dtype_r8, pbuf_get_index
      use units,                   only: getunit, freeunit 
#ifdef CLUBB_SGS
#ifdef SILHS
      use clubb_api_module,        only: iirrm, iiNrm, iirsm, iirim, &
                                         iirgm, iiNsm, &
                                         iiNim, iiNgm, &
                                         l_mix_rat_hm, l_frozen_hm, &

                                         hydromet_dim, &

                                         l_prescribed_avg_deltaz, &

                                         rr_tol, ri_tol, rs_tol, rg_tol, &
                                         Nr_tol, Ni_tol, Ns_tol, Ng_tol, Nc_tol, &

                                         l_diagnose_correlations, &
                                         l_calc_w_corr, l_use_cloud_cover, &
                                         l_fix_chi_eta_correlations, l_const_Nc_in_cloud, &

                                         d_variables, &
                                         setup_corr_varnce_array_api, &
                                         setup_pdf_indices_api, &

                                         hydromet_list, hydromet_tol, &
                                         hmp2_ip_on_hmm2_ip, &
                                         Ncnp2_on_Ncnm2, &

                                         set_clubb_debug_level_api

    use silhs_api_module, only :         l_lh_importance_sampling

#endif
#endif

      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

#ifdef CLUBB_SGS
#ifdef SILHS

      integer :: iunit = 501 ! Default value, will get iunit from CAM 
      !character(len=*), parameter :: default_corr_case = "arm_97"
      character(len=*), parameter :: &
            cloud_file_ext  = "_corr_array_cloud.in", & ! File extensions for corr files
            below_file_ext  = "_corr_array_below.in"
      character(len=256) :: corr_file_path_cloud, corr_file_path_below

      ! Set CLUBB's debug level
      ! This is called in module clubb_intr; no need to do it here.
!      call set_clubb_debug_level_api( 0 )

      !-------------------------------
      ! CLUBB-SILHS Parameters (global module variables)
      !-------------------------------

      l_fix_chi_eta_correlations = .true.
      l_lh_importance_sampling = .true.
      l_diagnose_correlations = .false.
      l_calc_w_corr = .false.
!      l_prescribed_avg_deltaz = .false.
      l_use_cloud_cover = .false.
      l_const_Nc_in_cloud = .true.

      ! Values from the namelist
      docldfracscaling = subcol_SILHS_use_clear_col

      ! Namelist "tuning" or set correlations
      ! KTC Todo: Move these to a tuning "in" file or into the namelist
      ! JHTODO: we might want these on CLUBB's API and ultimatively on a namelist for tuning
!      C6rt = subcol_SILHS_c6rt
!      C7 = subcol_SILHS_c7                                      ! to all ice clouds
!      C8 = subcol_SILHS_c8
!      C11 = subcol_SILHS_c11
!      C11b = subcol_SILHS_c11b
!      gamma_coef = subcol_SILHS_gamma_coef
!      mult_coef = subcol_SILHS_mult_coef
!      mu = subcol_SILHS_mu

      !call set_clubb_debug_level( 0 )  !#KTCtodo: Add a namelist variable to set debug level
     
      !-------------------------------
      ! Define physics buffer indexes
      !-------------------------------
      thlm_idx = pbuf_get_index('THLM')   
      num_subcols_idx = pbuf_get_index('num_subcols')
      pdf_params_idx = pbuf_get_index('PDF_PARAMS')
      rcm_idx = pbuf_get_index('RCM')
      rtm_idx = pbuf_get_index('RTM')
      cld_idx = pbuf_get_index('CLD')
      alst_idx = pbuf_get_index('ALST')  ! SILHS expects clubb's cloud_frac liq stratus fraction
      wp2_idx = pbuf_get_index('WP2_nadv')
      ztodt_idx = pbuf_get_index('ZTODT')
      ice_supersat_idx = pbuf_get_index('ISS_FRAC')
     
      !-------------------------------
      ! Set up SILHS hydrometeors #KTCtodo: move microphys specification to config time,
      !        Steve wants to set up a microphysics query so I can ask the microphysics
      !        scheme which hydrometeors to use. For the future.
      !-------------------------------
      iirrm = 1
      iirsm = 3
      iirim  = 5
      iirgm = -1

      iiNrm    = 2
      iiNsm = 4
      iiNim    = 6
      iiNgm = -1

      hydromet_dim = 6

 
      allocate( hydromet_list(hydromet_dim) )
      allocate( hydromet_tol(hydromet_dim) )
      allocate( l_mix_rat_hm(hydromet_dim) )
      allocate( l_frozen_hm(hydromet_dim) )
      allocate( hmp2_ip_on_hmm2_ip(hydromet_dim) )

      if ( iirrm > 0 ) then
         ! The microphysics scheme predicts rain water mixing ratio, rr.
         hydromet_list(iirrm) = "rrm"
         l_mix_rat_hm(iirrm)  = .true.
         l_frozen_hm(iirrm)   = .false.
         hydromet_tol(iirrm)  = rr_tol
         hmp2_ip_on_hmm2_ip(iirrm) = hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip 
      endif
      if ( iirim > 0 ) then
         ! The microphysics scheme predicts ice mixing ratio, ri.
         hydromet_list(iirim) = "rim"
         l_mix_rat_hm(iirim)  = .true.
         l_frozen_hm(iirim)   = .true.
         hydromet_tol(iirim)  = ri_tol
         hmp2_ip_on_hmm2_ip(iirim) = hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip
      endif
      if ( iirsm > 0 ) then
         ! The microphysics scheme predicts snow mixing ratio, rs.
         hydromet_list(iirsm) = "rsm"
         l_mix_rat_hm(iirsm)  = .true.
         l_frozen_hm(iirsm)   = .true.
         hydromet_tol(iirsm)  = rs_tol
         hmp2_ip_on_hmm2_ip(iirsm) = hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip
      endif
      if ( iirgm > 0 ) then
         ! The microphysics scheme predicts graupel mixing ratio, rg.
         hydromet_list(iirgm) = "rgm"
         l_mix_rat_hm(iirgm)  = .true.
         l_frozen_hm(iirgm)   = .true.
         hydromet_tol(iirgm)  = rg_tol
         hmp2_ip_on_hmm2_ip(iirgm) = hmp2_ip_on_hmm2_ip_ratios%rgp2_ip_on_rgm2_ip 
      endif
      if ( iiNrm > 0 ) then
         ! The microphysics scheme predicts rain drop concentration, Nr.
         hydromet_list(iiNrm) = "Nrm"
         l_frozen_hm(iiNrm)   = .false.
         l_mix_rat_hm(iiNrm)  = .false.
         hydromet_tol(iiNrm)  = Nr_tol
         hmp2_ip_on_hmm2_ip(iiNrm) = hmp2_ip_on_hmm2_ip_ratios%Nrp2_ip_on_Nrm2_ip
      endif
      if ( iiNim > 0 ) then
         ! The microphysics scheme predicts ice concentration, Ni.
         hydromet_list(iiNim) = "Nim"
         l_mix_rat_hm(iiNim)  = .false.
         l_frozen_hm(iiNim)   = .true.
         hydromet_tol(iiNim)  = Ni_tol
         hmp2_ip_on_hmm2_ip(iiNim) = hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip
      endif
      if ( iiNsm > 0 ) then
         ! The microphysics scheme predicts snowflake concentration, Ns.
         hydromet_list(iiNsm) = "Nsm"
         l_mix_rat_hm(iiNsm)  = .false.
         l_frozen_hm(iiNsm)   = .true.
         hydromet_tol(iiNsm)  = Ns_tol
         hmp2_ip_on_hmm2_ip(iiNsm) = hmp2_ip_on_hmm2_ip_ratios%Nsp2_ip_on_Nsm2_ip
      endif
      if ( iiNgm > 0 ) then
         ! The microphysics scheme predicts graupel concentration, Ng.
         hydromet_list(iiNgm) = "Ngm"
         l_mix_rat_hm(iiNgm)  = .false.
         l_frozen_hm(iiNgm)   = .true.
         hydromet_tol(iiNgm)  = Ng_tol
         hmp2_ip_on_hmm2_ip(iiNgm) = hmp2_ip_on_hmm2_ip_ratios%Ngp2_ip_on_Ngm2_ip
      endif

      Ncnp2_on_Ncnm2 = subcol_SILHS_ncnp2_on_ncnm2

      !-------------------------------
      ! Set up hydrometeors and correlation arrays for SILHS
      !-------------------------------
      corr_file_path_cloud = trim( subcol_SILHS_corr_file_path )//trim( subcol_SILHS_corr_file_name )//cloud_file_ext
      corr_file_path_below = trim( subcol_SILHS_corr_file_path )//trim( subcol_SILHS_corr_file_name )//below_file_ext

      iunit = getunit()


      call setup_pdf_indices_api( hydromet_dim, iirrm, iiNrm, iirim, &
                              iiNim, iirsm, iiNsm, iirgm, iiNgm )
      call setup_corr_varnce_array_api( corr_file_path_cloud, corr_file_path_below, &
                                    iunit )
      call freeunit(iunit) 

#endif
      !-------------------------------
      ! Register output fields from SILHS
      ! #KTCtodo: Remove these from the default output list
      !-------------------------------
      call addfld('SILHS_NCLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm^-3', &
           'Subcolumn Cloud Number Concentration', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_NRAIN_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm^-3', &
           'Subcolumn Number Concentration of Rain from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_OMEGA_SCOL', (/'psubcols', 'ilev    '/), 'I', 'Pa/s', &
           'Subcolumn vertical pressure velocity', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RCM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Liquid Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RICLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Ice Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_NICLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Ice Number Conc from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RRAIN_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Precipitating Liquid Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RT_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg ', &
           'Subcolumn Total Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_THLM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'K', &
           'Subcolumn liquid water pot temperature', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_WEIGHT_SCOL', (/'psubcols'/), 'I', 'frac', &
           'Weights for each subcolumn', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_WM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm/s', &
           'Subcolumn vertical velocity from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)

      call addfld('NR_IN_LH', (/ 'lev' /), 'I', 'm^-3', &
                  'Num Rain Conc as input to SILHS')
     call addfld('RTM_CLUBB', (/ 'ilev' /), 'I', 'kg/kg', &
                  'Input total water mixing ratio')
     call addfld('THLM_CLUBB', (/ 'ilev' /), 'I', 'K', &
                  'Input liquid water potential temperature')
     call addfld('SILHS_QC_IN', (/ 'lev' /), 'I', 'kg/kg', &
                  'Input cloud water mixing ratio')
     call addfld('SILHS_QI_IN', (/ 'lev' /), 'I', 'kg/kg', &
                  'Input cloud ice mixing ratio')
     call addfld('SILHS_NC_IN', (/ 'lev' /), 'I', '#/kg', &
                  'Input cloud water number concentration')
     call addfld('SILHS_NI_IN', (/ 'lev' /), 'I', '#/kg', &
                  'Input cloud ice number concentration')
     call addfld('AKM_CLUBB', (/ 'ilev' /), 'I', '(kg/kg)/s', &
                  'Exact Kessler autoconversion')
     call addfld('AKM_LH_CLUBB', (/ 'ilev' /), 'I', '(kg/kg)/s', &
                  'Monte Carlo estimate of Kessler autoconversion')
     call addfld('INVS_EXNER', (/ 'lev' /), 'I', 'none', &
                  'inverse EXNER function from state in subcol_SILHS')
     call addfld('SILHS_ZTODT', horiz_only, 'I', 's', & 
                  'Length of Physics timestep (for debugging)')
     call addfld('SILHS_MSC_CLDICE', (/ 'lev' /), 'A', 'kg/kg', &
                  'Mean Cloud Ice across subcolumns')
     call addfld('SILHS_STDSC_CLDICE', (/ 'lev' /), 'A', 'kg/kg', &
                  'Standard deviation of Ice across subcolumns')
     call addfld('SILHS_MSC_CLDLIQ', (/ 'lev' /), 'A', 'kg/kg', &
                  'Mean Cloud Liquid across subcolumns')
     call addfld('SILHS_STDSC_CLDLIQ', (/ 'lev' /), 'A', 'kg/kg', &
                  'Standard deviation of Liquid across subcolumns')
     call addfld('SILHS_MSC_Q', (/ 'lev' /), 'A', 'kg/kg', &
                  'Mean water vapor across subcolumns')
     call addfld('SILHS_STDSC_Q', (/ 'lev' /), 'A', 'kg/kg', &
                  'Standard deviation of water vapor across subcolumns')
     call addfld('SILHS_EFF_CLDFRAC', (/ 'lev' /), 'A', 'frac', &
                  'Calculated cloud fraction from subcolumn liq or ice') 

     call addfld('SILHS_CLUBB_PRECIP_FRAC', (/ 'lev' /), 'A', 'frac', &
                  'Precipitation fraction from CLUBB (setup_pdf_params)')
     call addfld('SILHS_CLUBB_ICE_SS_FRAC', (/ 'lev' /), 'A', 'frac', &
                  'Ice supersaturation fraction from CLUBB')

      !call add_default('SILHS_NCLD_SCOL', 1, ' ')
      !call add_default('SILHS_NRAIN_SCOL', 1, ' ')
      !call add_default('SILHS_OMEGA_SCOL', 1, ' ')
      !call add_default('SILHS_RCM_SCOL', 1, ' ')
      !call add_default('SILHS_RICLD_SCOL', 1, ' ')
      !call add_default('SILHS_NICLD_SCOL', 1, ' ')
      !call add_default('SILHS_RRAIN_SCOL', 1, ' ')
      !call add_default('SILHS_RT_SCOL', 1, ' ')
      !call add_default('SILHS_THLM_SCOL', 1, ' ')
      !call add_default('SILHS_WEIGHT_SCOL', 1, ' ')
      !call add_default('SILHS_WM_SCOL', 1, ' ')

      !call add_default('NR_IN_LH', 1, ' ')
      !call add_default('RTM_CLUBB', 1, ' ')
      !call add_default('THLM_CLUBB', 1, ' ')
      !call add_default('SILHS_QC_IN', 1, ' ')
      !call add_default('SILHS_QI_IN', 1, ' ')
      !call add_default('AKM_CLUBB', 1, ' ')
      !call add_default('AKM_LH_CLUBB', 1, ' ')
      !call add_default('INVS_EXNER', 1, ' ')

#endif
   end subroutine subcol_init_SILHS
   
   subroutine subcol_gen_SILHS(state, tend, state_sc, tend_sc, pbuf)
      !-------------------------------
      ! This is where the subcolumns are created, and the call to
      !      lh_subcolumn_generator_mod_api
      !    goes out. Variables needed to make this call are pulled from the 
      !    pbuf, from module data, and calculated based on the CAM state.
      !-------------------------------

      use physics_buffer,         only : physics_buffer_desc, pbuf_get_index, &
                                         pbuf_get_field
      use ppgrid,                 only : pver, pverp, pcols
      use time_manager,           only : get_nstep
      use subcol_utils,           only : subcol_set_subcols, subcol_set_weight
      use phys_control,           only : phys_getopts
      use spmd_utils,             only : masterproc
      use shr_const_mod,          only : SHR_CONST_PI, SHR_CONST_RHOFW

#ifdef CLUBB_SGS
#ifdef SILHS
      use clubb_api_module,       only : pdf_parameter, unpack_pdf_params_api, &
                                         num_pdf_params, &
                                         hydromet_dim, &

                                         Lscale, wp2_zt, &
                                         wphydrometp, &

                                         setup_pdf_parameters_api, &

                                         l_stats_samp, &

                                         hydromet_pdf_parameter, &

                                         zm2zt_api, setup_grid_heights_api, &

                                         iirrm, iiNrm, iirsm, iirim, &
                                         iirgm, iiNsm, &
                                         iiNim, iiNgm, &

                                         core_rknd, &

                                         w_tol_sqd, zero_threshold, cloud_frac_min, & ! rc_tol, &

                                         d_variables, &
                                         corr_array_n_cloud, &
                                         corr_array_n_below, &
                                         iiPDF_chi, iiPDF_rr, &
                                         iiPDF_w, iiPDF_Nr, &
                                         iiPDF_ri, iiPDF_Ni, &
                                         iiPDF_Ncn, iiPDF_rs, iiPDF_Ns
   
      use silhs_api_module, only :       lh_subcolumn_generator_api, & ! Ncn_to_Nc, &
                                         lh_clipped_variables_type, &
                                         clip_transform_silhs_output_api, &
                                         est_kessler_microphys_api
#endif
#endif
      
      ! CAM data structures
      type(physics_state), intent(inout) :: state
      type(physics_tend),  intent(inout) :: tend
      type(physics_state), intent(inout) :: state_sc        ! sub-column state
      type(physics_tend),  intent(inout) :: tend_sc         ! sub-column tend
      type(physics_buffer_desc), pointer :: pbuf(:)

#ifdef CLUBB_SGS
#ifdef SILHS
      !----------------
      ! Local variables
      !----------------
      logical, parameter :: &
                 l_implemented = .true.   ! Implemented in a host model
      logical, parameter :: rx_Nc = .false. ! Use NC calculated based on grid mean effective radius
      integer, parameter :: &
                 grid_type = 3            ! The 3 grid centered on momentum levels
      real(r8), parameter :: cldmin = 0.001_r8 ! To use when cld frac = 0.0 to be consistant with micro_mg
      real(r8), parameter :: min_num_conc = 1.0e-12_r8
      real(r8), parameter :: qsmall = 1.0e-18  ! Microphysics cut-off for cloud

      integer :: i, j, k, ngrdcol, ncol, lchnk, stncol
      integer :: ixcldice, ixnumice, ixq, ixcldliq, ixnumliq, ixrain, ixnumrain, ixsnow, ixnumsnow
      integer :: begin_height, end_height ! Output from setup_grid call
      real(r8) :: sfc_elevation  ! Surface elevation
      real(r8), dimension(pverp) :: zt_g, zi_g ! Thermo & Momentum grids for clubb
      real(r8), dimension(pverp) :: scfrac     ! cloud fraction based on sc distributions
      real(r8) :: msc, std, maxcldfrac, maxsccldfrac
      real(r8) :: scale = 1.0_r8

      !----------------
      ! Required for setup_pdf_params
      !----------------
      real(r8), dimension(pverp) :: cld_frac_in  ! Cloud fraction
      real(r8), dimension(pverp) :: wp2_flip     ! CLUBB vert vel var flipped from pbuf
      type(hydromet_pdf_parameter), dimension(pverp) :: &
                                    hydromet_pdf_params  ! Hydrometeor PDF parameters
      real(r8), dimension(:,:,:), allocatable :: &       ! Correlation matrix for pdf components
                                    corr_array_1, corr_array_2 
      real(r8), dimension(:,:), allocatable :: &
                                    mu_x_1, mu_x_2, &    ! Mean array for PDF components
                                    sigma_x_1, sigma_x_2 ! Std dev arr for PDF components
      real(r8), dimension(:,:,:), allocatable :: &       ! Transposed corr cholesky mtx
                                    corr_cholesky_mtx_1, corr_cholesky_mtx_2
      real(r8), dimension(pverp) :: Nc_in_cloud
      real(r8), dimension(pverp) :: ice_supersat_frac_in
      real(r8), dimension(pverp,hydromet_dim) :: hydrometp2


      !----------------
      ! Input to lh_subcolumn_generator
      !----------------
      integer :: iter                            ! CLUBB iteration 
      integer :: num_subcols                     ! Number of subcolumns
      integer, dimension(pcols) :: numsubcol_arr ! To set up the state struct
      integer, parameter :: sequence_length = 1  ! Number of timesteps btn subcol calls
      type(pdf_parameter), dimension(pverp) :: pdf_params     ! PDF parameters
      real(r8), dimension(pverp, num_pdf_params) :: pdf_params_packed
      real(r8), dimension(pverp) :: rho_ds_zt    ! Dry static density (kg/m^3) on thermo levs
      real(r8), dimension(pver)  :: dz_g         ! thickness of layer
      real(r8), dimension(pverp) :: delta_zm     ! Difference in u wind altitudes
      real(r8), dimension(pverp) :: invs_dzm     ! 1/delta_zm
      real(r8), dimension(pverp) :: rcm_in       ! Cld water mixing ratio on CLUBB levs
      real(r8), dimension(pverp,hydromet_dim) :: hydromet  ! Hydrometeor spcies
      real(r8), dimension(pverp)              :: Ncm ! Mean cloud droplet concentration, <N_c>
      

      !---------------
      !Output from lh_subcolumn_generator
      !--------------
      real(r8), allocatable, dimension(:,:,:) :: X_nl_all_levs ! Sample transformed to normal-lognormal
      real(r8), allocatable, dimension(:) :: LH_sample_point_weights ! Subcolumn weights
      integer, allocatable, dimension(:,:) :: X_mixt_comp_all_levs ! Which Mixture Component

      real(r8), allocatable, dimension(:,:) :: rc_all_points ! Calculate RCM from LH output
      real(r8), allocatable, dimension(:,:) :: rain_all_pts  ! Calculate Rain from LH output
      real(r8), allocatable, dimension(:,:) :: nrain_all_pts ! Calculate Rain Conc from LH
      real(r8), allocatable, dimension(:,:) :: snow_all_pts  ! Calculate Snow from LH output
      real(r8), allocatable, dimension(:,:) :: nsnow_all_pts ! Calculate Snow Conc from LH
      real(r8), allocatable, dimension(:,:) :: w_all_points  ! Calculate W from LH output
      ! real(r8), allocatable, dimension(:,:) :: RVM_LH_out    ! Vapor mixing ratio sent away
      real(r8), allocatable, dimension(:,:) :: ice_all_pts   ! Calculate Cld Ice from LH output
      real(r8), allocatable, dimension(:,:) :: nice_all_pts  ! Calculate Num cld ice from LH
      real(r8), allocatable, dimension(:,:) :: nclw_all_pts  ! Calculate Num cld wat from LH

      !----------------
      ! Output from clip_transform_silhs_output_api
      !----------------
      type(lh_clipped_variables_type), dimension(:,:), allocatable :: &
        lh_clipped_vars

      logical, parameter :: &
        l_use_Ncn_to_Nc = .true.  ! Whether to call Ncn_to_Nc (.true.) or not (.false.);
                                  ! Ncn_to_Nc might cause problems with the MG microphysics 
                                  ! since the changes made here (Nc-tendency) are not fed into 
                                  ! the microphysics
        

      !----------------
      ! Output to history
      !----------------
      real(r8), dimension(pcols*psubcols, pverp) :: RT_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: THL_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: OMEGA_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: WM_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: RVM_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: RCM_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: NCLW_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: ICE_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: NICE_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: RAIN_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: NRAIN_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: SNOW_LH_out
      real(r8), dimension(pcols*psubcols, pverp) :: NSNOW_LH_out

      real(r8), dimension(state_sc%psetcols) :: weights ! Subcol weights

      real(r8), dimension(pcols, pver) :: meansc_ice
      real(r8), dimension(pcols, pver) :: stdsc_ice

      real(r8), dimension(pcols, pver) :: meansc_liq
      real(r8), dimension(pcols, pver) :: stdsc_liq

      real(r8), dimension(pcols, pver) :: meansc_vap
      real(r8), dimension(pcols, pver) :: stdsc_vap
      real(r8), dimension(pcols, pver) :: grmn_eff_rad
      real(r8), dimension(pcols, pver) :: eff_cldfrac
      real(r8), dimension(pcols, pver) :: precip_frac_out

      real(r8) :: tmp_mean, diff_mean, rcubed

      !----------------
      ! Output from Est_Kessler_microphys
      !----------------
      real(r8), dimension(pverp) :: lh_Akm     ! Monte Carlo estimate of Kessler Autoconversion
      real(r8), dimension(pverp) :: AKm        ! Exact Kessler autoconversion
      real(r8), dimension(pverp) :: AKstd      ! Exact Stdev of gba Kessler
      real(r8), dimension(pverp) :: AKstd_cld  ! Exact w/in cloud stdev of gba Kessler
      real(r8), dimension(pverp) :: AKm_rcm    ! Exact local gba Kessler auto based on rcm
      real(r8), dimension(pverp) :: AKm_rcc    ! Exact local gba Kessler based on w/in cloud rc
      real(r8), dimension(pverp) :: LH_rcm_avg ! LH estimate of grid box avg liquid water
      real(r8), dimension(pcols,pverp) :: lh_AKm_out, AKm_out

      !----------------
      ! Needed to update State
      !----------------
      real(r8), dimension(pver)  :: Temp_prof  ! Subcolumn LWPT converted to Abs Temp
      real(r8), dimension(pver)  :: SE_prof    ! Static Energy calculated from Abs Temp
      real(r8), dimension(pver)  :: No_cloud = 0.0_r8     ! Clear air condensate profile
      real(r8), dimension(pcols, pver)  :: invs_exner  ! inverse exner sent to conversion codw
                                                       ! pcols for output to history
      real(r8) :: eff_rad_coef = 1.0_r8/(4.0_r8/3.0_r8*SHR_CONST_RHOFW*SHR_CONST_PI)
      real(r8), dimension(pver) :: eff_rad_prof ! r^3 as calculated from grid mean MR & NC
     
      !----------------
      ! Pointers
      !----------------
      real(r8), pointer, dimension(:) :: num_subcol_ptr
      real(r8), pointer, dimension(:) :: ztodt_ptr
      real(r8), pointer, dimension(:,:) :: thlm      ! Mean temperature
      real(r8), pointer, dimension(:,:,:) :: pdf_params_ptr  ! Packed PDF_Params
      real(r8), pointer, dimension(:,:) :: ice_supersat_frac ! ice cloud fraction
      real(r8), pointer, dimension(:,:) :: rcm       ! CLUBB cld water mr
      real(r8), pointer, dimension(:,:) :: rtm       ! mean moisture mixing ratio
      real(r8), pointer, dimension(:,:) :: cld       ! CAM cloud fraction
      real(r8), pointer, dimension(:,:) :: alst      ! CLUBB liq cloud fraction
      real(r8), pointer, dimension(:,:) :: wp2       ! CLUBB vert vel variance
      real(r8), pointer, dimension(:,:) :: qrain     ! micro_mg rain from previous step
      real(r8), pointer, dimension(:,:) :: qsnow     
      real(r8), pointer, dimension(:,:) :: nrain     ! micro_mg rain num conc 
      real(r8), pointer, dimension(:,:) :: nsnow


      if (.not. allocated(state_sc%lat)) then
         call endrun('subcol_gen error: state_sc must be allocated before calling subcol_gen')
      end if

      ! Determine num of columns and which chunk we're working on and what timestep
      ngrdcol = state%ngrdcol
      ncol = state%ncol
      lchnk = state%lchnk
      iter = get_nstep() ! #KTCtodo: The model iteration is passed into SILHS without taking
                         !           substepping into account. I may need to change this in 
                         !           the future. Also, why does SILHS need an iter, but CLUBB
                         !           does not?
                         ! #ERDBG:   The model iteration number is not used in SILHS unless
                         !           sequence_length > 1, but nobody runs with that option.
      !----------------
      ! Establish associations between pointers and physics buffer fields
      ! (Do this now so that num_subcol_ptr is available for the state copy below)
      !----------------
      call pbuf_get_field(pbuf, thlm_idx, thlm)
      call pbuf_get_field(pbuf, num_subcols_idx, num_subcol_ptr)
      call pbuf_get_field(pbuf, ztodt_idx, ztodt_ptr)
      call pbuf_get_field(pbuf, pdf_params_idx, pdf_params_ptr)
      call pbuf_get_field(pbuf, ice_supersat_idx, ice_supersat_frac)
      call pbuf_get_field(pbuf, rcm_idx, rcm)
      call pbuf_get_field(pbuf, rtm_idx, rtm)
      call pbuf_get_field(pbuf, alst_idx, alst)
      call pbuf_get_field(pbuf, cld_idx, cld)
      call pbuf_get_field(pbuf, wp2_idx, wp2)
      call pbuf_get_field(pbuf, qrain_idx, qrain)
      call pbuf_get_field(pbuf, qsnow_idx, qsnow)
      call pbuf_get_field(pbuf, nrain_idx, nrain)
      call pbuf_get_field(pbuf, nsnow_idx, nsnow)

      !----------------
      ! Copy state and populate numbers and values of sub-columns
      !----------------
      ztodt = ztodt_ptr(1)
      numsubcol_arr(:) = 0  ! Start over each chunk
      numsubcol_arr(:ngrdcol) = subcol_SILHS_numsubcol ! Only set for valid grid columns
      call subcol_set_subcols(state, tend, numsubcol_arr, state_sc, tend_sc)

      !----------------
      ! Get indices for ice mass and number
      ! This is the same code from clubb_intr.F90
      !----------------
      call cnst_get_ind('Q', ixq)
      call cnst_get_ind('CLDICE', ixcldice)
      call cnst_get_ind('NUMICE', ixnumice)
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('NUMLIQ', ixnumliq)
      call cnst_get_ind('RAINQM', ixrain, abort=.false.)
      call cnst_get_ind('NUMRAI', ixnumrain, abort=.false.)
      call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)
      call cnst_get_ind('NUMSNO', ixnumsnow, abort=.false.)

      !----------------
      ! Loop over all the active grid columns in the chunk
      !----------------
      do i=1,ngrdcol
      
         ! JHDBG: Big suspicion about that code
         num_subcols = numsubcol_arr(i)
         stncol = 0         ! Each grid column needs to know how many subcolumns have gone by
         do k=1,i-1
            ! stncol = stncol + numsubcol_arr(i-1)
            ! Eric Raut replaced i-1 with k in line immediately above.
            stncol = stncol + numsubcol_arr(k)
         enddo

         ! Setup the CLUBB vertical grid object. This must be done for each
         ! column as the z-distance between hybrid pressure levels can 
         ! change easily. This code is taken directly from clubb_intr.F90.
         sfc_elevation = state%zi(i,pver+1)
         ! Define the CLUBB momentum grid (in height, units of m)
         do k=1,pverp
            zi_g(k) = state%zi(i,pverp-k+1)-sfc_elevation
         enddo
         ! Define the CLUBB thermodynamic grid (in units of m)
         do k=1,pver
            zt_g(k+1) = state%zm(i,pver-k+1)-state%zi(i,pver+1)
            dz_g(k) = state%zi(i,k)-state%zi(i,k+1) ! compute thickness for SILHS
         enddo
         ! Thermodynamic ghost point is below surface
         zt_g(1) = -1._r8*zt_g(2)
         ! allocate grid object
         call setup_grid_heights_api(l_implemented, grid_type, &
                         zi_g(2), zi_g(1), zi_g(1:pverp), &
                         zt_g(1:pverp) )

         ! Pull pdf params out of 2-D real array
         call unpack_pdf_params_api(pdf_params_ptr(i,:,:), pverp, pdf_params)

         ! Inverse delta_zm is required for the 3-level L-scale averaging
         do k=1,pverp-1
            delta_zm(k+1) = state%zi(i,pverp-k)-state%zi(i,pverp-k+1)
            invs_dzm(k+1) = 1.0_r8/delta_zm(k+1)
         enddo
         ! Handle CLUBB sub-sfc ghost point as done in clubb grid_class.F90
         delta_zm(1) = delta_zm(2) 
         invs_dzm(1) = invs_dzm(2)

         ! Compute dry static density on CLUBB vertical grid
         do k=1, pver
            rho_ds_zt(k+1) = (1._r8/gravit)*state%pdel(i,pver-k+1)/dz_g(pver-k+1)
         enddo
         ! CLUBB ghost point under the surface
         rho_ds_zt(1) = rho_ds_zt(2)

         ! Set up hydromet array, flipped from CAM vert grid to CLUBB
         do k=1,pver  
            if ( iirrm > 0 ) then
              ! If ixrain and family are greater than zero, then MG2 is
              ! being used, and rain and snow are part of state. Otherwise,
              ! diagnostic rain and snow from MG1 are used in hydromet.
              if (ixrain > 0) then
                hydromet(k+1,iirrm) = state%q(i,pver-k+1,ixrain)
              else
                hydromet(k+1,iirrm) = qrain(i,pver-k+1)
              end if
            endif
            if ( iiNrm > 0 ) then
              if (ixnumrain > 0) then
                hydromet(k+1,iiNrm) = state%q(i,pver-k+1,ixnumrain)
              else
                hydromet(k+1,iiNrm) = nrain(i,pver-k+1)
              end if
            endif
            if ( iirsm > 0 ) then
              if (ixsnow > 0) then
                hydromet(k+1,iirsm) = state%q(i,pver-k+1,ixsnow)
              else
                hydromet(k+1,iirsm) = qsnow(i,pver-k+1)
              end if
            endif
            if ( iiNsm > 0 ) then
              if (ixnumsnow > 0) then
                hydromet(k+1,iiNsm) = state%q(i,pver-k+1,ixnumsnow)
              else
                hydromet(k+1,iiNsm) = nsnow(i,pver-k+1)
              end if
            endif
            if ( iirim > 0 ) then
               hydromet(k+1,iirim) = state%q(i,pver-k+1,ixcldice)
            endif
            if ( iiNim > 0 ) then
               hydromet(k+1,iiNim) = state%q(i,pver-k+1,ixnumice)
            endif
     
            Ncm(k+1)            = state%q(i,pver-k+1,ixnumliq)

            ! Calculate effective radius cubed, CAM-grid oriented for use in subcolumns
            eff_rad_prof(k) = eff_rad_coef*state%q(i,k,ixcldliq)/state%q(i,k,ixnumliq)
            ! Test a fixed effective radius
            ! eff_rad_prof(k) = 5.12e-16_r8 ! 8 microns
         enddo

         Ncm(1) = Ncm(2)

         do k=1,hydromet_dim ! ghost point below the surface
            hydromet(1,k) = hydromet(2,k)                  
         enddo

         ! Allocate arrays for setup_pdf_params
         allocate( corr_array_1(d_variables, d_variables, pverp) )
         allocate( corr_array_2(d_variables, d_variables, pverp) )
         allocate( mu_x_1(d_variables, pverp) )
         allocate( mu_x_2(d_variables, pverp) )
         allocate( sigma_x_1(d_variables, pverp) )
         allocate( sigma_x_2(d_variables, pverp) )
         allocate( corr_cholesky_mtx_1(d_variables, d_variables, pverp) )
         allocate( corr_cholesky_mtx_2(d_variables, d_variables, pverp) )
         ! Allocate arrays for SILHS output
         allocate( LH_sample_point_weights(num_subcols) )
         allocate( X_mixt_comp_all_levs(pverp,num_subcols) )
         allocate( X_nl_all_levs(pverp,num_subcols,d_variables) )
         allocate( lh_clipped_vars(pverp,num_subcols) )
         ! Allocate arrays for output to either history files or for updating state_sc
         allocate( rc_all_points(pverp, num_subcols) )
         allocate( rain_all_pts(pverp, num_subcols) )
         allocate( nrain_all_pts(pverp, num_subcols) )
         allocate( snow_all_pts(pverp, num_subcols) )
         allocate( nsnow_all_pts(pverp, num_subcols) )
         allocate( w_all_points(pverp, num_subcols) )
         ! allocate( RVM_LH_out(num_subcols, pverp) )  ! This one used only to update state
         allocate( ice_all_pts(pverp, num_subcols) )
         allocate( nice_all_pts(pverp, num_subcols) )
         allocate( nclw_all_pts(pverp, num_subcols) )
         
         ! Convert from CAM vertical grid to CLUBB
         do k=1,pverp 
            rcm_in(k)  = rcm(i,pverp-k+1)
            wp2_flip(k) = wp2(i,pverp-k+1)
            ice_supersat_frac_in(k) = ice_supersat_frac(i,pverp-k+1)
         enddo
         do k=1,pver
            cld_frac_in(k+1) = alst(i,pver-k+1)
         enddo
         cld_frac_in(1) = cld_frac_in(2) ! Ghost pt below surface
         ! Calculate a clubb-specific exner function
         ! (This is grid mean, as pressure levels do not change in 
         !  the subcolumn state)
         invs_exner(i,:) = ((state%pmid(i,:)/p0_clubb)**(rair/cpair))

         ! Call setup_pdf_parameters to get the CLUBB PDF ready for SILHS
         ! Compute Num concentration of cloud nuclei
         Nc_in_cloud = Ncm / max( cld_frac_in, cloud_frac_min )
     
         ! interpolate wp2 to midpoint levels
         wp2_zt  = max( zm2zt_api( wp2_flip ), w_tol_sqd)
         
         ! make the call
         call setup_pdf_parameters_api(pverp, d_variables, ztodt, &    ! In
                                   Nc_in_cloud, rcm_in, cld_frac_in, &            ! In
                                   ice_supersat_frac_in, hydromet, wphydrometp, & ! In
                                   corr_array_n_cloud, corr_array_n_below, &          ! In
                                   pdf_params, l_stats_samp, &                    ! In
                                   hydrometp2, &                                  ! Out
                                   mu_x_1, mu_x_2, &                              ! Out
                                   sigma_x_1, sigma_x_2, &                        ! Out
                                   corr_array_1, corr_array_2, &                  ! Out
                                   corr_cholesky_mtx_1, corr_cholesky_mtx_2, &    ! Out
                                   hydromet_pdf_params )                          ! Out

         ! Let's generate some subcolumns!!!!!
         call lh_subcolumn_generator_api &
              (iter, d_variables, num_subcols, sequence_length, & ! In
              pverp, pdf_params, delta_zm, rcm_in, Lscale, &      ! In
              rho_ds_zt, mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &  ! In 
              corr_cholesky_mtx_1, corr_cholesky_mtx_2, &         ! In
              hydromet_pdf_params, &                              ! In
              X_nl_all_levs, X_mixt_comp_all_levs, &              ! Out
              LH_sample_point_weights)                            ! Out

         ! Extract clipped variables from subcolumns
         call clip_transform_silhs_output_api( pverp, num_subcols, d_variables, &     ! In
                                               X_mixt_comp_all_levs, X_nl_all_levs, & ! In
                                               pdf_params, l_use_Ncn_to_Nc, &         ! In
                                               lh_clipped_vars )                      ! Out

         ! Test subcolumns by comparing to an estimate of kessler autoconversion
         call est_kessler_microphys_api &
              (pverp, num_subcols, d_variables, X_nl_all_levs, pdf_params, &
              rcm_in, cld_frac_in, X_mixt_comp_all_levs, LH_sample_point_weights, &
              lh_AKm, AKm, AKstd, AKstd_cld, AKm_rcm, AKm_rcc, LH_rcm_avg)

         ! Calc column liquid water for output (rcm)
         rc_all_points = lh_clipped_vars(:,:)%rc
         ! Calc subcolumn precipitating liq water for output (rrm)
         rain_all_pts = real( X_nl_all_levs(:,:,iiPDF_rr), kind=r8 )
         ! Calc subcolumn number rain conc for output (nrainm)
         nrain_all_pts = real( X_nl_all_levs(:,:,iiPDF_Nr), kind=r8 )
         ! Calc subcolumn precipitating snow      for output (rsm)
         snow_all_pts = real( X_nl_all_levs(:,:,iiPDF_rs), kind=r8 )
         ! Calc subcolumn precipitating snow conc for output (Nsm)
         nsnow_all_pts = real( X_nl_all_levs(:,:,iiPDF_Ns), kind=r8 )
         ! Calc subcolumn cloud ice mixing ratio
         ice_all_pts = real( X_nl_all_levs(:,:,iiPDF_ri), kind=r8)
         ! Calc subcolumn cloud ice number
         nice_all_pts = real( X_nl_all_levs(:,:,iiPDF_Ni), kind=r8)
         ! Calc subcolumn vert velocity for output (wm)
         w_all_points = real( X_nl_all_levs(:,:,iiPDF_w), kind=r8 )
         ! Calc cloud liq water number conc 
         nclw_all_pts = lh_clipped_vars(:,:)%Nc
         ! Calc mean liquid water potential temp for clear air
         !call THL_profile(pver, state%t(i,:), invs_exner(i,:), No_cloud, Temp_prof)

         ! Calc effective cloud fraction for testing
         eff_cldfrac(:,:) = 0.0_r8
         do k=1,pver
           do j=1, num_subcols

              if((rc_all_points(pverp-k+1,j).gt.qsmall) &
                  .or.(ice_all_pts(pverp-k+1,j).gt.qsmall)) then
                  eff_cldfrac(i,k) = eff_cldfrac(i,k)+LH_sample_point_weights(j)
              endif
           enddo 

           eff_cldfrac(i,k) = eff_cldfrac(i,k)/real(num_subcols, kind=r8)
         enddo

         ! Pack precip_frac for output
         do k=2, pverp
           precip_frac_out(i,pver-k+2) = hydromet_pdf_params(k)%precip_frac
         end do

         ! Pack up weights for output
         do j = 1, num_subcols      
               if (subcol_SILHS_weight) then 
                  weights(stncol+j) = LH_sample_point_weights(j)
               else
                  weights(stncol+j) = 1._r8
               endif
         enddo

         ! Convert from CLUBB vertical grid to CAM grid for history output and
         ! Updating state variables
         do k=1,pverp
            do j=1,num_subcols
               RT_LH_out(    stncol+j,k ) = lh_clipped_vars(pverp-k+1,j)%rt
               RCM_LH_out(   stncol+j,k ) = rc_all_points(pverp-k+1,j)
               NCLW_LH_out(  stncol+j,k ) = nclw_all_pts(pverp-k+1,j)
               ICE_LH_out(   stncol+j,k) = ice_all_pts(pverp-k+1,j)
               NICE_LH_out(  stncol+j,k) = nice_all_pts(pverp-k+1,j)
!               RVM_LH_out(j,k) = RT_LH_out(stncol+j,k)-RCM_LH_out(stncol+j,k)-ICE_LH_out(stncol+j,k)
               RVM_LH_out(stncol+j,k) = lh_clipped_vars(pverp-k+1,j)%rv
               THL_LH_out(   stncol+j,k ) = lh_clipped_vars(pverp-k+1,j)%thl
               RAIN_LH_out(  stncol+j,k ) = rain_all_pts(pverp-k+1,j)
               NRAIN_LH_out( stncol+j,k ) = nrain_all_pts(pverp-k+1,j)
               SNOW_LH_out(  stncol+j,k ) = snow_all_pts(pverp-k+1,j)
               NSNOW_LH_out( stncol+j,k ) = nsnow_all_pts(pverp-k+1,j)
               WM_LH_out(    stncol+j,k ) = w_all_points(pverp-k+1,j)
               OMEGA_LH_out( stncol+j,k ) = -1._r8*WM_LH_out(stncol+j,k)*rho_ds_zt(pverp-k+1)*gravit
               AKm_out(i,k) = AKm(pverp-k+1)
               lh_AKm_out(i,k) = lh_AKm(pverp-k+1)
            enddo
         enddo

         ! Constrain the sample distribution of cloud water and ice to the same mean
         ! as the grid to prevent negative condensate errors
         if(subcol_SILHS_constrainmn) then
            call subcol_constrainmn( num_subcols, ICE_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixcldice), meansc_ice(i,:), stdsc_ice(i,:) )
            if ( ixrain > 0 ) &
            call subcol_constrainmn( num_subcols, RAIN_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixrain) )
            if ( ixsnow > 0 ) &
            call subcol_constrainmn( num_subcols, SNOW_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixsnow) )
            call subcol_constrainmn( num_subcols, RCM_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixcldliq), meansc_liq(i,:), stdsc_liq(i,:) )
            call subcol_constrainmn( num_subcols, RVM_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixq), meansc_vap(i,:), stdsc_vap(i,:) )
            call subcol_constrainmn( num_subcols, NICE_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumice) )
            if ( ixnumrain > 0 ) &
            call subcol_constrainmn( num_subcols, NRAIN_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumrain) )
            if ( ixnumsnow > 0 ) &
            call subcol_constrainmn( num_subcols, NSNOW_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumsnow) )
            call subcol_constrainmn( num_subcols, NCLW_LH_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumliq) )
            do k=1,pver
               ! Look for exceptionally large values of condensate
               if(ANY(ICE_LH_out(stncol+1:stncol+num_subcols,k) .gt. 0.01_r8)) then
                  ! Clip the large values
                  where(ICE_LH_out(stncol+1:stncol+num_subcols,k) .gt. 0.01_r8)
                     ICE_LH_out(stncol+1:stncol+num_subcols,k) = 0.01_r8
                     NICE_LH_out(stncol+1:stncol+num_subcols,k) = 1.5e+7_r8
                  end where
                  ! Recalculate the weighted subcolumn mean
                  tmp_mean = meansc( ICE_LH_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  ! Calculate the difference between the weighted mean and grid mean
                  diff_mean = state%q(i,k,ixcldice)-tmp_mean
                  ! Add the difference to each subcolumn
                  ICE_LH_out(stncol+1:stncol+num_subcols,k) = &
                     ICE_LH_out(stncol+1:stncol+num_subcols,k)+diff_mean
                  ! Recalculate the weight subcolumn mean for ice num conc
                  tmp_mean = meansc( NICE_LH_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  ! Calculate the difference between the weighted mean and grid mean
                  diff_mean = state%q(i,k,ixnumice)-tmp_mean
                  ! Add the difference to each subcolumn
                  if(diff_mean.gt.0.0_r8) then
                     NICE_LH_out(stncol+1:stncol+num_subcols,k) = &
                         NICE_LH_out(stncol+1:stncol+num_subcols,k)+diff_mean
                  else ! just use the grid mean in each subcolumn
                     NICE_LH_out(stncol+1:stncol+num_subcols,k) = &
                         state%q(i,k,ixnumice)
                  end if
                  ! Test adjusted means for debugging
                  tmp_mean = meansc( ICE_LH_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  diff_mean = state%q(i,k,ixcldice)-tmp_mean
                  tmp_mean = meansc( NICE_LH_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  diff_mean = state%q(i,k,ixnumice)-tmp_mean
               endif
            enddo ! k = 1,pver
         endif ! subcol_silhs_constrainm

         ! Code to update the state variables for interactive runs
         ! Set state variables
         do j=1,numsubcol_arr(i)

            call Abs_Temp_profile(pver, THL_LH_out(stncol+j,1:pver), &
                              invs_exner(i,:), RCM_LH_out(stncol+j,1:pver), Temp_prof)
            state_sc%t(stncol+j,:) = Temp_prof(:)
            call StaticEng_profile(pver, Temp_prof, state%zm(i,:), state%phis(i), &
                                      SE_prof)
            state_sc%s(stncol+j,:) = SE_prof(:)

            ! Vertical Velocity is not part of the energy conservation checks, but
            ! we need to be careful here, because the SILHS output VV is noisy.
            state_sc%omega(stncol+j,:) = OMEGA_LH_out(stncol+j,1:pver)
            state_sc%q(stncol+j,:,ixq) = RVM_LH_out(stncol+j,1:pver) 

            if( rx_Nc ) then
                stop "rx_Nc not enabled in subcol_gen_SILHS"
!               ! Test calculating num const based on grid mean eff radius
!                      where(eff_rad_prof.gt.0.0)
!                state_sc%q(stncol+j,:,ixnumliq) = (RCM_ADJ_out(stncol+j,1:pver) &
!                                                         /eff_rad_prof(:))*eff_rad_coef
!                      elsewhere
!                        state_sc%q(stncol+j,:,ixnumliq) = NCLW_ADJ_out(stncol+j,1:pver)
!                      end where
!                      NCLW_LH_out(stncol+j,1:pver) = state_sc%q(stncol+j,:,ixnumliq) 
!           else
!                       state_sc%q(stncol+j,:,ixnumliq) = NCLW_ADJ_out(stncol+j,1:pver)
!                       NCLW_LH_out(stncol+j,:) = NCLW_ADJ_out(stncol+j,:)
            endif


            if(subcol_SILHS_meanice) then
                stop "subcol_SILHS_meanice = T not currently available"
                state_sc%q(stncol+j,:,ixcldice) = state%q(i,:,ixcldice)
                state_sc%q(stncol+j,:,ixnumice) = state%q(i,:,ixnumice)
                state_sc%q(stncol+j,:,ixcldliq) = RCM_LH_out(stncol+j,1:pver)
                state_sc%q(stncol+j,:,ixnumliq) = NCLW_LH_out(stncol+j,1:pver)
            else
               if(subcol_SILHS_q_to_micro) then ! Send SILHS predicted constituents to microp
                  state_sc%q(stncol+j,:,ixcldliq) = RCM_LH_out(stncol+j,1:pver)
                  state_sc%q(stncol+j,:,ixcldice) = ICE_LH_out(stncol+j,1:pver)
                  if (ixrain > 0) &
                    state_sc%q(stncol+j,:,ixrain) = RAIN_LH_out(stncol+j,1:pver)
                  if (ixsnow > 0) &
                    state_sc%q(stncol+j,:,ixsnow) = SNOW_LH_out(stncol+j,1:pver)
               else            
                  state_sc%q(stncol+j,:,ixcldliq) = state%q(i,:,ixcldliq)
                  state_sc%q(stncol+j,:,ixcldice) = state%q(i,:,ixcldice)
                  if (ixrain > 0) &
                    state_sc%q(stncol+j,:,ixrain) = state%q(i,:,ixrain)
                  if (ixsnow > 0) &
                    state_sc%q(stncol+j,:,ixsnow) = state%q(i,:,ixsnow)
               endif
               if(subcol_SILHS_n_to_micro) then ! Send SILHS predicted number conc to microp
                  state_sc%q(stncol+j,:,ixnumice) = NICE_LH_out(stncol+j,1:pver)
                  state_sc%q(stncol+j,:,ixnumliq) = NCLW_LH_out(stncol+j,1:pver)
                  if (ixnumrain > 0) &
                    state_sc%q(stncol+j,:,ixnumrain) = NRAIN_LH_out(stncol+j,1:pver)
                  if (ixnumsnow > 0) &
                    state_sc%q(stncol+j,:,ixnumsnow) = NSNOW_LH_out(stncol+j,1:pver)
               else            
                  state_sc%q(stncol+j,:,ixnumliq) = state%q(i,:,ixnumliq)
                  state_sc%q(stncol+j,:,ixnumice) = state%q(i,:,ixnumice)
                  if (ixnumrain > 0) &
                    state_sc%q(stncol+j,:,ixnumrain) = state%q(i,:,ixnumrain)
                  if (ixnumsnow > 0) &
                    state_sc%q(stncol+j,:,ixnumsnow) = state%q(i,:,ixnumsnow)
               endif
            endif ! meanice

            ! Change liq and ice (and rain and snow) num conc zeros to min values (1e-12)
            where(state_sc%q(stncol+j,:,ixnumliq).lt.min_num_conc) 
               state_sc%q(stncol+j,:,ixnumliq)=min_num_conc
            end where
            where(state_sc%q(stncol+j,:,ixnumice).lt.min_num_conc)
               state_sc%q(stncol+j,:,ixnumice)=min_num_conc
            end where
            if (ixnumrain > 0) then
               where(state_sc%q(stncol+j,:,ixnumrain).lt.min_num_conc)
                  state_sc%q(stncol+j,:,ixnumrain)=min_num_conc
               end where
            end if
            if (ixnumsnow > 0) then
               where(state_sc%q(stncol+j,:,ixnumsnow).lt.min_num_conc)
                  state_sc%q(stncol+j,:,ixnumsnow)=min_num_conc
               end where
            end if
               
         enddo

         ! Only use weights if namelist variable turned on
         if (subcol_SILHS_weight) call subcol_set_weight(state_sc%lchnk, weights)

     
         ! Deallocate the dynamic arrays used
         deallocate( LH_sample_point_weights, X_mixt_comp_all_levs, X_nl_all_levs, &
                     lh_clipped_vars, corr_array_1, corr_array_2, mu_x_1, mu_x_2, &
                     sigma_x_1, sigma_x_2, corr_cholesky_mtx_1, corr_cholesky_mtx_2 )
         ! deallocate( RVM_LH_out ) 
         deallocate( rc_all_points, rain_all_pts, nrain_all_pts, snow_all_pts, nsnow_all_pts, ice_all_pts, &
                     nice_all_pts, nclw_all_pts, w_all_points )
      enddo ! ngrdcol

      call outfld( 'SILHS_THLM_SCOL', THL_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RT_SCOL', RT_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_OMEGA_SCOL', OMEGA_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_WM_SCOL', WM_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RCM_SCOL', RCM_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RICLD_SCOL', ICE_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NICLD_SCOL', NICE_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NCLD_SCOL', NCLW_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RRAIN_SCOL', RAIN_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NRAIN_SCOL', NRAIN_LH_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_WEIGHT_SCOL', weights, pcols*psubcols, lchnk )
      call outfld( 'NR_IN_LH', nrain, pcols, lchnk )
      call outfld( 'RTM_CLUBB', rtm, pcols, lchnk )
      call outfld( 'THLM_CLUBB', thlm, pcols, lchnk )
      call outfld( 'SILHS_QC_IN', state%q(:,:,ixcldliq), pcols, lchnk )
      call outfld( 'SILHS_QI_IN', state%q(:,:,ixcldice), pcols, lchnk )
      call outfld( 'SILHS_NC_IN', state%q(:,:,ixnumliq), pcols, lchnk )
      call outfld( 'SILHS_NI_IN', state%q(:,:,ixnumice), pcols, lchnk )
      call outfld( 'AKM_CLUBB', AKm_out, pcols, lchnk )
      call outfld( 'AKM_LH_CLUBB', lh_AKm_out, pcols, lchnk )
      call outfld( 'INVS_EXNER', invs_exner, pcols, lchnk )
      call outfld( 'SILHS_ZTODT', ztodt_ptr, pcols, lchnk )
      call outfld( 'SILHS_MSC_CLDICE', meansc_ice, pcols, lchnk )
      call outfld( 'SILHS_STDSC_CLDICE', stdsc_ice, pcols, lchnk )
      call outfld( 'SILHS_MSC_CLDLIQ', meansc_liq, pcols, lchnk )
      call outfld( 'SILHS_STDSC_CLDLIQ', stdsc_liq, pcols, lchnk )
      call outfld( 'SILHS_MSC_Q', meansc_vap, pcols, lchnk )
      call outfld( 'SILHS_STDSC_Q', stdsc_vap, pcols, lchnk )
      call outfld( 'SILHS_EFF_CLDFRAC', eff_cldfrac, pcols, lchnk )
      call outfld( 'SILHS_CLUBB_PRECIP_FRAC', precip_frac_out, pcols, lchnk )
      call outfld( 'SILHS_CLUBB_ICE_SS_FRAC', ice_supersat_frac, pcols, lchnk )

#endif
#endif
   end subroutine subcol_gen_SILHS

   subroutine subcol_ptend_avg_SILHS(ptend_sc, ngrdcol, lchnk, ptend)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_ptend_get_firstsubcol, subcol_ptend_avg_shr, &
                                  subcol_get_weight, subcol_get_filter, &
                                  is_filter_set, is_weight_set

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols)
      !-----------------------------------

      type(physics_ptend), intent(in)             :: ptend_sc        ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      type(physics_ptend), intent(inout)          :: ptend
      ! Because we can't get a state passed in here, we might have to use values from the 
      ! subcolumn generation. This would make any conservation checks invalid if this
      ! function is called after another parameterization... hmm.

       call subcol_ptend_avg_shr(ptend_sc, ngrdcol, lchnk, ptend, is_filter_set(), is_weight_set())

   end subroutine subcol_ptend_avg_SILHS

#ifdef SILHS
   real(r8) function meansc(arr_in, w_in, ns) result(val)
      real(r8), intent(in) :: ns                         ! Length of Array
      real(r8), dimension(ns), intent(in) :: arr_in      ! Input array
      real(r8), dimension(ns), intent(in) :: w_in        ! Weights
      real(r8) :: acc  ! accumulator
      integer :: i
      acc = 0
      val = 0
      do i=1,ns
         acc = acc + arr_in(i)*w_in(i)
      enddo
      val = acc/ns
   end function

   real(r8) function stdsc(arr_in, w_in, mn_in, ns) result(val)
      real(r8), intent(in) :: ns  ! Number of elements (subcolumns)
      real(r8), dimension(ns), intent(in) :: arr_in, w_in  !Input array and weights
      real(r8), intent(in) :: mn_in   ! The mean of arr_in
      real(r8) :: accvar, var
      integer :: i
      accvar = 0
      do i=1,ns
         accvar = accvar + ((arr_in(i)-mn_in)**2)*w_in(i)
      enddo
      var = accvar/ns
      val = sqrt(var)
   end function

   subroutine Abs_Temp_profile(nz, LWPT_prof, ex_prof, rcm_prof, ABST_prof)

      use clubb_api_module,              only : thlm2T_in_K_api

      integer,                 intent(in)  :: nz         ! Num vert levels
      real(r8), dimension(nz), intent(in)  :: LWPT_prof  ! Temp prof in LWPT
      real(r8), dimension(nz), intent(in)  :: ex_prof    ! Profile of Exner func
      real(r8), dimension(nz), intent(in)  :: rcm_prof   ! Profile of Cld Wat MR
      real(r8), dimension(nz), intent(out) :: ABST_prof  ! Abs Temp prof
      integer :: i
 
      do i=1,nz
         ABST_prof(i) = thlm2T_in_K_api(LWPT_prof(i), ex_prof(i), rcm_prof(i))
      enddo
      
   end subroutine

   subroutine THL_profile(nz, ABST_prof, ex_prof, rcm_prof, THL_prof)

      use clubb_api_module,              only : T_in_K2thlm_api

      integer,                 intent(in)  :: nz         ! Num vert levels
      real(r8), dimension(nz), intent(in)  :: ABST_prof  ! Abs Temp prof
      real(r8), dimension(nz), intent(in)  :: ex_prof    ! Profile of Exner func
      real(r8), dimension(nz), intent(in)  :: rcm_prof   ! Profile of Cld Wat MR
      real(r8), dimension(nz), intent(out) :: THL_prof  ! LWPT prof
      integer :: i
 
      do i=1,nz
         THL_prof(i) = T_in_K2thlm_api(ABST_prof(i), ex_prof(i), rcm_prof(i))
      enddo
      
   end subroutine

   subroutine StaticEng_profile(nz, ABST_prof, zm_prof, zsfc, s_prof)
      integer,                 intent(in) :: nz
      real(r8), dimension(nz), intent(in) :: ABST_prof
      real(r8), dimension(nz), intent(in) :: zm_prof
      real(r8),                intent(in) :: zsfc
      real(r8), dimension(nz), intent(out) :: s_prof
      integer :: i

      do i=1,nz
         s_prof(i) = cpair*(ABST_prof(i)) + gravit*zm_prof(i)+zsfc
      enddo

   end subroutine

   subroutine subcol_constrainmn( num_subcols, samples, weights, grid_mean, mean_sc, std_sc )

      ! Input/Output Variables
      integer, intent(in) :: num_subcols
      real(r8), dimension(num_subcols, pverp), intent(inout) :: samples
      real(r8), dimension(num_subcols), intent(in) :: weights
      real(r8), dimension(pverp), intent(in) :: grid_mean
      real(r8), dimension(pver), intent(out), optional :: mean_sc, std_sc

      ! Local Variables
      real(r8) :: meansc_loc, adj_rat
      integer :: k
   !------------------------------------------------------------------
      !----- Begin Code -----
      do k=1, pver
         meansc_loc = meansc( samples(:,k), weights(:), real(num_subcols, r8) )

         if (present(mean_sc)) &
            mean_sc(k) = meansc_loc
         if (present(std_sc)) &
            std_sc(k) = stdsc( samples(:,k), weights(:), meansc_loc, &
                               real(num_subcols, r8) )

         if ( meansc_loc > 0.0_r8 ) then
            adj_rat = grid_mean(k)/meansc_loc
         else 
            ! If the mean is zero, then zero out all subcolumns to avoid
            ! negative samples
            adj_rat = 0.0_r8
         end if
         samples(:,k) = samples(:,k) * adj_rat
      end do
   end subroutine subcol_constrainmn

   ! =============================================================================== !
   !                                                                                 !
   ! =============================================================================== !
   function clubb_flip_grid ( profile ) result( profile_flipped )

     ! Description:
     !   Swaps the elements in profile so they are in reverse order. CAM and
     !   CLUBB's grids are flipped with respect to each other.
     !
     !   Usage:
     !     clubb_var = clubb_flip_grid( cam_var )
     !     cam_var   = clubb_flip_grid( clubb_var )

     implicit none

     ! Input Variable
     real(r8), dimension(pverp), intent(in) :: profile

     ! Output Variable
     real(r8), dimension(pverp) :: profile_flipped

     ! Local Variable
     integer :: k

     do k=1, pverp
       profile_flipped(k) = profile(pverp-k+1)
     end do ! k=1, pverp

     return
   end function clubb_flip_grid

   ! =============================================================================== !
   !                                                                                 !
   ! =============================================================================== !

   subroutine subcol_SILHS_var_covar_driver &
              ( ztodt, state_sc, ptend_sc, &
                pbuf )

     ! This subroutine calculates microphysical effects on five variances and
     ! covariances: rtp2, thlp2, wprtp, wpthlp, and rtpthlp.
     !
     ! This code is experimental!!

     use physics_buffer,          only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
     use subcol_utils,            only: subcol_unpack, subcol_get_nsubcol, subcol_get_weight
     use clubb_api_module,        only: T_in_K2thlm_api
     use silhs_api_module,        only: lh_microphys_var_covar_driver_api

     implicit none

     ! Parameters
     !  This fill value is set to catch errors; it should not be read.
     real(r8), parameter                   :: fillvalue = -999._r8

     ! Input Variables
     real(r8), intent(in)                  :: ztodt        ! model time increment
     type(physics_state), intent(in)       :: state_sc     ! state for sub-columns
     type(physics_ptend), intent(in)       :: ptend_sc     ! ptend for sub-columns

     ! Pointers
     type(physics_buffer_desc), pointer    :: pbuf(:)

     ! Local Variables
     integer :: lchnk, ngrdcol, igrdcol, isubcol, ns, k
     integer, dimension(pcols) :: nsubcol
     real(r8), dimension(pcols*psubcols)       :: weights_packed
     real(r8), dimension(pcols,psubcols)       :: weights
     real(r8), dimension(pcols,psubcols,pverp) :: rc_all, rv_all, rt_all, w_all, thl_all
     real(r8), dimension(pcols,psubcols,pver ) :: s_all, t_all, zm_all, omega_all, pmid_all
     real(r8), dimension(pcols,psubcols)       :: phis_all
     real(r8), dimension(pcols,psubcols,pver ) :: stend, ttend
     real(r8), dimension(pcols,psubcols,pverp) :: thltend, qctend, qvtend

     real(r8), dimension(pcols,psubcols,pver)  :: dz_g, pdel_all, rho
     real(r8), dimension(pcols,psubcols,pverp) :: zi_all

     integer :: ixcldliq, ixq

     real(r8), dimension(pcols,psubcols,pver ) :: exner

     ! Inputs to lh_microphys_var_covar_driver
     real(r8), dimension(pcols,pverp,psubcols) :: rt_all_clubb, thl_all_clubb, w_all_clubb, &
                                                  qctend_clubb, qvtend_clubb, thltend_clubb
     ! Outputs from lh_microphys_var_covar_driver
     real(r8), dimension(:,:), pointer :: rtp2_mc_zt, thlp2_mc_zt, wprtp_mc_zt, &
                                          wpthlp_mc_zt, rtpthlp_mc_zt

     ! pbuf indices
     integer :: &
       rtp2_mc_zt_idx,   &
       thlp2_mc_zt_idx,  &
       wprtp_mc_zt_idx,  &
       wpthlp_mc_zt_idx, &
       rtpthlp_mc_zt_idx

     !----- Begin Code -----

     ! Don't do anything if this option isn't enabled.
     if ( .not. subcol_SILHS_var_covar_src ) return

     lchnk = state_sc%lchnk
     ngrdcol  = state_sc%ngrdcol

     ! Obtain indices
     call cnst_get_ind('Q', ixq)
     call cnst_get_ind('CLDLIQ', ixcldliq)
     rtp2_mc_zt_idx = pbuf_get_index('rtp2_mc_zt')
     thlp2_mc_zt_idx = pbuf_get_index('thlp2_mc_zt')
     wprtp_mc_zt_idx = pbuf_get_index('wprtp_mc_zt')
     wpthlp_mc_zt_idx = pbuf_get_index('wpthlp_mc_zt')
     rtpthlp_mc_zt_idx = pbuf_get_index('rtpthlp_mc_zt')

     ! Obtain pbuf fields for output
     call pbuf_get_field(pbuf, rtp2_mc_zt_idx, rtp2_mc_zt)
     call pbuf_get_field(pbuf, thlp2_mc_zt_idx, thlp2_mc_zt)
     call pbuf_get_field(pbuf, wprtp_mc_zt_idx, wprtp_mc_zt)
     call pbuf_get_field(pbuf, wpthlp_mc_zt_idx, wpthlp_mc_zt)
     call pbuf_get_field(pbuf, rtpthlp_mc_zt_idx, rtpthlp_mc_zt)

     ! Unpack needed tendencies from subcolumn ptends
     call subcol_unpack(lchnk, ptend_sc%s(:,:), stend, fillvalue)
     call subcol_unpack(lchnk, ptend_sc%q(:,:,ixcldliq), qctend(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, ptend_sc%q(:,:,ixq), qvtend(:,:,1:pver), fillvalue)

     ! Unpack sample point values from subcolumn states
     call subcol_unpack(lchnk, state_sc%q(:,:,ixcldliq), rc_all(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, state_sc%q(:,:,ixq),      rv_all(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, state_sc%omega(:,:),      omega_all (:,:,:),  fillvalue)
     call subcol_unpack(lchnk, state_sc%s(:,:),          s_all,              fillvalue)
     call subcol_unpack(lchnk, state_sc%zm,              zm_all,             fillvalue)
     call subcol_unpack(lchnk, state_sc%phis,            phis_all,           fillvalue)
     call subcol_unpack(lchnk, state_sc%zi,              zi_all,             fillvalue)
     call subcol_unpack(lchnk, state_sc%pdel,            pdel_all,           fillvalue)
     call subcol_unpack(lchnk, state_sc%pmid,            pmid_all,           fillvalue)
     rt_all(:,:,1:pver) = rc_all(:,:,1:pver) + rv_all(:,:,1:pver)

     ! Compute dry static density on CLUBB vertical grid
     do k=1, pver
       dz_g(:,:,k) = zi_all(:,:,k) - zi_all(:,:,k+1) ! thickness
       rho(:,:,k) = (1._r8/gravit)*pdel_all(:,:,k)/dz_g(:,:,k)
     enddo
     ! Compute w from omega
     w_all(:,:,1:pver) = -omega_all / ( rho * gravit )

     ! Convert stend and s_all to ttend and t_all
     !  Note 1: With subcolumns, cpair is truly a constant (I think).
     !  Note 2: For tendencies, the extra terns zm and phis should
     !          not be included in the calculation.
     ttend = stend / cpair
     do k=1, pver
       t_all(:,:,k) = (s_all(:,:,k) - gravit*zm_all(:,:,k) - phis_all(:,:)) / cpair
     end do ! k=1, pver

     ! This formula is taken from earlier in this file.
     exner(:,:,:) = (pmid_all(:,:,:)/p0_clubb)**(rair/cpair)

     ! Note: all tendencies or all means should be used in the call to
     !       T_in_K2thlm_api (with the exception of exner)
     thltend(:,:,1:pver) = T_in_K2thlm_api(ttend(:,:,:), exner(:,:,:), qctend(:,:,1:pver))
     thl_all(:,:,1:pver) = T_in_K2thlm_api(t_all(:,:,:), exner(:,:,:), rc_all(:,:,1:pver))

     ! Add ghost points
     rt_all (:,:,pverp)   = rt_all (:,:,pver)
     thl_all(:,:,pverp)   = thl_all(:,:,pver)
     w_all  (:,:,pverp)   = w_all  (:,:,pver)
     qctend (:,:,pverp)   = qctend (:,:,pver)
     qvtend (:,:,pverp)   = qvtend (:,:,pver)
     thltend(:,:,pverp)   = thltend(:,:,pver)

     ! How many subcolumns in each column?
     call subcol_get_nsubcol(lchnk, nsubcol)

     ! Flip inputs to CLUBB's grid. Note the dimension ordering change.
     do igrdcol=1, ngrdcol
       do isubcol=1, nsubcol(igrdcol)
         rt_all_clubb(igrdcol,:,isubcol) = clubb_flip_grid( rt_all(igrdcol,isubcol,:) )
         thl_all_clubb(igrdcol,:,isubcol) = clubb_flip_grid( thl_all(igrdcol,isubcol,:) )
         w_all_clubb(igrdcol,:,isubcol) = clubb_flip_grid( w_all(igrdcol,isubcol,:) )
         qctend_clubb(igrdcol,:,isubcol) = clubb_flip_grid( qctend(igrdcol,isubcol,:) )
         qvtend_clubb(igrdcol,:,isubcol) = clubb_flip_grid( qvtend(igrdcol,isubcol,:) )
         thltend_clubb(igrdcol,:,isubcol) = clubb_flip_grid( thltend(igrdcol,isubcol,:) )
       end do
     end do

     ! Obtain weights
     call subcol_get_weight(lchnk, weights_packed)
     call subcol_unpack(lchnk, weights_packed, weights, fillvalue)

     ! Call lh_microphys_var_covar_driver for each column
     do igrdcol=1, ngrdcol
       ns = nsubcol(igrdcol)
       ! Make the call!!!!!
       call lh_microphys_var_covar_driver_api &
            ( pverp, ns, ztodt, weights(igrdcol,1:ns), &
              rt_all_clubb(igrdcol,:,1:ns), thl_all_clubb(igrdcol,:,1:ns), w_all_clubb(igrdcol,:,1:ns), &
              qctend_clubb(igrdcol,:,1:ns), qvtend_clubb(igrdcol,:,1:ns), thltend_clubb(igrdcol,:,1:ns), &
              rtp2_mc_zt(igrdcol,:), thlp2_mc_zt(igrdcol,:), wprtp_mc_zt(igrdcol,:), &
              wpthlp_mc_zt(igrdcol,:), rtpthlp_mc_zt(igrdcol,:) )
     end do ! igrdcol=1, ngrdcol

     return
   end subroutine subcol_SILHS_var_covar_driver

   ! =============================================================================== !
   !                                                                                 !
   ! =============================================================================== !

   subroutine subcol_SILHS_massless_droplet_destroyer &
              ( ztodt, state, &
                ptend )

     ! This subroutine eradicates cloud droplets in grid boxes with no cloud
     ! mass. This code is actually not at all related to SILHS. It should be
     ! moved out of this cage they call "subcol_SILHS.F90" someday.
     !
     ! This code is experimental!!

     use micro_mg_utils, only: qsmall

     implicit none

     ! Input Variables
     real(r8), intent(in)                  :: ztodt     ! model time increment
     type(physics_state), intent(in)       :: state     ! state for columns

     ! Input/Output Variables
     type(physics_ptend), intent(inout)    :: ptend     ! ptend for columns

     ! Local Variables
     integer :: icol, k

     integer :: ixcldliq, ixnumliq

     !----- Begin Code -----

     ! Don't do anything if this option isn't enabled.
     if ( .not. subcol_SILHS_destroy_massless_droplets ) return

     ! Indices!
     call cnst_get_ind('CLDLIQ', ixcldliq)
     call cnst_get_ind('NUMLIQ', ixnumliq)

     ! These "labels" in loops are really cool. We should start doing this in
     ! CLUBB.
     col_loop: do icol=1, state%ncol
       ! If updated qc (after microphysics) is zero, then ensure updated nc is also zero!!
       vert_loop: do k=1, pver
         if ( state%q(icol,k,ixcldliq) + (ztodt*ptend%q(icol,k,ixcldliq)) < qsmall ) then
           ptend%lq(ixnumliq) = .true. ! This is probably already true, but it doesn't
                                       ! hurt to set it.
           ptend%q(icol,k,ixnumliq) = -(state%q(icol,k,ixnumliq) / ztodt)
         end if
       end do vert_loop
     end do col_loop

     return
   end subroutine subcol_SILHS_massless_droplet_destroyer

#endif
   
end module subcol_SILHS 
