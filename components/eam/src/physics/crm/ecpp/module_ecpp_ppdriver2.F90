module module_ecpp_ppdriver2

!-------------------------------------------------------------------------------------
!      Purpose: 
!      Provide the CAM interface to the Explicit-Cloud Parameterized-Pollutant hygrid
!      approach for aerosol-cloud interactions in the MMF models.  
! 
!      This module was adopted from the one written for the WRF-chem by Dick Easter. 
!
!      Minghuai Wang (Minghuai.Wang@pnl.gov), 2009-11 
!---------------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use ppgrid,         only: pcols, pver, pverp 
  use constituents,   only: pcnst, cnst_name
  use cam_abortutils, only: endrun
  use crmdims,        only: crm_nz
  use ecppvars,       only: nupdraft_in, ndndraft_in, ncls_ecpp_in, ncc_in, nprcp_in 
  use crmclouds_camaerosols, only: ecpp_mixnuc_tend => crmclouds_mixnuc_tend 
  use module_data_ecpp1 
  use module_data_mosaic_asect

  implicit none

  public :: parampollu_driver2
  public :: papampollu_init
  public :: ecpp_mixnuc_tend

  !+++mhwang follow what is done in ndrop.F90 for qqcw
  ! ptr2d_t is used to create arrays of pointers to 2D fields
  type ptr2d_t
     real(r8), pointer :: fldcw(:,:)
  end type ptr2d_t

  contains

  !-----------------------------------------------------------------------------------------------
  !
  ! rce 2005-mar-10 - created
  !
  !-----------------------------------------------------------------------------------------------

  !==================================================================================================
  !==================================================================================================
  !==================================================================================================
  subroutine   papampollu_init ( )
    !------------------------------------------------------
    ! Description:
    ! initialize some data used in ECPP
    ! map aerosol inforation in cam4 into mosaic.
    ! - Minghuai Wang, 2009-11
    !------------------------------------------------------
    use cam_history, only: addfld,  add_default,horiz_only
    use modal_aero_data
    use module_ecpp_td2clm, only:  set_of_aerosol_stuff 
    use module_ecpp_util,   only:  parampollu_1clm_set_opts

    ! Local variables
    integer  :: n, ll
    integer  :: ichem, ichem2
    real(r8) :: pi
    real(r8) :: tmpa

    pi = 4._r8*atan(1._r8)

    !
    ! set pp options (should this be done from driver?)
    !

    num_moist_ecpp   = 9 
    num_moist        = 9
    num_chem_ecpp    = 2 * pcnst      ! 2x for cloud-borne and interstitial aerosol
    num_chem         = num_chem_ecpp
    param_first_ecpp = num_moist+1    ! the first index for non-water species
    p_qv = 1
    p_qc = 2

    allocate (is_aerosol(1:num_chem_ecpp))
    allocate (iphase_of_aerosol(1:num_chem_ecpp))
    allocate (isize_of_aerosol(1:num_chem_ecpp))
    allocate (itype_of_aerosol(1:num_chem_ecpp))
    allocate (inmw_of_aerosol(1:num_chem_ecpp))
    allocate (laicwpair_of_aerosol(1:num_chem_ecpp))

    !------------------------------------------------------------------------------------------------------------
    ! Map the modal aerosol information in modal_aero_data.F90 to module_data_mosaic_asect.F90 
    ! In the ECPP written for the WRF-chem, it used the MOSAIC aerosol data. MOSAIC have different
    ! classifications, and use aeroso types, aerosol size bins, chemical components, and aerosol phases
    ! to describe aerosols. In the CAM4's modal aerosol treatment, it use aerosol modes, and chemical
    ! components to describe aerosols, and interstial and cloud-borne aerosols are separately tracked.  
    ! When the ECPP codes are ported from the WRF-chem into the MMF model (CAM4.0_SAM), 
    ! the MOSAIC's description of the aerosols are kept, in order to minimize
    ! the codes changes, but the aerosol information in CAM4.0 is mapped into the MOSAIC one in the
    ! following way: aeroso type is equivalent to aerosol modes in CAM4, and aerosol size is one for each aerosol type, 
    ! and the aerosol chemical composition is just the same as that in CAM4. Interstitial aerosols in CAM4 is put into
    ! the phase 1, and cloud-borne aerosol in CAM4 is put into the pase 2. -Minghuai Wang (minghuai.wang@pnl.gov)
    !------------------------------------------------------------------------------------------------------------

    maxd_atype = ntot_amode
    maxd_asize = 1
    maxd_acomp = maxd_aspectype
    maxd_aphase = 2

    ai_phase = 1      ! index for interstial aerosols
    cw_phase = 2      ! index for cloud-borne aerosols

    ntype_aer = ntot_amode
    nphase_aer = 2

    allocate (nsize_aer( 1:maxd_atype ))
    allocate (ncomp_aer( 1:maxd_atype ))
    allocate (massptr_aer( 1:maxd_acomp, 1:maxd_asize, 1:maxd_atype, 1:maxd_aphase ))
    allocate (numptr_aer( 1:maxd_asize, 1:maxd_atype, 1:maxd_aphase ))
    allocate (dens_aer( 1:maxd_acomp, 1:maxd_atype ))
    allocate (hygro_aer( 1:maxd_acomp, 1:maxd_atype ))
    allocate (volumhi_sect( 1:maxd_asize, 1:maxd_atype ))
    allocate (volumlo_sect( 1:maxd_asize, 1:maxd_atype ))
    allocate (sigmag_aer( 1:maxd_asize, 1:maxd_atype ))
    allocate (dcen_sect(1:maxd_asize, 1:maxd_atype ))
    allocate (dlo_sect(1:maxd_asize, 1:maxd_atype ))
    allocate (dhi_sect(1:maxd_asize, 1:maxd_atype ))


    nsize_aer(1:maxd_atype) = 1
    ncomp_aer(1:maxd_atype) = nspec_amode(1:ntot_amode)

    massptr_aer(1:maxd_acomp, 1, 1:maxd_atype, 1) = lmassptr_amode(1:maxd_aspectype, 1:ntot_amode)
    massptr_aer(1:maxd_acomp, 1, 1:maxd_atype, 2) = lmassptrcw_amode(1:maxd_aspectype, 1:ntot_amode) + pcnst

    numptr_aer(1, 1:maxd_atype, 1) = numptr_amode(1:ntot_amode)
    numptr_aer(1, 1:maxd_atype, 2) = numptrcw_amode(1:ntot_amode) + pcnst

    do n=1, ntype_aer
      do ll=1, ncomp_aer(n)
        dens_aer(ll, n) = specdens_amode(lspectype_amode(ll, n))  
        hygro_aer(ll, n) = spechygro(lspectype_amode(ll, n))
      end do

      sigmag_aer(1, n) = sigmag_amode(n)
      !-----------------------------------------------------------------------
      ! Notes: 
      !        the tmpa factor is because
      !        dcen_sect, dlo_sect, dhi_sect are used as,
      !            and are compared to, volume-mean diameters
      !        dgnum_amode, dgnumlo_amode, dgnumhi_amode are used as,
      !            and are compared to, number-distribution geometric-mean diameters
      !        volume_mixing_ratio/(number_mixing_ratio*pi/6)
      !                = volume_mean_diameter**3
      !                = (number_geometric_mean_diameter*tmpa)**3
      !-----------------------------------------------------------------------

      tmpa = exp( 1.5 * log(sigmag_amode(n))**2 )
      dcen_sect(1, n) =   dgnum_amode(n)*tmpa
      dlo_sect( 1, n) = dgnumlo_amode(n)*tmpa
      dhi_sect( 1, n) = dgnumhi_amode(n)*tmpa

      volumlo_sect(1, n) = pi/6 * (dgnumlo_amode(n)*tmpa)**3
      volumhi_sect(1, n) = pi/6 * (dgnumhi_amode(n)*tmpa)**3
    end do

    afrac_cut = aw_draft_cut/w_draft_max
    afrac_cut_bb  = afrac_cut*0.5_r8
    afrac_cut_0p5 = afrac_cut*0.5_r8
    afrac_cut_0p2 = afrac_cut*0.2_r8
    afrac_cut_0p1 = afrac_cut*0.1_r8

    ! set flags 
    activat_onoff_ecpp = 1          ! droplet activation; 1 turns on activation
    cldchem_onoff_ecpp = 1          ! cloud chemistry 
    rename_onoff_ecpp  = 1          ! renaming (modal merging)

    wetscav_onoff_ecpp = 400        ! wet removable 400 turn on wet scaving
          
    ! set convection lifetime
    draft_lifetime = 7200           ! seconds, 2 hours lifetime for the momement

    ! set flag for a/c partition
    iflag_ecpp_startup_acw_partition = 1  ! 1 to turn on a/c parition

    ! set flag for whether update changs from host codes
    iflag_ecpp_startup_host_chemtend = 0

    ! set other flags
    iflag_ecpp_test_bypass_1 = 0
    iflag_ecpp_test_fixed_fcloud = 0

    parampollu_opt = 2223  !  method flag for parameterized-pollutants module

    !
    !   set pp options (should this be done from driver?)
    !
    call parampollu_1clm_set_opts(ppopt_updn_prof_aa_wfull,  &
                                  ppopt_quiescn_mf_byppmx1,  &
                                  ppopt_quiescn_sosi_x1,     &
                                  ppopt_chemtend_wq_wfullx1, &
                                  ppopt_chemtend_dtsub_x1,   &
                                  ppopt_chemtend_updnfreq_x1 )

    !
    !   some other initialization
    !
    call set_of_aerosol_stuff(is_aerosol, iphase_of_aerosol,        &
                              isize_of_aerosol, itype_of_aerosol,   &
                              inmw_of_aerosol, laicwpair_of_aerosol )

    ! add fields into history file
    do ichem=param_first_ecpp, pcnst
      if(trim(cnst_name(ichem))//'EP' == 'EP') then
        write(0, *) ichem, trim(cnst_name(ichem))//'EP'
        call endrun('ecpp init1')
      end if 

      !==Guangxing Lin, modify the addfld arguments
      call addfld(trim(cnst_name(ichem))//'EP',(/'lev'/),'A', 'kg/kg/s', trim(cnst_name(ichem))//' tendency from ECPP' )
      call addfld(trim(cnst_name(ichem))//'ACHEM_EP',(/'lev'/),'A','kg/kg/s',  trim(cnst_name(ichem))//' tendency from aqueous chemistry in ECPP' )
      call addfld(trim(cnst_name(ichem))//'RENM_EP',(/'lev'/),'A','kg/kg/s',  trim(cnst_name(ichem))//' tendency from renaming in ECPP')
      call addfld(trim(cnst_name(ichem))//'ACT_EP', (/'lev'/),'A','kg/kg/s',  trim(cnst_name(ichem))//' tendency from activation/resuspension in ECPP' )
      call addfld(trim(cnst_name(ichem))//'WET_EP', (/'lev'/),'A', 'kg/kg/s',trim(cnst_name(ichem))//' tendency from wet removable in ECPP')
      call addfld(trim(cnst_name(ichem))//'WRESU_EP', (/'lev'/),'A','kg/kg/s',  trim(cnst_name(ichem))//' tendency from resuspension in wet removable in ECPP' )
      call addfld(trim(cnst_name(ichem))//'CONV_EP', (/'lev'/),'A','kg/kg/s', trim(cnst_name(ichem))//' tendency from convective tansport in ECPP' )

      call addfld(trim(cnst_name(ichem))//'SFEP', horiz_only, 'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from ECPP' )
      call addfld(trim(cnst_name(ichem))//'SFACHEM_EP',horiz_only,'A', 'kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from aqueus chemistry in ECPP' )
      call addfld(trim(cnst_name(ichem))//'SFRENM_EP',horiz_only,'A', 'kg/m2/s', trim(cnst_name(ichem))//' column-integrated tendency from renaming in ECPP' )
      call addfld(trim(cnst_name(ichem))//'SFACT_EP', horiz_only,'A', 'kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from activation/resuspension ECPP') 
      call addfld(trim(cnst_name(ichem))//'SFWET_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from wet removable in ECPP')
      call addfld(trim(cnst_name(ichem))//'SFWRESU_EP',horiz_only,'A', 'kg/m2/s', trim(cnst_name(ichem))//' column-integrated tendency from resupspension in wet removable in ECPP')
      call addfld(trim(cnst_name(ichem))//'SFCONV_EP',horiz_only,'A', 'kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from convective transport in ECPP')
      ! call addfld(trim(cnst_name(ichem))//'EP', 'kg/kg/s', pver, 'A', trim(cnst_name(ichem))//' tendency from ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'ACHEM_EP', 'kg/kg/s', pver, 'A', trim(cnst_name(ichem))//' tendency from aqueous chemistry in ECPP' , phys_decomp )
      !call addfld(trim(cnst_name(ichem))//'RENM_EP', 'kg/kg/s', pver, 'A', trim(cnst_name(ichem))//' tendency from renaming in ECPP' , phys_decomp )
      !call addfld(trim(cnst_name(ichem))//'ACT_EP', 'kg/kg/s', pver, 'A', trim(cnst_name(ichem))//' tendency from activation/resuspension in ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'WET_EP', 'kg/kg/s', pver, 'A', trim(cnst_name(ichem))//' tendency from wet removable in ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'WRESU_EP', 'kg/kg/s', pver, 'A', trim(cnst_name(ichem))//' tendency from resuspension in wet removable in ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'CONV_EP', 'kg/kg/s', pver, 'A', trim(cnst_name(ichem))//' tendency from convective tansport in ECPP' , phys_decomp)

      !call addfld(trim(cnst_name(ichem))//'SFEP', 'kg/m2/s', 1, 'A', trim(cnst_name(ichem))//' column-integrated tendency from ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'SFACHEM_EP', 'kg/m2/s', 1, 'A', trim(cnst_name(ichem))//' column-integrated tendency from aqueus chemistry in ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'SFRENM_EP', 'kg/m2/s', 1, 'A', trim(cnst_name(ichem))//' column-integrated tendency from renaming in ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'SFACT_EP', 'kg/m2/s', 1, 'A', trim(cnst_name(ichem))//' column-integrated tendency from activation/resuspension ECPP' , phys_decomp) 
      !call addfld(trim(cnst_name(ichem))//'SFWET_EP', 'kg/m2/s', 1, 'A', trim(cnst_name(ichem))//' column-integrated tendency from wet removable in ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'SFWRESU_EP', 'kg/m2/s', 1, 'A', trim(cnst_name(ichem))//' column-integrated tendency from resupspension in wet removable in ECPP' , phys_decomp)
      !call addfld(trim(cnst_name(ichem))//'SFCONV_EP', 'kg/m2/s', 1, 'A', trim(cnst_name(ichem))//' column-integrated tendency from convective transport in ECPP' , phys_decomp)

      ! Quiescent class
      call addfld(trim(cnst_name(ichem))//'SFACHQU_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from aqueus chemistry in ECPP (quiescent)' )
      call addfld(trim(cnst_name(ichem))//'SFREMQU_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from renaming in ECPP (quiescent)' )
      call addfld(trim(cnst_name(ichem))//'SFACTQU_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from activation/resuspension ECPP (quiescent)' )
      call addfld(trim(cnst_name(ichem))//'SFWETQU_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from wet removable in ECPP (quiescent)' )
      call addfld(trim(cnst_name(ichem))//'SFRESQU_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from resupspension in wet removable in ECPP (quiescent)' )
      call addfld(trim(cnst_name(ichem))//'SFCONQU_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from convective transport in ECPP (quiescent)' )

      ! Updraft class
      call addfld(trim(cnst_name(ichem))//'SFACHUP_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from aqueus chemistry in ECPP (updraft)' )
      call addfld(trim(cnst_name(ichem))//'SFREMUP_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from renaming in ECPP (updraft)' )
      call addfld(trim(cnst_name(ichem))//'SFACTUP_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from activation/resuspension ECPP (updraft)' )
      call addfld(trim(cnst_name(ichem))//'SFWETUP_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from wet removable in ECPP (updraft)' )
      call addfld(trim(cnst_name(ichem))//'SFRESUP_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from resupspension in wet removable in ECPP (updraft)' )
      call addfld(trim(cnst_name(ichem))//'SFCONUP_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from convective transport in ECPP (updraft)' )

      ! Downdraft class
      call addfld(trim(cnst_name(ichem))//'SFACHDN_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from aqueus chemistry in ECPP (downdraft)' )
      call addfld(trim(cnst_name(ichem))//'SFREMDN_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from renaming in ECPP (downdraft)' )
      call addfld(trim(cnst_name(ichem))//'SFACTDN_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from activation/resuspension ECPP (downdraft)' )
      call addfld(trim(cnst_name(ichem))//'SFWETDN_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from wet removable in ECPP (downdraft)' )
      call addfld(trim(cnst_name(ichem))//'SFRESDN_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from resupspension in wet removable in ECPP (downdraft)' )           
      call addfld(trim(cnst_name(ichem))//'SFCONDN_EP', horiz_only,'A','kg/m2/s',  trim(cnst_name(ichem))//' column-integrated tendency from convective transport in ECPP (downdraft)' )

    end do
    do ichem=param_first_ecpp, pcnst
      if(.not. (cnst_name_cw(ichem) == ' ')) then
        call addfld(trim(cnst_name_cw(ichem))//'EP', (/ 'lev' /), 'A', 'kg/kg/s', &
          trim(cnst_name_cw(ichem))//' tendency from ECPP' )
        call addfld(trim(cnst_name_cw(ichem))//'ACHEM_EP', (/ 'lev' /), 'A', 'kg/kg/s', &
          trim(cnst_name_cw(ichem))//' tendency from aqueous chemistry in ECPP' )
        call addfld(trim(cnst_name_cw(ichem))//'RENM_EP', (/ 'lev' /), 'A', 'kg/kg/s', &
          trim(cnst_name_cw(ichem))//' tendency from renaming in ECPP' )
        call addfld(trim(cnst_name_cw(ichem))//'ACT_EP', (/ 'lev' /), 'A', 'kg/kg/s', &
          trim(cnst_name_cw(ichem))//' tendency from activation/resuspension in ECPP' )
        call addfld(trim(cnst_name_cw(ichem))//'WET_EP', (/ 'lev' /), 'A', 'kg/kg/s', &
          trim(cnst_name_cw(ichem))//' tendency from wet removable in ECPP' )
        call addfld(trim(cnst_name_cw(ichem))//'WRESU_EP', (/ 'lev' /), 'A', 'kg/kg/s', &
          trim(cnst_name_cw(ichem))//' tendency from resuspension in wet removable in ECPP' )
        call addfld(trim(cnst_name_cw(ichem))//'CONV_EP', (/ 'lev' /), 'A', 'kg/kg/s', &
          trim(cnst_name_cw(ichem))//' tendency from convective tansport in ECPP' )
       
         call addfld(trim(cnst_name_cw(ichem))//'SFEP'      ,horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from ECPP' )
         call addfld(trim(cnst_name_cw(ichem))//'SFACHEM_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from aqueus chemistry in ECPP' )
         call addfld(trim(cnst_name_cw(ichem))//'SFRENM_EP' ,horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from renaming chemistry in ECPP' )
         call addfld(trim(cnst_name_cw(ichem))//'SFACT_EP'  ,horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from activation/resuspension ECPP' ) 
         call addfld(trim(cnst_name_cw(ichem))//'SFWET_EP'  ,horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from wet removable in ECPP' )
         call addfld(trim(cnst_name_cw(ichem))//'SFWRESU_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from resuspension in wet removable in ECPP' )
         call addfld(trim(cnst_name_cw(ichem))//'SFCONV_EP' ,horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from convective transport in ECPP' )

        ! Quiescent class
        call addfld(trim(cnst_name_cw(ichem))//'SFACHQU_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from aqueus chemistry in ECPP (quiescent)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFREMQU_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from renaming in ECPP (quiescent)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFACTQU_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from activation/resuspension ECPP (quiescent)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFWETQU_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from wet removable in ECPP (quiescent)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFRESQU_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from resupspension in wet removable in ECPP (quiescent)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFCONQU_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from convective transport in ECPP (quiescent)' )

        ! Updraft class
        call addfld(trim(cnst_name_cw(ichem))//'SFACHUP_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from aqueus chemistry in ECPP (updraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFREMUP_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from renaming in ECPP (updraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFACTUP_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from activation/resuspension ECPP (updraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFWETUP_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from wet removable in ECPP (updraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFRESUP_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from resupspension in wet removable in ECPP (updraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFCONUP_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from convective transport in ECPP (updraft)' )

        ! Downdraft class
        call addfld(trim(cnst_name_cw(ichem))//'SFACHDN_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from aqueus chemistry in ECPP (downdraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFREMDN_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from renaming in ECPP (downdraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFACTDN_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from activation/resuspension ECPP (downdraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFWETDN_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from wet removable in ECPP (downdraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFRESDN_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from resupspension in wet removable in ECPP (downdraft)' )
        call addfld(trim(cnst_name_cw(ichem))//'SFCONDN_EP',horiz_only,'A','kg/m2/s',trim(cnst_name_cw(ichem))//' column-integrated tendency from convective transport in ECPP (downdraft)' )

      end if
    end do

    call addfld('AQSO4_H2O2_EP', horiz_only,'A','kg/m2/s',  'SO4  aqueous phase chemistry due to H2O2 (kg/m2/s) in ECPP' )
    call addfld('AQSO4_O3_EP', horiz_only,'A','kg/m2/s', 'SO4  aqueous phase chemistry due to O3 (kg/m2/s) in ECPP' )
    call addfld('XPH_LWC_EP', (/ 'lev' /), 'A', ' ', 'pH value multiplied by lwc in ECPP')
    call add_default('AQSO4_H2O2_EP', 1, ' ')
    call add_default('AQSO4_O3_EP', 1, ' ')
    call add_default('XPH_LWC_EP', 1, ' ')

    do ichem=param_first_ecpp, pcnst 

      if(.not. (cnst_name_cw(ichem) == ' ')) then
        call add_default(trim(cnst_name_cw(ichem))//'SFEP'      , 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFACHEM_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFRENM_EP' , 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFACT_EP'  , 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFWET_EP'  , 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFWRESU_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFCONV_EP' , 1, ' ')

        call add_default(trim(cnst_name_cw(ichem))//'SFACHQU_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFREMQU_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFACTQU_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFWETQU_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFRESQU_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFCONQU_EP', 1, ' ')

        call add_default(trim(cnst_name_cw(ichem))//'SFACHUP_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFREMUP_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFACTUP_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFWETUP_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFRESUP_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFCONUP_EP', 1, ' ')

        call add_default(trim(cnst_name_cw(ichem))//'SFACHDN_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFREMDN_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFACTDN_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFWETDN_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFRESDN_EP', 1, ' ')
        call add_default(trim(cnst_name_cw(ichem))//'SFCONDN_EP', 1, ' ')
      end if

      call add_default(trim(cnst_name(ichem))//'SFEP'      , 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFACHEM_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFRENM_EP' , 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFACT_EP'  , 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFWET_EP'  , 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFWRESU_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFCONV_EP' , 1, ' ')

      call add_default(trim(cnst_name(ichem))//'SFACHQU_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFREMQU_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFACTQU_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFWETQU_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFRESQU_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFCONQU_EP', 1, ' ')

      call add_default(trim(cnst_name(ichem))//'SFACHUP_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFREMUP_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFACTUP_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFWETUP_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFRESUP_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFCONUP_EP', 1, ' ')

      call add_default(trim(cnst_name(ichem))//'SFACHDN_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFREMDN_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFACTDN_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFWETDN_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFRESDN_EP', 1, ' ')
      call add_default(trim(cnst_name(ichem))//'SFCONDN_EP', 1, ' ')

    end do

    ! for test purpose, additional 3D tendency 
    do ichem=param_first_ecpp, pcnst
      if(trim(cnst_name(ichem)) == 'DMS' .or. trim(cnst_name(ichem)) == 'SO2' .or. &
        trim(cnst_name(ichem)) == 'so4_a1') then
        if(.not. (cnst_name_cw(ichem) == ' ')) then
          call add_default(trim(cnst_name_cw(ichem))//'EP'      , 1, ' ')
          call add_default(trim(cnst_name_cw(ichem))//'ACHEM_EP', 1, ' ')
          call add_default(trim(cnst_name_cw(ichem))//'RENM_EP' , 1, ' ')
          call add_default(trim(cnst_name_cw(ichem))//'ACT_EP'  , 1, ' ')
          call add_default(trim(cnst_name_cw(ichem))//'WET_EP'  , 1, ' ')
          call add_default(trim(cnst_name_cw(ichem))//'CONV_EP' , 1, ' ')
        end if
        call add_default(trim(cnst_name(ichem))//'EP'      , 1, ' ')
        call add_default(trim(cnst_name(ichem))//'ACHEM_EP', 1, ' ')
        call add_default(trim(cnst_name(ichem))//'RENM_EP' , 1, ' ')
        call add_default(trim(cnst_name(ichem))//'ACT_EP'  , 1, ' ')
        call add_default(trim(cnst_name(ichem))//'WET_EP'  , 1, ' ')
        call add_default(trim(cnst_name(ichem))//'CONV_EP' , 1, ' ')
      end if 
    end do

  end subroutine papampollu_init
  !==================================================================================================
  !==================================================================================================
  !==================================================================================================
  subroutine parampollu_driver2(state, ptend,  pbuf,                              &
                                dtstep_in, dtstep_pp_in,                          &
                                acen_3d, abnd_3d,                                 &
                                acen_tf_3d, abnd_tf_3d,                           &
                                massflxbnd_3d,                                    &
                                rhcen_3d, qcloudcen_3d, qlsinkcen_3d,             &
                                precrcen_3d, precsolidcen_3d,                     &
                                acldy_cen_tbeg_3d                                 &
                               )

    ! modules from CAM
    use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
    use physics_buffer, only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
    use physconst,      only: gravit 
    use time_manager,   only: get_nstep, is_first_step
    use constituents,   only: cnst_name
    use cam_history,    only: outfld
#ifdef MODAL_AERO
    use modal_aero_data, only: ntot_amode, cnst_name_cw,  qqcw_get_field
#endif

    ! modules from ECPP
    use module_ecpp_td2clm, only:  parampollu_td240clm

    implicit none

    !-----------------------------------------------------------------------
    ! DESCRIPTION
    !
    ! parampollu_driver2 is the interface between wrf-chem and the
    ! parameterized pollutants "1 column" routine
    !
    ! main inputs : 
    !    aerosol and trace gas mixing ratios for a subset of the
    ! host-code domain : 
    !    ecpp (sub-grid) cloud statistics for the same subset of domain
    ! main outputs :
    !    updated aerosol and trace gas mixing ratios, with changes due
    !    to sub-grid vertical transport, activation/resuspension,
    !    cloud chemistry, and wet removal
    !
    !-----------------------------------------------------------------------

    !!! Interface Arguments

    real(r8), intent(in) :: dtstep_in, dtstep_pp_in
    ! dtstep_in - main model time step (s)
    ! dtstep_pp_in - time step (s) for "parameterized pollutants" calculations

    type(physics_state),       intent(in)    :: state      ! Physics state variables
    type(physics_ptend),       intent(inout) :: ptend      ! individual parameterization tendencies
    type(physics_buffer_desc), pointer       :: pbuf(:)    ! physics buffer 

    real(r8), intent(in), dimension(pcols,pverp,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: abnd_3d
    real(r8), intent(in), dimension(pcols,pverp,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: abnd_tf_3d
    real(r8), intent(in), dimension(pcols,pverp,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: massflxbnd_3d
    real(r8), intent(in), dimension(pcols,pver ,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: acen_3d
    real(r8), intent(in), dimension(pcols,pver ,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: acen_tf_3d
    real(r8), intent(in), dimension(pcols,pver ,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: rhcen_3d
    real(r8), intent(in), dimension(pcols,pver ,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: qcloudcen_3d
    real(r8), intent(in), dimension(pcols,pver ,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: qlsinkcen_3d
    real(r8), intent(in), dimension(pcols,pver ,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: precrcen_3d
    real(r8), intent(in), dimension(pcols,pver ,1:ncc_in,1:ncls_ecpp_in,1:nprcp_in ) :: precsolidcen_3d
    !-----------------------------------------------------------------------
    ! *** note - these are not "3d" now but probably could be in the mmf code
    !   abnd_3d and abnd_tf_3d  - sub-class frac area (--) at layer bottom boundary
    !       abnd_3d             - average for full time period (=dtstep_pp_in)
    !       abnd_tf_3d          - average for end-portion of time period
    !   acen_3d and acen_tf_3d  - sub-class frac area (--) at layer center
    !       acen_3d             - average for full time period (=dtstep_pp_in)
    !       acen_tf_3d          - average for end-portion of time period
    !   massflxbnd_3d - sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
    !-----------------------------------------------------------------------
    ! *** note - These are calculated using wfull, not wprime.
    !   rhcen_3d         - relative humidity (0-1) at layer center
    !   qcloudcen_3d     - cloud water mixing ratio (kg/kg) at layer center
    !   qlsinkcen_3d     - cloud-water first-order loss rate to precipitation (/s) at layer center
    !   precrcen_3d      - liquid (rain) precipitation rate (kg/m2/s) at layer center
    !   precsolidcen_3d  - solid (snow,graupel,...) precipitation rate (kg/m2/s) at layer center
    !-----------------------------------------------------------------------

    real(r8), intent(inout), dimension( pcols, pver)  :: acldy_cen_tbeg_3d
    !-----------------------------------------------------------------------
    ! acldy_cen_tbeg_3d = total (all sub-classes) cloudy fractional area
    !   on input,  = value from end of the previous time step
    !   on output, = value from end of the current  time step

    !-----------------------------------------------------------------------
    ! local variables
    integer :: ncol, lchnk
    integer :: mbuf
    integer :: id
    integer :: i, icc, ipass, ipp, itmpa, it, ichem, ichem2
    integer :: j, jclrcld, jcls, jclsaa, jclsbb, jt
    integer :: nstep, nstep_pp
    integer :: k, ka, kb, lk
    integer :: l, ll, levdbg_err, levdbg_info
    integer :: lun, lun60, lun61, lun131, lun132, lun133, lun134, lun135
    integer :: n, ncls_ecpp, nupdraft, ndndraft
    integer :: itmpcnt(pver+1,4)
    integer :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

    integer, dimension( 1:2, 1:maxcls_ecpp ) :: kdraft_bot_ecpp
    integer, dimension( 1:2, 1:maxcls_ecpp ) :: kdraft_top_ecpp
    integer, dimension( 1:2, 1:maxcls_ecpp ) :: mtype_updnenv_ecpp
    
    real(r8) :: dtstep, dtstep_pp
    real(r8) :: tmpa, tmpb, tmpc, tmpd
    real(r8) :: za, zb, zc

    integer, dimension( 1:nupdraft_in ) :: kupdraftbase
    integer, dimension( 1:nupdraft_in ) :: kupdrafttop
    integer, dimension( 1:ndndraft_in ) :: kdndraftbase
    integer, dimension( 1:ndndraft_in ) :: kdndrafttop

    !-----------------------------------------------------------------------
    ! kupdraftbase, kupdrafttop - lower-most and upper-most level for each updraft class
    ! *** note1- these refer to layer centers, not layer boundaries.  Thus
    !     acen > 0 for kupdraftbase:kupdrafttop and = 0 at other k
    !     abnd > 0 for kupdraftbase+1:kupdrafttop and = 0 at other k
    !     massflxbnd > 0 for kupdraftbase+1:kupdrafttop and = 0 at other k
    ! kdndraftbase, kdndrafttop - lower-most and upper-most level for each downdraft class
    ! *** note2- these get checked/adjusted later, so simply setting k--draftbase = kts
    !     and k--drafttop = ktecen is OK
    !-----------------------------------------------------------------------

    real(r8)  ::  tcen_bar   (pver)                 ! temperature at layer centers (K)  
    real(r8)  ::  pcen_bar   (pver)                 ! pressure at layer centers (K)
    real(r8)  ::  rhocen_bar (pver)                 ! air density at layer centers (kg/m3)
    real(r8)  ::  dzcen      (pver)                 ! layer depth (m)
    real(r8)  ::  wcen_bar   (pver)                 ! vertical velocity at layer centers (m/s)
    real(r8)  ::  rhobnd_bar (pverp)                ! air density at layer boundaries (kg/m3)
    real(r8)  ::  zbnd       (pverp)                ! elevation at layer boundaries         (m) 
    real(r8)  ::  wbnd_bar   (pverp)                ! vertical velocity at layer boundaries (m/s)  
    real(r8)  ::  chem_bar (pver, 1:num_chem_ecpp)  ! mixing ratios of trace gase (ppm) and aerosol species 
                                                    ! ( ug/kg for mass species, #/kg for number species )
#ifdef MODAL_AERO
    ! real(r8), pointer, dimension(:, :, :) :: qqcw  ! cloud-borne aerosol
    type(ptr2d_t) :: qqcw(pcnst)
    ! real(r8) :: qqcwold(pcols, pver, pcnst)
#endif
    real(r8), dimension(pverp,0:2,0:maxcls_ecpp)      :: abnd_tavg
    real(r8), dimension(pverp,0:2,0:maxcls_ecpp)      :: abnd_tfin
    real(r8), dimension(pverp,0:2,0:maxcls_ecpp)      :: mfbnd
    real(r8), dimension(pver ,0:2,0:maxcls_ecpp)      :: acen_tavg
    real(r8), dimension(pver ,0:2,0:maxcls_ecpp)      :: acen_tfin
    real(r8), dimension(pver ,0:2,0:maxcls_ecpp)      :: acen_tbeg
    real(r8), dimension(pver ,0:2,0:maxcls_ecpp)      :: acen_prec
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2)  :: rh_sub2
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2)  :: qcloud_sub2
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2)  :: qlsink_sub2
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2)  :: precr_sub2
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2)  :: precs_sub2
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2,1:num_chem_ecpp )  :: del_cldchem       ! tendency of chem_sub from aqueous chemistry
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2,1:num_chem_ecpp )  :: del_rename        ! tendency of chem_sub from renaming.
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2,1:num_chem_ecpp )  :: del_wetscav       ! tendency of chem_sub from wet deposition
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,1:2,1:num_chem_ecpp )  :: del_wetresu       !
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,    1:num_chem_ecpp )  :: del_activate      ! tendency of chem_sub from activation/resuspension
    real(r8), dimension(pver ,1:2,1:maxcls_ecpp,    1:num_chem_ecpp )  :: del_conv          ! tendency of chem_sub from convective transport

    real(r8), dimension(pcols,pver,1:2,1:maxcls_ecpp,1:2,1:num_chem_ecpp) :: del_cldchem3d   ! tendency of chem_sub from aqueous chemistry
    real(r8), dimension(pcols,pver,1:2,1:maxcls_ecpp,1:2,1:num_chem_ecpp) :: del_rename3d    ! tendency of chem_sub from renaming.
    real(r8), dimension(pcols,pver,1:2,1:maxcls_ecpp,1:2,1:num_chem_ecpp) :: del_wetscav3d   ! tendency of chem_sub from wet deposition
    real(r8), dimension(pcols,pver,1:2,1:maxcls_ecpp,1:2,1:num_chem_ecpp) :: del_wetresu3d   ! tendency of chem_sub from resuspension in wet deposition
    real(r8), dimension(pcols,pver,1:2,1:maxcls_ecpp,    1:num_chem_ecpp) :: del_activate3d  ! tendency of chem_sub from activation/resuspension
    real(r8), dimension(pcols,pver,1:2,1:maxcls_ecpp,    1:num_chem_ecpp) :: del_conv3d      ! tendency of chem_sub from convective transport

    real(r8), dimension(pcols)  :: aqso4_h2o2  ! SO4 aqueous phase chemistry due to H2O2 (kg/m2/s)
    real(r8), dimension(pcols)  :: aqso4_o3    ! SO4 aqueous phase chemistry due to O3 (kg/m2/s)

    real(r8), dimension(pver,      1:2,1:maxcls_ecpp,1:2) :: xphlwc         ! pH value multiplied by lwc
    real(r8), dimension(pcols,pver,1:2,1:maxcls_ecpp,1:2) :: xphlwc3d
    real(r8), dimension(pcols,pver)                       :: xphlwc_gcm

    real(r8), dimension(pcols,pver,1:num_chem_ecpp) :: ptend_cldchem      ! tendency at GCM grid
    real(r8), dimension(pcols,pver,1:num_chem_ecpp) :: ptend_rename       ! tendency at GCM grid
    real(r8), dimension(pcols,pver,1:num_chem_ecpp) :: ptend_wetscav      ! tendency at GCM grid
    real(r8), dimension(pcols,pver,1:num_chem_ecpp) :: ptend_wetresu      ! tendency at GCM grid
    real(r8), dimension(pcols,pver,1:num_chem_ecpp) :: ptend_activate     ! tendency at GCM grid
    real(r8), dimension(pcols,pver,1:num_chem_ecpp) :: ptend_conv         ! tendency at GCM grid

    real(r8), dimension(pcols,pver,1:maxcls_ecpp,1:num_chem_ecpp) :: ptend_activate_cls   ! activation tendency for sub transport class
    real(r8), dimension(pcols,pver,1:maxcls_ecpp,1:num_chem_ecpp) :: ptend_cldchem_cls    ! aqueous chemistry
    real(r8), dimension(pcols,pver,1:maxcls_ecpp,1:num_chem_ecpp) :: ptend_rename_cls     ! renaming
    real(r8), dimension(pcols,pver,1:maxcls_ecpp,1:num_chem_ecpp) :: ptend_wetscav_cls    ! wet deposition
    real(r8), dimension(pcols,pver,1:maxcls_ecpp,1:num_chem_ecpp) :: ptend_wetresu_cls    ! resuspension
    real(r8), dimension(pcols,pver,1:maxcls_ecpp,1:num_chem_ecpp) :: ptend_conv_cls       ! convective transport
                                

    real(r8), dimension(pcols, 1:maxcls_ecpp, 1:num_chem_ecpp) :: ptend_activate_cls_col  ! column-integrated activation tendency for sub transport class
    real(r8), dimension(pcols, 1:maxcls_ecpp, 1:num_chem_ecpp) :: ptend_cldchem_cls_col   ! aqueous chemistry
    real(r8), dimension(pcols, 1:maxcls_ecpp, 1:num_chem_ecpp) :: ptend_rename_cls_col    ! renaming
    real(r8), dimension(pcols, 1:maxcls_ecpp, 1:num_chem_ecpp) :: ptend_wetscav_cls_col   ! wet deposition
    real(r8), dimension(pcols, 1:maxcls_ecpp, 1:num_chem_ecpp) :: ptend_wetresu_cls_col   ! resuspension
    real(r8), dimension(pcols, 1:maxcls_ecpp, 1:num_chem_ecpp) :: ptend_conv_cls_col      ! convective transport


    real(r8), dimension(pcols, 1:num_chem_ecpp)  :: ptend_cldchem_col     ! column-integrated tendency 
    real(r8), dimension(pcols, 1:num_chem_ecpp)  :: ptend_rename_col      ! column-integrated tendency 
    real(r8), dimension(pcols, 1:num_chem_ecpp)  :: ptend_wetscav_col     ! column-integrated tendency 
    real(r8), dimension(pcols, 1:num_chem_ecpp)  :: ptend_wetresu_col     ! column-integrated tendency 
    real(r8), dimension(pcols, 1:num_chem_ecpp)  :: ptend_activate_col    ! column-integrated tendency 
    real(r8), dimension(pcols, 1:num_chem_ecpp)  :: ptend_conv_col        ! column-integrated tendency 
    real(r8), dimension(pcols, 1:num_chem_ecpp)  :: ptendq_col            ! column-integrated tendency 

    real(r8), dimension(pcols, pver, 1:pcnst)  :: ptend_qqcw              ! tendency for cloud-borne aerosols

    real(r8), dimension(pcols, 1:num_chem_ecpp) :: del_chem_col_cldchem   ! column tendency calcuated in ECPP
    real(r8), dimension(pcols, 1:num_chem_ecpp) :: del_chem_col_rename    ! column tendency calcuated in ECPP
    real(r8), dimension(pcols, 1:num_chem_ecpp) :: del_chem_col_wetscav   ! column tendency calcuated in ECPP

    character(len=100) :: msg
    logical ::  lq(pcnst)

    !-----------------------------------------------------------------------
    !   set flags that turn diagnostic output on/off
    !
    !   for a specific output to be "on", both the 
    ! idiagaa_ecpp(--) and ldiagaa_ecpp(--) be positive
    !   the ldiagaa_ecpp(--) is the output unit number
    !
    !    60 -  from subr parampollu_driver2
    ! short messages on entry and exit
    !    61 -  from subr parampollu_driver2
    ! "rcetestpp diagnostics" block
    !    62 - from subr parampollu_td240clm
    ! short messages on entry and exit, and showing sub-time-step
    !    63 - from subr parampollu_check_adjust_inputs 
    ! shows some summary statistics about the check/adjust process
    !   115, 116, 117 - from subr parampollu_1clm_dumpaa
    ! shows various statistics on transport class and subarea
    ! fractional areas and mass fluxes
    ! 116 is before    call to parampollu_check_adjust_inputs
    ! 117 is after 1st call to parampollu_check_adjust_inputs
    ! 115 is after 2nd call to parampollu_check_adjust_inputs
    !   118 - from subr parampollu_tdx_main_integ and parampollu_tdx_area_change
    ! diagnostics involving changes to species 9 in those subrs
    !   119 - from subr parampollu_tdx_cleanup
    ! diagnostics involving changes to species 9 in that subr
    !   121 - from subr parampollu_tdx_cleanup
    ! diagnostics involving mass conservation
    !   122 - from subr parampollu_tdx_entdet_sub1 and parampollu_tdx_entdet_diag01
    ! diagnostics involving entrainment/detrainment and area changes
    !   123 - from subr parampollu_tdx_entdet_sub1
    ! diagnostics involving entrainment/detrainment and area changes
    !   124 - from subr parampollu_tdx_main_integ
    ! diagnostics involving sub-time-step for "main integration", 
    ! related to stability and courant number
    !   125 - from subr parampollu_tdx_activate_intface
    ! diagnostics involving aerosol activation and associated vertical velocities
    !   131-135 -  from subr parampollu_driver2
    ! shows various statistics on transport class and subarea
    ! fractional areas and mass fluxes
    !   141-143 -  from subr parampollu_tdx_wetscav_2
    !       diagnostics for the "new" wetscav code designed for the mmf-with-ecpp
    !   155 - from subr parampollu_check_adjust_inputs 
    ! shows "history" of acen_tavg_use thru the check/adjust process
    !   161, 162, 164  - from subr parampollu_tdx_startup & parampollu_tdx_partition_acw
    ! involves partitioning of cloudborne/interstitial aerosol between clear
    ! and cloudy subareas


    idiagaa_ecpp(:) = 0
    ! idiagaa_ecpp(60:63) = 1
    idiagaa_ecpp(60:63) = -1 
    idiagaa_ecpp(115:119) = 1 ; idiagaa_ecpp(118) = 111
    idiagaa_ecpp(121:125) = 1
    idiagaa_ecpp(131:135) = 1
    idiagaa_ecpp(141:143) = 1
    !==Guangxing Lin
    !idiagaa_ecpp(155) = 1
    idiagaa_ecpp(155) = -1
    !==Guangxing Lin
    idiagaa_ecpp(161) = 1 ; idiagaa_ecpp(162) = 1 ; idiagaa_ecpp(164) = 1

    idiagaa_ecpp(131:135) = -1  ! not output in the MMF model
    idiagaa_ecpp(115:117) = -1  ! not dump the original field in parampollu_td240clm
    idiagaa_ecpp(118:119) = -1 
    idiagaa_ecpp(121:125) = -1
    idiagaa_ecpp(141:143) = -1
    idiagaa_ecpp(165:167) = -1
    idiagaa_ecpp(164) = -1
    idiagaa_ecpp(161) = -1 
    idiagaa_ecpp(162) = -1

    idiagaa_ecpp(121) = -1

    do i = 1, 199
      ldiagaa_ecpp(i) = i
    end do
    ldiagaa_ecpp(60:69) = 6
    ldiagaa_ecpp(62) = 62

    !-----------------------------------------------------------------------

    lun60 = -1
    lun61 = -1
    lun131 = -1
    lun132 = -1
    lun133 = -1
    lun134 = -1
    lun135 = -1
    
    if (idiagaa_ecpp(60)  > 0) lun60  = ldiagaa_ecpp(60)
    if (idiagaa_ecpp(61)  > 0) lun61  = ldiagaa_ecpp(61)
    if (idiagaa_ecpp(131) > 0) lun131 = ldiagaa_ecpp(131)
    if (idiagaa_ecpp(132) > 0) lun132 = ldiagaa_ecpp(132)
    if (idiagaa_ecpp(133) > 0) lun133 = ldiagaa_ecpp(133)
    if (idiagaa_ecpp(134) > 0) lun134 = ldiagaa_ecpp(134)
    if (idiagaa_ecpp(135) > 0) lun135 = ldiagaa_ecpp(135)

        
    ncol = state%ncol
    lchnk = state%lchnk

    !------------------------------------------------------
    ! Initialize ptend
    !------------------------------------------------------
    lq(:) = .true.
    call physics_ptend_init(ptend, state%psetcols,'ecpp',lq=lq)
    ptend%lq(:) = .true.
    ptend%q(:,:,:) = 0.0_r8
    !------------------------------------------------------
    !------------------------------------------------------

    dtstep = dtstep_in
    dtstep_pp = dtstep_pp_in

    !rcetestpp diagnostics --------------------------------------------------
    if (lun61 > 0) then
      write(lun61,93010) ' '
      write(lun61,93010) 'rcetestpp diagnostics from parampollu_driver2'
      write(lun61,93020) 'dtstep, dtstep_pp              ',   &
      dtstep, dtstep_pp
    end if ! (lun61 > 0)
    !rcetestpp diagnostics --------------------------------------------------
93010   format( a, 8(1x,i6) )
93020   format( a, 8(1p,e14.6) )


    if (num_chem_ecpptmp < num_chem_ecpp)  then 
      msg = '*** parampollu_driver -- bad num_chem_ecpptmp'
      call endrun(msg)
    end if

    ! check for valid ncls_ecpptmp
    nupdraft = nupdraft_in
    ndndraft = ndndraft_in
    ncls_ecpp = (nupdraft + ndndraft + 1)
    if (ncls_ecpp > maxcls_ecpp) then
      write(msg,'(a,2(1x,i6))')   &
          '*** parampollu_driver - ncls_ecpp > maxcls_ecpp, values =',   &
          ncls_ecpp, maxcls_ecpp
      call endrun( msg )
    end if
    if (ncls_ecpp /= ncls_ecpp_in) then
      write(msg,'(a,2(1x,i8))') &
        '*** parampollu_driver -- bad ncls_ecpp_in', &
        ncls_ecpp_in, ncls_ecpp
      call endrun( msg )
    end if

    ! on very first time step, initialize acldy_cen_tbeg
    !
    ! *** this code should probably go into parampollu_init0 (or somewhere else)
    nstep = get_nstep()
    nstep_pp = nstep
    if (is_first_step()) then
      acldy_cen_tbeg_3d(:,:) = 0.0

      do k=1,pver
        do i=1,ncol 
          tmpa = 0.0 ; tmpb = 0.0
          do ipp = 1, nprcp_in
          do jcls = 1, ncls_ecpp
              tmpa = tmpa + max( 0.0_r8, acen_3d(i,k,1,jcls,ipp) )
              tmpb = tmpb + max( 0.0_r8, acen_3d(i,k,2,jcls,ipp) )
          end do
          end do

          if (abs(tmpa+tmpb-1.0_r8) > 1.0e-3_r8) then
              write(msg,'(a,3i5,1pe15.7)') &
            '*** parampollu_driver -- bad acen_tbeg - i,j,k,acen', &
            i, j, k, (tmpa+tmpb)
              call endrun(msg)
          end if
          tmpa = tmpa/(tmpa+tmpb)

          tmpa = 1.0_r8   ! force to initially clear -- might want to change this

          ! when iflag_ecpp_test_fixed_fcloud = 2/3/4/5, 
          ! force acen_tbeg 100%/0%/70%/30% clear 
          if ((iflag_ecpp_test_fixed_fcloud >= 2) .and. &
              (iflag_ecpp_test_fixed_fcloud <= 5)) then
              if      (iflag_ecpp_test_fixed_fcloud == 2) then
            tmpa = 1.0_r8
              else if (iflag_ecpp_test_fixed_fcloud == 3) then
            tmpa = 0.0_r8
              else if (iflag_ecpp_test_fixed_fcloud == 4) then
            tmpa = 0.7_r8
              else
            tmpa = 0.3_r8
              end if
          end if

          acldy_cen_tbeg_3d(i,k) = 1.0_r8 - tmpa
        end do ! i=1,ncol 
      end do ! k=1,pver
    end if ! is_first_step()


    ! set some variables to their wrf-chem "standard" values
    levdbg_err = 0
    levdbg_info = 15

#ifdef MODAL_AERO
    ! mbuf = pbuf_get_fld_idx( 'QQCW' )
    ! if ( associated(pbuf(mbuf)%fld_ptr) ) then
    !   qqcw => pbuf(mbuf)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:pcnst )
    ! else
    !   call endrun( 'pbuf for QQCW not allocated in aerosol_wet_intr' )
    ! end if
    !+++mhwang 2012-02-22
    ! qqcw_get_field is no longer used in ndrop.F90. Make sure
    ! it is still valid !!!!
    do i=1,pcnst
      qqcw(i)%fldcw   =>  qqcw_get_field(pbuf, i,lchnk,.true.)
    end do
#endif
    
    !---------------------------------------------------------------
    ! Begin loop over columns
    !---------------------------------------------------------------
    do i=1,ncol
      !
      ! load column arrays
      !
      zbnd(1) = 0.0_r8
      wbnd_bar(1) = 0.0_r8
      do k=pver,1,-1
        tcen_bar(pver-k+1)   = state%t(i,k)
        pcen_bar(pver-k+1)   = state%pmid(i,k)
        rhocen_bar(pver-k+1) = state%pmiddry(i,k)/(287.0*state%t(i,k))    ! dry air density is calcualted, because tracer mixing ratios
                                                                          ! are defined with respect to dry air in CAM.  
        wbnd_bar(pver-k+2)   = -1*state%omega(i,k)/(rhocen_bar(pver-k+1)*gravit)   ! pressure vertical velocity (Pa/s) to height vertical velocity (m/s)
        dzcen(pver-k+1)      = state%pdeldry(i,k)/gravit/rhocen_bar(pver-k+1)
        zbnd(pver-k+2)       = zbnd(pver-k+1) + dzcen(pver-k+1)
      end do ! k=pver,1,-1

      do k = 1, pver+1
        ka = max( 1, min(pver-1, k-1 ) )
        kb = ka + 1
        za = 0.5*(zbnd(ka) + zbnd(ka+1))
        zb = 0.5*(zbnd(kb) + zbnd(kb+1))
        rhobnd_bar(k) = rhocen_bar(ka)    &
        + (rhocen_bar(kb)-rhocen_bar(ka))*(zbnd(k)-za)/(zb-za)
      end do ! k = 1, pver+1

      chem_bar(:,:) = 0.0
      ! Load chem
      do k=pver,1,-1
        do ichem=1,num_chem_ecpp
          if(ichem.le.pcnst) then
            chem_bar(pver-k+1, ichem) = state%q(i, k, ichem)
#ifdef MODAL_AERO
          else
            ! chem_bar(pver-k+1, ichem) = qqcw(i, k, ichem-pcnst)
            if(associated(qqcw(ichem-pcnst)%fldcw)) then
              chem_bar(pver-k+1, ichem) = qqcw(ichem-pcnst)%fldcw(i, k)
            else
              chem_bar(pver-k+1, ichem) = 0.0
            end if
#endif
          end if
        end do ! ichem=1,num_chem_ecpp
      end do ! k=pver,1,-1

      !
      ! load transport-class arrays
      !

      ! load other/quiescent
      jcls = 1

      kupdraftbase = 1
      kupdrafttop  = pver  
      kdndraftbase = 1
      kdndrafttop  = pver

      kdraft_bot_ecpp(   1:2,jcls) = 1
      kdraft_top_ecpp(   1:2,jcls) = pver 
      mtype_updnenv_ecpp(1:2,jcls) = mtype_quiescn_ecpp

      ! load updrafts
      do n=1,nupdraft
        jcls = jcls + 1

        kdraft_bot_ecpp(   1:2,jcls) = max( kupdraftbase(n), 1 )
#ifdef ECPP_LEV_MOD
        kdraft_top_ecpp(   1:2,jcls) = min( kupdrafttop(n), crm_nz )
#else
        kdraft_top_ecpp(   1:2,jcls) = min( kupdrafttop(n), pver )
#endif
        mtype_updnenv_ecpp(1:2,jcls) = mtype_updraft_ecpp
      end do ! n=1,nupdraft

      ! load downdrafts
      do n=1,ndndraft
        jcls = jcls + 1

        kdraft_bot_ecpp(   1:2,jcls) = max( kdndraftbase(n), 1 )
        kdraft_top_ecpp(   1:2,jcls) = min( kdndrafttop(n), pver )
        mtype_updnenv_ecpp(1:2,jcls) = mtype_dndraft_ecpp
      end do ! n=1,ndndraft

      ! load mfbnd and "area" arrays for all classes 
      mfbnd(    :,:,:) = 0.0
      abnd_tavg(:,:,:) = 0.0
      abnd_tfin(:,:,:) = 0.0
      acen_tavg(:,:,:) = 0.0
      acen_tfin(:,:,:) = 0.0

      do jcls = 1, ncls_ecpp
      do icc = 1, 2
      do k = 1, pver+1 
                lk=pver+1-k+1
          mfbnd(    lk,icc,jcls) = massflxbnd_3d(i, k,icc,jcls,1) &
                                + massflxbnd_3d(i, k,icc,jcls,2)
          abnd_tavg(lk,icc,jcls) = abnd_3d(i, k,icc,jcls,1) &
                                + abnd_3d(i, k,icc,jcls,2)
          abnd_tfin(lk,icc,jcls) = abnd_tf_3d(i, k,icc,jcls,1) &
                                + abnd_tf_3d(i, k,icc,jcls,2)
      end do ! k
      end do ! icc
      end do ! jcls

      ! load these arrays
      acen_prec(  :,:,:  ) = 0.0
      qcloud_sub2(:,:,:,:) = 0.0
      qlsink_sub2(:,:,:,:) = 0.0
      precr_sub2( :,:,:,:) = 0.0
      precs_sub2( :,:,:,:) = 0.0
      rh_sub2(    :,:,:,:) = 0.0
      do k=1,pver
        lk=pver-k+1
        acen_tavg(  lk,1:2,1:ncls_ecpp    ) = acen_3d(i, k,1:2,1:ncls_ecpp,1)+ &
                                              acen_3d(i, k,1:2,1:ncls_ecpp,2)
        acen_tfin(  lk,1:2,1:ncls_ecpp    ) = acen_tf_3d(i, k,1:2,1:ncls_ecpp,1)+ &
                                              acen_tf_3d(i, k,1:2,1:ncls_ecpp,2)
        acen_prec(  lk,1:2,1:ncls_ecpp    ) = acen_3d(i, k,1:2,1:ncls_ecpp,2)
        qcloud_sub2(lk,1:2,1:ncls_ecpp,1:2) = qcloudcen_3d(i, k,1:2,1:ncls_ecpp,1:2)
        qlsink_sub2(lk,1:2,1:ncls_ecpp,1:2) = qlsinkcen_3d(i, k,1:2,1:ncls_ecpp,1:2)
        precr_sub2( lk,1:2,1:ncls_ecpp,1:2) = precrcen_3d(i,  k,1:2,1:ncls_ecpp,1:2)
        precs_sub2( lk,1:2,1:ncls_ecpp,1:2) = precsolidcen_3d(i, k,1:2,1:ncls_ecpp,1:2)
        rh_sub2(    lk,1:2,1:ncls_ecpp,1:2) = rhcen_3d(i, k,1:2,1:ncls_ecpp,1:2)
        if( sum(acen_tfin(  lk,1:2,jcls_qu)).lt.0.05) then
          write(0, *) 'test acen_tfin < 0.40', sum(acen_tfin(  lk,1:2,jcls_qu)), pcen_bar(lk), i,lk  !+++mhwang
        end if
      end do ! k=1,pver

      ! force kdraft_top > kdraft_bot
      ! (note:  need to change the wrf3d post-processor so this is not needed)
      do jcls=1,ncls_ecpp
        do jclrcld=1,2
          kdraft_top_ecpp(jclrcld,jcls) = max( kdraft_top_ecpp(jclrcld,jcls),   &
                                               kdraft_bot_ecpp(jclrcld,jcls)+1 )
        end do ! jclrcld=1,2
      end do ! jcls=1,ncls_ecpp

      ! load acen_tbeg from 3d saved values
      acen_tbeg(:,:,:) = 0.0
      jcls = 1
      do k=1,pver
        lk=pver-k+1
        acen_tbeg(lk,2,jcls) = acldy_cen_tbeg_3d(i,k)
        acen_tbeg(lk,1,jcls) = 1.0_r8 - acen_tbeg(lk,2,jcls)
      end do

      !----------------------------------------------------------------
      !   start of temporary diagnostics ------------------------------
      !----------------------------------------------------------------
      do ipass = 1, 3

        do ll = 131, 133
          lun = -1
          if (ll == 131) lun = lun131
          if (ll == 132) lun = lun132
          if (ll == 133) lun = lun133
          if (lun <= 0) cycle

          write(lun,*)
          if (ipass .eq. 1) then
            n = nupdraft
            write(lun,'(a,3i5)') 'updrafts, nup, ktau', n, nstep,  nstep_pp
          else if (ipass .eq. 2) then
            n = ndndraft
            write(lun,'(a,3i5)') 'dndrafts, nup, ktau', n, nstep, nstep_pp
          else
            n = ncls_ecpp
            write(lun,'(a,3i5)') 'quiescents, ncls_ecpp, ktau', n, nstep, nstep_pp 
          end if
        end do

        do ka = (2*((pver+1)/2)-1), 1, -2
          tmpa = 0.0
          tmpb = 0.0
          tmpc = 0.0
          tmpd = 0.0
          kb = ka+1
          ! kb = ka

          if (ipass .eq. 1) then
            jclsaa = 1 + 1
            jclsbb = 1 + nupdraft
          else if (ipass .eq. 2) then
            jclsaa = 1 + nupdraft + 1
            jclsbb = 1 + nupdraft + ndndraft
          else
            jclsaa = 1
            jclsbb = 1
          end if
          do ipp = 1, 2
            do jcls = jclsaa, jclsbb
              tmpa = tmpa + abnd_3d(i,ka,1,jcls,ipp) + abnd_3d(i,kb,1,jcls,ipp) 
              tmpb = tmpb + abnd_3d(i,ka,2,jcls,ipp) + abnd_3d(i,kb,2,jcls,ipp) 
              tmpc = tmpc + massflxbnd_3d(i,ka,1,jcls,ipp) + massflxbnd_3d(i,kb,1,jcls,ipp) 
              tmpd = tmpd + massflxbnd_3d(i,ka,2,jcls,ipp) + massflxbnd_3d(i,kb,2,jcls,ipp) 
            end do
          end do

          tmpa = tmpa*0.5 ; tmpb = tmpb*0.5 ;
          tmpc = tmpc*0.5 ; tmpd = tmpd*0.5
          if (lun131 > 0) &
          write(lun131,'(i3,2(3x,1p,3e10.2))') ka,   &
          tmpa, tmpb, (tmpa+tmpb), tmpc, tmpd, (tmpc+tmpd)

          tmpa = tmpa*100.0 ; tmpb = tmpb*100.0
          tmpc = tmpc*100.0 ; tmpd = tmpd*100.0
          if (lun132 > 0) &
          write(lun132,'(i3,2(2x,    3f8.3))') ka,   &
          tmpa, tmpb, (tmpa+tmpb), tmpc, tmpd, (tmpc+tmpd)

          if (lun133 > 0) &
          write(lun133,'(i3,2(2x,    3f7.2))') ka,   &
          tmpa, tmpb, (tmpa+tmpb), tmpc, tmpd, (tmpc+tmpd)
        end do ! ka
      end do ! ipass

      if (lun134 > 0) then
        do n = 1, nupdraft
          write(lun134,'(/a,5i5)') 'updraft -- n, kbase, ktop, ktaus',   &
                n, kupdraftbase(n), kupdrafttop(n), nstep, nstep_pp
          do k = pver+1, 1, -1
            jcls = 1 + n
            write(lun134,'(i3,2(2x,2f10.5))') k,   &
                  sum(abnd_3d(i,k,1,jcls,1:2))*100.0, sum(abnd_3d(i,k,2,jcls,1:2))*100.0,   &
                  sum(massflxbnd_3d(i,k,1,jcls,1:2))*100.0,   &
                  sum(massflxbnd_3d(i,k,2,jcls,1:2))*100.0
          end do
        end do

        do n = 1, ndndraft
          write(lun134,'(/a,5i5)') 'dndraft -- n, kbase, ktop, ktaus',   &
                n, kdndraftbase(n), kdndrafttop(n), nstep, nstep_pp
          do k = pver+1, 1, -1
            jcls = 1 + nupdraft + n
            write(lun134,'(i3,2(2x,2f10.5))') k,   &
                  sum(abnd_3d(i,k,1,jcls,1:2))*100.0, sum(abnd_3d(i,k,2,jcls,1:2))*100.0,   &
                  sum(massflxbnd_3d(i,k,1,jcls,1:2))*100.0,   &
                  sum(massflxbnd_3d(i,k,2,jcls,1:2))*100.0
          end do
        end do
      end if ! (lun134 > 0)


      if (lun135 > 0) then
        itmpcnt(:,:) = 0
        do n = 1, nupdraft
          write(lun135,'(/a,5i5)') 'updraft -- n, kbase, ktop, ktaus',   &
                n, kupdraftbase(n), kupdrafttop(n), nstep, nstep_pp
          do k = pver+1, 1, -1
            jcls = 1 + n
            tmpa = sum(abnd_3d(i,k,1,jcls,1:2))
            tmpb = sum(abnd_3d(i,k,2,jcls,1:2))
            tmpc = sum(massflxbnd_3d(i,k,1,jcls,1:2))
            tmpd = sum(massflxbnd_3d(i,k,2,jcls,1:2))
            write(lun135,'(i3,2(2x,1p,2e10.2))') k, tmpa, tmpb, tmpc, tmpd
            if (tmpa .gt. 0.0) itmpcnt(k,1) = itmpcnt(k,1) + 1
            if (tmpb .gt. 0.0) itmpcnt(k,2) = itmpcnt(k,2) + 1
            if (tmpc .gt. 0.0) itmpcnt(k,3) = itmpcnt(k,3) + 1
            if (tmpd .gt. 0.0) itmpcnt(k,4) = itmpcnt(k,4) + 1
          end do
        end do
        write(lun135,'(/a,5i5)') 'updraft non-zero counts -- ktaus',   &
              nstep, nstep_pp
        do k = pver+1, 1, -1
          write(lun135,'(i3,2(5x,2i5))') k, itmpcnt(k,1:4)
        end do

        itmpcnt(:,:) = 0
        do n = 1, ndndraft
          write(lun135,'(/a,5i5)') 'dndraft -- n, kbase, ktop, ktaus',   &
                n, kdndraftbase(n), kdndrafttop(n), nstep, nstep_pp
          do k = pver+1, 1, -1
            jcls = 1 + nupdraft + n
            tmpa = sum(abnd_3d(i,k,1,jcls,1:2))
            tmpb = sum(abnd_3d(i,k,2,jcls,1:2))
            tmpc = sum(massflxbnd_3d(i,k,1,jcls,1:2))
            tmpd = sum(massflxbnd_3d(i,k,2,jcls,1:2))
            write(lun135,'(i3,2(2x,1p,2e10.2))') k, tmpa, tmpb, tmpc, tmpd
            if (tmpa .gt. 0.0) itmpcnt(k,1) = itmpcnt(k,1) + 1
            if (tmpb .gt. 0.0) itmpcnt(k,2) = itmpcnt(k,2) + 1
            if (tmpc .lt. 0.0) itmpcnt(k,3) = itmpcnt(k,3) + 1
            if (tmpd .lt. 0.0) itmpcnt(k,4) = itmpcnt(k,4) + 1
          end do
        end do
        write(lun135,'(/a,5i5)') 'dndraft non-zero counts -- ktaus',   &
              nstep, nstep_pp
        do k = pver+1, 1, -1
          write(lun135,'(i3,2(5x,2i5))') k, itmpcnt(k,1:4)
        end do
      end if ! (lun135 > 0)
      !----------------------------------------------------------------
      !   end   of temporary diagnostics ------------------------------
      !----------------------------------------------------------------

      !
      ! do parameterized pollutant calculations on current column
      !
      itmpa = parampollu_opt

      if ((itmpa == 2220) .or.   &
          (itmpa == 2223)) then
        if (lun60 > 0) write(lun60,93010) &
        'calling parampollu_td240clm - i=', i
        ! write (0, *) i, lchnk, 'before parampollu_td240clm', nstep
        call parampollu_td240clm(nstep, dtstep, nstep_pp, dtstep_pp,        &
                                 idiagaa_ecpp, ldiagaa_ecpp,                &
                                 tcen_bar, pcen_bar, rhocen_bar, dzcen,     &
                                 rhobnd_bar, zbnd, wbnd_bar,                &
                                 chem_bar, ncls_ecpp,                       &
                                 kdraft_bot_ecpp, kdraft_top_ecpp,          &
                                 mtype_updnenv_ecpp, mfbnd,                 &
                                 abnd_tavg, acen_tavg,                      & 
                                 acen_tfin, acen_tbeg,                      &
                                 acen_prec, rh_sub2,                        &
                                 qcloud_sub2, qlsink_sub2,                  & 
                                 precr_sub2, precs_sub2,                    &
                                 del_cldchem,  del_rename,                  & 
                                 del_wetscav, del_wetresu,                  &
                                 del_activate, del_conv,                    &
                                 del_chem_col_cldchem(i,:),                 &
                                 del_chem_col_rename(i, :),                 &
                                 del_chem_col_wetscav(i, :),                &
                                 aqso4_h2o2(i), aqso4_o3(i), xphlwc,        &
                                 i, lchnk, 1, pver+1, pver, pbuf            & 
                                )
        ! write (0, *) i, lchnk, 'after parampollu_td240clm', nstep

        aqso4_h2o2(i) = aqso4_h2o2(i)/dtstep
        aqso4_o3(i) = aqso4_o3(i)/dtstep 

      else 

      end if


      !
      ! put selected arrays back into 3d arrays
      !
      if (itmpa > 0) then

        do k = 1, pver 
          lk=pver-k+1
          acldy_cen_tbeg_3d(i,k) = sum( acen_tfin(lk,2,1:ncls_ecpp) )
        end do

        ! Interstial species 
        ptend_qqcw(i,:,:) = 0.0
        do k=1, pver
          lk=pver-k+1 
          do ichem=param_first_ecpp, pcnst 
            ptend%q(i,k,ichem)= (chem_bar(lk, ichem)-state%q(i,k,ichem))/dtstep
            ! ptend_qqcw(i,k,ichem)=(chem_bar(lk, ichem+pcnst)-qqcw(i,k,ichem))/dtstep
            ! qqcw(i,k,ichem) = chem_bar(lk, ichem+pcnst)
            if(associated(qqcw(ichem)%fldcw)) then
              ptend_qqcw(i,k,ichem)=(chem_bar(lk, ichem+pcnst)-qqcw(ichem)%fldcw(i,k))/dtstep
              qqcw(ichem)%fldcw(i,k) = chem_bar(lk, ichem+pcnst)
            else 
              ptend_qqcw(i,k,ichem)= 0.0
            endif 
          end do
          del_cldchem3d(i,k,:,:,:,:) = del_cldchem(lk,:,:,:,:)/dtstep
          del_rename3d(i,k,:,:,:,:) = del_rename(lk,:,:,:,:)/dtstep
          del_wetscav3d(i,k,:,:,:,:) = del_wetscav(lk,:,:,:,:)/dtstep
          del_wetresu3d(i,k,:,:,:,:) = del_wetresu(lk,:,:,:,:)/dtstep
          del_activate3d(i,k,:,:,:) = del_activate(lk,:,:,:)/dtstep
          del_conv3d(i,k,:,:,:) = del_conv(lk,:,:,:)/dtstep
          xphlwc3d(i,k,:,:,:) = xphlwc(lk,:,:,:)
        end do
        ! cloud borne species  
      end if ! (itmpa > 0)

    end do ! i=1,ncol
    !---------------------------------------------------------------
    ! End column loop
    !---------------------------------------------------------------


    ptend_cldchem  = 0.0
    ptend_rename   = 0.0
    ptend_wetscav  = 0.0
    ptend_wetresu  = 0.0
    ptend_activate = 0.0
    ptend_conv     = 0.0
    xphlwc_gcm     = 0.0

    ptend_cldchem_cls  = 0.0
    ptend_rename_cls   = 0.0
    ptend_wetscav_cls  = 0.0
    ptend_wetresu_cls  = 0.0
    ptend_activate_cls = 0.0
    ptend_conv_cls     = 0.0

    ptend_cldchem_col  = 0.0
    ptend_rename_col   = 0.0
    ptend_wetscav_col  = 0.0
    ptend_wetresu_col  = 0.0
    ptend_activate_col = 0.0
    ptend_conv_col     = 0.0
    ptendq_col         = 0.0

    ptend_cldchem_cls_col  = 0.0
    ptend_rename_cls_col   = 0.0
    ptend_wetscav_cls_col  = 0.0
    ptend_wetresu_cls_col  = 0.0
    ptend_activate_cls_col = 0.0
    ptend_conv_cls_col     = 0.0

    do i=1,ncol
      do k=1,pver
        do jcls = 1,ncls_ecpp
          do icc = 1,2
            ! tendency at GCM grids
            do ipp=1, 2
              ptend_cldchem(i,k,:) = ptend_cldchem(i,k,:)+del_cldchem3d(i,k,icc,jcls,ipp,:)
              ptend_rename(i,k,:)  = ptend_rename(i,k,:) +del_rename3d(i,k,icc,jcls,ipp,:)
              ptend_wetscav(i,k,:) = ptend_wetscav(i,k,:)+del_wetscav3d(i,k,icc,jcls,ipp,:)
              ptend_wetresu(i,k,:) = ptend_wetresu(i,k,:)+del_wetresu3d(i,k,icc,jcls,ipp,:)
              xphlwc_gcm(i,k)      = xphlwc_gcm(i,k) + xphlwc3d(i,k,icc,jcls,ipp)
              ! tendency at each transport class:
              ptend_cldchem_cls(i,k,jcls,:) = ptend_cldchem_cls(i,k,jcls,:)+del_cldchem3d(i,k,icc,jcls,ipp,:)
              ptend_rename_cls(i,k,jcls,:)  = ptend_rename_cls(i,k,jcls,:) +del_rename3d(i,k,icc,jcls,ipp,:)
              ptend_wetscav_cls(i,k,jcls,:) = ptend_wetscav_cls(i,k,jcls,:)+del_wetscav3d(i,k,icc,jcls,ipp,:)
              ptend_wetresu_cls(i,k,jcls,:) = ptend_wetresu_cls(i,k,jcls,:)+del_wetresu3d(i,k,icc,jcls,ipp,:)
            end do

            ptend_activate(i,k,:)  = ptend_activate(i,k,:)+del_activate3d(i,k,icc,jcls,:)
            ptend_activate_cls(i,k,jcls, :) = ptend_activate_cls(i,k,jcls, :) + del_activate3d(i,k,icc,jcls,:) 
            ptend_conv(i,k,:)      = ptend_conv(i,k,:)+del_conv3d(i,k,icc,jcls,:)
            ptend_conv_cls(i,k,jcls,:)      = ptend_conv_cls(i,k,jcls,:)+del_conv3d(i,k,icc,jcls,:)
          end do ! icc = 1,2
        end do ! jcls = 1,ncls_ecpp

        ! column-integrated tendency
        ptend_cldchem_col(i,:)  = ptend_cldchem_col(i,:)  + ptend_cldchem(i,k,:)  *state%pdeldry(i,k)/gravit
        ptend_rename_col(i,:)   = ptend_rename_col(i,:)   + ptend_rename(i,k,:)   *state%pdeldry(i,k)/gravit
        ptend_wetscav_col(i,:)  = ptend_wetscav_col(i,:)  + ptend_wetscav(i,k,:)  *state%pdeldry(i,k)/gravit
        ptend_wetresu_col(i,:)  = ptend_wetresu_col(i,:)  + ptend_wetresu(i,k,:)  *state%pdeldry(i,k)/gravit
        ptend_activate_col(i,:) = ptend_activate_col(i,:) + ptend_activate(i,k,:) *state%pdeldry(i,k)/gravit
        ptend_conv_col(i,:)     = ptend_conv_col(i,:)     + ptend_conv(i,k,:)     *state%pdeldry(i,k)/gravit

        ptend_cldchem_cls_col(i,:,:)  = ptend_cldchem_cls_col(i,:,:)  + ptend_cldchem_cls(i,k,:,:)  *state%pdeldry(i,k)/gravit
        ptend_rename_cls_col(i,:,:)   = ptend_rename_cls_col(i,:,:)   + ptend_rename_cls(i,k,:,:)   *state%pdeldry(i,k)/gravit
        ptend_wetscav_cls_col(i,:,:)  = ptend_wetscav_cls_col(i,:,:)  + ptend_wetscav_cls(i,k,:,:)  *state%pdeldry(i,k)/gravit
        ptend_wetresu_cls_col(i,:,:)  = ptend_wetresu_cls_col(i,:,:)  + ptend_wetresu_cls(i,k,:,:)  *state%pdeldry(i,k)/gravit
        ptend_activate_cls_col(i,:,:) = ptend_activate_cls_col(i,:,:) + ptend_activate_cls(i,k,:,:) *state%pdeldry(i,k)/gravit
        ptend_conv_cls_col(i,:,:)     = ptend_conv_cls_col(i,:,:)     + ptend_conv_cls(i,k,:,:)     *state%pdeldry(i,k)/gravit


        ptendq_col(i,param_first_ecpp:pcnst) = ptendq_col(i,param_first_ecpp:pcnst)+      &
                                               ptend%q(i,k,param_first_ecpp:pcnst)*state%pdeldry(i,k)/gravit
        ptendq_col(i,param_first_ecpp+pcnst:pcnst+pcnst) = ptendq_col(i,  param_first_ecpp+pcnst:pcnst+pcnst)+      &
                                                           ptend_qqcw(i,k,param_first_ecpp:pcnst)*state%pdeldry(i,k)/gravit 
      end do ! k=1,pver
    end do ! i=1,ncol

    !---------------------------------------------------------------
    ! History output
    !---------------------------------------------------------------
    do ichem=param_first_ecpp, pcnst 

      call outfld(trim(cnst_name(ichem))//'EP'      , ptend%q(:,:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'ACHEM_EP', ptend_cldchem(:,:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'RENM_EP' , ptend_rename(:,:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'ACT_EP'  , ptend_activate(:,:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'WET_EP'  , ptend_wetscav(:,:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'WRESU_EP', ptend_wetresu(:,:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'CONV_EP' , ptend_conv(:,:,ichem), pcols, lchnk)

      call outfld(trim(cnst_name(ichem))//'SFEP'      , ptendq_col(:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFACHEM_EP', ptend_cldchem_col(:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFRENM_EP' , ptend_rename_col(:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFACT_EP'  , ptend_activate_col(:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFWET_EP'  , ptend_wetscav_col(:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFWRESU_EP', ptend_wetresu_col(:,ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFCONV_EP' , ptend_conv_col(:,ichem), pcols, lchnk)

      call outfld(trim(cnst_name(ichem))//'SFACHQU_EP', ptend_cldchem_cls_col(:,1, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFACHUP_EP', ptend_cldchem_cls_col(:,2, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFACHDN_EP', ptend_cldchem_cls_col(:,3, ichem), pcols, lchnk)

      call outfld(trim(cnst_name(ichem))//'SFREMQU_EP', ptend_rename_cls_col(:,1, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFREMUP_EP', ptend_rename_cls_col(:,2, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFREMDN_EP', ptend_rename_cls_col(:,3, ichem), pcols, lchnk)

      call outfld(trim(cnst_name(ichem))//'SFACTQU_EP', ptend_activate_cls_col(:,1, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFACTUP_EP', ptend_activate_cls_col(:,2, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFACTDN_EP', ptend_activate_cls_col(:,3, ichem), pcols, lchnk)

      call outfld(trim(cnst_name(ichem))//'SFWETQU_EP', ptend_wetscav_cls_col(:,1, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFWETUP_EP', ptend_wetscav_cls_col(:,2, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFWETDN_EP', ptend_wetscav_cls_col(:,3, ichem), pcols, lchnk)

      call outfld(trim(cnst_name(ichem))//'SFRESQU_EP', ptend_wetresu_cls_col(:,1, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFRESUP_EP', ptend_wetresu_cls_col(:,2, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFRESDN_EP', ptend_wetresu_cls_col(:,3, ichem), pcols, lchnk)

      call outfld(trim(cnst_name(ichem))//'SFCONQU_EP', ptend_conv_cls_col(:,1, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFCONUP_EP', ptend_conv_cls_col(:,2, ichem), pcols, lchnk)
      call outfld(trim(cnst_name(ichem))//'SFCONDN_EP', ptend_conv_cls_col(:,3, ichem), pcols, lchnk)

    end do ! ichem=param_first_ecpp, pcnst 
    !---------------------------------------------------------------
    ! More History output
    !---------------------------------------------------------------
    do ichem=param_first_ecpp, pcnst 
      ichem2=ichem+pcnst
      if(.not. (cnst_name_cw(ichem) == ' ')) then
        call outfld(trim(cnst_name_cw(ichem))//'EP', ptend_qqcw(:,:,ichem), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'ACHEM_EP', ptend_cldchem(:,:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'RENM_EP', ptend_rename(:,:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'ACT_EP', ptend_activate(:,:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'WET_EP', ptend_wetscav(:,:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'WRESU_EP', ptend_wetresu(:,:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'CONV_EP', ptend_conv(:,:,ichem2), pcols, lchnk)

        call outfld(trim(cnst_name_cw(ichem))//'SFEP', ptendq_col(:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFACHEM_EP', ptend_cldchem_col(:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFRENM_EP', ptend_rename_col(:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFACT_EP', ptend_activate_col(:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFWET_EP', ptend_wetscav_col(:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFWRESU_EP', ptend_wetresu_col(:,ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFCONV_EP', ptend_conv_col(:,ichem2), pcols, lchnk)

        call outfld(trim(cnst_name_cw(ichem))//'SFACTQU_EP', ptend_activate_cls_col(:,1, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFACTUP_EP', ptend_activate_cls_col(:,2, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFACTDN_EP', ptend_activate_cls_col(:,3, ichem2), pcols, lchnk)

        call outfld(trim(cnst_name_cw(ichem))//'SFACHQU_EP', ptend_cldchem_cls_col(:,1, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFACHUP_EP', ptend_cldchem_cls_col(:,2, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFACHDN_EP', ptend_cldchem_cls_col(:,3, ichem2), pcols, lchnk)

        call outfld(trim(cnst_name_cw(ichem))//'SFREMQU_EP', ptend_rename_cls_col(:,1, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFREMUP_EP', ptend_rename_cls_col(:,2, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFREMDN_EP', ptend_rename_cls_col(:,3, ichem2), pcols, lchnk)

        call outfld(trim(cnst_name_cw(ichem))//'SFWETQU_EP', ptend_wetscav_cls_col(:,1, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFWETUP_EP', ptend_wetscav_cls_col(:,2, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFWETDN_EP', ptend_wetscav_cls_col(:,3, ichem2), pcols, lchnk)

        call outfld(trim(cnst_name_cw(ichem))//'SFRESQU_EP', ptend_wetresu_cls_col(:,1, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFRESUP_EP', ptend_wetresu_cls_col(:,2, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFRESDN_EP', ptend_wetresu_cls_col(:,3, ichem2), pcols, lchnk)

        call outfld(trim(cnst_name_cw(ichem))//'SFCONQU_EP', ptend_conv_cls_col(:,1, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFCONUP_EP', ptend_conv_cls_col(:,2, ichem2), pcols, lchnk)
        call outfld(trim(cnst_name_cw(ichem))//'SFCONDN_EP', ptend_conv_cls_col(:,3, ichem2), pcols, lchnk)

        ! do i=1,ncol
        !  do k=1,pver
        !  ! if(cnst_name_cw(ichem) == 'bc_c1') then
        !  !   if(abs(ptend_wetscav(i, k, ichem2)).gt.1.0e-16 .and. qqcwold(i, k, ichem).gt. 1.0e-13) then
        !  !     if(abs(ptend_conv(i, k, ichem2)).lt.1.0e-20 .and. abs(ptend_activate(i, k, ichem2)).lt.1.0e-20) then
        !  !      write(0, *) 'nstep, ecpp wet, qqcw',  nstep, qqcwold(i, k, ichem), qqcw(i,k,ichem), state%q(i, k, ichem), ptend_wetscav(i, k, ichem2)*1800, ptend_wetscav(i, k, ichem2)*86400/qqcwold(i, k, ichem)
        !  !      write(0, *) 'ecpp acen', acen_3d(i, k,2,1:ncls_ecpp,1), acen_3d(i, k,2,1:ncls_ecpp,2)
        !  !      write(0, *) 'ecpp qlsink' , qlsinkcen_3d(i, k,2,1:ncls_ecpp,1)*86400, qlsinkcen_3d(i, k,2,1:ncls_ecpp,2)*86400
        !  !      write(0, *) 'ecpp wetscav', del_wetscav3d(i,k,2,1:ncls_ecpp,1, ichem2)*1800, del_wetscav3d(i,k,2,1:ncls_ecpp,2, ichem2)*1800
        !  !     call endrun('ptend_conv error')
        !  !    end if
        !  !   end if
        !  !    if(abs(ptend_conv_col(i, ichem2)).gt.1.0e-15) then
        !  !      write(0, *) 'ptend_conv error', ptend_wetresu_col(i,ichem2)+ptend_wetscav_col(i,ichem2),  &
        !  !                 ptend_cldchem_col(i,ichem2), ptend_activate_col(i,ichem2), ptend_conv_col(i,ichem2), ptendq_col(i,ichem2)
        !  !      write(0, *) 'ptend_conv error2' , del_chem_col_wetscav(i, ichem2)/dtstep, del_chem_col_cldchem(i,ichem2)/dtstep
        !  !      write(0, *) 'ptend_conv error3' , ptendq_col(i,ichem2), ptend_wetresu_col(i,ichem2)+ptend_wetscav_col(i,ichem2)  &
        !  !               +ptend_cldchem_col(i,ichem2)+ptend_activate_col(i,ichem2), del_chem_col_wetscav(i, ichem2)/dtstep+ptend_cldchem_col(i,ichem2)+ptend_activate_col(i,ichem2)
        !  !      call endrun('ptend_conv error')
        !  !    end if
        !  ! end if 
        !  end do ! k=1,pver
        ! end do ! i=1,ncol
      end if
    end do ! ichem=param_first_ecpp, pcnst 

    call outfld('AQSO4_H2O2_EP', aqso4_h2o2, pcols, lchnk)
    call outfld('AQSO4_O3_EP', aqso4_o3, pcols, lchnk)
    call outfld('XPH_LWC_EP', xphlwc_gcm, pcols, lchnk)
    !---------------------------------------------------------------
    !---------------------------------------------------------------

    !
    ! qqcw is updated above, and q is upated in tphysbc
    !

    return
  end subroutine parampollu_driver2
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
end module module_ecpp_ppdriver2
