module module_mosaic_cam_init
  
  use shr_kind_mod,            only: r8 => shr_kind_r8

  !---------------------------------------------------------------------------------------!
  !BSINGH: This module initilizes Mosaic chemistry variables.
  !---------------------------------------------------------------------------------------!

  implicit none
  private

  public:: mosaic_cam_init

  integer, parameter :: bigint = huge(1)

contains

#if ( defined MOSAIC_SPECIES )
  subroutine mosaic_cam_init
    !---------------------------------------------------------------------------------------!
    !BSINGH: This subroutine initialzies some Mosaic conastans and inpput parameters
    ! Called by: modal_aero_initialize_data.F90
    !---------------------------------------------------------------------------------------!
    use spmd_utils,                   only: masterproc
    use cam_logfile,                  only: iulog

    use modal_aero_amicphys,          only: max_mode, mam_amicphys_check_mosaic_mw

    use modal_aero_data,              only: &
         modeptr_pcarbon, ntot_amode, &
         lptr_so4_a_amode, lptr_no3_a_amode, lptr_cl_a_amode, &
         lptr_msa_a_amode, lptr_nh4_a_amode, lptr2_soa_a_amode

    use module_data_mosaic_aero,      only: nbin_a_max, nbin_a, mhyst_method, mhyst_force_up, &
         mGAS_AER_XFER, mDYNAMIC_SOLVER, msize_framework, mmodal, alpha_ASTEM, rtol_eqb_ASTEM,      &
         overall_massbal_atoler, overall_massbal_rtoler,      &
         ptol_mol_ASTEM, method_bcfrac, method_kappa, maersize_init_flag1, mcoag_flag1,     &
         ifreq_coag, mmovesect_flag1, mnewnuc_flag1, msectional_flag1, &
         use_cam5mam_soa_params, use_cam5mam_accom_coefs, &
         i_gas2bin_uptk_flag, m_gas2bin_uptk_flag, ngas_aerchtot, &
         ih2so4_g, ihno3_g, ihcl_g, inh3_g, imsa_g, ilim2_g

    use module_data_mosaic_main,      only: ipmcmos,    &
         mgas, maer, mcld, maeroptic, mshellcore, msolar, mphoto

    use module_data_mosaic_asecthp,     only: lunerr, lunout, ntype_md1_aer, ntype_md2_aer

    use module_mosaic_init_aerpar,    only: mosaic_init_aer_params

    ! local variables
    integer :: ibin, iv
    logical :: lgtmpa

    !Initialize Mosaic constants with values from CAM constants
    nbin_a_max = max_mode !*BALLI* Ask Dick about it
    nbin_a     = max_mode !Maximum # of modes is equal to # of bins in Mosaic
    if(masterproc) then
       write(iulog,*) 'mosaic_cam_init: nbin_a_max=', nbin_a_max
    endif

    lunerr = iulog
    lunout = iulog

    use_cam5mam_soa_params  = 1  ! use cam5-mam soa/soag parameter values
    use_cam5mam_accom_coefs = 1  ! use cam5-mam accomodation coefficient values

    !BSINGH - Initialize other constants which sit in the input file of Mosaic
    !and are used in the present code(**BALLI Ask Dick about it)
    mhyst_method    = mhyst_force_up  !rceaster !mhyst_method (1=uporlo_jhyst, 2=uporlo_waterhyst, 3=force_up, 4=force_low)
    mGAS_AER_XFER   = 1    !mGAS_AER_XFER: 1=do gas-aerosol partitioning 0=do not partition
    mDYNAMIC_SOLVER = 1    !mDYNAMIC_SOLVER: 1=astem  2=lsodes
    msize_framework = mmodal  ! rceaster (1=modal, 2=unstructured, 3=sectional)
    alpha_ASTEM     = 0.5  !Solver parameter. range: 0.01 - 1.0
    rtol_eqb_ASTEM  = 0.01 !Relative eqb tolerance. range: 0.01 - 0.03
    ptol_mol_ASTEM  = 0.01 !Percent mol tolerance.  range: 0.01 - 1.0
    overall_massbal_rtoler = 1.0e-4_r8   ! relative tolerance (0-1 frac) for massbal warning message
    overall_massbal_atoler = 1.0e-8_r8  ! absolute tolerance (ug/m^3  ) for massbal warning message
    ipmcmos         = 0    !Additional inputs needed when ipmcmos > 0


    !BSINGH - Initialize constants to 'bigint' which sit in the input file of Mosaic
    !and are NOT used in the present code(**BALLI Ask Dick about it)
    !'bigint' initialized variables will cause the code to halt on their first use

    ntype_md1_aer       = bigint !(number of aerosol types)
    ntype_md2_aer       = bigint !(number of aerosol types)
    method_bcfrac       = bigint !(only used for sectional and ntype>1)
    method_kappa        = bigint !(only used for sectional and ntype>1)
    maersize_init_flag1 = bigint !(only used for sectional and ntype>1)

    mcoag_flag1         = bigint !(only used for sectional)
    ifreq_coag          = bigint !(only used for sectional)
    mmovesect_flag1     = bigint !(only used for sectional)
    mnewnuc_flag1       = bigint !(only used for sectional)
    msectional_flag1    = bigint !(currently not used)

! these variables are now in module_data_mosaic_boxmod, and can be ignored by cam and cambox codes
!   iprint              = bigint !freq of output. Every iprint*dt_min mins.
!   iwrite_gas          = bigint
!   iwrite_aer_bin      = bigint
!   iwrite_aer_dist     = bigint
!   iwrite_aer_species  = bigint
!   mmode               = bigint !: 1=time integration 2=parametric analysis

    mgas                = bigint !: 1=gas chem on,  0=gas chem off**
    maer                = bigint !: 1=aer chem on,  0=aer chem off**
    mcld                = bigint !: 1=cld chem on,  0=cld chem off**
    maeroptic           = bigint !: 1=aer_optical on,  0=aer_optical off **
    mshellcore          = bigint !: 0=no shellcore,  1=core is BC only,  2=core is BC and DUST **
    msolar              = bigint !: 1=diurnally varying phot, 2=fixed phot**
    mphoto              = bigint !: 1=Rick's param 2=Yang's param**

    call mosaic_init_aer_params

    call mam_amicphys_check_mosaic_mw

! set the m_gas2bin_uptk_flag=0, which indicates that some gases cannot condense onto every bin/mode
    m_gas2bin_uptk_flag = 0
! set the i_gas2bin_uptk_flag, which control which semi-volatile gases
!    are allowed to condense on a bin
! h2so4, hno3, hcl, and nh3 can condense onto all modes
! soag currently only condenses onto modes 1-3
    if ( .not. allocated( i_gas2bin_uptk_flag ) ) then
       allocate( i_gas2bin_uptk_flag(ngas_aerchtot,nbin_a) )
    endif
    i_gas2bin_uptk_flag(1:ngas_aerchtot,1:nbin_a) = 0
    do ibin = 1, ntot_amode
       lgtmpa = .false. ; if (ibin == modeptr_pcarbon) lgtmpa = .true.
       if (lptr_so4_a_amode(ibin)    > 0 .or. lgtmpa) i_gas2bin_uptk_flag(ih2so4_g,ibin) = 1
       if (lptr_no3_a_amode(ibin)    > 0 .or. lgtmpa) i_gas2bin_uptk_flag(ihno3_g, ibin) = 1
       if (lptr_cl_a_amode( ibin)    > 0 .or. lgtmpa) i_gas2bin_uptk_flag(ihcl_g,  ibin) = 1
       if (lptr_nh4_a_amode(ibin)    > 0 .or. lgtmpa) i_gas2bin_uptk_flag(inh3_g,  ibin) = 1
       if (lptr2_soa_a_amode(ibin,1) > 0 .or. lgtmpa) i_gas2bin_uptk_flag(ilim2_g, ibin) = 1
    end do
    if ( masterproc ) then
        write(iulog,'(a,10i5)') 'mosaic_cam_init m/i_gas2bin_uptk_flag', m_gas2bin_uptk_flag
        write(iulog,'(a,10i5)') 'h2so4 ', i_gas2bin_uptk_flag(ih2so4_g,1:ntot_amode)
        write(iulog,'(a,10i5)') 'hno3  ', i_gas2bin_uptk_flag(ihno3_g, 1:ntot_amode)
        write(iulog,'(a,10i5)') 'hcl   ', i_gas2bin_uptk_flag(ihcl_g,  1:ntot_amode)
        write(iulog,'(a,10i5)') 'nh3   ', i_gas2bin_uptk_flag(inh3_g,  1:ntot_amode)
        write(iulog,'(a,10i5)') 'lim2  ', i_gas2bin_uptk_flag(ilim2_g, 1:ntot_amode)
        do iv = 1, ngas_aerchtot
        write(iulog,'(i3.3,3x,10i5)') mod(iv,1000), i_gas2bin_uptk_flag(iv,1:ntot_amode)
        end do
    end if

  end subroutine mosaic_cam_init

#else
  subroutine mosaic_cam_init
  use cam_abortutils,  only:  endrun
  call endrun( '*** error -- mosaic_cam_init should not have been called' )
  end subroutine mosaic_cam_init

#endif

end module module_mosaic_cam_init
