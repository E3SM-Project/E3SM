
module mo_usrrxt

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use ppgrid,       only : pver, pcols
  use abortutils,   only : endrun
#ifdef MODAL_AERO
  use modal_aero_data,only : ntot_amode, nspec_amode, modename_amode
#endif

  implicit none

  private
  public :: usrrxt, usrrxt_inti, usrrxt_hrates

  save

  integer :: usr_O_O2_ndx
  integer :: usr_HO2_HO2_ndx
  integer :: usr_N2O5_M_ndx
  integer :: usr_HNO3_OH_ndx
  integer :: usr_HO2NO2_M_ndx
  integer :: usr_N2O5_aer_ndx
  integer :: usr_NO3_aer_ndx
  integer :: usr_NO2_aer_ndx
  integer :: usr_CO_OH_a_ndx
  integer :: usr_CO_OH_b_ndx
  integer :: usr_PAN_M_ndx
  integer :: usr_CH3COCH3_OH_ndx
  integer :: usr_MCO3_NO2_ndx
  integer :: usr_MPAN_M_ndx
  integer :: usr_XOOH_OH_ndx
  integer :: usr_SO2_OH_ndx
  integer :: usr_DMS_OH_ndx
  integer :: usr_HO2_aer_ndx
  
  integer :: tag_NO2_NO3_ndx
  integer :: tag_NO2_OH_ndx
  integer :: tag_NO2_HO2_ndx
  integer :: tag_C2H4_OH_ndx
  integer :: tag_C3H6_OH_ndx
  integer :: tag_CH3CO3_NO2_ndx
  
  integer :: usr_OA_O2_NDX
  integer :: usr_XNO2NO3_M_ndx
  integer :: usr_NO2XNO3_M_ndx
  integer :: usr_XHNO3_OH_ndx
  integer :: usr_XHO2NO2_M_ndx
  integer :: usr_XNO2NO3_aer_ndx
  integer :: usr_NO2XNO3_aer_ndx
  integer :: usr_XNO3_aer_ndx
  integer :: usr_XNO2_aer_ndx
  integer :: usr_XPAN_M_ndx
  integer :: usr_XMPAN_M_ndx
  integer :: usr_MCO3_XNO2_ndx
  
  integer :: usr_C2O3_NO2_ndx
  integer :: usr_C2H4_OH_ndx
  integer :: usr_XO2N_HO2_ndx
  integer :: usr_C2O3_XNO2_ndx
  
  integer :: tag_XO2N_NO_ndx
  integer :: tag_XO2_HO2_ndx
  integer :: tag_XO2_NO_ndx
  
  integer :: usr_O_O_ndx
  integer :: usr_CL2O2_M_ndx
  integer :: usr_SO3_H2O_ndx
  integer :: tag_CLO_CLO_ndx
  
  integer :: ion1_ndx, ion2_ndx, ion3_ndx, ion11_ndx
  integer :: elec1_ndx, elec2_ndx, elec3_ndx
  integer :: het1_ndx

  integer :: usr_oh_co_ndx, het_no2_h2o_ndx, usr_oh_dms_ndx, aq_so2_h2o2_ndx, aq_so2_o3_ndx

  integer :: h2o_ndx, so4_ndx, cb2_ndx, oc2_ndx, soa_ndx, nit_ndx

!lke++
  integer :: usr_COhc_OH_ndx
  integer :: usr_COme_OH_ndx
  integer :: usr_CO01_OH_ndx
  integer :: usr_CO02_OH_ndx
  integer :: usr_CO03_OH_ndx
  integer :: usr_CO04_OH_ndx
  integer :: usr_CO05_OH_ndx
  integer :: usr_CO06_OH_ndx
  integer :: usr_CO07_OH_ndx
  integer :: usr_CO08_OH_ndx
  integer :: usr_CO09_OH_ndx
  integer :: usr_CO10_OH_ndx
  integer :: usr_CO11_OH_ndx
  integer :: usr_CO12_OH_ndx
  integer :: usr_CO13_OH_ndx
  integer :: usr_CO14_OH_ndx
  integer :: usr_CO15_OH_ndx
  integer :: usr_CO16_OH_ndx
  integer :: usr_CO17_OH_ndx
  integer :: usr_CO18_OH_ndx
  integer :: usr_CO19_OH_ndx
  integer :: usr_CO20_OH_ndx
  integer :: usr_CO21_OH_ndx
  integer :: usr_CO22_OH_ndx
  integer :: usr_CO23_OH_ndx
  integer :: usr_CO24_OH_ndx
  integer :: usr_CO25_OH_ndx
  integer :: usr_CO26_OH_ndx
  integer :: usr_CO27_OH_ndx
  integer :: usr_CO28_OH_ndx
  integer :: usr_CO29_OH_ndx
  integer :: usr_CO30_OH_ndx
  integer :: usr_CO31_OH_ndx
  integer :: usr_CO32_OH_ndx
  integer :: usr_CO33_OH_ndx
  integer :: usr_CO34_OH_ndx
  integer :: usr_CO35_OH_ndx
  integer :: usr_CO36_OH_ndx
  integer :: usr_CO37_OH_ndx
  integer :: usr_CO38_OH_ndx
  integer :: usr_CO39_OH_ndx
  integer :: usr_CO40_OH_ndx
  integer :: usr_CO41_OH_ndx
  integer :: usr_CO42_OH_ndx
!lke--

  logical :: has_aerosols

  real(r8), parameter :: t0     = 300._r8                ! K
  real(r8), parameter :: trlim2 = 17._r8/3._r8           ! K
  real(r8), parameter :: trlim3 = 15._r8/3._r8           ! K

  logical :: has_ion_rxts

#ifdef MODAL_AERO
  integer :: dgnumwet_idx = -1
  integer :: aitken_idx = -1
  integer, dimension(ntot_amode) :: num_idx = -1
  integer :: index_tot_mass(ntot_amode,10) = -1
  integer :: index_chm_mass(ntot_amode,10) = -1
#endif

contains

  subroutine usrrxt_inti
    !-----------------------------------------------------------------
    !        ... intialize the user reaction constants module
    !-----------------------------------------------------------------

    use mo_chem_utls,   only : get_rxt_ndx, get_spc_ndx
    use spmd_utils,     only : masterproc
    use physics_buffer, only : pbuf_get_index

    implicit none
#ifdef MODAL_AERO
    integer :: l
    character(len=6) :: test_name
    character(len=64) :: errmes
#endif
!
! full tropospheric chemistry
!
    usr_O_O2_ndx         = get_rxt_ndx( 'usr_O_O2' )
    usr_HO2_HO2_ndx      = get_rxt_ndx( 'usr_HO2_HO2' )
    usr_N2O5_M_ndx       = get_rxt_ndx( 'usr_N2O5_M' )
    usr_HNO3_OH_ndx      = get_rxt_ndx( 'usr_HNO3_OH' )
    usr_HO2NO2_M_ndx     = get_rxt_ndx( 'usr_HO2NO2_M' )
    usr_N2O5_aer_ndx     = get_rxt_ndx( 'usr_N2O5_aer' )
    usr_NO3_aer_ndx      = get_rxt_ndx( 'usr_NO3_aer' )
    usr_NO2_aer_ndx      = get_rxt_ndx( 'usr_NO2_aer' )
    usr_CO_OH_a_ndx      = get_rxt_ndx( 'usr_CO_OH_a' )
    usr_CO_OH_b_ndx      = get_rxt_ndx( 'usr_CO_OH_b' )
    usr_PAN_M_ndx        = get_rxt_ndx( 'usr_PAN_M' )
    usr_CH3COCH3_OH_ndx  = get_rxt_ndx( 'usr_CH3COCH3_OH' )
    usr_MCO3_NO2_ndx     = get_rxt_ndx( 'usr_MCO3_NO2' )
    usr_MPAN_M_ndx       = get_rxt_ndx( 'usr_MPAN_M' )
    usr_XOOH_OH_ndx      = get_rxt_ndx( 'usr_XOOH_OH' )
    usr_SO2_OH_ndx       = get_rxt_ndx( 'usr_SO2_OH' )
    usr_DMS_OH_ndx       = get_rxt_ndx( 'usr_DMS_OH' )
    usr_HO2_aer_ndx      = get_rxt_ndx( 'usr_HO2_aer' )
 !
    tag_NO2_NO3_ndx      = get_rxt_ndx( 'tag_NO2_NO3' )
    tag_NO2_OH_ndx       = get_rxt_ndx( 'tag_NO2_OH' )
    tag_NO2_HO2_ndx      = get_rxt_ndx( 'tag_NO2_HO2' )
    tag_C2H4_OH_ndx      = get_rxt_ndx( 'tag_C2H4_OH' )
    tag_C3H6_OH_ndx      = get_rxt_ndx( 'tag_C3H6_OH' )
    tag_CH3CO3_NO2_ndx   = get_rxt_ndx( 'tag_CH3CO3_NO2' )     
 !
 ! additional reactions for O3A/XNO
 !
    usr_OA_O2_ndx        = get_rxt_ndx( 'usr_OA_O2' )
    usr_XNO2NO3_M_ndx    = get_rxt_ndx( 'usr_XNO2NO3_M' )
    usr_NO2XNO3_M_ndx    = get_rxt_ndx( 'usr_NO2XNO3_M' )
    usr_XNO2NO3_aer_ndx  = get_rxt_ndx( 'usr_XNO2NO3_aer' )
    usr_NO2XNO3_aer_ndx  = get_rxt_ndx( 'usr_NO2XNO3_aer' )
    usr_XHNO3_OH_ndx     = get_rxt_ndx( 'usr_XHNO3_OH' )
    usr_XNO3_aer_ndx     = get_rxt_ndx( 'usr_XNO3_aer' )
    usr_XNO2_aer_ndx     = get_rxt_ndx( 'usr_XNO2_aer' )
    usr_MCO3_XNO2_ndx    = get_rxt_ndx( 'usr_MCO3_XNO2' )
    usr_XPAN_M_ndx       = get_rxt_ndx( 'usr_XPAN_M' )
    usr_XMPAN_M_ndx      = get_rxt_ndx( 'usr_XMPAN_M' )
    usr_XHO2NO2_M_ndx    = get_rxt_ndx( 'usr_XHO2NO2_M' )
!
! reduced hydrocarbon chemistry
!
    usr_C2O3_NO2_ndx     = get_rxt_ndx( 'usr_C2O3_NO2' )
    usr_C2H4_OH_ndx      = get_rxt_ndx( 'usr_C2H4_OH' )
    usr_XO2N_HO2_ndx     = get_rxt_ndx( 'usr_XO2N_HO2' )
    usr_C2O3_XNO2_ndx    = get_rxt_ndx( 'usr_C2O3_XNO2' )
!
    tag_XO2N_NO_ndx      = get_rxt_ndx( 'tag_XO2N_NO' )
    tag_XO2_HO2_ndx      = get_rxt_ndx( 'tag_XO2_HO2' )
    tag_XO2_NO_ndx       = get_rxt_ndx( 'tag_XO2_NO' )
!
! stratospheric chemistry
!
    usr_O_O_ndx          = get_rxt_ndx( 'usr_O_O' )
    usr_CL2O2_M_ndx      = get_rxt_ndx( 'usr_CL2O2_M' )
    usr_SO3_H2O_ndx      = get_rxt_ndx( 'usr_SO3_H2O' )
!    
    tag_CLO_CLO_ndx      = get_rxt_ndx( 'tag_CLO_CLO' )
!
! stratospheric aerosol chemistry
!
    het1_ndx             = get_rxt_ndx( 'het1' )
!
! ion chemistry
!   
    ion1_ndx  = get_rxt_ndx( 'ion_Op_O2' )
    ion2_ndx  = get_rxt_ndx( 'ion_Op_N2' )
    ion3_ndx  = get_rxt_ndx( 'ion_N2p_Oa' )
    ion11_ndx = get_rxt_ndx( 'ion_N2p_Ob' )

    elec1_ndx  = get_rxt_ndx( 'elec1' )
    elec2_ndx  = get_rxt_ndx( 'elec2' )
    elec3_ndx  = get_rxt_ndx( 'elec3' )

    has_ion_rxts = ion1_ndx>0 .and. ion2_ndx>0 .and. ion3_ndx>0 .and. elec1_ndx>0 &
                 .and. elec2_ndx>0 .and. elec3_ndx>0

    so4_ndx    = get_spc_ndx( 'SO4' )
    cb2_ndx    = get_spc_ndx( 'CB2' )
    oc2_ndx    = get_spc_ndx( 'OC2' )
    soa_ndx    = get_spc_ndx( 'SOA' )
    nit_ndx    = get_spc_ndx( 'NH4NO3' )
    h2o_ndx    = get_spc_ndx( 'H2O' )

    !
    ! llnl super fast
    !
    usr_oh_co_ndx  = get_rxt_ndx( 'usr_oh_co' )
    het_no2_h2o_ndx  = get_rxt_ndx( 'het_no2_h2o' )
    usr_oh_dms_ndx  = get_rxt_ndx( 'usr_oh_dms' )
    aq_so2_h2o2_ndx  = get_rxt_ndx( 'aq_so2_h2o2' )
    aq_so2_o3_ndx  = get_rxt_ndx( 'aq_so2_o3' )
    
!lke++
! CO tags
!
    usr_COhc_OH_ndx      = get_rxt_ndx( 'usr_COhc_OH' )
    usr_COme_OH_ndx      = get_rxt_ndx( 'usr_COme_OH' )
    usr_CO01_OH_ndx      = get_rxt_ndx( 'usr_CO01_OH' )
    usr_CO02_OH_ndx      = get_rxt_ndx( 'usr_CO02_OH' )
    usr_CO03_OH_ndx      = get_rxt_ndx( 'usr_CO03_OH' )
    usr_CO04_OH_ndx      = get_rxt_ndx( 'usr_CO04_OH' )
    usr_CO05_OH_ndx      = get_rxt_ndx( 'usr_CO05_OH' )
    usr_CO06_OH_ndx      = get_rxt_ndx( 'usr_CO06_OH' )
    usr_CO07_OH_ndx      = get_rxt_ndx( 'usr_CO07_OH' )
    usr_CO08_OH_ndx      = get_rxt_ndx( 'usr_CO08_OH' )
    usr_CO09_OH_ndx      = get_rxt_ndx( 'usr_CO09_OH' )
    usr_CO10_OH_ndx      = get_rxt_ndx( 'usr_CO10_OH' )
    usr_CO11_OH_ndx      = get_rxt_ndx( 'usr_CO11_OH' )
    usr_CO12_OH_ndx      = get_rxt_ndx( 'usr_CO12_OH' )
    usr_CO13_OH_ndx      = get_rxt_ndx( 'usr_CO13_OH' )
    usr_CO14_OH_ndx      = get_rxt_ndx( 'usr_CO14_OH' )
    usr_CO15_OH_ndx      = get_rxt_ndx( 'usr_CO15_OH' )
    usr_CO16_OH_ndx      = get_rxt_ndx( 'usr_CO16_OH' )
    usr_CO17_OH_ndx      = get_rxt_ndx( 'usr_CO17_OH' )
    usr_CO18_OH_ndx      = get_rxt_ndx( 'usr_CO18_OH' )
    usr_CO19_OH_ndx      = get_rxt_ndx( 'usr_CO19_OH' )
    usr_CO20_OH_ndx      = get_rxt_ndx( 'usr_CO20_OH' )
    usr_CO21_OH_ndx      = get_rxt_ndx( 'usr_CO21_OH' )
    usr_CO22_OH_ndx      = get_rxt_ndx( 'usr_CO22_OH' )
    usr_CO23_OH_ndx      = get_rxt_ndx( 'usr_CO23_OH' )
    usr_CO24_OH_ndx      = get_rxt_ndx( 'usr_CO24_OH' )
    usr_CO25_OH_ndx      = get_rxt_ndx( 'usr_CO25_OH' )
    usr_CO26_OH_ndx      = get_rxt_ndx( 'usr_CO26_OH' )
    usr_CO27_OH_ndx      = get_rxt_ndx( 'usr_CO27_OH' )
    usr_CO28_OH_ndx      = get_rxt_ndx( 'usr_CO28_OH' )
    usr_CO29_OH_ndx      = get_rxt_ndx( 'usr_CO29_OH' )
    usr_CO30_OH_ndx      = get_rxt_ndx( 'usr_CO30_OH' )
    usr_CO31_OH_ndx      = get_rxt_ndx( 'usr_CO31_OH' )
    usr_CO32_OH_ndx      = get_rxt_ndx( 'usr_CO32_OH' )
    usr_CO33_OH_ndx      = get_rxt_ndx( 'usr_CO33_OH' )
    usr_CO34_OH_ndx      = get_rxt_ndx( 'usr_CO34_OH' )
    usr_CO35_OH_ndx      = get_rxt_ndx( 'usr_CO35_OH' )
    usr_CO36_OH_ndx      = get_rxt_ndx( 'usr_CO36_OH' )
    usr_CO37_OH_ndx      = get_rxt_ndx( 'usr_CO37_OH' )
    usr_CO38_OH_ndx      = get_rxt_ndx( 'usr_CO38_OH' )
    usr_CO39_OH_ndx      = get_rxt_ndx( 'usr_CO39_OH' )
    usr_CO40_OH_ndx      = get_rxt_ndx( 'usr_CO40_OH' )
    usr_CO41_OH_ndx      = get_rxt_ndx( 'usr_CO41_OH' )
    usr_CO42_OH_ndx      = get_rxt_ndx( 'usr_CO42_OH' )
!lke--

#ifdef MODAL_AERO
    if( usr_NO2_aer_ndx > 0 .or. usr_NO3_aer_ndx > 0 .or. usr_N2O5_aer_ndx > 0 .or. usr_HO2_aer_ndx > 0 ) then
       do l=1,ntot_amode
          if ( trim(modename_amode(l)) == 'aitken' ) then
             aitken_idx = l
          end if
          test_name = ' '
          write(test_name,fmt='(a5,i1)') 'num_a',l
          num_idx(l) = get_spc_ndx( trim(test_name) )
          if (num_idx(l) < 0) then
             write(errmes,fmt='(a,i1)') 'usrrxt_inti: cannot find MAM num_idx ',l
             write(iulog,*) errmes
             call endrun(errmes)
          endif
       end do
       dgnumwet_idx = pbuf_get_index('DGNUMWET')
       if ( aitken_idx < 0 ) then
          errmes = 'usrrxt_inti: cannot find aitken_idx'
          call endrun(errmes)
       end if
!
! define indeces associated with the various names (defined in
! chemistry/modal_aero/modal_aero_initialize_data.F90)
!
#if ( defined MODAL_AERO_3MODE )
!
! accumulation mode #1
!
       index_tot_mass(1,1) = get_spc_ndx('so4_a1')
       index_tot_mass(1,2) = get_spc_ndx('pom_a1')
       index_tot_mass(1,3) = get_spc_ndx('soa_a1')
       index_tot_mass(1,4) = get_spc_ndx('bc_a1' )
       index_tot_mass(1,5) = get_spc_ndx('dst_a1')
       index_tot_mass(1,6) = get_spc_ndx('ncl_a1')
       index_chm_mass(1,1) = get_spc_ndx('so4_a1')
       index_chm_mass(1,2) = get_spc_ndx('soa_a1')
       index_chm_mass(1,3) = get_spc_ndx('bc_a1' )
!
! aitken mode
!
       index_tot_mass(2,1) = get_spc_ndx('so4_a2')
       index_tot_mass(2,2) = get_spc_ndx('soa_a2')
       index_tot_mass(2,3) = get_spc_ndx('ncl_a2')
       index_chm_mass(2,1) = get_spc_ndx('so4_a2')
       index_chm_mass(2,2) = get_spc_ndx('soa_a2')
!
! coarse mode
!
       index_tot_mass(3,1) = get_spc_ndx('dst_a3')
       index_tot_mass(3,2) = get_spc_ndx('ncl_a3')
       index_tot_mass(3,3) = get_spc_ndx('so4_a3')
       index_chm_mass(3,1) = get_spc_ndx('so4_a3')
!
#endif
#if ( defined MODAL_AERO_7MODE )
!
! accumulation mode #1
!
       index_tot_mass(1,1) = get_spc_ndx('so4_a1')
       index_tot_mass(1,2) = get_spc_ndx('nh4_a1')
       index_tot_mass(1,3) = get_spc_ndx('pom_a1')
       index_tot_mass(1,4) = get_spc_ndx('soa_a1')
       index_tot_mass(1,5) = get_spc_ndx('bc_a1' )
       index_tot_mass(1,6) = get_spc_ndx('ncl_a1')
       index_chm_mass(1,1) = get_spc_ndx('so4_a1')
       index_chm_mass(1,2) = get_spc_ndx('nh4_a1')
       index_chm_mass(1,3) = get_spc_ndx('soa_a1')
       index_chm_mass(1,4) = get_spc_ndx('bc_a1' )
!
! aitken mode
!
       index_tot_mass(2,1) = get_spc_ndx('so4_a2')
       index_tot_mass(2,2) = get_spc_ndx('nh4_a2')
       index_tot_mass(2,3) = get_spc_ndx('soa_a2')
       index_tot_mass(2,4) = get_spc_ndx('ncl_a2')
       index_chm_mass(2,1) = get_spc_ndx('so4_a2')
       index_chm_mass(2,2) = get_spc_ndx('nh4_a2')
       index_chm_mass(2,3) = get_spc_ndx('soa_a2')
!
! primary carbon mode not added 
!
! fine sea salt 
!
       index_tot_mass(4,1) = get_spc_ndx('so4_a4')
       index_tot_mass(4,2) = get_spc_ndx('nh4_a4')
       index_tot_mass(4,3) = get_spc_ndx('ncl_a4')
       index_chm_mass(4,1) = get_spc_ndx('so4_a4')
       index_chm_mass(4,2) = get_spc_ndx('nh4_a4')
!
! fine soil dust 
!
       index_tot_mass(5,1) = get_spc_ndx('so4_a5')
       index_tot_mass(5,2) = get_spc_ndx('nh4_a5')
       index_tot_mass(5,3) = get_spc_ndx('dst_a5')
       index_chm_mass(5,1) = get_spc_ndx('so4_a5')
       index_chm_mass(5,2) = get_spc_ndx('nh4_a5')
!
! coarse sea salt 
!
       index_tot_mass(6,1) = get_spc_ndx('so4_a6')
       index_tot_mass(6,2) = get_spc_ndx('nh4_a6')
       index_tot_mass(6,3) = get_spc_ndx('ncl_a6')
       index_chm_mass(6,1) = get_spc_ndx('so4_a6')
       index_chm_mass(6,2) = get_spc_ndx('nh4_a6')
!
! coarse soil dust 
!
       index_tot_mass(7,1) = get_spc_ndx('so4_a7')
       index_tot_mass(7,2) = get_spc_ndx('nh4_a7')
       index_tot_mass(7,3) = get_spc_ndx('dst_a7')
       index_chm_mass(7,1) = get_spc_ndx('so4_a7')
       index_chm_mass(7,2) = get_spc_ndx('nh4_a7')
!
#endif


    endif
#endif

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'usrrxt_inti: diagnostics '
       write(iulog,'(10i5)') usr_O_O2_ndx,usr_HO2_HO2_ndx,tag_NO2_NO3_ndx,usr_N2O5_M_ndx,tag_NO2_OH_ndx,usr_HNO3_OH_ndx &
                            ,tag_NO2_HO2_ndx,usr_HO2NO2_M_ndx,usr_N2O5_aer_ndx,usr_NO3_aer_ndx,usr_NO2_aer_ndx &
                            ,usr_CO_OH_b_ndx,tag_C2H4_OH_ndx,tag_C3H6_OH_ndx,tag_CH3CO3_NO2_ndx,usr_PAN_M_ndx,usr_CH3COCH3_OH_ndx &
                            ,usr_MCO3_NO2_ndx,usr_MPAN_M_ndx,usr_XOOH_OH_ndx,usr_SO2_OH_ndx,usr_DMS_OH_ndx,usr_HO2_aer_ndx
    end if

  end subroutine usrrxt_inti

  subroutine usrrxt( rxt, temp, tempi, tempe, invariants, h2ovmr,  ps, &
                     pmid, m, sulfate, mmr, relhum, strato_sad, &
                     ltrop, ncol, sad_total, cwat, mbar, pbuf )

!-----------------------------------------------------------------
!        ... set the user specified reaction rates
!-----------------------------------------------------------------
    
    use mo_constants,  only : pi, avo => avogadro, boltz_cgs, rgas
    use chem_mods,     only : nfs, rxntot, gas_pcnst, inv_m_ndx=>indexm
    use mo_chem_utls,  only : get_rxt_ndx, get_spc_ndx
    use mo_setinv,     only : inv_o2_ndx=>o2_ndx, inv_h2o_ndx=>h2o_ndx
    use physics_buffer,only : physics_buffer_desc
    use carma_flags_mod,  only : carma_do_hetchem

    implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    integer, intent(in)     :: ltrop(pcols)               ! tropopause vertical index
    real(r8), intent(in)    :: temp(pcols,pver)           ! temperature (K); neutral temperature
    real(r8), intent(in)    :: tempi(pcols,pver)          ! ionic temperature (K); only used if ion chemistry
    real(r8), intent(in)    :: tempe(pcols,pver)          ! electronic temperature (K); only used if ion chemistry
    real(r8), intent(in)    :: m(ncol,pver)               ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(ncol,pver)         ! sulfate aerosol (mol/mol)
    real(r8), intent(in)    :: strato_sad(pcols,pver)     ! stratospheric aerosol sad (1/cm)
    real(r8), intent(in)    :: h2ovmr(ncol,pver)          ! water vapor (mol/mol)
    real(r8), intent(in)    :: relhum(ncol,pver)          ! relative humidity
    real(r8), intent(in)    :: pmid(pcols,pver)           ! midpoint pressure (Pa)
    real(r8), intent(in)    :: ps(pcols)                  ! surface pressure (Pa)
    real(r8), intent(in)    :: invariants(ncol,pver,nfs)  ! invariants density (/cm^3)
    real(r8), intent(in)    :: mmr(pcols,pver,gas_pcnst)  ! species concentrations (kg/kg)
    real(r8), intent(in)    :: cwat(ncol,pver) !PJC Condensed Water (liquid+ice) (kg/kg)
    real(r8), intent(in)    :: mbar(ncol,pver) !PJC Molar mass of air (g/mol)
    real(r8), intent(inout) :: rxt(ncol,pver,rxntot)      ! gas phase rates
    real(r8), intent(inout) :: sad_total(pcols,pver)      ! total surface area density (cm2/cm3)
    type(physics_buffer_desc), pointer :: pbuf(:)
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------
    
    real(r8), parameter :: dg = 0.1_r8            ! mole diffusion =0.1 cm2/s (Dentener, 1993)

!-----------------------------------------------------------------
! 	... reaction probabilities for heterogeneous reactions
!-----------------------------------------------------------------
    real(r8), parameter :: gamma_n2o5 = 0.10_r8         ! from Jacob, Atm Env, 34, 2131, 2000
    real(r8), parameter :: gamma_ho2  = 0.20_r8         ! 
    real(r8), parameter :: gamma_no2  = 0.0001_r8       ! 
    real(r8), parameter :: gamma_no3  = 0.001_r8        ! 

    integer  ::  i, k
    real(r8) ::  tp(ncol)                       ! 300/t
    real(r8) ::  tinv(ncol)                     ! 1/t
    real(r8) ::  ko(ncol)   
    real(r8) ::  term1(ncol)
    real(r8) ::  term2(ncol)
    real(r8) ::  kinf(ncol)   
    real(r8) ::  fc(ncol)   
    real(r8) ::  xr(ncol)   
    real(r8) ::  sur(ncol)   
    real(r8) ::  sqrt_t(ncol)                   ! sqrt( temp )
    real(r8) ::  exp_fac(ncol)                  ! vector exponential
    real(r8) ::  lwc(ncol)   
    real(r8) ::  ko_m(ncol)   
    real(r8) ::  k0(ncol)   
    real(r8) ::  kinf_m(ncol)   
    real(r8) ::  o2(ncol)
    real(r8) ::  c_n2o5, c_ho2, c_no2, c_no3
    real(r8) ::  amas
    !-----------------------------------------------------------------
    !	... density of sulfate aerosol
    !-----------------------------------------------------------------
    real(r8), parameter :: gam1 = 0.04_r8                 ! N2O5+SUL ->2HNO3
    real(r8), parameter :: wso4 = 98._r8
    real(r8), parameter :: den  = 1.15_r8                 ! each molecule of SO4(aer) density g/cm3
    !-------------------------------------------------
    ! 	... volume of sulfate particles
    !           assuming mean rm 
    !           continient 0.05um  0.07um  0.09um
    !           ocean      0.09um  0.25um  0.37um
    !                      0.16um                  Blake JGR,7195, 1995
    !-------------------------------------------------
    real(r8), parameter :: rm1  = 0.16_r8*1.e-4_r8             ! mean radii in cm
    real(r8), parameter :: fare = 4._r8*pi*rm1*rm1             ! each mean particle(r=0.1u) area   cm2/cm3

    !-----------------------------------------------------------------------
    !        ... Aqueous phase sulfur quantities for SO2 + H2O2 and SO2 + O3
    !-----------------------------------------------------------------------
    real(r8), parameter  :: HENRY298_H2O2 =  7.45e+04_r8
    real(r8), parameter  :: H298_H2O2     = -1.45e+04_r8
    real(r8), parameter  :: HENRY298_SO2  =  1.23e+00_r8
    real(r8), parameter  :: H298_SO2      = -6.25e+03_r8
    real(r8), parameter  :: K298_SO2_HSO3 =  1.3e-02_r8
    real(r8), parameter  :: H298_SO2_HSO3 = -4.16e+03_r8
    real(r8), parameter  :: R_CONC        =  82.05e+00_r8 / avo
    real(r8), parameter  :: R_CAL         =  rgas * 0.239006e+00_r8
    real(r8), parameter  :: K_AQ          =  7.57e+07_r8
    real(r8), parameter  :: ER_AQ         =  4.43e+03_r8

    real(r8), parameter  :: HENRY298_O3   =  1.13e-02_r8
    real(r8), parameter  :: H298_O3       = -5.04e+03_r8
    real(r8), parameter  :: K298_HSO3_SO3 =  6.6e-08_r8
    real(r8), parameter  :: H298_HSO3_SO3 = -2.23e+03_r8
    real(r8), parameter  :: K0_AQ         =  2.4e+04_r8
    real(r8), parameter  :: ER0_AQ        =  0.0e+00_r8
    real(r8), parameter  :: K1_AQ         =  3.7e+05_r8
    real(r8), parameter  :: ER1_AQ        =  5.53e+03_r8
    real(r8), parameter  :: K2_AQ         =  1.5e+09_r8
    real(r8), parameter  :: ER2_AQ        =  5.28e+03_r8

    real(r8), parameter  :: pH            =  4.5e+00_r8
    
    real(r8), pointer :: sfc(:), dm_aer(:)

#ifdef MODAL_AERO
    real(r8), target  :: sfc_array(pcols,pver,ntot_amode), dm_array(pcols,pver,ntot_amode)
#else
    real(r8), target  :: sfc_array(pcols,pver,4), dm_array(pcols,pver,4)
#endif
    
    if( usr_NO2_aer_ndx > 0 .or. usr_NO3_aer_ndx > 0 .or. usr_N2O5_aer_ndx > 0 .or. usr_HO2_aer_ndx > 0 ) then
       if( carma_do_hetchem ) then
          sad_total(:ncol,:pver)=strato_sad(:ncol,:pver)
       else
#if (defined MODAL_AERO)
         call surfarea( mmr, pmid, temp, pbuf, ncol, sfc_array, dm_array, sad_total )
#else
         call surfarea( mmr, relhum, pmid, temp, strato_sad, sulfate, m, ltrop, ncol, sfc_array, dm_array, sad_total )
#endif
       endif
    endif

    level_loop : do k = 1,pver
       tinv(:)           = 1._r8 / temp(:ncol,k)
       tp(:)             = 300._r8 * tinv(:)
       sqrt_t(:)         = sqrt( temp(:ncol,k) )

!-----------------------------------------------------------------
!	... o + o2 + m --> o3 + m
!-----------------------------------------------------------------
       if( usr_O_O2_ndx > 0 ) then
          rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8
       end if
       if( usr_OA_O2_ndx > 0 ) then
          rxt(:,k,usr_OA_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8
       end if

!-----------------------------------------------------------------
!	... o + o + m -> o2 + m
!-----------------------------------------------------------------
       if ( usr_O_O_ndx > 0 ) then
          rxt(:,k,usr_O_O_ndx) = 2.76e-34_r8 * exp( 720.0_r8*tinv(:) )
       end if
         
!-----------------------------------------------------------------
! 	... cl2o2 + m -> 2*clo + m
!-----------------------------------------------------------------
       if ( usr_CL2O2_M_ndx > 0 ) then
          if ( tag_CLO_CLO_ndx > 0 ) then
             ko(:)            = 9.3e-28_r8 * exp( 8835.0_r8* tinv(:) )
             rxt(:,k,usr_CL2O2_M_ndx) = rxt(:,k,tag_CLO_CLO_ndx)/ko(:)         
          else
             rxt(:,k,usr_CL2O2_M_ndx) = 0._r8
          end if
       end if
       
!-----------------------------------------------------------------
!       ... so3 + 2*h2o --> h2so4 + h2o
!       Note: this reaction proceeds by the 2 intermediate steps below
!           so3 + h2o --> adduct
!           adduct + h2o --> h2so4 + h2o  
!               (Lovejoy et al., JCP, pp. 19911-19916, 1996)
!	The first order rate constant used here is recommended by JPL 2011.
!	This rate involves the water vapor number density.
!-----------------------------------------------------------------

       if ( usr_SO3_H2O_ndx > 0 ) then
          call comp_exp( exp_fac, 6540.0_r8*tinv(:), ncol )
          if( h2o_ndx > 0 ) then
             fc(:) = 8.5e-21_r8 * m(:,k) * h2ovmr(:,k) * exp_fac(:)
          else
             fc(:) = 8.5e-21_r8 * invariants(:,k,inv_h2o_ndx) * exp_fac(:)
          end if
          rxt(:,k,usr_SO3_H2O_ndx) = 1.0e-20_r8 * fc(:)
       end if
         
!-----------------------------------------------------------------
!	... n2o5 + m --> no2 + no3 + m
!-----------------------------------------------------------------
       if( usr_N2O5_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -11000.0_r8*tinv, ncol )
             rxt(:,k,usr_N2O5_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) * 3.703704e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_N2O5_M_ndx) = 0._r8
          end if
       end if
       if( usr_XNO2NO3_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -11000._r8*tinv, ncol )
             rxt(:,k,usr_XNO2NO3_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) * 3.703704e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XNO2NO3_M_ndx) = 0._r8
          end if
       end if
       if( usr_NO2XNO3_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -11000._r8*tinv, ncol )
             rxt(:,k,usr_NO2XNO3_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) * 3.703704e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_NO2XNO3_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!	set rates for:
! 	... hno3 + oh --> no3 + h2o
!           ho2no2 + m --> ho2 + no2 + m
!-----------------------------------------------------------------
       if( usr_HNO3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 1335._r8*tinv, ncol )
          ko(:) = m(:,k) * 6.5e-34_r8 * exp_fac(:)
          call comp_exp( exp_fac, 2199._r8*tinv, ncol )
          ko(:) = ko(:) / (1._r8 + ko(:)/(2.7e-17_r8*exp_fac(:)))
          call comp_exp( exp_fac, 460._r8*tinv, ncol )
          rxt(:,k,usr_HNO3_OH_ndx) = ko(:) + 2.4e-14_r8*exp_fac(:)
       end if
       if( usr_XHNO3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 1335._r8*tinv, ncol )
          ko(:) = m(:,k) * 6.5e-34_r8 * exp_fac(:)
          call comp_exp( exp_fac, 2199._r8*tinv, ncol )
          ko(:) = ko(:) / (1._r8 + ko(:)/(2.7e-17_r8*exp_fac(:)))
          call comp_exp( exp_fac, 460._r8*tinv, ncol )
          rxt(:,k,usr_XHNO3_OH_ndx) = ko(:) + 2.4e-14_r8*exp_fac(:)
       end if
       if( usr_HO2NO2_M_ndx > 0 ) then
          if( tag_NO2_HO2_ndx > 0 ) then
             call comp_exp( exp_fac, -10900._r8*tinv, ncol )
             rxt(:,k,usr_HO2NO2_M_ndx) = rxt(:,k,tag_NO2_HO2_ndx) * exp_fac(:) / 2.1e-27_r8
          else
             rxt(:,k,usr_HO2NO2_M_ndx) = 0._r8
          end if
       end if
       if( usr_XHO2NO2_M_ndx > 0 ) then
          if( tag_NO2_HO2_ndx > 0 ) then
             call comp_exp( exp_fac, -10900._r8*tinv, ncol )
             rxt(:,k,usr_XHO2NO2_M_ndx) = rxt(:,k,tag_NO2_HO2_ndx) * exp_fac(:) / 2.1e-27_r8
          else
             rxt(:,k,usr_XHO2NO2_M_ndx) = 0._r8
          end if
       end if
!-----------------------------------------------------------------
!           co + oh --> co2 + ho2     CAM-Chem
!-----------------------------------------------------------------
       if( usr_CO_OH_a_ndx > 0 ) then
          rxt(:,k,usr_CO_OH_a_ndx) = 1.5e-13_r8 * &
               (1._r8 + 6.e-7_r8*boltz_cgs*m(:,k)*temp(:ncol,k))
       end if
!-----------------------------------------------------------------
! 	... co + oh --> co2 + h (second branch JPL06; pg2.2; 2.10) WACCM
!-----------------------------------------------------------------
       if( usr_CO_OH_b_ndx > 0 ) then
         kinf(:)  = 2.1e+09_r8 * (temp(:ncol,k)/ t0)**(6.1_r8)
         ko  (:)  = 1.5e-13_r8 * (temp(:ncol,k)/ t0)**(0.6_r8)

         term1(:) = ko(:) / ( (kinf(:) / m(:,k)) )
         term2(:) = ko(:) / (1._r8 + term1(:))

         term1(:) = log10( term1(:) )
         term1(:) = 1.0_r8 / (1.0_r8 + term1(:)*term1(:))

         rxt(:ncol,k,usr_CO_OH_b_ndx) = term2(:) * (0.6_r8)**term1(:)
       end if

!-----------------------------------------------------------------
!	... ho2 + ho2 --> h2o2
!	note: this rate involves the water vapor number density
!-----------------------------------------------------------------
       if( usr_HO2_HO2_ndx > 0 ) then

          call comp_exp( exp_fac, 430._r8*tinv, ncol )
          ko(:)   = 3.5e-13_r8 * exp_fac(:)
          call comp_exp( exp_fac, 1000._r8*tinv, ncol )
          kinf(:) = 1.7e-33_r8 * m(:,k) * exp_fac(:)
          call comp_exp( exp_fac, 2200._r8*tinv, ncol )

          if( h2o_ndx > 0 ) then
             fc(:) = 1._r8 + 1.4e-21_r8 * m(:,k) * h2ovmr(:,k) * exp_fac(:)
          else
             fc(:) = 1._r8 + 1.4e-21_r8 * invariants(:,k,inv_h2o_ndx) * exp_fac(:)
          end if
          rxt(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)

       end if

!-----------------------------------------------------------------
!    	... mco3 + no2 -> mpan
!-----------------------------------------------------------------
       if( usr_MCO3_NO2_ndx > 0 ) then
          rxt(:,k,usr_MCO3_NO2_ndx) = 1.1e-11_r8 * tp(:) / m(:,k)
       end if
       if( usr_MCO3_XNO2_ndx > 0 ) then
          rxt(:,k,usr_MCO3_XNO2_ndx) = 1.1e-11_r8 * tp(:) / m(:,k)
       end if

!-----------------------------------------------------------------
!	... pan + m --> ch3co3 + no2 + m
!-----------------------------------------------------------------
       call comp_exp( exp_fac, -14000._r8*tinv, ncol )
       if( usr_PAN_M_ndx > 0 ) then
          if( tag_CH3CO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_PAN_M_ndx) = rxt(:,k,tag_CH3CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_PAN_M_ndx) = 0._r8
          end if
       end if
       if( usr_XPAN_M_ndx > 0 ) then
          if( tag_CH3CO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_XPAN_M_ndx) = rxt(:,k,tag_CH3CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XPAN_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!	... mpan + m --> mco3 + no2 + m
!-----------------------------------------------------------------
       if( usr_MPAN_M_ndx > 0 ) then
          if( usr_MCO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_MPAN_M_ndx) = rxt(:,k,usr_MCO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_MPAN_M_ndx) = 0._r8
          end if
       end if
       if( usr_XMPAN_M_ndx > 0 ) then
          if( usr_MCO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_XMPAN_M_ndx) = rxt(:,k,usr_MCO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XMPAN_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!       ... xooh + oh -> h2o + oh
!-----------------------------------------------------------------
       if( usr_XOOH_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 253._r8*tinv, ncol )
          rxt(:,k,usr_XOOH_OH_ndx) = temp(:ncol,k)**2._r8 * 7.69e-17_r8 * exp_fac(:)
       end if

!-----------------------------------------------------------------
!       ... ch3coch3 + oh -> ro2 + h2o
!-----------------------------------------------------------------
       if( usr_CH3COCH3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, -2000._r8*tinv, ncol )
          rxt(:,k,usr_CH3COCH3_OH_ndx) = 3.82e-11_r8 * exp_fac(:) + 1.33e-13_r8
       end if

!-----------------------------------------------------------------
!       ... DMS + OH  --> .5 * SO2
!-----------------------------------------------------------------
       if( usr_DMS_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 7460._r8*tinv, ncol )
          ko(:) = 1._r8 + 5.5e-31_r8 * exp_fac * m(:,k) * 0.21_r8
          call comp_exp( exp_fac, 7810._r8*tinv, ncol )
          rxt(:,k,usr_DMS_OH_ndx) = 1.7e-42_r8 * exp_fac * m(:,k) * 0.21_r8 / ko(:)
       end if

!-----------------------------------------------------------------
!       ... SO2 + OH  --> SO4  (REFERENCE?? - not Liao)
!-----------------------------------------------------------------
       if( usr_SO2_OH_ndx > 0 ) then
          fc(:) = 3.0e-31_r8 *(300._r8*tinv(:))**3.3_r8
          ko(:) = fc(:)*m(:,k)/(1._r8 + fc(:)*m(:,k)/1.5e-12_r8) 
          rxt(:,k,usr_SO2_OH_ndx) = ko(:)*.6_r8**(1._r8 + (log10(fc(:)*m(:,k)/1.5e-12_r8))**2._r8)**(-1._r8)
       end if
!
! reduced hydrocarbon scheme
!
       if ( usr_C2O3_NO2_ndx > 0 ) then
          ko(:)   = 2.6e-28_r8 * m(:,k)
          kinf(:) = 1.2e-11_r8
          rxt(:,k,usr_C2O3_NO2_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log10(ko/kinf))**2))
       end if
       if ( usr_C2O3_XNO2_ndx > 0 ) then
          ko(:)   = 2.6e-28_r8 * m(:,k)
          kinf(:) = 1.2e-11_r8
          rxt(:,k,usr_C2O3_XNO2_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log10(ko/kinf))**2))
       end if
       if ( usr_C2H4_OH_ndx > 0 ) then
          ko(:)   = 1.0e-28_r8 * m(:,k)
          kinf(:) = 8.8e-12_r8
          rxt(:,k,usr_C2H4_OH_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log(ko/kinf))**2))
       end if
       if ( usr_XO2N_HO2_ndx > 0 ) then
          rxt(:,k,usr_XO2N_HO2_ndx) = rxt(:,k,tag_XO2N_NO_ndx)*rxt(:,k,tag_XO2_HO2_ndx)/(rxt(:,k,tag_XO2_NO_ndx)+1.e-36_r8)
       end if
       
!
! hydrolysis reactions on wetted aerosols
!      
       if( usr_NO2_aer_ndx > 0 .or. usr_NO3_aer_ndx > 0 .or. usr_N2O5_aer_ndx > 0 .or. usr_HO2_aer_ndx > 0 ) then

          long_loop : do i = 1,ncol

             sfc    => sfc_array(i,k,:)
             dm_aer => dm_array(i,k,:)

             c_n2o5 = 1.40e3_r8 * sqrt_t(i)         ! mean molecular speed of n2o5
             c_no3  = 1.85e3_r8 * sqrt_t(i)         ! mean molecular speed of no3
             c_no2  = 2.15e3_r8 * sqrt_t(i)         ! mean molecular speed of no2
             c_ho2  = 2.53e3_r8 * sqrt_t(i)         ! mean molecular speed of ho2

             !-------------------------------------------------------------------------
             !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
             !    rxt = sfc / ( (rad_aer/Dg_gas) + (4/(c_gas*gamma_gas)))
             !-------------------------------------------------------------------------
             !-------------------------------------------------------------------------
             ! 	... n2o5 -> 2 hno3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_N2O5_aer_ndx > 0 ) then
                rxt(i,k,usr_N2O5_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             if( usr_XNO2NO3_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO2NO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             if( usr_NO2XNO3_aer_ndx > 0 ) then
                rxt(i,k,usr_NO2XNO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             !-------------------------------------------------------------------------
             ! 	... no3 -> hno3  (on sulfate, nh4no3, oc, soa)
             !-------------------------------------------------------------------------
             if( usr_NO3_aer_ndx > 0 ) then
                rxt(i,k,usr_NO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3, gamma_no3 )
             end if
             if( usr_XNO3_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3, gamma_no3 )
             end if
             !-------------------------------------------------------------------------
             ! 	... no2 -> 0.5 * (ho+no+hno3)  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NO2_aer_ndx > 0 ) then
                rxt(i,k,usr_NO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no2, gamma_no2 )
             end if
             if( usr_XNO2_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no2, gamma_no2 )
            end if
             !-------------------------------------------------------------------------
             ! 	... ho2 -> 0.5 * h2o2  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_HO2_aer_ndx > 0 ) then
                rxt(i,k,usr_HO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_ho2, gamma_ho2 )
             end if
          end do long_loop
       end if

       ! LLNL super fast chem reaction rates

       !-----------------------------------------------------------------------
       !     ... CO + OH --> CO2 + HO2
       !-----------------------------------------------------------------------
       if ( usr_oh_co_ndx > 0 ) then
          ko(:)     = 5.9e-33_r8 * tp(:)**1.4_r8
          kinf(:)   = 1.1e-12_r8 * (temp(:ncol,k) / 300._r8)**1.3_r8
          ko_m(:)   = ko(:) * m(:,k)
          k0(:)     = 1.5e-13_r8 * (temp(:ncol,k) / 300._r8)**0.6_r8
          kinf_m(:) = (2.1e+09_r8 * (temp(:ncol,k) / 300._r8)**6.1_r8) / m(:,k)
          rxt(:,k,usr_oh_co_ndx) = (ko_m(:)/(1._r8+(ko_m(:)/kinf(:)))) * &
               0.6_r8**(1._r8/(1._r8+(log10(ko_m(:)/kinf(:)))**2._r8)) + &
               (k0(:)/(1._r8+(k0(:)/kinf_m(:)))) * &
               0.6_r8**(1._r8/(1._r8+(log10(k0(:)/kinf_m(:)))**2._r8))
       endif
       !-----------------------------------------------------------------------
       !     ... NO2 + H2O --> 0.5 HONO + 0.5 HNO3
       !-----------------------------------------------------------------------
       if ( het_no2_h2o_ndx > 0 ) then
          rxt(:,k,het_no2_h2o_ndx) = 4.0e-24_r8
       endif
       !-----------------------------------------------------------------------
       !     ... DMS + OH --> 0.75 SO2 + 0.25 MSA
       !-----------------------------------------------------------------------
       if ( usr_oh_dms_ndx > 0 ) then
          o2(:ncol) = invariants(:ncol,k,inv_o2_ndx)
          rxt(:,k,usr_oh_dms_ndx) = 2.000e-10_r8 * exp(5820.0_r8 * tinv(:)) / &
               ((2.000e29_r8 / o2(:)) + exp(6280.0_r8 * tinv(:)))
       endif
       if ( aq_so2_h2o2_ndx > 0 .or. aq_so2_o3_ndx > 0 ) then
          lwc(:) = cwat(:ncol,k) * invariants(:ncol,k,inv_m_ndx) * mbar(:ncol,k) /avo !PJC convert kg/kg to g/cm3
          !-----------------------------------------------------------------------
          !     ... SO2 + H2O2 --> S(VI)
          !-----------------------------------------------------------------------
          if ( aq_so2_h2o2_ndx > 0 ) then
             rxt(:,k,aq_so2_h2o2_ndx) = lwc(:) * 1.0e-03_r8 * avo * &
                  K_AQ * &
                  exp(ER_AQ * ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  HENRY298_SO2 * &
                  K298_SO2_HSO3 * &
                  HENRY298_H2O2 * &
                  exp(((H298_SO2 + H298_SO2_HSO3 + H298_H2O2) / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (R_CONC * temp(:ncol,k))**2.0e+00_r8 / &
                  (1.0e+00_r8 + 13.0e+00_r8 * 10.0e+00_r8**(-pH))
          endif
          !-----------------------------------------------------------------------
          !     ... SO2 + O3 --> S(VI)
          !-----------------------------------------------------------------------
          if (aq_so2_o3_ndx >0) then
             rxt(:,k,aq_so2_o3_ndx) = lwc(:) * 1.0e-03_r8 * avo * &
                  HENRY298_SO2 * exp((H298_SO2 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (K0_AQ * exp(ER0_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) + &
                  K298_SO2_HSO3 * exp((H298_SO2_HSO3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (K1_AQ * exp(ER1_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) / &
                  10.0e+00_r8**(-pH) + K2_AQ * exp(ER2_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  K298_HSO3_SO3 * exp((H298_HSO3_SO3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) / &
                  (10.0e+00_r8**(-pH))**2.0e+00_r8) ) * &
                  HENRY298_O3 * exp((H298_O3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (R_CONC * temp(:ncol,k))**2.0e+00_r8
          endif
       endif

    end do level_loop
    
!-----------------------------------------------------------------
! 	... the ionic rates
!-----------------------------------------------------------------
    if ( has_ion_rxts ) then
       level_loop2 : do k = 1,pver
 	   tp(:ncol)         = (2._r8*tempi(:ncol,k) + temp(:ncol,k)) / ( 3._r8 * t0 )
	   tp(:)             = max( min( tp(:),20._r8 ),1._r8 )
	   rxt(:,k,ion1_ndx) = 2.82e-11_r8 + tp(:)*(-7.74e-12_r8 + tp(:)*(1.073e-12_r8  &
			 + tp(:)*(-5.17e-14_r8 + 9.65e-16_r8*tp(:))))
	   tp(:ncol)         = (.6363_r8*tempi(:ncol,k) + .3637_r8*temp(:ncol,k)) / t0
	   tp(:)             = max( min( tp(:),trlim2 ),1._r8 )
	   rxt(:,k,ion2_ndx) = 1.533e-12_r8 + tp(:)*(-5.92e-13_r8 + tp(:)*8.6e-14_r8)
	   tp(:ncol)         = 2._r8 * t0 /(tempi(:ncol,k) + temp(:ncol,k))
	   where( tp(:ncol) < trlim3 )
		  rxt(:,k,ion3_ndx)  = 1.4e-10_r8 * tp(:)**.44_r8
		  rxt(:,k,ion11_ndx) = 1.e-11_r8 * tp(:)**.23_r8
       elsewhere
		  rxt(:,k,ion3_ndx)  = 5.2e-11_r8 / tp(:)**.2_r8
	      rxt(:,k,ion11_ndx) = 3.6e-12_r8 / tp(:)**.41_r8
	   end where
	   tp(:ncol)          = t0 / tempe(:ncol,k)
	   rxt(:,k,elec1_ndx) = 4.e-7_r8 * tp(:)**.85_r8
	   rxt(:,k,elec3_ndx) = 1.8e-7_r8 * tp(:)**.39_r8
	   where( tp(:ncol) < 4._r8 )
	      rxt(:,k,elec2_ndx) = 2.7e-7_r8 * tp(:)**.7_r8
	   elsewhere
	      rxt(:,k,elec2_ndx) = 1.6e-7_r8 * tp(:)**.55_r8
	   end where
	end do level_loop2
     endif

!-----------------------------------------------------------------
!	... tropospheric "aerosol" rate constants
!-----------------------------------------------------------------
     if ( het1_ndx > 0 .AND. (.NOT. usr_N2O5_aer_ndx > 0) ) then
         amas = 4._r8*pi*rm1**3*den/3._r8            ! each mean particle(r=0.1u) mass (g)
         do k = 1,pver
!-------------------------------------------------------------------------
! 	... estimate humidity effect on aerosols (from Shettle and Fenn, 1979)
!           xr is a factor of the increase aerosol radii with hum (hum=0., factor=1)
!-------------------------------------------------------------------------
            xr(:)     = .999151_r8 + relhum(:ncol,k)*(1.90445_r8 + relhum(:ncol,k)*(-6.35204_r8 + relhum(:ncol,k)*5.32061_r8))
!-------------------------------------------------------------------------
! 	... estimate sulfate particles surface area (cm2/cm3) in each grid
!-------------------------------------------------------------------------
            if ( carma_do_hetchem ) then
               sur(:ncol) = strato_sad(:ncol,k)
            else
              sur(:)  = sulfate(:,k)*m(:,k)/avo*wso4 &              ! xform mixing ratio to g/cm3
                        / amas &                                    ! xform g/cm3 to num particels/cm3
                        * fare &                                    ! xform num particels/cm3 to cm2/cm3
                        * xr(:)*xr(:)                               ! humidity factor
            endif
!-----------------------------------------------------------------
!	... compute the "aerosol" reaction rates
!-----------------------------------------------------------------
!             k = gam * A * velo/4
!
!       where velo = sqrt[ 8*bk*T/pi/(w/av) ]
!             bk = 1.381e-16
!             av = 6.02e23
!             w  = 108 (n2o5)  HO2(33)  CH2O (30)  NH3(15)  
!
!       so that velo = 1.40e3*sqrt(T)  (n2o5)   gama=0.1
!       so that velo = 2.53e3*sqrt(T)  (HO2)    gama>0.2
!       so that velo = 2.65e3*sqrt(T)  (CH2O)   gama>0.022
!       so that velo = 3.75e3*sqrt(T)  (NH3)    gama=0.4
!--------------------------------------------------------
!-----------------------------------------------------------------
!	... use this n2o5 -> 2*hno3 only in tropopause
!-----------------------------------------------------------------
	    rxt(:,k,het1_ndx) = rxt(:,k,het1_ndx) &
                                +.25_r8 * gam1 * sur(:) * 1.40e3_r8 * sqrt( temp(:ncol,k) )
         end do
      end if

!lke++
!-----------------------------------------------------------------
!      ... CO tags
!-----------------------------------------------------------------
      if( usr_CO_OH_b_ndx > 0 ) then
         if( usr_COhc_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_COhc_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_COme_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_COme_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO01_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO01_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO02_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO02_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO03_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO03_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO04_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO04_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO05_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO05_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO06_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO06_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO07_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO07_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO08_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO08_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO09_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO09_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO10_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO10_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO11_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO11_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO12_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO12_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO13_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO13_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO14_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO14_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO15_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO15_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO16_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO16_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO17_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO17_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO18_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO18_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO19_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO19_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO20_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO20_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO21_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO21_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO22_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO22_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO23_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO23_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO24_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO24_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO25_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO25_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO26_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO26_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO27_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO27_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO28_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO28_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO29_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO29_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO30_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO30_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO31_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO31_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO32_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO32_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO33_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO33_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO34_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO34_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO35_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO35_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO36_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO36_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO37_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO37_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO38_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO38_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO39_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO39_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO40_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO40_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO41_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO41_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO42_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO42_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
      end if
!lke--

  end subroutine usrrxt

      subroutine usrrxt_hrates( rxt, tempn, tempi, tempe, invariants, &
				h2ovmr, pmid, m, ncol, kbot )
!-----------------------------------------------------------------
!        ... set the user specified reaction rates for heating
!-----------------------------------------------------------------

      use shr_kind_mod,  only : r8 => shr_kind_r8
      use chem_mods,     only : nfs, rxntot
      use ppgrid,        only : pver, pcols

      implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
      integer, intent(in)     :: ncol                         ! number columns in chunk
      integer, intent(in)     :: kbot                         ! heating levels
      real(r8), intent(in)    :: tempn(pcols,pver)            ! neutral temperature (K)
      real(r8), intent(in)    :: tempi(pcols,pver)            ! ion temperature (K)
      real(r8), intent(in)    :: tempe(pcols,pver)            ! electron temperature (K)
      real(r8), intent(in)    :: m(ncol,pver)                 ! total atm density (1/cm^3)
      real(r8), intent(in)    :: h2ovmr(ncol,pver)            ! water vapor (vmr)
      real(r8), intent(in)    :: pmid(pcols,pver)             ! midpoint pressure (Pa)
      real(r8), intent(in)    :: invariants(ncol,pver,nfs)    ! invariants density (1/cm^3)
      real(r8), intent(inout) :: rxt(ncol,pver,rxntot)        ! gas phase rates
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------

      integer  ::  k
      real(r8), dimension(ncol) :: &
                   tp, &
                   tinv, &
                   ko, &
                   kinf, &
                   fc, &
                   xr                       ! factor to increase particle radii depending on rel hum

!-----------------------------------------------------------------
!	... o + o2 + m --> o3 + m
!-----------------------------------------------------------------
      do k = 1,kbot
         tinv(:ncol)       = 1._r8 / tempn(:ncol,k)
         tp(:)             = 300._r8 * tinv(:)
         rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8

!-----------------------------------------------------------------
!	... o + o + m -> o2 + m
!-----------------------------------------------------------------
         rxt(:,k,usr_O_O_ndx) = 2.76e-34_r8 * exp( 720.0_r8*tinv(:) )

!-----------------------------------------------------------------
!	... ho2 + ho2 --> h2o2
!	Note: this rate involves the water vapor number density
!-----------------------------------------------------------------
         ko(:)   = 3.5e-13_r8 * exp( 430._r8*tinv(:) )
         kinf(:) = 1.7e-33_r8 * m(:,k) * exp( 1000._r8*tinv(:) )
         fc(:)   = 1._r8 + 1.4e-21_r8 * m(:,k) * h2ovmr(:,k) * exp( 2200._r8*tinv(:) )
         rxt(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)

      end do

!-----------------------------------------------------------------
! 	... the ionic rates
!-----------------------------------------------------------------
      if ( has_ion_rxts ) then
         level_loop2 :  do k = 1,kbot
            tp(:ncol)         = (2._r8*tempi(:ncol,k) + tempn(:ncol,k)) / ( 3._r8 * t0 )
            tp(:)             = max( min( tp(:),20._r8 ),1._r8 )
            rxt(:,k,ion1_ndx) = 2.82e-11_r8 + tp(:)*(-7.74e-12_r8 + tp(:)*(1.073e-12_r8  &
                 + tp(:)*(-5.17e-14_r8 + 9.65e-16_r8*tp(:))))
            tp(:ncol)         = (.6363_r8*tempi(:ncol,k) + .3637_r8*tempn(:ncol,k)) / t0
            tp(:)             = max( min( tp(:),trlim2 ),1._r8 )
            rxt(:,k,ion2_ndx) = 1.533e-12_r8 + tp(:)*(-5.92e-13_r8 + tp(:)*8.6e-14_r8)
            tp(:ncol)         = 2._r8 * t0 /(tempi(:ncol,k) + tempn(:ncol,k))
            where( tp(:ncol) < trlim3 )
               rxt(:,k,ion3_ndx)  = 1.4e-10_r8 * tp(:)**.44_r8
            elsewhere
               rxt(:,k,ion3_ndx)  = 5.2e-11_r8 / tp(:)**.2_r8
            endwhere
            tp(:ncol)          = t0 / tempe(:ncol,k)
            rxt(:,k,elec1_ndx) = 4.e-7_r8 * tp(:)**.85_r8
            rxt(:,k,elec3_ndx) = 1.8e-7_r8 * tp(:)**.39_r8
            where( tp(:ncol) < 4._r8 )
               rxt(:,k,elec2_ndx) = 2.7e-7_r8 * tp(:)**.7_r8
            elsewhere
               rxt(:,k,elec2_ndx) = 1.6e-7_r8 * tp(:)**.55_r8
            endwhere
         end do level_loop2
      endif
      end subroutine usrrxt_hrates

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine comp_exp( x, y, n )

    implicit none

    real(r8), intent(out) :: x(:)
    real(r8), intent(in)  :: y(:)
    integer,  intent(in)  :: n
    
#ifdef IBM
    call vexp( x, y, n )
#else
    x(:n) = exp( y(:n) )
#endif

  end subroutine comp_exp

  !-------------------------------------------------------------------------
  !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
  !-------------------------------------------------------------------------
  function hetrxtrate( sfc, dm_aer, dg_gas, c_gas, gamma_gas ) result(rate)

    real(r8), intent(in) :: sfc(:)
    real(r8), intent(in) :: dm_aer(:)
    real(r8), intent(in) :: dg_gas
    real(r8), intent(in) :: c_gas
    real(r8), intent(in) :: gamma_gas
    real(r8) :: rate

    real(r8),allocatable :: rxt(:)
    integer :: n, i

    n = size(sfc)

    allocate(rxt(n))
    do i=1,n
       rxt(i) = sfc(i) / (0.5_r8*dm_aer(i)/dg_gas + (4._r8/(c_gas*gamma_gas)))
    enddo
    if (n==4) then
       ! b4b kludge
       rate = rxt(1) + rxt(2) + rxt(3) + rxt(4)
    else
       rate = sum(rxt)
    endif
    deallocate(rxt)

  endfunction hetrxtrate

#ifdef MODAL_AERO
  !-------------------------------------------------------------------------
  ! provides aerosol surface area info for modal aerosols
  !-------------------------------------------------------------------------
  subroutine surfarea( mmr, pmid, temp, pbuf, ncol, sfc, dm_aer, sad_total )
    use modal_aero_data, only : ntot_amode,nspec_amode,alnsg_amode
    use physics_buffer,  only : physics_buffer_desc, pbuf_get_field, pbuf_get_index
    use mo_constants,    only : pi
    use ref_pres,        only : top_lev=>trop_cloud_top_lev

    ! arguments
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer,  intent(in)    :: ncol
    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(out)   :: sad_total(:,:)

    ! local vars

    real(r8), target :: sad_mode(pcols,pver,ntot_amode)
    real(r8) :: rho_air
    real(r8), pointer, dimension(:,:,:) :: dgnumwet
    integer :: l,m
    integer :: i,k
!
    real(r8) :: chm_mass,tot_mass
!

    call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet )
    !
    ! compute surface aero for each mode; however, at this point we only use Aitken mode (mode 2 in MAM3; how
    ! can we move from hard-wiring this?) as the surface area for chemical reactions.
    !
    sad_mode = 0._r8
    sad_total = 0._r8
    do k = top_lev,pver
       do i = 1,ncol
          rho_air = pmid(i,k)/(temp(i,k)*287.04_r8)
          do l=1,ntot_amode
!
! compute a mass weighting of the number
!
             tot_mass = 0._r8
             chm_mass = 0._r8
             do m=1,nspec_amode(l)
               if ( index_tot_mass(l,m) > 0 ) &
                    tot_mass = tot_mass + mmr(i,k,index_tot_mass(l,m))
               if ( index_chm_mass(l,m) > 0 ) &
                    chm_mass = chm_mass + mmr(i,k,index_chm_mass(l,m))
             end do
             if ( tot_mass > 0._r8 ) then
               sad_mode(i,k,l) = chm_mass/tot_mass * &
                    mmr(i,k,num_idx(l))*rho_air*pi*dgnumwet(i,k,l)**2*&
                    exp(2*alnsg_amode(l)**2)  ! m^2/m^3
               sad_mode(i,k,l) = 1.e-2_r8 * sad_mode(i,k,l) ! cm^2/cm^3
             else
               sad_mode(i,k,l) = 0._r8
             end if
          end do
!
! old code
!
!          sad_total(i,k) = sad_mode(i,k,aitken_idx)
!
! new code
!
          sad_total(i,k) = sum(sad_mode(i,k,:))
!
       enddo
    enddo

    sfc(:,:,:) = sad_mode(:,:,:)     ! aitken_idx:aitken_idx) 
    dm_aer(:,:,:)  = dgnumwet(:,:,:) ! aitken_idx:aitken_idx) 
    dm_aer(:,1:top_lev-1,:)  = 0._r8
    
  end subroutine surfarea
#else
  !-------------------------------------------------------------------------
  ! provides aerosol surface area info
  !-------------------------------------------------------------------------
  subroutine surfarea( mmr,relhum, pmid,temp, strato_sad, sulfate, m, ltrop, ncol, sfc, dm_aer, sad_total )
    use mo_constants, only : pi, avo => avogadro

    ! arguments
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: strato_sad(:,:)
    real(r8), intent(in)    :: sulfate(:,:)
    real(r8), intent(in)    :: m(:,:)
    integer,  intent(in)    :: ltrop(:)
    integer,  intent(in)    :: ncol
    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(out)   :: sad_total(:,:)

    ! local vars

    integer  :: i,k
    real(r8) :: rho_air
    real(r8) :: v, n, n_exp, r_rd, r_sd
    real(r8) :: dm_sulf, dm_sulf_wet, log_sd_sulf, sfc_sulf, sfc_nit
    real(r8) :: dm_orgc, dm_orgc_wet, log_sd_orgc, sfc_oc, sfc_soa
    real(r8) :: dm_bc, dm_bc_wet, log_sd_bc, sfc_bc
    real(r8) :: rxt_sulf, rxt_nit, rxt_oc, rxt_soa
    real(r8) :: c_n2o5, c_ho2, c_no2, c_no3
    real(r8) :: s_exp

    !-----------------------------------------------------------------
    ! 	... parameters for log-normal distribution by number
    ! references:
    !   Chin et al., JAS, 59, 461, 2003
    !   Liao et al., JGR, 108(D1), 4001, 2003
    !   Martin et al., JGR, 108(D3), 4097, 2003
    !-----------------------------------------------------------------
    real(r8), parameter :: rm_sulf  = 6.95e-6_r8        ! mean radius of sulfate particles (cm) (Chin)
    real(r8), parameter :: sd_sulf  = 2.03_r8           ! standard deviation of radius for sulfate (Chin)
    real(r8), parameter :: rho_sulf = 1.7e3_r8          ! density of sulfate aerosols (kg/m3) (Chin) 

    real(r8), parameter :: rm_orgc  = 2.12e-6_r8        ! mean radius of organic carbon particles (cm) (Chin)
    real(r8), parameter :: sd_orgc  = 2.20_r8           ! standard deviation of radius for OC (Chin)
    real(r8), parameter :: rho_orgc = 1.8e3_r8          ! density of OC aerosols (kg/m3) (Chin)

    real(r8), parameter :: rm_bc    = 1.18e-6_r8        ! mean radius of soot/BC particles (cm) (Chin)
    real(r8), parameter :: sd_bc    = 2.00_r8           ! standard deviation of radius for BC (Chin)
    real(r8), parameter :: rho_bc   = 1.0e3_r8          ! density of BC aerosols (kg/m3) (Chin)

    real(r8), parameter :: mw_so4 = 98.e-3_r8     ! so4 molecular wt (kg/mole)

    integer  ::  irh, rh_l, rh_u
    real(r8) ::  factor, rfac_sulf, rfac_oc, rfac_bc, rfac_ss

    !-----------------------------------------------------------------
    ! 	... table for hygroscopic growth effect on radius (Chin et al)
    !           (no growth effect for mineral dust)
    !-----------------------------------------------------------------
    real(r8), dimension(7) :: table_rh, table_rfac_sulf, table_rfac_bc, table_rfac_oc, table_rfac_ss

    data table_rh(1:7)        / 0.0_r8, 0.5_r8, 0.7_r8, 0.8_r8, 0.9_r8, 0.95_r8, 0.99_r8/
    data table_rfac_sulf(1:7) / 1.0_r8, 1.4_r8, 1.5_r8, 1.6_r8, 1.8_r8, 1.9_r8,  2.2_r8/
    data table_rfac_oc(1:7)   / 1.0_r8, 1.2_r8, 1.4_r8, 1.5_r8, 1.6_r8, 1.8_r8,  2.2_r8/
    data table_rfac_bc(1:7)   / 1.0_r8, 1.0_r8, 1.0_r8, 1.2_r8, 1.4_r8, 1.5_r8,  1.9_r8/
    data table_rfac_ss(1:7)   / 1.0_r8, 1.6_r8, 1.8_r8, 2.0_r8, 2.4_r8, 2.9_r8,  4.8_r8/

    !-----------------------------------------------------------------
    ! 	... exponent for calculating number density
    !-----------------------------------------------------------------
    n_exp = exp( -4.5_r8*log(sd_sulf)*log(sd_sulf) )

    dm_sulf = 2._r8 * rm_sulf
    dm_orgc = 2._r8 * rm_orgc
    dm_bc   = 2._r8 * rm_bc

    log_sd_sulf = log(sd_sulf)
    log_sd_orgc = log(sd_orgc)
    log_sd_bc   = log(sd_bc)

    ver_loop: do k = 1,pver
       col_loop: do i = 1,ncol
          !-------------------------------------------------------------------------
          ! 	... air density (kg/m3)
          !-------------------------------------------------------------------------
          rho_air = pmid(i,k)/(temp(i,k)*287.04_r8)
          !-------------------------------------------------------------------------
          !       ... aerosol growth interpolated from M.Chin's table
          !-------------------------------------------------------------------------
          if (relhum(i,k) >= table_rh(7)) then
             rfac_sulf = table_rfac_sulf(7)
             rfac_oc = table_rfac_oc(7)
             rfac_bc = table_rfac_bc(7)
          else
             do irh = 2,7
                if (relhum(i,k) <= table_rh(irh)) then
                   exit
                end if
             end do
             rh_l = irh-1
             rh_u = irh

             factor = (relhum(i,k) - table_rh(rh_l))/(table_rh(rh_u) - table_rh(rh_l))

             rfac_sulf = table_rfac_sulf(rh_l) + factor*(table_rfac_sulf(rh_u) - table_rfac_sulf(rh_l))
             rfac_oc = table_rfac_oc(rh_u) + factor*(table_rfac_oc(rh_u) - table_rfac_oc(rh_l))
             rfac_bc = table_rfac_bc(rh_u) + factor*(table_rfac_bc(rh_u) - table_rfac_bc(rh_l))
          end if

          dm_sulf_wet = dm_sulf * rfac_sulf
          dm_orgc_wet = dm_orgc * rfac_oc
          dm_bc_wet = dm_bc * rfac_bc

          dm_bc_wet   = min(dm_bc_wet  ,50.e-6_r8) ! maximum size is 0.5 micron (Chin)
          dm_orgc_wet = min(dm_orgc_wet,50.e-6_r8) ! maximum size is 0.5 micron (Chin)


          !-------------------------------------------------------------------------
          ! 	... sulfate aerosols
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !       ... use ubvals climatology for stratospheric sulfate surface area density
          !-------------------------------------------------------------------------
          if( k < ltrop(i) ) then
             sfc_sulf = strato_sad(i,k)
             if ( het1_ndx > 0 ) then
                sfc_sulf = 0._r8        ! reaction already taken into account in mo_strato_rates.F90
             end if
          else

             if( so4_ndx > 0 ) then
                !-------------------------------------------------------------------------
                ! convert mass mixing ratio of aerosol to cm3/cm3 (cm^3_aerosol/cm^3_air)
                ! v=volume density (m^3/m^3)
                ! rho_aer=density of aerosol (kg/m^3)
                ! v=m*rho_air/rho_aer   [kg/kg * (kg/m3)_air/(kg/m3)_aer]
                !-------------------------------------------------------------------------
                v = mmr(i,k,so4_ndx) * rho_air/rho_sulf
                !-------------------------------------------------------------------------
                ! calculate the number density of aerosol (aerosols/cm3)
                ! assuming a lognormal distribution
                ! n  = (aerosols/cm3)
                ! dm = geometric mean diameter
                !
                ! because only the dry mass of the aerosols is known, we
                ! use the mean dry radius
                !-------------------------------------------------------------------------
                n  = v * (6._r8/pi)*(1._r8/(dm_sulf**3._r8))*n_exp
                !-------------------------------------------------------------------------
                ! find surface area of aerosols using dm_wet, log_sd 
                !  (increase of sd due to RH is negligible)
                ! and number density calculated above as distribution
                ! parameters
                ! sfc = surface area of wet aerosols (cm^2/cm^3)
                !-------------------------------------------------------------------------
                s_exp    = exp(2._r8*log_sd_sulf*log_sd_sulf)
                sfc_sulf = n * pi * (dm_sulf_wet**2._r8) * s_exp

             else
                !-------------------------------------------------------------------------
                !  if so4 not simulated, use off-line sulfate and calculate as above
                !  convert sulfate vmr to volume density of aerosol (cm^3_aerosol/cm^3_air)           
                !-------------------------------------------------------------------------
                v = sulfate(i,k) * m(i,k) * mw_so4 / (avo * rho_sulf) *1.e6_r8
                n  = v * (6._r8/pi)*(1._r8/(dm_sulf**3._r8))*n_exp
                s_exp    = exp(2._r8*log_sd_sulf*log_sd_sulf)
                sfc_sulf = n * pi * (dm_sulf_wet**2._r8) * s_exp

             end if
          end if

          !-------------------------------------------------------------------------
          ! ammonium nitrate (follow same procedure as sulfate, using size and density of sulfate)
          !-------------------------------------------------------------------------
          if( nit_ndx > 0 ) then
             v = mmr(i,k,nit_ndx) * rho_air/rho_sulf
             n  = v * (6._r8/pi)*(1._r8/(dm_sulf**3._r8))*n_exp
             s_exp   = exp(2._r8*log_sd_sulf*log_sd_sulf)
             sfc_nit = n * pi * (dm_sulf_wet**2._r8) * s_exp
          else
             sfc_nit = 0._r8
          end if

          !-------------------------------------------------------------------------
          ! hydrophylic organic carbon (follow same procedure as sulfate)
          !-------------------------------------------------------------------------
          if( oc2_ndx > 0 ) then
             v = mmr(i,k,oc2_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3))*n_exp
             s_exp    = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_oc   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_oc = 0._r8
          end if

          !-------------------------------------------------------------------------
          ! secondary organic carbon (follow same procedure as sulfate)
          !-------------------------------------------------------------------------
          if( soa_ndx > 0 ) then
             v = mmr(i,k,soa_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_soa   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_soa = 0._r8
          end if

          !-------------------------------------------------------------------------
          ! black carbon (follow same procedure as sulfate)
          !-------------------------------------------------------------------------
          if( cb2_ndx > 0 ) then
             v = mmr(i,k,cb2_ndx) * rho_air/rho_bc
             n  = v * (6._r8/pi)*(1._r8/(dm_bc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_bc*log_sd_bc)
             sfc_bc   = n * pi * (dm_bc_wet**2._r8) * s_exp
          else
             sfc_bc = 0._r8
          end if

          sfc(i,k,:) = (/ sfc_sulf, sfc_nit, sfc_oc, sfc_soa /)
          dm_aer(i,k,:) = (/ dm_sulf_wet,dm_sulf_wet,dm_orgc_wet,dm_orgc_wet /)

          !-------------------------------------------------------------------------
          !  	... add up total surface area density for output
          !-------------------------------------------------------------------------
          sad_total(i,k) = sfc_sulf + sfc_nit + sfc_oc + sfc_soa + sfc_bc

       enddo col_loop
    enddo ver_loop

  endsubroutine surfarea
#endif

end module mo_usrrxt
