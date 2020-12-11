    module module_ecpp_td2clm

        use ecpp_modal_aero_activate, only: parampollu_tdx_activate1
        use ecpp_modal_cloudchem,  only: parampollu_tdx_cldchem
        use ecpp_modal_wetscav,    only: parampollu_tdx_wetscav_2
        use perf_mod
!==Guangxing Lin
        !use abortutils, only: endrun
        use cam_abortutils, only: endrun
!==Guangxing Lin
        use physics_buffer, only : physics_buffer_desc

    implicit none


    integer, parameter :: jgrp_up=2, jgrp_dn=3


    contains

!-----------------------------------------------------------------------
!
! rce 2005-mar-10 - created
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine parampollu_td240clm(                           &
        ktau, dtstep, ktau_pp_in, dtstep_pp,          &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, zbnd, wbnd_bar,                       &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        kdraft_bot_ecpp, kdraft_top_ecpp,                 &
        mtype_updnenv_ecpp,                               &
        mfbnd_ecpp,                                       &
        abnd_tavg_ecpp, acen_tavg_ecpp,                   &
        acen_tfin_ecpp, acen_tbeg_ecpp, acen_prec_ecpp,   &
        rh_sub2, qcloud_sub2, qlsink_sub2,                &
        precr_sub2, precs_sub2,                           &
                del_cldchem3d,  del_rename3d,                     &
                del_wetscav3d, del_wetresu3d,     & 
                del_activate3d, del_conv3d,                       &
                del_chem_clm_cldchem, del_chem_clm_rename, del_chem_clm_wetscav,       &
                aqso4_h2o2, aqso4_o3, xphlwc3d,                   &
        it,      jt,      kts,ktebnd,ktecen, pbuf         )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_td240clm is a top level routine for doing
!    ecpp parameterized pollutants calculations on a single column
!    of the host-code grid
!
! this version uses the hybrid time-dependent up/dndraft formulation
!   the up and dndrafts are time-dependent, rather than steady state, 
!   with a lifetime equal "draft_lifetime"
!   in the hybrid formulation, the host-code column is conceptually
!   divided into ntstep_hybrid == (draft_lifetime/dtstep_pp) pieces
!   time integrations over dtstep_pp are done for each piece, sequentially
!   the up and downdrafts start "fresh" in the first piece
!   at the end of each "piece integration", the up and downdrafts are
!   shifted into the next piece
!   the the drafts evolve over time = draft_lifetime, but different 
!   pieces of the environment are affected by different aged drafts
!   the hybrid approach avoids two problems of the original time-dependent
!   up/dndraft formulation:
!   (a) having to store draft information (specifically aerosol mixing
!   ratios in the drafts sub-classes) from one host-code time-step to
!   the next
!   (b) having to determine when drafts should be re-initialized
!
!-----------------------------------------------------------------------


    use module_data_mosaic_asect, only:  ai_phase, cw_phase, nphase_aer

    use module_data_ecpp1

    use module_data_mosaic_asect, only:  is_aerosol, iphase_of_aerosol, isize_of_aerosol, itype_of_aerosol,   &
                inmw_of_aerosol, laicwpair_of_aerosol

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message, &
                                     parampollu_1clm_set_opts

!==Guangxing Lin                                    
        !use abortutils, only: endrun
        use cam_abortutils, only: endrun
!==Guangxing Lin                                    

!   arguments
    integer, intent(in) ::                  &
        ktau, ktau_pp_in,           &
        it, jt, kts, ktebnd, ktecen
!   ktau - time step number
!   ktau_pp_in - time step number for "parameterized pollutants" calculations
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!   chem_driver and routines under it do calculations
!   over these spatial indices.

    integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)
!   these control diagnostic output
    
    real(r8), intent(in) :: dtstep, dtstep_pp
!   dtstep - main model time step (s)
!   dtstep_pp - time step (s) for "parameterized pollutants" calculations

    real(r8), intent(in), dimension( kts:ktecen ) ::   &
        tcen_bar, pcen_bar, rhocen_bar, dzcen
    real(r8), intent(in), dimension( kts:ktebnd ) ::   &
        rhobnd_bar, wbnd_bar, zbnd
!   tcen_bar - temperature (K) at layer centers
!   rhocen_bar, rhobnd_bar - dry air density (kg/m^3) at layer centers and boundaries
!   pcen_bar - air pressure (Pa) at layer centers
!   wbnd_bar - vertical velocity (m/s) at layer boundaries
!   zbnd - elevation (m) at layer boundaries
!   dzcen - layer thicknesses (m)

    real(r8), intent(inout), dimension( kts:ktecen, 1:num_chem_ecpp ) :: &
        chem_bar
!   chem_bar - mixing ratios of trace gase (ppm) and aerosol species
!   (ug/kg for mass species, #/kg for number species)

!   NOTE - tcen_bar through chem_bar are all grid-cell averages
!       (on the host-code grid)

    integer, intent(in) :: ncls_ecpp
!   ncls_ecpp - number of ecpp transport classes in the grid column

    integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        kdraft_bot_ecpp, kdraft_top_ecpp,   &
        mtype_updnenv_ecpp
!   kdraft_bot_ecpp = lowest  layer in/thru which sub-area transport occurs
!                   = lowest  layer for which massflux != 0 at layer upper boundary
!                                          OR areafrac != 0 at layer center
!                   >= kts
!   kdraft_top_ecpp = highest layer in/thru which sub-area transport occurs
!                   = highest layer for which massflux != 0 at layer lower boundary
!                                          OR areafrac != 0 at layer center
!                   <= kte-1
!   mtype_updnenv_ecpp - transport-class (updraft, downdraft, or quiescent)

    real(r8), intent(in), dimension( kts:ktebnd, 0:2, 0:maxcls_ecpp ) ::   &
        abnd_tavg_ecpp, mfbnd_ecpp 
!   real(r8), intent(in), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
!   acen_tavg_ecpp, acen_tbeg_ecpp, acen_prec_ecpp
       real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
        acen_tavg_ecpp, acen_tbeg_ecpp, acen_prec_ecpp
    real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
        acen_tfin_ecpp
!   abnd_tavg_ecpp - sub-class fractional area (--) at layer bottom boundary
!   acen_tavg_ecpp, acen_tbeg_ecpp, acen_tfin_ecpp - sub-class fractional area (--) 
!   at layer centers
!   _tavg_ is average for full time period (=dtstep_pp_in)
!       _tbeg_ is average at beginning of time period
!       _tfin_ is average for end-portion of time period
!   acen_prec_ecpp - fractional area (---) of the portion of a sub-class that
!   has precipitation
!   0 <= acen_prec_ecpp(:,:,:)/acen_tavg_ecpp(:,:,:) <= 1
!   mfbnd_ecpp - sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
!
!   NOTE 1 - these 6 xxx_ecpp arrays contain statistics from the crm 
!   post-processor or interface.
!   Each array has a xxx_use array that contains "checked and adjusted values",
!   and those values are the ones that are used.
!   NOTE 2 - indexing for these arrays
!   the first index is vertical layer
!   the second index (0:2):  1=clear, 2=cloudy, and 0=clear+cloudy combined
!   the third index is transport class

    real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
        rh_sub2, qcloud_sub2, qlsink_sub2, precr_sub2, precs_sub2
!   rh_sub2 - relative humidity (0-1) at layer center
!   qcloud_sub2 - cloud water mixing ratio (kg/kg) at layer center
!   qlsink_sub2 - cloud-water first-order loss rate to precipitation (kg/kg/s) at layer center
!   precr_sub2 - liquid (rain) precipitation rate (kg/m2/s) at layer center
!   precsolid_sub2 - solid (snow,graupel,...) precipitation rate (kg/m2/s) at layer center
!
!   NOTE - indexing for these arrays
!   the first index is vertical layer
!   the second index (0:2) is:  1=clear, 2=cloudy
!   the third index is transport class
!   the fourth index (0:2) is:  1=non-precipitating, 2=precipitating

      real(r8), intent(out), dimension(kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:num_chem_ecpp )  ::       &
                                  del_cldchem3d,      &            ! 3D change in chem_sub from aqueous chemistry
                                  del_rename3d,       &            ! 3D change in chem_sub from renaming (modal merging)
                                  del_wetscav3d,      &            ! 3D change in chem_sub from wet deposition
                                  del_wetresu3d

      real(r8), intent(out), dimension(kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp )  ::       &
                                  del_activate3d            ! 3D change in chem_sub from activation/resuspension

      real(r8), intent(out), dimension(kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp )  ::       &
                                  del_conv3d            ! 3D change in chem_sub from convective transport 

        real(r8), intent(out) :: aqso4_h2o2,            &         ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
                                   aqso4_o3                          ! SO4 aqueous phase chemistry due to O3 (kg/m2)

        real(r8), intent(out), dimension(kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2)  ::  &
                                   xphlwc3d                          ! pH value multiplied by lwc

        real(r8), intent(out), dimension( 1:num_chem_ecpp ) :: del_chem_clm_cldchem, del_chem_clm_rename, del_chem_clm_wetscav
        type(physics_buffer_desc), pointer :: pbuf(:)

!   local variables
    integer :: activate_onoff_use
    integer :: icc, iccy, idiag,   &
        ipass_area_change, ipass_check_adjust_inputs,   &
        itstep_hybrid
    integer :: jcls, jclsbb, jgrp, jgrpbb
    integer :: k, ktau_pp
    integer :: l, laa, lbb, ll, lun, lun62
    integer :: ncls_use, ntstep_hybrid

    integer, dimension( 1:2, 1:maxcls_ecpp ) ::   &
        kdraft_bot_use, kdraft_top_use,   &
        mtype_updnenv_use

    real(r8) :: draft_area_fudge, draft_area_fudge_1m
    real(r8) :: tmpa 
        real(r8) :: tmpd, tmpe, tmpf, tmpg, tmph
        real(r8) :: tmpveca(100)
        real(r8), save :: tmpvecsva(100), tmpvecsvb(100), tmpvecsvc(100)

    real(r8), dimension( kts:ktebnd ) :: wbnd_bar_use

    real(r8), dimension( kts:ktecen ) :: rhodz_cen

    real(r8), dimension( kts:ktebnd, 0:2, 0:maxcls_ecpp ) ::   &
        abnd_tavg_use, mfbnd_use,   &
        abnd_tavg_usex1, mfbnd_usex1,   &
        ar_bnd_tavg

    real(r8), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
        acen_tavg_usex1, acen_tbeg_usex1, acen_tfin_usex1,   &
        acen_tavg_use, acen_tbeg_use, acen_tfin_use, acen_prec_use,   &
        ardz_cen_tbeg, ardz_cen_tfin,   &
        ardz_cen_tavg,   &
        ardz_cen_old, ardz_cen_new

    real(r8), dimension( kts:ktebnd, 0:2, 0:2 ) ::   &
        mfbnd_quiescn_up, mfbnd_quiescn_dn

    real(r8), dimension( kts:ktecen, 1:maxcls_ecpp, 1:num_chem_ecpp ) :: &
        chem_cls
 
    real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
        chem_sub_new, chem_sub_beg, chem_sub_ac1sv, chem_sub_hysum

    real(r8), dimension( 1:2, num_chem_ecpp ) :: chem_bar_iccfactor

        real(r8), dimension(kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp )  ::       &
                                  del_activate3da            ! 3D change in chem_sub from activation/resuspension

        

    character(len=120) :: msg


    ktau_pp = 10

        lun62 = -1
        if (idiagaa_ecpp(62) > 0) lun62 = ldiagaa_ecpp(62)

    activate_onoff_use = 0
    if ( (nphase_aer >= 2) .and.   &
         (ai_phase > 0) .and. (cw_phase > 0) )   &
        activate_onoff_use = activat_onoff_ecpp

!   in sub-classes with area ~= 0, chem_sub is set to chem_bar
!   EXCEPT for aerosol species, where activated=0 in clear,
!   and activated=interstitial=0.5*chem_bar in cloudy
    chem_bar_iccfactor(:,:) = 1.0
    if (activate_onoff_use > 0) then
        do l = param_first_ecpp, num_chem_ecpp
        if ( is_aerosol(l) ) then
            if (iphase_of_aerosol(l) == ai_phase) then
            chem_bar_iccfactor(2,l) = 1.0
            else if (iphase_of_aerosol(l) == cw_phase) then
            chem_bar_iccfactor(2,l) = 1.0
            chem_bar_iccfactor(1,l) = 1.0
            end if
        end if
        end do
    end if

!
!   output the original fields with same format as ppboxmakeinp01
!
       ll = 116
       lun = ldiagaa_ecpp(ll)
       if ((idiagaa_ecpp(ll) > 0) .and. (lun > 0)) then
           call parampollu_1clm_dumpaa(                          &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        tcen_bar,   pcen_bar, rhocen_bar, dzcen,          &
        rhobnd_bar, zbnd, wbnd_bar,                       &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        kdraft_bot_ecpp, kdraft_top_ecpp,                 &
        mtype_updnenv_ecpp,                               &
        mfbnd_ecpp, abnd_tavg_ecpp,                       &
        acen_tavg_ecpp, acen_tbeg_ecpp, acen_tfin_ecpp,   &
        it,      jt,      kts,ktebnd,ktecen,              &
        lun                                               )
    end if


!
!   check and adjust input information
!   and do startup calcs (for this parampollu timestep)
!
    do ipass_check_adjust_inputs = 1, 2

    call parampollu_check_adjust_inputs(                      &
        ipass_check_adjust_inputs,                        &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar,   pcen_bar, rhocen_bar, dzcen,          &
        rhobnd_bar, zbnd, wbnd_bar,                       &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        kdraft_bot_ecpp, kdraft_top_ecpp,                 &
        mtype_updnenv_ecpp,                               &
        mfbnd_ecpp, abnd_tavg_ecpp,                       &
        acen_tavg_ecpp, acen_tfin_ecpp, acen_prec_ecpp,   &
        wbnd_bar_use,                                     &
        ncls_use,                                         &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        mfbnd_use, mfbnd_quiescn_up, mfbnd_quiescn_dn,    &
        abnd_tavg_use,                                    &
        acen_tavg_use, acen_tfin_use, acen_prec_use,      &
        rhodz_cen,                                        &
        it,      jt,      kts,ktebnd,ktecen               )
 
!   do startup calcs (for this parampollu timestep)
    if (ipass_check_adjust_inputs == 1) then
        acen_tbeg_use(:,:,:) = acen_tbeg_ecpp(:,:,:)
    else
        call parampollu_tdx_startup(                          &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        rhocen_bar, dzcen,                                &
        chem_bar, chem_cls,                               &
        ncls_ecpp,                                        &
        acen_tbeg_ecpp,                                   &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use,                                         &
        chem_sub_beg,                                     &
        acen_tbeg_use, ardz_cen_tbeg, rhodz_cen,          &
        activate_onoff_use,                               &
        iphase_of_aerosol, laicwpair_of_aerosol           )
    end if

!   output the adjusted fields with same format as ppboxmakeinp01
    if (ipass_check_adjust_inputs == 1) then
        acen_tavg_usex1(:,:,:) = acen_tavg_use(:,:,:)
        acen_tfin_usex1(:,:,:) = acen_tfin_use(:,:,:)
        abnd_tavg_usex1(:,:,:) = abnd_tavg_use(:,:,:)
        mfbnd_usex1(    :,:,:) = mfbnd_use(    :,:,:)
        ll = 117
    else
        ll = 115
    end if

    lun = ldiagaa_ecpp(ll)
    if ((idiagaa_ecpp(ll) > 0) .and. (lun > 0)) then
        call parampollu_1clm_dumpaa(                          &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        tcen_bar,   pcen_bar, rhocen_bar, dzcen,          &
        rhobnd_bar, zbnd, wbnd_bar_use,                   &
        chem_bar,                                         &
        ncls_use,                                         &
        kdraft_bot_use,  kdraft_top_use,                  &
        mtype_updnenv_use,                                &
        mfbnd_use, abnd_tavg_use,                         &
        acen_tavg_use, acen_tbeg_use, acen_tfin_use,      &
        it,      jt,      kts,ktebnd,ktecen,              &
        lun                                               )
    end if

    end do   ! ipass_check_adjust_inputs


 
!   *** temporary exit
    if (iflag_ecpp_test_bypass_1 > 0) return
 

!   save values in these arrays
    acen_tbeg_usex1(:,:,:) = acen_tbeg_use(:,:,:)
    chem_sub_new(:,:,:,:) = chem_sub_beg(:,:,:,:)

        del_activate3d(:,:,:,:) = 0.0

!   calc "area*rho*dz" and "area*rho" arrays
    ardz_cen_tavg(:,:,:) = 0.0
    ardz_cen_tfin(:,:,:) = 0.0
    ar_bnd_tavg(:,:,:) = 0.0
    do k = kts, ktebnd
        do icc = 0, 2
        ar_bnd_tavg(  k,icc,0:ncls_use) = abnd_tavg_use(k,icc,0:ncls_use)*rhobnd_bar(k)
        if (k > ktecen) cycle
        ardz_cen_tavg(k,icc,0:ncls_use) = acen_tavg_use(k,icc,0:ncls_use)*rhodz_cen(k)
        ardz_cen_tfin(k,icc,0:ncls_use) = acen_tfin_use(k,icc,0:ncls_use)*rhodz_cen(k)
        end do
    end do


!
!   apply area changes (acen_tbeg_use --> ... --> acen_tfin_use) here
!   parampollu_opt == 2220
!   apply area changes in one step, before 15000 loop
!   parampollu_opt == 2223
!   apply area changes in two steps, before and after 15000 loop
!
    ardz_cen_old(:,:,:) = ardz_cen_tbeg(:,:,:)
    if      (parampollu_opt == 2220) then
        ardz_cen_new(:,:,:) = ardz_cen_tfin(:,:,:)
    else if (parampollu_opt == 2223) then
        ardz_cen_new(:,:,:) = ardz_cen_tavg(:,:,:)
    else
        stop
    end if

!   note about parampollu_tdx_area_change and parampollu_tdx_main_integ
!     initial values are taken from chem_sub_new
!     final   values are put   into chem_sub_new
    ipass_area_change = 1
    call parampollu_tdx_area_change(                          &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, wbnd_bar,                             &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use, ipass_area_change,                      &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        chem_sub_new,                                     &
                del_activate3d,                                   &
        mfbnd_use, ar_bnd_tavg,                           &
        ardz_cen_old, ardz_cen_new, rhodz_cen,            &
        chem_bar_iccfactor, activate_onoff_use,           &
        iphase_of_aerosol, isize_of_aerosol,              &
        itype_of_aerosol, inmw_of_aerosol,                &
        laicwpair_of_aerosol                              )



! save current chem_sub values
    chem_sub_ac1sv(:,:,:,:) = 0.0
    chem_sub_ac1sv(kts:ktecen,1:2,1:ncls_use,1:num_chem_ecpp) = &
      chem_sub_new(kts:ktecen,1:2,1:ncls_use,1:num_chem_ecpp) 
! initialize chem_sub hybrid-sum
    chem_sub_hysum(:,:,:,:) = 0.0

    ntstep_hybrid = nint( draft_lifetime / dtstep )
    ntstep_hybrid = max( 1, ntstep_hybrid ) 
    if (lun62 > 0) write(lun62,'(a,2i10)') &
        'parampollu_td240clm - ktau, ntstep_hybrid', &
        ktau, ntstep_hybrid


    del_chem_clm_cldchem(:) = 0.0
        del_chem_clm_rename(:) = 0.0
        del_cldchem3d(:,:,:,:,:) = 0.0
        del_rename3d(:,:,:,:,:) = 0.0
    del_chem_clm_wetscav(:) = 0.0
        del_wetscav3d(:,:,:,:,:) = 0.0
        del_wetresu3d(:,:,:,:,:) = 0.0
        del_activate3da(:,:,:,:) = 0.0

        aqso4_h2o2  = 0.0_r8
        aqso4_o3    = 0.0_r8
        xphlwc3d(:,:,:,:) = 0.0_r8

itstep_hybrid_loop:   &
    do itstep_hybrid = 1, ntstep_hybrid
    ktau_pp = itstep_hybrid + 100

!
!   main integration
!
    ardz_cen_old(:,:,:) = ardz_cen_new(:,:,:)

    call parampollu_tdx_main_integ(                           &
        ktau, dtstep, ktau_pp, dtstep_pp,                 &
                itstep_hybrid, ntstep_hybrid,                     &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, zbnd, wbnd_bar,                       &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use,                                         &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        chem_sub_new,                                     &
        del_chem_clm_cldchem, del_chem_clm_rename, del_chem_clm_wetscav,       &
                del_cldchem3d, del_rename3d,                      &
                del_wetscav3d, del_wetresu3d,      &
                del_activate3da,                                  &
                aqso4_h2o2, aqso4_o3, xphlwc3d,                   &
        mfbnd_use, mfbnd_quiescn_up, mfbnd_quiescn_dn,    &
        ar_bnd_tavg,                                      &
        ardz_cen_old, ardz_cen_new, rhodz_cen,            &
        acen_tavg_use, acen_prec_use,                     &
        rh_sub2, qcloud_sub2, qlsink_sub2,                &
        precr_sub2, precs_sub2,                           &
        chem_bar_iccfactor, activate_onoff_use,           &
        iphase_of_aerosol, isize_of_aerosol,              &
        itype_of_aerosol, inmw_of_aerosol,                &
        laicwpair_of_aerosol, pbuf                        )


    do l = param_first_ecpp, num_chem_ecpp
    do jcls = 1, ncls_use
!   increment chem_sub_hysum
        if ((jcls == jcls_qu) .or. (itstep_hybrid == ntstep_hybrid)) then
!   for quiescent (all steps) or up/dndrafts (final step), use chem_sub_new
        chem_sub_hysum(kts:ktecen,1:2,jcls,l) = &
        chem_sub_hysum(kts:ktecen,1:2,jcls,l) + &
          chem_sub_new(kts:ktecen,1:2,jcls,l)
        else
!   for up/dndrafts (all but final step), use chem_sub_ac1sv
        chem_sub_hysum(kts:ktecen,1:2,jcls,l) = &
        chem_sub_hysum(kts:ktecen,1:2,jcls,l) + &
        chem_sub_ac1sv(kts:ktecen,1:2,jcls,l)
        end if

!   on all but final step, prepare for next main_integ by
!   restoring jcls_qu to chem_sub_ac1sv values
        if ((jcls == jcls_qu) .and. (itstep_hybrid < ntstep_hybrid)) then
        chem_sub_new(  kts:ktecen,1:2,jcls,l) = &
        chem_sub_ac1sv(kts:ktecen,1:2,jcls,l)
        end if

!   on (after) final step, convert chem_sub_hysum to an average
!   and load into chem_sub_new
        if (itstep_hybrid == ntstep_hybrid) then
        tmpa = 1.0_r8/ntstep_hybrid
        chem_sub_new(  kts:ktecen,1:2,jcls,l) = &
        chem_sub_hysum(kts:ktecen,1:2,jcls,l)*tmpa
        end if

    end do ! jcls
    end do ! l


    end do itstep_hybrid_loop

    tmpa = ntstep_hybrid ; tmpa = 1.0_r8/tmpa
    del_chem_clm_cldchem(:) = del_chem_clm_cldchem(:)*tmpa
        del_chem_clm_rename(:) = del_chem_clm_rename(:)*tmpa
        del_cldchem3d(:,:,:,:,:) = del_cldchem3d(:,:,:,:,:) * tmpa
        del_rename3d(:,:,:,:,:) = del_rename3d(:,:,:,:,:) * tmpa
    del_chem_clm_wetscav(:) = del_chem_clm_wetscav(:)*tmpa
        del_wetscav3d(:,:,:,:,:) = del_wetscav3d(:,:,:,:,:)*tmpa
        del_wetresu3d(:,:,:,:,:) = del_wetresu3d(:,:,:,:,:)*tmpa
        del_activate3d(:,:,:,:) = del_activate3d(:,:,:,:) + del_activate3da(:,:,:,:) * tmpa

        aqso4_h2o2 = aqso4_h2o2 * tmpa
        aqso4_o3   = aqso4_o3 * tmpa
        xphlwc3d(:,:,:,:) = xphlwc3d(:,:,:,:) * tmpa


    ktau_pp = 20


!   when parampollu_opt == 2223, do 2nd half of area change here
    if (parampollu_opt == 2223) then
        ipass_area_change = 2
        ardz_cen_old(:,:,:) = ardz_cen_new(:,:,:)
        ardz_cen_new(:,:,:) = ardz_cen_tfin(:,:,:)

        call parampollu_tdx_area_change(                      &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, wbnd_bar,                             &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use, ipass_area_change,                      &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        chem_sub_new,                                     &
                del_activate3d,                                   &
        mfbnd_use, ar_bnd_tavg,                           &
        ardz_cen_old, ardz_cen_new, rhodz_cen,            &
        chem_bar_iccfactor, activate_onoff_use,           &
        iphase_of_aerosol, isize_of_aerosol,              &
        itype_of_aerosol, inmw_of_aerosol,                &
        laicwpair_of_aerosol                              )

    end if


!   do "cleanup"
    call parampollu_tdx_cleanup(                              &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        chem_bar, chem_cls,                               &
        ncls_ecpp,                                        &
        acen_tfin_ecpp,                                   &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use,                                         &
        chem_sub_beg, chem_sub_new,                       &
        del_chem_clm_cldchem, del_chem_clm_wetscav,       &
                del_cldchem3d, del_rename3d,                      &
                del_wetscav3d, del_wetresu3d,                     &
                del_activate3d, del_conv3d,                       &
        acen_tbeg_use, acen_tfin_use, rhodz_cen,          &
        activate_onoff_use,                               &
        iphase_of_aerosol, isize_of_aerosol,              &
        itype_of_aerosol, inmw_of_aerosol,                &
        laicwpair_of_aerosol                              )


!   output precip info
!
        if (ktau <= 1) then
            tmpvecsva(:) = 0.0 ; tmpvecsvb(:) = 0.0 ; tmpvecsvc(:) = 0.0
        end if
        tmpveca(:) = 0.0
        do jcls = 1, ncls_use
        do icc = 1, 2
            tmpe = max( 0.0_r8, acen_prec_use(kts,icc,jcls) )
            tmpf = max( 0.0_r8, acen_tavg_use(kts,icc,jcls) - tmpe )
            tmpg = max( 0.0_r8, precr_sub2(kts,icc,jcls,2) ) + &
                   max( 0.0_r8, precs_sub2(kts,icc,jcls,2) )
            tmph = max( 0.0_r8, precr_sub2(kts,icc,jcls,1) ) + &
                   max( 0.0_r8, precs_sub2(kts,icc,jcls,1) )
            tmpveca(1) = tmpveca(1) + tmpg
            tmpveca(2) = tmpveca(2) + tmph
            tmpveca(3) = tmpveca(3) + tmpg*tmpe
            tmpveca(4) = tmpveca(4) + tmph*tmpf
            do k = kts, ktecen
                tmpe = max( 0.0_r8, acen_prec_use(k,icc,jcls) )
                tmpf = max( 0.0_r8, acen_tavg_use(k,icc,jcls) - tmpe )
                tmpg = max( 0.0_r8, precr_sub2(k,icc,jcls,2) ) + &
                       max( 0.0_r8, precs_sub2(k,icc,jcls,2) )
                tmph = max( 0.0_r8, precr_sub2(k,icc,jcls,1) ) + &
                       max( 0.0_r8, precs_sub2(k,icc,jcls,1) )
                tmpa = tmpg*tmpe + tmph*tmpf
                if (icc == 1) then
                    tmpvecsvb(k) = tmpvecsvb(k) + tmpa
                else
                    tmpvecsvc(k) = tmpvecsvc(k) + tmpa
                end if
            end do
        end do
        end do
        tmpvecsva(1:4) = tmpvecsva(1:4) + tmpveca(1:4)
!   3600 factor converts from kg/m2/s=mm/s to mm/h
!        write(165,'(a,f7.1,4f10.3,3x,4f10.3)') 'precp at sfc', (ktau*dtstep/3600.0), &
!            tmpveca(1:4)*3600.0, tmpvecsva(1:4)*3600.0/max(ktau,1)

        if (mod(ktau,18) == 0 .and. ktau.ge.1) then
!            write(165,'(/a,i5)') 'precp for clear, cldy, both at ktau =', ktau
            tmpa = 3600.0/ktau  ! converts accumulated precip to time avg and mm/h
!            do k = ktecen, kts, -1
!                write(165,'(i5,4f10.3)') k, tmpvecsvb(k)*tmpa, tmpvecsvc(k)*tmpa, &
!                        (tmpvecsvb(k)+tmpvecsvc(k))*tmpa
!            end do
!            write(165,'(a)')
        end if


        if (mod(ktau,18) == 0 .and. ktau.ge.1) then
!            if (ktau == 18) write(166,'(a)') &
!                'ktau, k, sum(rh), sum(precr), sum(precall), sum(qcloud)'
!            write(166,'(a)')
!            do k = kts, ktecen, 5
!                write(166,'(2i4,1p,5e11.3)') &
!                    ktau, k, sum( rh_sub2(k,1:2,1:ncls_use,1:2) ), &
!                    sum( precr_sub2(k,1:2,1:ncls_use,1:2) ), &
!                    sum( precr_sub2(k,1:2,1:ncls_use,1:2)    &
!                       + precs_sub2(k,1:2,1:ncls_use,1:2) ), &
!                    sum( qcloud_sub2(k,1:2,1:ncls_use,1:2) )
!            end do
        end if

        if (mod(ktau,18) == 0 .and. ktau.ge.1) then
!            if (ktau == 18) write(167,'(a)') &
!                'ktau, k, sum(a*rh), sum(a*precr), sum(a*precall), sum(a*qcloud)'
!            write(167,'(a)')
            do k = kts, ktecen, 5
                tmpveca(:) = 0.0
                do jcls = 1, ncls_use
                do icc = 1, 2
                    tmpe = max( 0.0_r8, acen_prec_use(k,icc,jcls) )
                    tmpf = max( 0.0_r8, acen_tavg_use(k,icc,jcls) - tmpe )
                    tmpveca(1) = tmpveca(1) + &
                        tmpe*max( 0.0_r8, rh_sub2(k,icc,jcls,2) ) + &
                        tmpf*max( 0.0_r8, rh_sub2(k,icc,jcls,1) )
                    tmpveca(2) = tmpveca(2) + &
                        tmpe*max( 0.0_r8, precr_sub2(k,icc,jcls,2) ) + &
                        tmpf*max( 0.0_r8, precr_sub2(k,icc,jcls,1) )
                    tmpveca(3) = tmpveca(3) + &
                        tmpe*max( 0.0_r8, precs_sub2(k,icc,jcls,2) ) + &
                        tmpf*max( 0.0_r8, precs_sub2(k,icc,jcls,1) )
                    tmpveca(4) = tmpveca(4) + &
                        tmpe*max( 0.0_r8, qcloud_sub2(k,icc,jcls,2) ) + &
                        tmpf*max( 0.0_r8, qcloud_sub2(k,icc,jcls,1) )
                end do
                end do
                tmpveca(3) = tmpveca(3) + tmpveca(2)
!                write(167,'(2i4,1p,5e11.3)') &
!                    ktau, k, tmpveca(1:4)
            end do
        end if

!
!  all done
!
    if (lun62 > 0) write(lun62,*) '*** leaving parampollu_td240clm'
    return
    end subroutine parampollu_td240clm



!-----------------------------------------------------------------------
    subroutine parampollu_tdx_main_integ(                     &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
                itstep_hybrid, ntstep_hybrid,                     &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, zbnd, wbnd_bar,                       &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use,                                         &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        chem_sub_new,                                     &
        del_chem_clm_cldchem, del_chem_clm_rename, del_chem_clm_wetscav,       &
                del_cldchem3d,  del_rename3d,                     &
                del_wetscav3d, del_wetresu3d,                     &
                del_activate3d,                                   &
                aqso4_h2o2, aqso4_o3, xphlwc3d,                   &
        mfbnd_use, mfbnd_quiescn_up, mfbnd_quiescn_dn,    &
        ar_bnd_tavg,                                      &
        ardz_cen_old, ardz_cen_new, rhodz_cen,            &
        acen_tavg_use, acen_prec_use,                     &
        rh_sub2, qcloud_sub2, qlsink_sub2,                &
        precr_sub2, precs_sub2,                           &
        chem_bar_iccfactor, activate_onoff_use,           &
        iphase_of_aerosol, isize_of_aerosol,              &
        itype_of_aerosol, inmw_of_aerosol,                &
        laicwpair_of_aerosol, pbuf                        )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_main_integ does the "main integration"
!    of the trace-species conservation equations over time-step dtstep_pp
!
! incoming chem_sub_new holds current sub-class mixing ratios
! outgoing chem_sub_new holds updated sub-class mixing ratios
!
! treats 
!    sub-grid vertical transport and associated horizontal exchange 
!       (entrainment and detrainment)
!    activation/resuspension
!    cloud chemistry and wet removal
!
! does not treat
!    horizontal exchange associated with sub-class area changes
!
!-----------------------------------------------------------------------

    use module_data_radm2, only:  epsilc

    use module_data_mosaic_asect, only:  ai_phase, cw_phase,   &
        massptr_aer, maxd_asize, maxd_atype,   &
        ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer

    use module_data_ecpp1

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message

!   arguments
    integer, intent(in) ::                  &
        ktau, ktau_pp,              &
                itstep_hybrid, ntstep_hybrid,   &
        it, jt, kts, ktebnd, ktecen
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!   chem_driver and routines under it do calculations
!   over these spatial indices.

    integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

    real(r8), intent(in) :: dtstep, dtstep_pp
!   dtstep - main model time step (s)
!   dtstep_pp - time step (s) for "parameterized pollutants" calculations

    real(r8), intent(in), dimension( kts:ktecen ) ::   &
        tcen_bar, pcen_bar, rhocen_bar, dzcen
    real(r8), intent(in), dimension( kts:ktebnd ) ::   &
        rhobnd_bar, wbnd_bar, zbnd
!   tcen_bar - temperature (K) at layer centers
!   rhocen_bar, rhobnd_bar - dry air density (kg/m^3) at layer centers and boundaries
!   pcen_bar - air pressure (Pa) at layer centers
!   wbnd_bar - vertical velocity (m/s) at layer boundaries
!   zbnd - elevation (m) at layer boundaries
!   dzcen - layer thicknesses (m)

    real(r8), intent(in), dimension( kts:ktecen, 1:num_chem_ecpp ) :: &
        chem_bar
!   chem_bar - mixing ratios of trace gase (ppm) and aerosol species
!   (ug/kg for mass species, #/kg for number species)

    integer, intent(in) :: ncls_ecpp
!   ncls_ecpp - number of ecpp transport classes in the grid column

    integer, intent(in)    :: ncls_use

    integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        kdraft_bot_use, kdraft_top_use,   &
        mtype_updnenv_use

    real(r8), intent(inout), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
        chem_sub_new

    real(r8), intent(inout), dimension( 1:num_chem_ecpp ) :: del_chem_clm_cldchem, del_chem_clm_rename, del_chem_clm_wetscav

        real(r8), intent(inout), dimension(  kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:num_chem_ecpp ) ::    &
                  del_cldchem3d,              &              ! 3D change from aqueous chemistry
                  del_rename3d,              &              ! 3D change from renaming (modal merging)
                  del_wetscav3d,              &              ! 3D change from wet deposition
                  del_wetresu3d

        real(r8), intent(inout), dimension(  kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::    &
                  del_activate3d             ! 3D change from activation/resuspension 

        real(r8), intent(inout) :: aqso4_h2o2,            &         ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
                                   aqso4_o3                          ! SO4 aqueous phase chemistry due to O3 (kg/m2)

        real(r8), intent(inout), dimension(kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2)  ::  &
                                   xphlwc3d                          ! pH value multiplied by lwc

    real(r8), intent(inout), dimension( kts:ktebnd, 0:2, 0:maxcls_ecpp ) :: &
        mfbnd_use, ar_bnd_tavg

    real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) :: &
        ardz_cen_old, ardz_cen_new, acen_tavg_use, acen_prec_use

    real(r8), intent(inout), dimension( kts:ktebnd, 0:2, 0:2 ) ::   &
        mfbnd_quiescn_up, mfbnd_quiescn_dn

    real(r8), intent(inout), dimension( kts:ktecen ) :: rhodz_cen

    real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
        rh_sub2, qcloud_sub2, qlsink_sub2, precr_sub2, precs_sub2

    real(r8), intent(in), dimension( 1:2, num_chem_ecpp ) :: chem_bar_iccfactor

    integer, intent(in) :: activate_onoff_use

    integer, intent(in), dimension( 1:num_chem_ecpp ) ::   &
        iphase_of_aerosol, isize_of_aerosol, itype_of_aerosol,   &
        inmw_of_aerosol, laicwpair_of_aerosol
        type(physics_buffer_desc), pointer :: pbuf(:)



!   local variables
    integer, parameter :: activate_onoff_testaa = 1
    integer :: icc, iccb, iccy, ido_actres_tmp, ifrom_where,   &
        itstep_sub, itmpa, iupdn
    integer :: idiag118_pt1, idiag118_pt2, idiag118_pt3
        integer :: idiagbb_wetscav
    integer :: jcls, jclsy
    integer :: k, kb, l, la, laa, lbb, lc, lun118, lun124
    integer :: m, n, ntstep_sub
    integer, save :: ntstep_sub_sum = 0
    integer :: p1st

    integer, dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp ) ::   &
        ido_actres_horz

    logical :: not_aicw

    real(r8) :: ardz_cut
    real(r8) :: dtstep_sub
    real(r8) :: tmpa, tmpb, tmpc, tmpd
    real(r8) :: tmpcourout, tmpcourmax
    real(r8) :: tmp_ardz,   tmp_del_ardz
    real(r8) :: tmp_ardzqa, tmp_del_ardzqa
    real(r8) :: tmp_ardzqc, tmp_del_ardzqc
        real(r8) :: tmp_del_ardzqa_act, tmp_del_ardzqc_act
    real(r8) :: tmp_fmnact
    real(r8) :: tmp_qyla, tmp_qylc
    real(r8) :: tmp2dxa(0:2,0:maxcls_ecpp), tmp2dxb(0:2,0:maxcls_ecpp)
    real(r8) :: xntstep_sub_inv

    real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
        chem_sub_old
    real(r8), dimension( 1:2, 1:maxcls_ecpp, kts:ktecen ) :: &
        ent_airamt_tot, det_airamt_tot
    real(r8), dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
        ent_airamt, det_airamt
    real(r8), dimension( 1:maxd_asize, 1:maxd_atype, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
        fmact_horz, fnact_horz
    real(r8), dimension( 1:maxd_asize, 1:maxd_atype, kts:ktecen ) ::   &
        fmact_vert, fnact_vert
    real(r8), dimension( kts:ktebnd, 0:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
        tmpverta, tmphoriz

        real(r8) :: frc_ent_act    ! the fraction of updraft entrainment that may experince activation   +++mhwang
        real(r8) :: frc_tmp
        real(r8) :: abnd_up        ! cloud fraction in the upper boundary
        real(r8) :: abnd_dn        ! cloud fraction in the lower boundary

      !  call t_startf('ecpp_mainintegr') !==Guangxing Lin

    p1st = param_first_ecpp

    idiag118_pt1 = 10 * mod( max(idiagaa_ecpp(118),0)/1,   10 )
    idiag118_pt2 = 10 * mod( max(idiagaa_ecpp(118),0)/10,  10 )
    idiag118_pt3 = 10 * mod( max(idiagaa_ecpp(118),0)/100, 10 )

    lun124 = -1
    if (idiagaa_ecpp(124) > 0) lun124 = ldiagaa_ecpp(124)

        idiagbb_wetscav = 0

!
!   calc entrain/detrain amounts
!
!   first calc net (entrainment-detrainment) amount = area change
    ent_airamt_tot(:,:,:) = 0.0
    det_airamt_tot(:,:,:) = 0.0
    do jcls = 1, ncls_use
    do icc = 1, 2
    do k = kts, ktecen
        ardz_cut = afrac_cut*rhodz_cen(k)*0.3
        tmpa = max( ardz_cen_new(k,icc,jcls), ardz_cen_old(k,icc,jcls) )
        if (tmpa < ardz_cut) cycle   ! k loop

        if (jcls /= jcls_qu) then
!   this is for area change
!       tmpb = ardz_cen_new(k,icc,jcls) - ardz_cen_old(k,icc,jcls)
!   this is for vertical mass flux divergence/convergence
        tmpb = (mfbnd_use(k+1,icc,jcls) - mfbnd_use(k,icc,jcls))*dtstep_pp
        if (tmpb > 0.0) then
            ent_airamt_tot(icc,jcls,k) = tmpb
        else if (tmpb < 0.0) then
            det_airamt_tot(icc,jcls,k) = -tmpb
        end if

        else
        ! +mfbnd_quiescn_up(k+1,icc,0  ) is upwards outflow from sub-class 
        !                                at top of layer (and is >= 0)
        ! +mfbnd_quiescn_dn(k+1,0  ,icc) is dnwards inflow  to   sub-class 
        !                                at top of layer (and is <= 0)
        ! -mfbnd_quiescn_up(k  ,0, ,icc) is upwards inflow  to   sub-class 
        !                                at bottom of layer (and is <= 0)
        ! -mfbnd_quiescn_dn(k  ,icc,0  ) is dnwards outflow from sub-class 
        !                                at bottom of layer (and is >= 0)
        ! tmpb = net vertical in/outflows 
        !        (positive if net outflow, negative if net inflow)
        tmpb = ( mfbnd_quiescn_up(k+1,icc,0  )   &
               + mfbnd_quiescn_dn(k+1,0  ,icc)   &
               - mfbnd_quiescn_up(k  ,0  ,icc)   &
               - mfbnd_quiescn_dn(k  ,icc,0  ) )*dtstep_pp
        if (tmpb > 0.0) then
            ent_airamt_tot(icc,jcls,k) = tmpb
        else if (tmpb < 0.0) then
            det_airamt_tot(icc,jcls,k) = -tmpb
        end if

        end if
    end do
    end do
    end do

!   next calc detailed ent/det amounts
        call t_startf('ecpp_entdet')
    ifrom_where = 10
    call parampollu_tdx_entdet_sub1(                          &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use, ifrom_where,                            &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        ardz_cen_old, ardz_cen_new, rhodz_cen,            &
        ent_airamt_tot, det_airamt_tot,                   &
        ent_airamt, det_airamt                            )
         call t_stopf('ecpp_entdet')


!
!   calc activation/resuspension fractions associated with ent/det
!   and vertical transport
!   
    if (activate_onoff_use > 0) then
        call t_startf('ecpp_activate')
    ifrom_where = 10
    call parampollu_tdx_activate1(                            &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, wbnd_bar,                             &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use, ifrom_where, activate_onoff_use,        &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        chem_sub_new,                                     &
        mfbnd_use,                                        &
        ar_bnd_tavg,                                      &
        ent_airamt,                                       &
        ido_actres_horz, fmact_horz, fnact_horz,          &
        fmact_vert, fnact_vert, mfbnd_quiescn_up          )
        call t_stopf('ecpp_activate')
    end if


!
!   determine number of integration sub-steps
!   calc "outflow" courant number for each sub-class
!       = (sum of outflow air-mass fluxes) * dt / ardz_cen
!   calc tmpcourmax = maximum outflow courant number 
!       for all layers and sub-classes
!   select ntstep_sub (number of integration sub-steps) so that
!       (tmpcourmax/ntstep_sub) <= 1.0
!
    if (lun124 > 0) &
        write( lun124, '(/a,2i5/a)' ) 'new courout stuff -- ktau, ktau_pp', ktau, ktau_pp, &
           'k, tmpcouroutc(qu), tmpcouroutb(up), tmpcouroutb(dn)' 
    tmp2dxb(:,:) = -1.0
    tmpcourmax = 0.0
    do k = ktecen, kts, -1
        ardz_cut = afrac_cut*rhodz_cen(k)*0.3
        do jcls = 1, ncls_use
        do icc = 1, 2

!   tmpa = (air-mass leaving sub-class over dtstep_pp by vertical mass flux)
        if (jcls == jcls_qu) then
            tmpa = mfbnd_quiescn_up(k+1,icc,0) - mfbnd_quiescn_dn(k,icc,0)
        else
            tmpa = max(0.0_r8,mfbnd_use(k+1,icc,jcls)) + max(0.0_r8,-mfbnd_use(k,icc,jcls))
        end if
        tmpa = tmpa*dtstep_pp
!   tmpb = tmpa + (air-mass leaving sub-class over dtstep_pp by horizontal detrainment)
        tmpb = tmpa + max(0.0_r8,det_airamt_tot(icc,jcls,k))

!   (area*rho*dz) is fixed at ardz_cen_new during the integration loop
        tmp_ardz = ardz_cen_new(k,icc,jcls)

        if (tmp_ardz < ardz_cut) then
            tmpcourout = 0.0
        else if (tmpb > 1.0e3*tmp_ardz) then
            tmpcourout = 1.0e3
        else
            tmpcourout = tmpb/tmp_ardz
        end if

        tmpcourmax = max( tmpcourmax, tmpcourout )
        tmp2dxa(icc,jcls) = tmpcourout
        tmp2dxb(icc,jcls) = max( tmp2dxb(icc,jcls), tmpcourout )
        end do   ! icc
        end do   ! jcls
        if (lun124 > 0) &
        write( lun124, '(i3,1p,3e12.4,2x,3e12.4)' ) k, (tmp2dxa(iccy,1:3), iccy=1,2)
    end do   ! k
    if (lun124 > 0) &
        write( lun124, '( a,1p,3e12.4,2x,3e12.4)' ) 'max', (tmp2dxb(iccy,1:3), iccy=1,2)

    if (tmpcourmax > 1.0_r8) then
        tmpa = max( 0.0_r8, tmpcourmax-1.0e-7_r8 )
        ntstep_sub = 1 + int( tmpa )
    else
        ntstep_sub = 1
    end if
    ntstep_sub_sum = ntstep_sub_sum + ntstep_sub
    dtstep_sub = dtstep_pp/ntstep_sub
    xntstep_sub_inv = 1.0_r8/ntstep_sub

    lun118 = -1
    if (idiag118_pt2 > 0) lun118 = ldiagaa_ecpp(118)
    if (lun118 > 0) then
    write(lun118,'(a,1p,2e12.4,2i12)')   &
        ' tmpcourmax, dtstep_sub, nstep_sub =', tmpcourmax, dtstep_sub,   &
        ntstep_sub, ntstep_sub_sum
    end if
    if (lun124 > 0) &
        write( lun124, '(a,1p,2e12.4,2i12)' )   &
        ' tmpcourmax, dtstep_sub, nstep_sub =', tmpcourmax, dtstep_sub,   &
        ntstep_sub, ntstep_sub_sum



!
!   do multiple integration sub-steps
!   apply vertical transport and balancing entrainment/detrainment
!
!   area change is done elsewhere, so area is fixed at ardz_cen_new 
!       during the integration loop
!
main_itstep_sub_loop:   &
    do itstep_sub = 1, ntstep_sub

         call t_startf('ecpp_vertical')

!   copy "current" chem_sub values to chem_sub_old
    chem_sub_old(:,:,:,:) = chem_sub_new(:,:,:,:)


    tmpverta(:,:,:) = 0.0
    tmphoriz(  :,:,:) = 0.0

!   calculate "transport" changes to chem_sub over one time sub-step
!   (vertical transport and horizontal exchange, including activation/resuspension)
main_trans_jcls_loop:   &
    do jcls = 1, ncls_use
main_trans_icc_loop:   &
    do icc = 1, 2
main_trans_k_loop:   &
    do k = kts, ktecen


!   if area ~= 0, then just set chem_sub_new to chem_bar
    ardz_cut = afrac_cut*rhodz_cen(k)*0.3
    if (ardz_cen_new(k,icc,jcls) < ardz_cut) then
        do l = param_first_ecpp, num_chem_ecpp
        chem_sub_new(k,icc,jcls,l) =   &
            chem_bar(k,l)*chem_bar_iccfactor(icc,l)
        end do
        cycle main_trans_k_loop
    end if


!   la loop goes over all species
!      for la = non-aerosol species, loop is executed with lc=0
!      for la = interstitial aerosol species, loop is excecuted with 
!         lc=activated counterpart
!      for la = activated aerosol species, loop is skipped
main_trans_la_loop:   &
    do la = p1st, num_chem_ecpp

            tmp_del_ardzqa_act = 0.0
            tmp_del_ardzqc_act = 0.0

        lc = 0
        l = -999888777
        not_aicw = .true.
!       if (activate_onoff_use > 999888777) then
        if (activate_onoff_use > 0) then
        if (iphase_of_aerosol(la) == ai_phase) then
            lc = laicwpair_of_aerosol(la)
            not_aicw = .false.
        else if (iphase_of_aerosol(la) == cw_phase) then
            cycle main_trans_la_loop
        end if
        end if
        if ((lc < p1st) .or. (lc > num_chem_ecpp)) lc = -999888777
        m = isize_of_aerosol(la) ; if (m <= 0) m = -999888777
        n = itype_of_aerosol(la) ; if (n <= 0) n = -999888777

        tmp_ardz   = ardz_cen_old(k,icc,jcls)
        tmp_ardzqa = chem_sub_old(k,icc,jcls,la)*tmp_ardz
        tmp_ardzqc = 0.0
        if (lc > 0)   &
        tmp_ardzqc = chem_sub_old(k,icc,jcls,lc)*tmp_ardz

!   subtract detrainment loss (no activation/resuspension here)
        tmp_del_ardz = -det_airamt_tot(icc,jcls,k)*xntstep_sub_inv
        if (tmp_del_ardz < 0.0) then
        tmp_ardz   = tmp_ardz   + tmp_del_ardz
        tmp_ardzqa = tmp_ardzqa + chem_sub_old(k,icc,jcls,la)*tmp_del_ardz
        tmphoriz(k,jcls,la) = tmphoriz(k,jcls,la)    &
                      + chem_sub_old(k,icc,jcls,la)*tmp_del_ardz
        if (lc > 0) then
        tmp_ardzqc = tmp_ardzqc + chem_sub_old(k,icc,jcls,lc)*tmp_del_ardz
        tmphoriz(k,jcls,lc) = tmphoriz(k,jcls,lc)    &
                      + chem_sub_old(k,icc,jcls,lc)*tmp_del_ardz
        end if
        end if

!   add entrainment contributions (need activation/resuspension here)
!
!+++mhwang
!   Calculate the fraction of entrainment that may expericence activations. 
!     (we assume only the new cloudy updraft may experience activation, and 
!      old updraft do not experience activation
!    Minghuai Wang, 2010-05
            frc_ent_act = 1.0_r8
            if(mtype_updnenv_use(icc, jcls) == mtype_updraft_ecpp) then
              abnd_up = 0.0_r8
              abnd_dn = 0.0_r8
              if(rhobnd_bar(k+1).gt.1.0e-10) then
                abnd_up = ar_bnd_tavg(k+1, icc, jcls)/rhobnd_bar(k+1)
              end if
              if(rhobnd_bar(k).gt.1.0e-10) then
                abnd_dn = ar_bnd_tavg(k, icc, jcls)/rhobnd_bar(k)
              end if
              if(k.eq.kts) then
                frc_ent_act = 1.0_r8
              else if(abnd_up.gt.1.0e-5) then
                frc_ent_act = 1.0_r8 - min(1.0_r8, abnd_dn/abnd_up)

                if(mfbnd_use(k+1, icc, jcls).gt.1.0e-20) then
                  frc_tmp =  max(1.0e-5_r8, 1.0_r8-mfbnd_use(k, icc, jcls)/mfbnd_use(k+1, icc, jcls))
                  frc_ent_act = min(1.0_r8, frc_ent_act / frc_tmp)
                endif
              end if
            end if    ! end mtype_updnenv_use
!---mhwang


entrain_jclsy_loop:   &
        do jclsy = 1, ncls_use
entrain_iccy_loop:   &
        do iccy = 1, 2
        tmp_del_ardz = ent_airamt(icc,jcls,iccy,jclsy,k)*xntstep_sub_inv
        if (tmp_del_ardz <= 0.0) cycle entrain_iccy_loop

        if ( not_aicw ) then
            ido_actres_tmp = 0
        else
            ido_actres_tmp = ido_actres_horz(icc,jcls,iccy,jclsy)
        end if

        tmp_qyla = chem_sub_old(k,iccy,jclsy,la)
        if (lc > 0) then
            tmp_qylc = chem_sub_old(k,iccy,jclsy,lc)
        else
            tmp_qylc = 0.0
        end if
        tmp_ardz = tmp_ardz + tmp_del_ardz

        if (activate_onoff_testaa <= 0) ido_actres_tmp = 0   ! for testing
!+++mhwangtest
! turn activation in entrainment off
!                ido_actres_tmp = 0   ! +++mhwangtest
                
        if (ido_actres_tmp == 0) then
            ! non aicw-aerosol species OR no activation or resuspension
            tmp_del_ardzqa = tmp_qyla*tmp_del_ardz
            tmp_del_ardzqc = tmp_qylc*tmp_del_ardz

        else if (ido_actres_tmp > 0) then
            ! activation of (la+lc)
            if (inmw_of_aerosol(la) == 1) then
!           tmp_fmnact = fnact_horz(m,n,jcls,iccy,jclsy,k)
                       tmp_fmnact = fnact_horz(m,n,jcls,iccy,jclsy,k) * frc_ent_act    !+++mhwang
            else
!           tmp_fmnact = fmact_horz(m,n,jcls,iccy,jclsy,k)
                        tmp_fmnact = fmact_horz(m,n,jcls,iccy,jclsy,k) * frc_ent_act    ! +++mhwang
            end if
            if (ido_actres_tmp == 2) then
            tmp_del_ardzqa = (tmp_qyla+tmp_qylc)*(1.0-tmp_fmnact)*tmp_del_ardz
            tmp_del_ardzqc = (tmp_qyla+tmp_qylc)*(tmp_fmnact    )*tmp_del_ardz
            else
            tmp_del_ardzqa = (tmp_qyla*(1.0-tmp_fmnact)     )*tmp_del_ardz
            tmp_del_ardzqc = (tmp_qyla*tmp_fmnact + tmp_qylc)*tmp_del_ardz
            end if

        else 
            ! resuspension of lc
            tmp_del_ardzqa = (tmp_qyla+tmp_qylc)*tmp_del_ardz
            tmp_del_ardzqc = 0.0

        end if

        tmp_ardzqa = tmp_ardzqa + tmp_del_ardzqa
        tmp_ardzqc = tmp_ardzqc + tmp_del_ardzqc
        tmphoriz(k,jcls,la) = tmphoriz(k,jcls,la) + tmp_del_ardzqa
        if (lc > 0)   &
        tmphoriz(k,jcls,lc) = tmphoriz(k,jcls,lc) + tmp_del_ardzqc

! change from activation/resuspension
                tmp_del_ardzqa_act = tmp_del_ardzqa_act + (tmp_del_ardzqa - tmp_qyla*tmp_del_ardz)
                if (lc > 0)    &
                tmp_del_ardzqc_act = tmp_del_ardzqc_act + (tmp_del_ardzqc - tmp_qylc*tmp_del_ardz)

        end do entrain_iccy_loop
        end do entrain_jclsy_loop


        if (jcls == jcls_qu) then
!   quiescent class -- calc change to layer k mixrat due to vertical transport at lower boundary
!   mfbnd_quiescn_up(k,icc1,icc2) is upwards mass flux from icc1 to icc2 
!       at bottom of layer k
!   mfbnd_quiescn_dn(k,icc1,icc2) is downwards ...
!   activation/resuspension calcs
!   k-1,clear  to k,cloudy - do activation
!   k-1,cloudy to k,clear  - do resuspension
!   k,either   to k-1,either - are just calculating loss to k here, so no act/res needed
vert_botqu_iupdn_loop:   &
        do iupdn = 1, 2
            if (k <= kts) cycle vert_botqu_iupdn_loop   ! skip k=kts
vert_botqu_iccy_loop:   &
        do iccy = 1, 2
            ! kb & iccy refer to the layer and sub-class from which 
            !    air and tracer mass are leaving
            ido_actres_tmp = 0
            if (iupdn == 1) then
            ! air is going from kb=k-1,iccb=iccy=1:2 to k,icc
            tmp_del_ardz = mfbnd_quiescn_up(k,iccy,icc)*dtstep_sub
            kb = k - 1
            iccb = iccy
            if (not_aicw .eqv. .false.) then
                if      ((iccy == 1) .and. (icc == 2)) then
                ido_actres_tmp = 1
                else if ((iccy == 2) .and. (icc == 1)) then
                ido_actres_tmp = -1
                end if
            end if
            else
            ! air is going from kb=k,iccb=icc to k-1,iccy=1:2
            ! since this is a loss from k, we can calc iccy=1&2 
            !    together using mfbnd_quiescn_dn(k,icc,0)
            if (iccy > 1) cycle vert_botqu_iccy_loop
            tmp_del_ardz = mfbnd_quiescn_dn(k,icc,0)*dtstep_sub
            kb = k
            iccb = icc
            end if

            if (tmp_del_ardz == 0.0) cycle vert_botqu_iccy_loop

            tmp_qyla = chem_sub_old(kb,iccb,jcls,la)
            if (lc > 0) then
            tmp_qylc = chem_sub_old(kb,iccb,jcls,lc)
            else
            tmp_qylc = 0.0
            end if

            tmp_ardz = tmp_ardz + tmp_del_ardz

            if (activate_onoff_testaa <= 0) ido_actres_tmp = 0   ! for testing
!+++mhwangtest
! turn activation in entrainment off
!                ido_actres_tmp = 0   ! +++mhwangtest
            if (ido_actres_tmp == 0) then
            ! non aicw-aerosol species OR no activation or resuspension
            tmp_del_ardzqa = tmp_qyla*tmp_del_ardz
            tmp_del_ardzqc = tmp_qylc*tmp_del_ardz

            else if (ido_actres_tmp > 0) then
            ! activation of (la+lc)
            if (inmw_of_aerosol(la) == 1) then
                tmp_fmnact = fnact_vert(m,n,k)
            else
                tmp_fmnact = fmact_vert(m,n,k)
            end if
            tmp_del_ardzqa = (tmp_qyla*(1.0-tmp_fmnact)     )*tmp_del_ardz
            tmp_del_ardzqc = (tmp_qyla*tmp_fmnact + tmp_qylc)*tmp_del_ardz

            else 
            ! resuspension of lc
            tmp_del_ardzqa = (tmp_qyla+tmp_qylc)*tmp_del_ardz
            tmp_del_ardzqc = 0.0

            end if

            tmp_ardzqa = tmp_ardzqa + tmp_del_ardzqa
            tmp_ardzqc = tmp_ardzqc + tmp_del_ardzqc
            if (icc == 1) then
            tmpverta(k,jcls,la) = tmpverta(k,jcls,la) + tmp_del_ardzqa
            if (lc > 0)   &
            tmpverta(k,jcls,lc) = tmpverta(k,jcls,lc) + tmp_del_ardzqc
            end if

! change from activation/resuspension
                   tmp_del_ardzqa_act = tmp_del_ardzqa_act + (tmp_del_ardzqa - tmp_qyla*tmp_del_ardz)
                   if (lc > 0)    &
                   tmp_del_ardzqc_act = tmp_del_ardzqc_act + (tmp_del_ardzqc - tmp_qylc*tmp_del_ardz)

            ! with "pgf90 -O2", code seg-faulted until following statement 
            !    was added.  (note that it is do-nothing, since la>0 always)
            if (la < 0) write(*,*)   &
            'vert_botqu gggg - icc,iupdn,ido', iccy, iupdn, ido_actres_tmp
        end do vert_botqu_iccy_loop
        end do vert_botqu_iupdn_loop

!   quiescent class -- calc change to layer k mixrat due to vertical transport at upper boundary
!   mfbnd_quiescn_up(k+1,icc1,icc2) is upwards mass flux from icc1 to icc2 
!       at top of layer k
!   mfbnd_quiescn_dn(k+1,icc1,icc2) is downwards ...
!   activation/resuspension calcs
!   k+1,clear  to k,cloudy - downwards motion so skip activation ???
!   k+1,cloudy to k,clear  - do resuspension
!   k,either   to k+1,either - are just calculating loss to k here, so no act/res needed
vert_topqu_iupdn_loop:   &
        do iupdn = 1, 2
            if (k >= ktebnd-1) cycle vert_topqu_iupdn_loop   ! skip k=ktebnd-1,ktebnd
vert_topqu_iccy_loop:   &
        do iccy = 1, 2
            ido_actres_tmp = 0
            if (iupdn == 1) then
            ! air is going from kb=k,iccb=icc to k+1,iccy=1:2
            ! since this is a loss from k, we can calc iccy=1&2 
            !    together using mfbnd_quiescn_up(k+1,icc,0)
            if (iccy > 1) cycle vert_topqu_iccy_loop
            tmp_del_ardz = -mfbnd_quiescn_up(k+1,icc,0)*dtstep_sub
            kb = k
            iccb = icc
            else
            ! air is going from kb=k+1,iccb=iccy=1:2 to k,icc
            tmp_del_ardz = -mfbnd_quiescn_dn(k+1,iccy,icc)*dtstep_sub
            kb = k+1
            iccb = iccy
            if (not_aicw .eqv. .false.) then
                if ((iccy == 2) .and. (icc == 1)) then
                ido_actres_tmp = -1
                end if
            end if
            end if

            if (tmp_del_ardz == 0.0) cycle vert_topqu_iccy_loop

            tmp_qyla = chem_sub_old(kb,iccb,jcls,la)
            if (lc > 0) then
            tmp_qylc = chem_sub_old(kb,iccb,jcls,lc)
            else
            tmp_qylc = 0.0
            end if

            tmp_ardz = tmp_ardz + tmp_del_ardz

            if (activate_onoff_testaa <= 0) ido_actres_tmp = 0   ! for testing
!+++mhwangtest
! turn activation in entrainment off
!                ido_actres_tmp = 0   ! +++mhwangtest
            if (ido_actres_tmp == 0) then
            ! non aicw-aerosol species OR no activation or resuspension
            tmp_del_ardzqa = tmp_qyla*tmp_del_ardz
            tmp_del_ardzqc = tmp_qylc*tmp_del_ardz

            else if (ido_actres_tmp > 0) then
            ! activation of (la+lc)
            if (inmw_of_aerosol(la) == 1) then
                tmp_fmnact = fnact_vert(m,n,k)
            else
                tmp_fmnact = fmact_vert(m,n,k)
            end if
            tmp_del_ardzqa = (tmp_qyla*(1.0-tmp_fmnact)     )*tmp_del_ardz
            tmp_del_ardzqc = (tmp_qyla*tmp_fmnact + tmp_qylc)*tmp_del_ardz

            else 
            ! resuspension of lc
            tmp_del_ardzqa = (tmp_qyla+tmp_qylc)*tmp_del_ardz
            tmp_del_ardzqc = 0.0

            end if

            tmp_ardzqa = tmp_ardzqa + tmp_del_ardzqa
            tmp_ardzqc = tmp_ardzqc + tmp_del_ardzqc

! change from activation/resuspension
                    tmp_del_ardzqa_act = tmp_del_ardzqa_act + (tmp_del_ardzqa - tmp_qyla*tmp_del_ardz)
                    if (lc > 0)    &
                    tmp_del_ardzqc_act = tmp_del_ardzqc_act + (tmp_del_ardzqc - tmp_qylc*tmp_del_ardz)

            ! with "pgf90 -O2", code seg-faulted until following statement 
            !    was added.  (note that it is do-nothing, since la>0 always)
            if (la < 0) write(*,*)   &
            'vert_topqu gggg - icc,iupdn,ido', iccy, iupdn, ido_actres_tmp
        end do vert_topqu_iccy_loop
        end do vert_topqu_iupdn_loop


        else
!   up/dndraft class -- add/subtract vertical transport at lower boundary 
!   no activation/resuspension here as the vertical transport within up/dndrafts
!      is clear-->clear or cloudy-->cloudy.   (The within up/dndraft
!      clear<-->cloudy is done by ent/detrainment.)
        if (k > kts) then
        tmp_del_ardz = mfbnd_use(k,icc,jcls)*dtstep_sub
        if (abs(tmp_del_ardz) > 0.0) then
            if (tmp_del_ardz > 0.0) then
            kb = k - 1
            else
            kb = k
            end if
            tmp_ardz   = tmp_ardz   + tmp_del_ardz
            tmp_ardzqa = tmp_ardzqa + chem_sub_old(kb,icc,jcls,la)*tmp_del_ardz
            if (lc > 0)   &
            tmp_ardzqc = tmp_ardzqc + chem_sub_old(kb,icc,jcls,lc)*tmp_del_ardz
            if (icc == 1) then
            tmpverta(k,jcls,la) = chem_sub_old(kb,icc,jcls,la)*tmp_del_ardz
            if (lc > 0)   &
            tmpverta(k,jcls,lc) = chem_sub_old(kb,icc,jcls,lc)*tmp_del_ardz
            end if
        end if
        end if   ! (k > kts)
!   up/dndraft class -- add/subtract vertical transport at upper boundary 
        if (k < ktebnd-1) then
        tmp_del_ardz = -mfbnd_use(k+1,icc,jcls)*dtstep_sub
        if (abs(tmp_del_ardz) > 0.0) then
            if (tmp_del_ardz > 0.0) then
            kb = k + 1
            else
            kb = k
            end if
            tmp_ardz   = tmp_ardz   + tmp_del_ardz
            tmp_ardzqa = tmp_ardzqa + chem_sub_old(kb,icc,jcls,la)*tmp_del_ardz
            if (lc > 0)   &
            tmp_ardzqc = tmp_ardzqc + chem_sub_old(kb,icc,jcls,lc)*tmp_del_ardz
        end if
        end if   ! (k < ktebnd-1)

        end if   ! (jcls == jcls_qu)


!   new mixing ratio
        chem_sub_new(k,icc,jcls,la) = tmp_ardzqa/ardz_cen_new(k,icc,jcls)
        if (lc > 0)   &
        chem_sub_new(k,icc,jcls,lc) = tmp_ardzqc/ardz_cen_new(k,icc,jcls)

!   change in mixing ratio (*fraction) from activation/resuspension
            del_activate3d(k,icc,jcls,la) =  del_activate3d(k,icc,jcls,la)+tmp_del_ardzqa_act/rhodz_cen(k)
            if (lc > 0)   &
            del_activate3d(k,icc,jcls,lc) =  del_activate3d(k,icc,jcls,lc)+tmp_del_ardzqc_act/rhodz_cen(k)

    end do main_trans_la_loop

    end do main_trans_k_loop
    end do main_trans_icc_loop
    end do main_trans_jcls_loop


!   fort.118 diagnostics
    lun118 = -1
    if (idiag118_pt3 > 0) then
        if (idiag118_pt3 >= 10) lun118 = ldiagaa_ecpp(118)
        if (itstep_sub == ntstep_sub) lun118 = ldiagaa_ecpp(118)
    end if
    if (lun118 > 0) then
    do l = param_first_ecpp, num_chem_ecpp
    if ((l == 9) .or. (l == 9)) then

    write(lun118,'(/a,3i5)') 'new_main_integ pt3 ktau_pp, istep_sub, l =', ktau_pp, itstep_sub, l
    write(lun118,'(2a)')   &
        '(chem_sub_old(k,icc,jcls,l), chem_sub_new(k,icc,jcls,l), jcls=2,1,-1); ',   &
        'updr ardz_cen_new and w; dumverta/b, dumhoriz for updr then env'

    icc = 1
    tmpc = 1.0_r8/dtstep_sub
    do k = ktecen, kts, -1
    tmpa = 0.0
    if (ar_bnd_tavg(k,icc,jgrp_up) > 0.0)   &
        tmpa = mfbnd_use(k,icc,jgrp_up)/ar_bnd_tavg(k,icc,jgrp_up)
    write(lun118,'(i3,1p,3(1x,2e10.3),2(1x,3e10.3))') k,   &
        ( chem_sub_old(k,icc,jcls,l), chem_sub_new(k,icc,jcls,l), jcls=2,1,-1 ),   &
        ardz_cen_new(k,icc,jgrp_up), tmpa,   &
        ( tmpverta(k,jcls,l)*tmpc,   &
         (tmpverta(k,jcls,l)-tmpverta(k+1,jcls,l))*tmpc,   &
          tmphoriz(k,jcls,l)*tmpc, jcls=2,1,-1 )
    end do   ! k

    end if   ! (l == ...)
    end do   ! l
    end if   ! (lun118 > 0)

        call t_stopf('ecpp_vertical')

        end do main_itstep_sub_loop

!
! +++mhwang
! move cloud chemistry and wetscavenging outside of istep_sub_loop
! inside of the itstep_sub_loop is too expanseive
! Minghuai Wang, 2010-04-28
! 
        itstep_sub = 1
        dtstep_sub = dtstep_pp

!   calculate cloud chemistry changes to chem_sub over one time sub-step
!        call t_startf('ecpp_cldchem')
!   call parampollu_tdx_cldchem(                     &
!       ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
!                itstep_hybrid,                                     &
!       idiagaa_ecpp, ldiagaa_ecpp,                        &
!       tcen_bar, pcen_bar, rhocen_bar, dzcen,             &
!       rhobnd_bar, zbnd, wbnd_bar,                        &
!       chem_bar,                                          &
!       ncls_ecpp,                                         &
!       it,      jt,      kts,ktebnd,ktecen,               &
!       ncls_use,                                          &
!       kdraft_bot_use, kdraft_top_use,                    &
!       mtype_updnenv_use,                                 &
!       chem_sub_new,                                      &
!       del_chem_clm_cldchem, del_chem_clm_rename, del_cldchem3d, del_rename3d, &
!                aqso4_h2o2, aqso4_o3, xphlwc3d,                    &
!       ardz_cen_old, ardz_cen_new, rhodz_cen,             &
!       acen_tavg_use, acen_prec_use,                      &
!       rh_sub2, qcloud_sub2, qlsink_sub2,                 &
!       precr_sub2, precs_sub2,                            &
!       chem_bar_iccfactor, activate_onoff_use,            &
!       iphase_of_aerosol, isize_of_aerosol,               &
!       itype_of_aerosol, inmw_of_aerosol,                 &
!       laicwpair_of_aerosol                               )
!         call t_stopf('ecpp_cldchem')


!   calculate wet removal changes to chem_sub over one time sub-step

        if (wetscav_onoff_ecpp >= 100) then
        call t_startf('ecpp_wetscav')
!        write(*,'(a,3i8)') 'main integ calling wetscav_2', ktau, ktau_pp, itstep_sub
        call parampollu_tdx_wetscav_2(                             &
                ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub,     &
                itstep_hybrid,                                     &
                idiagaa_ecpp, ldiagaa_ecpp, idiagbb_wetscav,       &
                tcen_bar, pcen_bar, rhocen_bar, dzcen,             &
!               rhobnd_bar, zbnd, wbnd_bar,                        &   not needed ?
!               chem_bar,                                          &   not needed ?
!               ncls_ecpp,                                         &
                it,      jt,      kts,ktebnd,ktecen,               &
                ncls_use,                                          &
!               kdraft_bot_use, kdraft_top_use,                    &   not needed ?
!               mtype_updnenv_use,                                 &   not needed ?
                chem_sub_new,                                      &
                del_chem_clm_wetscav,                              &
                del_wetscav3d, del_wetresu3d,                      &
!               ardz_cen_old, ardz_cen_new,                        &   not needed ?
                rhodz_cen,             &
                acen_tavg_use, acen_prec_use,                      &
                rh_sub2, qcloud_sub2, qlsink_sub2,                 &
                precr_sub2, precs_sub2,                            &
!               chem_bar_iccfactor,                                &   not needed ?
                activate_onoff_use,                                &
                iphase_of_aerosol, isize_of_aerosol,               &
                itype_of_aerosol, inmw_of_aerosol,                 &
                laicwpair_of_aerosol                               )
!        write(*,'(a,3i8)') 'main integ backfrm wetscav_2', ktau, ktau_pp, itstep_sub
        call t_stopf('ecpp_wetscav')
        end if ! (wetscav_onoff_ecpp >= 100)

!   calculate cloud chemistry changes to chem_sub over one time sub-step
        call t_startf('ecpp_cldchem')
!==Guangxing Lin==test
       if(1.eq.1) then        
        call parampollu_tdx_cldchem(                     &
                ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
                itstep_hybrid,                                     &
                idiagaa_ecpp, ldiagaa_ecpp,                        &
                tcen_bar, pcen_bar, rhocen_bar, dzcen,             &
                rhobnd_bar, zbnd, wbnd_bar,                        &
                chem_bar,                                          &
                ncls_ecpp,                                         &
                it,      jt,      kts,ktebnd,ktecen,               &
                ncls_use,                                          &
                kdraft_bot_use, kdraft_top_use,                    &
                mtype_updnenv_use,                                 &
                chem_sub_new,                                      &
                del_chem_clm_cldchem, del_chem_clm_rename, del_cldchem3d, del_rename3d, &
                aqso4_h2o2, aqso4_o3, xphlwc3d,                    &
                ardz_cen_old, ardz_cen_new, rhodz_cen,             &
                acen_tavg_use, acen_prec_use,                      &
                rh_sub2, qcloud_sub2, qlsink_sub2,                 &
                precr_sub2, precs_sub2,                            &
                chem_bar_iccfactor, activate_onoff_use,            &
                iphase_of_aerosol, isize_of_aerosol,               &
                itype_of_aerosol, inmw_of_aerosol,                 &
                laicwpair_of_aerosol, pbuf                         )
       end if  
       call t_stopf('ecpp_cldchem')

!   end do main_itstep_sub_loop

        !call  t_stopf('ecpp_mainintegr') !==Guangxing Lin


    return
    end subroutine parampollu_tdx_main_integ


!-----------------------------------------------------------------------
    subroutine parampollu_tdx_area_change(                    &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, wbnd_bar,                             &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use, ipass_area_change,                      &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        chem_sub_new,                                     &
                del_activate3d,                                   &
        mfbnd_use, ar_bnd_tavg,                           &
        ardz_cen_old, ardz_cen_new, rhodz_cen,            &
        chem_bar_iccfactor, activate_onoff_use,           &
        iphase_of_aerosol, isize_of_aerosol,              &
        itype_of_aerosol, inmw_of_aerosol,                &
        laicwpair_of_aerosol                              )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_area_change does
!    horizontal exchange associated with sub-class area changes
!
! incoming chem_sub_new holds current sub-class mixing ratios
! outgoing chem_sub_new holds updated sub-class mixing ratios
!
!-----------------------------------------------------------------------

    use module_data_radm2, only:  epsilc

    use module_data_mosaic_asect, only:  ai_phase, cw_phase,   &
                maxd_asize, maxd_atype   

    use module_data_ecpp1

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message   

!   arguments
    integer, intent(in) ::                  &
        ktau, ktau_pp,              &
        it, jt, kts, ktebnd, ktecen
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!   chem_driver and routines under it do calculations
!   over these spatial indices.

    integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

    real(r8), intent(in) :: dtstep, dtstep_pp
!   dtstep - main model time step (s)
!   dtstep_pp - time step (s) for "parameterized pollutants" calculations

    real(r8), intent(in), dimension( kts:ktecen ) ::   &
        tcen_bar, pcen_bar, rhocen_bar, dzcen
    real(r8), intent(in), dimension( kts:ktebnd ) ::   &
        rhobnd_bar, wbnd_bar
!   tcen_bar - temperature (K) at layer centers
!   rhocen_bar, rhobnd_bar - dry air density (kg/m^3) at layer centers and boundaries
!   pcen_bar - air pressure (Pa) at layer centers
!   wbnd_bar - vertical velocity (m/s) at layer boundaries
!   dzcen - layer thicknesses (m)

    real(r8), intent(in), dimension( kts:ktecen, 1:num_chem_ecpp ) :: &
        chem_bar
!   chem_bar - mixing ratios of trace gase (ppm) and aerosol species
!   (ug/kg for mass species, #/kg for number species)

    integer, intent(in) :: ncls_ecpp
!   ncls_ecpp - number of ecpp transport classes in the grid column
    
    integer, intent(inout) :: ipass_area_change
    integer, intent(in)    :: ncls_use

    integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        kdraft_bot_use, kdraft_top_use,   &
        mtype_updnenv_use

    real(r8), intent(inout), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
        chem_sub_new

        real(r8), intent(inout), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
                del_activate3d

    real(r8), intent(inout), dimension( kts:ktebnd, 0:2, 0:maxcls_ecpp ) :: &
        mfbnd_use, ar_bnd_tavg

    real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) :: &
        ardz_cen_old, ardz_cen_new

    real(r8), intent(inout), dimension( kts:ktecen ) :: rhodz_cen

    real(r8), intent(in), dimension( 1:2, num_chem_ecpp ) :: chem_bar_iccfactor

    integer, intent(in) :: activate_onoff_use

    integer, intent(in), dimension( 1:num_chem_ecpp ) ::   &
        iphase_of_aerosol, isize_of_aerosol, itype_of_aerosol,   &
        inmw_of_aerosol, laicwpair_of_aerosol


!   local variables
    integer :: icc, iccy, ido_actres_tmp, ifrom_where, itmpa
    integer :: idiag118_pt3
    integer :: jcls, jclsy
    integer :: k
    integer :: l, la, laa, lbb, lc, lun118
    integer :: m, n
    integer :: p1st

    integer, dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp ) ::   &
        ido_actres_horz

    logical :: not_aicw

    real(r8) :: ardz_cut
    real(r8) :: tmpa, tmpb, tmpc, tmpd
    real(r8) :: tmp_fmnact, tmp_qyla, tmp_qylc
    real(r8) :: tmpvecd(0:maxcls_ecpp), tmpvece(0:maxcls_ecpp)
        real(r8) :: tmp_del_ardzqa, tmp_del_ardzqc
        real(r8) :: tmp_del_ardzqa_act, tmp_del_ardzqc_act

    real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
        chem_sub_old
    real(r8), dimension( 1:2, 1:maxcls_ecpp, kts:ktecen ) :: &
        ent_airamt_tot, det_airamt_tot
    real(r8), dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
        ent_airamt, det_airamt
    real(r8), dimension( 1:maxd_asize, 1:maxd_atype, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
        fmact_horz, fnact_horz



    p1st = param_first_ecpp
    idiag118_pt3 = 10 * mod( max(idiagaa_ecpp(118),0)/100, 10 )

!
!   calc entrain/detrain amounts
!
!   first calc net (entrainment-detrainment) amount = area change
    ent_airamt_tot(:,:,:) = 0.0
    det_airamt_tot(:,:,:) = 0.0
    do jcls = 1, ncls_use
    do icc = 1, 2
    do k = kts, ktecen
        ardz_cut = afrac_cut*rhodz_cen(k)*0.3
        tmpa = max( ardz_cen_new(k,icc,jcls), ardz_cen_old(k,icc,jcls) )
        if (tmpa >= ardz_cut) then
        tmpb = ardz_cen_new(k,icc,jcls) - ardz_cen_old(k,icc,jcls)
        if (tmpb > 0.0) then
            ent_airamt_tot(icc,jcls,k) = tmpb
        else if (tmpb < 0.0) then
            det_airamt_tot(icc,jcls,k) = -tmpb
        end if
        end if
    end do
    end do
    end do

!   next calc detailed ent/det amounts
    ifrom_where = ipass_area_change
    call parampollu_tdx_entdet_sub1(                          &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use, ifrom_where,                            &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        ardz_cen_old, ardz_cen_new, rhodz_cen,            &
        ent_airamt_tot, det_airamt_tot,                   &
        ent_airamt, det_airamt                            )


!
!   calc activation/resuspension fractions associated with ent/det
!   
    if (activate_onoff_use > 0) then
    ifrom_where = ipass_area_change
    call parampollu_tdx_activate1(                            &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, wbnd_bar,                             &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use, ifrom_where, activate_onoff_use,        &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        chem_sub_new,                                     &
        mfbnd_use,                                        &
        ar_bnd_tavg,                                      &
        ent_airamt,                                       &
        ido_actres_horz, fmact_horz, fnact_horz           )
    end if


!   copy chem_sub_new (= incoming current chem_sub values) into chem_sub_old
    chem_sub_old(:,:,:,:) = chem_sub_new(:,:,:,:)

!   calculate new chem_sub
main_jcls_loop:   &
    do jcls = 1, ncls_use
main_icc_loop:   &
    do icc = 1, 2
main_k_loop:   &
    do k = kts, ktecen

!   if entrainment and detrainment) both ~= 0, then no change
    if ( (ent_airamt_tot(icc,jcls,k) < 1.0e-30_r8) .and.   &
         (det_airamt_tot(icc,jcls,k) < 1.0e-30_r8) ) cycle

!   if new area ~= 0, then just set chem_sub_new to chem_bar
    ardz_cut = afrac_cut*rhodz_cen(k)*0.3
    if (ardz_cen_new(k,icc,jcls) < ardz_cut) then
        do l = p1st, num_chem_ecpp
        chem_sub_new(k,icc,jcls,l) =   &
            chem_bar(k,l)*chem_bar_iccfactor(icc,l)
        end do
        cycle main_k_loop
    end if

!   la loop goes over all species
!      for la = non-aerosol species, loop is executed with lc=0
!      for la = interstitial aerosol species, loop is excecuted with 
!         lc=activated counterpart
!      for la = activated aerosol species, loop is skipped
main_la_loop:   &
    do la = p1st, num_chem_ecpp

            tmp_del_ardzqa_act = 0.0
            tmp_del_ardzqc_act = 0.0

        lc = 0
        not_aicw = .true.
        if (activate_onoff_use > 0) then
        if (iphase_of_aerosol(la) == ai_phase) then
            lc = laicwpair_of_aerosol(la)
            not_aicw = .false.
        else if (iphase_of_aerosol(la) == cw_phase) then
            cycle main_la_loop
        end if
        end if
        if ((lc < p1st) .or. (lc > num_chem_ecpp)) lc = -999888777

!   tmpd = (original area) - (detrainment to all others)
        tmpd = ardz_cen_old(k,icc,jcls) - det_airamt_tot(icc,jcls,k)
        tmpd = max( tmpd, 0.0_r8 )

!   tmpa holds sum_of( mix_ratio * area ) for interstitial
!   tmpc holds sum_of( mix_ratio * area ) for activated
        tmpa = chem_sub_old(k,icc,jcls,la)*tmpd
        if (lc > 0)   &
        tmpc = chem_sub_old(k,icc,jcls,lc)*tmpd

!   add entrainment contributions
        do jclsy = 1, ncls_use
        do iccy = 1, 2
        tmpd = ent_airamt(icc,jcls,iccy,jclsy,k)
        if (tmpd <= 0.0) cycle

        if ( not_aicw ) then
            ido_actres_tmp = 0
        else
            ido_actres_tmp = ido_actres_horz(icc,jcls,iccy,jclsy)
        end if

        tmp_qyla = chem_sub_old(k,iccy,jclsy,la)
        if (lc > 0) then
            tmp_qylc = chem_sub_old(k,iccy,jclsy,lc)
        else
            tmp_qylc = 0.0
        end if

        if (ido_actres_tmp == 0) then
            ! non aicw-aerosol species OR no activation or resuspension
!           tmpa = tmpa + tmp_qyla*tmpd
!           tmpc = tmpc + tmp_qylc*tmpd
                   tmp_del_ardzqa = tmp_qyla*tmpd
                   tmp_del_ardzqc = tmp_qylc*tmpd

        else if (ido_actres_tmp > 0) then
            ! activation of (la+lc)
            m = isize_of_aerosol(la)
            n = itype_of_aerosol(la)
            if (inmw_of_aerosol(la) == 1) then
            tmp_fmnact = fnact_horz(m,n,jcls,iccy,jclsy,k)
            else
            tmp_fmnact = fmact_horz(m,n,jcls,iccy,jclsy,k)
            end if
            if (ido_actres_tmp == 2) then
!           tmpa = tmpa + (tmp_qyla+tmp_qylc)*(1.0-tmp_fmnact)*tmpd
!           tmpc = tmpc + (tmp_qyla+tmp_qylc)*(tmp_fmnact    )*tmpd
                       tmp_del_ardzqa = (tmp_qyla+tmp_qylc)*(1.0-tmp_fmnact)*tmpd
                       tmp_del_ardzqc = (tmp_qyla+tmp_qylc)*(tmp_fmnact    )*tmpd
            else
!           tmpa = tmpa + (tmp_qyla*(1.0-tmp_fmnact)     )*tmpd
!           tmpc = tmpc + (tmp_qyla*tmp_fmnact + tmp_qylc)*tmpd
                       tmp_del_ardzqa = (tmp_qyla*(1.0-tmp_fmnact)     )*tmpd 
                       tmp_del_ardzqc = (tmp_qyla*tmp_fmnact + tmp_qylc)*tmpd
            end if

        else 
            ! resuspension of lc
!           tmpa = tmpa + (tmp_qyla+tmp_qylc)*tmpd
                   tmp_del_ardzqa = (tmp_qyla+tmp_qylc)*tmpd
                   tmp_del_ardzqc = 0.0 

        end if
                tmpa = tmpa + tmp_del_ardzqa
                if (lc > 0)   &
                tmpc = tmpc + tmp_del_ardzqc

! change from activation/resuspension
                tmp_del_ardzqa_act = tmp_del_ardzqa_act + (tmp_del_ardzqa - tmp_qyla*tmpd)
                if (lc > 0)    &
                tmp_del_ardzqc_act = tmp_del_ardzqc_act + (tmp_del_ardzqc - tmp_qylc*tmpd)
        end do   ! iccy
        end do   ! jclsy
        chem_sub_new(k,icc,jcls,la) = tmpa/ardz_cen_new(k,icc,jcls)
        if (lc > 0)   &
        chem_sub_new(k,icc,jcls,lc) = tmpc/ardz_cen_new(k,icc,jcls)

!   change in mixing ratio (*fraction) from activation/resuspension
            del_activate3d(k,icc,jcls,la) =  del_activate3d(k,icc,jcls,la)+tmp_del_ardzqa_act/rhodz_cen(k)
            if (lc > 0)   &
            del_activate3d(k,icc,jcls,lc) =  del_activate3d(k,icc,jcls,lc)+tmp_del_ardzqc_act/rhodz_cen(k)

    end do main_la_loop

    end do main_k_loop
    end do main_icc_loop
    end do main_jcls_loop


!   diagnostics
    lun118 = -1
    if (idiag118_pt3 >= 10) lun118 = ldiagaa_ecpp(118)
    if (lun118 > 0) then
    l = 9
    icc = 1
    write(lun118,'(/a,2i5,a,3i5)') 'pt3 ppopt, ipass', parampollu_opt,   &
        ipass_area_change, '    ktau_pp, istep_sub, l =', ktau_pp, -1, l
    write(lun118,'(2a)') '(chem_sub_old(k,icc,jcls,l), ',   &
        'chem_sub_new(k,icc,jcls,l), jcls=1,3); up,dn,env a_cen_tmpa/tmpb'
    do k = ktecen, kts, -1
        write(lun118,'(i3,1p,7(1x,2e10.3))') k,   &
            (chem_sub_old(k,icc,jcls,l), chem_sub_new(k,icc,jcls,l), jcls=1,3),   &
            (ardz_cen_old(k,icc,jcls)/rhodz_cen(k), ardz_cen_new(k,icc,jcls)/rhodz_cen(k), jcls=1,3)
    end do
    end if   ! (lun118 > 0)



    return
    end subroutine parampollu_tdx_area_change



!-----------------------------------------------------------------------
    subroutine parampollu_tdx_entdet_sub1(                    &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        ncls_ecpp,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use, ifrom_where,                            &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        ardz_cen_old, ardz_cen_new, rhodz_cen,            &
        ent_airamt_tot, det_airamt_tot,                   &
        ent_airamt, det_airamt                            )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_entdet_sub1 calculates 
!    the "horizontal exchange coefficients" associated with
!    area changes or vertical mass fluxes
!
! the net (entrainment-detrainment) for each sub-class is
!    obtained trivially
! determining where the entrainment comes from, and where
!    the detrainment goes to, is much more involved
!
!-----------------------------------------------------------------------

    use module_data_radm2, only:  epsilc

    use module_data_ecpp1

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message

!   arguments
    integer, intent(in) ::                  &
        ktau, ktau_pp,              &
        it, jt, kts, ktebnd, ktecen

    integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

    real(r8), intent(in) :: dtstep, dtstep_pp
!   dtstep - main model time step (s)
!   dtstep_pp - time step (s) for "parameterized pollutants" calculations

    integer, intent(in) :: ncls_ecpp
!   ncls_ecpp - number of ecpp transport classes in the grid column
    
    integer, intent(in)    :: ifrom_where
    integer, intent(in)    :: ncls_use

    integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        kdraft_bot_use, kdraft_top_use,   &
        mtype_updnenv_use

    real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) :: &
        ardz_cen_old, ardz_cen_new

    real(r8), intent(inout), dimension( kts:ktecen ) :: rhodz_cen

    real(r8), intent(inout), dimension( 1:2, 1:maxcls_ecpp, kts:ktecen  ) :: &
        ent_airamt_tot, det_airamt_tot
!   ent_airamt_tot(icc,jcls,k) is the total detrainment into layer k,
!      sub-class icc, class jcls from all other sub-classes 
!   det_airamt_tot(icc,jcls,k) is the total detrainment from layer k,
!      sub-class icc, class jcls to all other sub-classes 
!   units are (kg/m2)
!
!   define entdet_net == ent_airamt_tot - det_airamt_tot
!   for "area-change" ent/det, entdet_net = rho*dz*d(area) where
!      d(area) is the fractional area change over the time-step
!   for "vertical-transport" ent/det, entdet_net = d(mfbnd)*dtstep where
!      d(mfbnd) is the change in vertical mass flux across a layer
!      (mfbnd at layer top minus mfbnd at layer bottom)
!
!   up and dndrafts
!      in the current formulation, each draft either entrains or detrains
!         at a given level, but not both simultaneously
!      for incoming ent/det_airamt_tot, one will be >= 0 and the other will be =0
!      the outgoing ent/det_airamt_tot will be unchanged
!   quiescent class
!      the quiescent class can entrain and detrain simultaneously at a given level
!      for incoming ent/det_airamt_tot, one will be >= 0 and will hold the
!         net (entrainment-detrainment)
!      the outgoing ent/det_airamt_tot can both be >0
!

    real(r8), intent(out),   &
        dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
        ent_airamt, det_airamt
!   ent_airamt(iccaa,jclsaa,iccbb,jclsbb,k) is (positive) the entrainment amount
!      into sub-class (iccaa,jclsaa,k) from sub-class (iccbb,jclsbb,k)
!   det_airamt(iccaa,jclsaa,iccbb,jclsbb,k) is (positive) the detrainment amount
!      from sub-class (iccaa,jclsaa,k) into sub-class (iccbb,jclsbb,k) 
!   units for both are (kg/m2)


!   local variables
    integer :: icc, iccy, itmpa
    integer :: jcls, jclsy
    integer :: jgrp, jgrpy, jgrp_of_jcls(1:maxcls_ecpp)
    integer :: k
    integer :: l, laa, lbb, lunaa, lunbb
    integer :: m

    logical, dimension( 1:2, 1:maxcls_ecpp ) :: &
        empty_old, empty_new, empty_oldnew

    real(r8) :: tmpa4, tmpb4

    real(r8), dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
        ent_airamt_sv1, det_airamt_sv1
    real(r8), dimension( 1:2, 1:maxcls_ecpp, kts:ktecen  ) :: &
        ent_airamt_tot_sv0, det_airamt_tot_sv0, &
        ent_airamt_tot_sv1, det_airamt_tot_sv1
    real(r8), dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp ) ::   &
        ecls_aa, dcls_aa
    real(r8), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        ecls_aaunasi, dcls_aaunasi
    real(r8), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        dcls_aalimit
    real(r8), dimension( 1:2, 1:3, 1:2, 1:3 ) ::   &
        egrp_aa, dgrp_aa
    real(r8), dimension( 1:2, 1:3 ) ::   &
        egrp_aaunasi, dgrp_aaunasi

    real(r8), dimension( 1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
        ecls_aaunasi_sv2, dcls_aaunasi_sv2
    real(r8), dimension( 1:2, 1:3, kts:ktecen ) ::   &
        egrp_aaunasi_sv2, dgrp_aaunasi_sv2

    integer, dimension(3), save ::   &
        ecls_aaunasi_worst_i=0, dcls_aaunasi_worst_i=0,   &
        ecls_aaunasi_worst_j=0, dcls_aaunasi_worst_j=0,   &
        ecls_aaunasi_worst_k=0, dcls_aaunasi_worst_k=0,   &
        ecls_aaunasi_worst_ktau=0, dcls_aaunasi_worst_ktau=0,   &
        egrp_aaunasi_worst_i=0, dgrp_aaunasi_worst_i=0,   &
        egrp_aaunasi_worst_j=0, dgrp_aaunasi_worst_j=0,   &
        egrp_aaunasi_worst_k=0, dgrp_aaunasi_worst_k=0,   &
        egrp_aaunasi_worst_ktau=0, dgrp_aaunasi_worst_ktau=0 
    real(r8), dimension(3), save ::   &
        ecls_aaunasi_worst=0.0, dcls_aaunasi_worst=0.0,   &
        egrp_aaunasi_worst=0.0, dgrp_aaunasi_worst=0.0 

    real(r8) :: ardz_cut
    real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf
    real(r8) :: tmpmatbb(0:2,0:2)
    real(r8) :: tmpmatff(1:2,1:2)
    real(r8) :: tmpvecbb(0:maxcls_ecpp), tmpvecgg(0:maxcls_ecpp)


!   diagnostics to fort.122 at selected timesteps
    lunaa = -1
!   if ( (ktau <=  10) .or.   &
!        (ktau == 581) .or.   &
!        (ktau == 818) ) lunaa = 122
    if ( (ktau <=  10) .or.   &
         (ktau == 210) .or.   &
         (ktau == 682) ) lunaa = ldiagaa_ecpp(122)
    if (idiagaa_ecpp(122) <= 0) lunaa = -1

!   save the incoming values of ent/det_airamt_tot
    ent_airamt_tot_sv0(:,:,:) = ent_airamt_tot(:,:,:)
    det_airamt_tot_sv0(:,:,:) = det_airamt_tot(:,:,:)

!
!   do a very simple calculation that mimics previous code
!   up and dndrafts entrain-from and detrain-too
!      the quiescent class with the same icc
!   (currently the simple calculation results are only used for diagnostic 
!      purposes, but turning them off would mess up the diagnostics.)
!
    ent_airamt(:,:,:,:,:) = 0.0
    det_airamt(:,:,:,:,:) = 0.0

entdet_main_kloop_bb: &
    do k = kts, ktecen

    do jcls = 1, ncls_use
    do icc = 1, 2
        if (jcls == jcls_qu) cycle   ! skip quiescent

        tmpa4 = ent_airamt_tot(icc,jcls,k) 
        if (tmpa4 > 0.0) then
        ent_airamt( icc,jcls,    icc,jcls_qu, k) = tmpa4
        det_airamt( icc,jcls_qu, icc,jcls,    k) = tmpa4
        end if

        tmpa4 = det_airamt_tot(icc,jcls,k) 
        if (tmpa4 > 0.0) then
        det_airamt( icc,jcls,    icc,jcls_qu, k) = tmpa4
        ent_airamt( icc,jcls_qu, icc,jcls,    k) = tmpa4
        end if

    end do   ! icc
    end do   ! jcls

    end do entdet_main_kloop_bb


    do k = kts, ktecen

    do jcls = 1, ncls_use
    do icc = 1, 2
        tmpa4 = 0.0
        tmpb4 = 0.0
        if (k < ktebnd) then
        do jclsy = 1, ncls_use
        do iccy = 1, 2
            tmpa4 = tmpa4 + ent_airamt( icc,jcls, iccy,jclsy, k)
            tmpb4 = tmpb4 + det_airamt( icc,jcls, iccy,jclsy, k)
        end do
        end do
        end if
        ent_airamt_tot(icc,jcls,k) = tmpa4
        det_airamt_tot(icc,jcls,k) = tmpb4
    end do   ! icc
    end do   ! jcls

    end do   ! k

    ent_airamt_sv1(:,:,:,:,:) = ent_airamt(:,:,:,:,:)
    det_airamt_sv1(:,:,:,:,:) = det_airamt(:,:,:,:,:)
    ent_airamt_tot_sv1(:,:,:) = ent_airamt_tot(:,:,:)
    det_airamt_tot_sv1(:,:,:) = det_airamt_tot(:,:,:)
!   end of simple calculation



!
!
!   do the full calculation of horizontal exchanges
!
!

!   reload the incoming values of ent/det_airamt_tot
    ent_airamt_tot(:,:,:) = ent_airamt_tot_sv0(:,:,:)
    det_airamt_tot(:,:,:) = det_airamt_tot_sv0(:,:,:)

!   calc the jgrp_of_jcls array
    icc = 1
    do jcls = 1, ncls_use
        if (mtype_updnenv_use(icc,jcls) == mtype_quiescn_ecpp) then
        jgrp_of_jcls(jcls) = 1
        else if (mtype_updnenv_use(icc,jcls) == mtype_updraft_ecpp) then
        jgrp_of_jcls(jcls) = 2
        else
        jgrp_of_jcls(jcls) = 3
        end if
    end do
    if (lunaa > 0) write(lunaa,'(a,10(2x,2i3))')  & 
                       'jcls  and jgrp_of_cls', (jcls, jgrp_of_jcls(jcls), jcls=1,ncls_use)

    ent_airamt(:,:,:,:,:) = 0.0
    det_airamt(:,:,:,:,:) = 0.0

    ecls_aaunasi_sv2(:,:,:) = 0.0
    egrp_aaunasi_sv2(:,:,:) = 0.0
    dcls_aaunasi_sv2(:,:,:) = 0.0
    dgrp_aaunasi_sv2(:,:,:) = 0.0


entdet_main_kloop_aa: &
    do k = kts, ktecen

    ardz_cut = afrac_cut*rhodz_cen(k)*0.3

    empty_old(:,:) = .false.
    empty_new(:,:) = .false.
    empty_oldnew(:,:) = .false.
    if (lunaa > 0) write(lunaa,'(/a)') 'k, jcls, emptyold/new/oldnew for icc=1 then icc=2'
    do jcls = 1, ncls_use
    do icc = 1, 2
        if (ardz_cen_old(k,icc,jcls) < ardz_cut) empty_old(icc,jcls) = .true.
        if (ardz_cen_new(k,icc,jcls) < ardz_cut) empty_new(icc,jcls) = .true.
        empty_oldnew(icc,jcls) = empty_old(icc,jcls) .and. empty_new(icc,jcls)
    end do
    if (lunaa > 0) write(lunaa,'(2i3,2(3x,3l3))') k, jcls,   &
        (empty_old(icc,jcls), empty_new(icc,jcls), empty_oldnew(icc,jcls), icc=1,2)
    end do

    if (lunaa > 0) then
        write(lunaa,'(/a,1p,10e16.8)') 'ardz_cut,rdz', ardz_cut, rhodz_cen(k)
        write(lunaa,'( a,1p,10e16.8)') 'ardz_cen_old', ardz_cen_old(k,0,0), ardz_cen_old(k,1:2,1:3) 
        write(lunaa,'( a,1p,10e16.8)') 'ardz_cen_new', ardz_cen_new(k,0,0), ardz_cen_new(k,1:2,1:3) 
        write(lunaa,'( a,1p,10e16.8)') 'new-old     ', (ardz_cen_new(k,0,0)-ardz_cen_new(k,0,0)), &
                                                       (ardz_cen_new(k,1:2,1:3)-ardz_cen_old(k,1:2,1:3))
        tmpa = 1.0/rhodz_cen(k)
        tmpb = sum( ardz_cen_old(k,1:2,1:3) )
        tmpc = sum( ardz_cen_new(k,1:2,1:3) )
        write(lunaa,'( a,1p,10e16.8)') 'area_cen_old', tmpa*tmpb, tmpa*ardz_cen_old(k,1:2,1:3) 
        write(lunaa,'( a,1p,10e16.8)') 'area_cen_new', tmpa*tmpc, tmpa*ardz_cen_new(k,1:2,1:3) 
        write(lunaa,'( a,1p,10e16.8)') 'new-old     ', tmpa*(tmpc-tmpb), &
                                                       tmpa*(ardz_cen_new(k,1:2,1:3)-ardz_cen_old(k,1:2,1:3))
        write(lunaa,'( a/1p,4(1x,3e11.3))') 'ardz_cen_old(0:2,0:3)', ardz_cen_old(k,0:2,0:3) 
        write(lunaa,'( a/1p,4(1x,3e11.3))') 'ardz_cen_new(0:2,0:3)', ardz_cen_new(k,0:2,0:3) 
    end if


!   step 1
!   initialize class and group "assigned" ent/det arrays to zero
!   initialize class "unassigned" ent/det arrays to ent/det_airamt_tot
!   calc group "unassigned" arrays by summing over classes
!
!   ***************************************************************
!   should check here that total ent = total det (with very small error allowed)
!      then adjust them to be even closer
!   ***************************************************************
    ecls_aa(:,:,:,:) = 0.0
    dcls_aa(:,:,:,:) = 0.0
    egrp_aa(:,:,:,:) = 0.0
    dgrp_aa(:,:,:,:) = 0.0
    egrp_aaunasi( :,:) = 0.0
    dgrp_aaunasi( :,:) = 0.0
    do jcls = 1, ncls_use
    do icc = 1, 2
        ecls_aaunasi(icc,jcls) = ent_airamt_tot(icc,jcls,k)
        dcls_aaunasi(icc,jcls) = det_airamt_tot(icc,jcls,k)
        jgrp = jgrp_of_jcls(jcls)
        egrp_aaunasi(icc,jgrp) = egrp_aaunasi(icc,jgrp) + ecls_aaunasi(icc,jcls)
        dgrp_aaunasi(icc,jgrp) = dgrp_aaunasi(icc,jgrp) + dcls_aaunasi(icc,jcls)
        if (ifrom_where < 10) then
        ! for area-change, detrainment is limited to initial subarea mass
        dcls_aalimit(icc,jcls) = ardz_cen_old(k,icc,jcls)
        else
        dcls_aalimit(icc,jcls) = 1.0e30
        end if
    end do
    end do
    call parampollu_tdx_entdet_diag01(   &
        1, lunaa,   &
        ifrom_where, ktau, k, kts, ktebnd, ktecen, ncls_use,   &
        ent_airamt_tot_sv1, ecls_aa, ecls_aaunasi, egrp_aa, egrp_aaunasi,   &
        det_airamt_tot_sv1, dcls_aa, dcls_aaunasi, dgrp_aa, dgrp_aaunasi,   &
        dcls_aalimit )


!   step 2
!   for up and dndrafts, if cloudy is entraining and clear is detraining
!       (or vice-versa), then assign as much as possible of the ent/det
!       as "clear up/dndraft" <--> "cloudy up/dndraft"
    do jcls = 1, ncls_use
        if (jcls == jcls_qu) cycle
        jgrp  = jgrp_of_jcls(jcls)
        jclsy = jcls
        jgrpy  = jgrp_of_jcls(jclsy)
        do icc = 1, 2
        iccy = 3 - icc
        if ( empty_old(icc ,jcls )  ) cycle
        if ( empty_new(iccy,jclsy)  ) cycle
        tmpa = min( dcls_aaunasi(icc,jcls), ecls_aaunasi(iccy,jcls) )
        if (tmpa > 0.0) then
            dcls_aaunasi(icc ,jcls ) = dcls_aaunasi(icc ,jcls ) - tmpa
            ecls_aaunasi(iccy,jclsy) = ecls_aaunasi(iccy,jclsy) - tmpa
            dcls_aalimit(icc ,jcls ) = dcls_aalimit(icc ,jcls ) - tmpa
            dcls_aa(icc ,jcls ,iccy,jclsy) = dcls_aa(icc ,jcls ,iccy,jclsy) + tmpa
            ecls_aa(iccy,jclsy,icc ,jcls ) = ecls_aa(iccy,jclsy,icc ,jcls ) + tmpa

            dgrp_aaunasi(icc ,jgrp ) = dgrp_aaunasi(icc ,jgrp ) - tmpa
            egrp_aaunasi(iccy,jgrpy) = egrp_aaunasi(iccy,jgrpy) - tmpa
            dgrp_aa(icc ,jgrp ,iccy,jgrpy) = dgrp_aa(icc ,jgrp ,iccy,jgrpy) + tmpa
            egrp_aa(iccy,jgrpy,icc ,jgrp ) = egrp_aa(iccy,jgrpy,icc ,jgrp ) + tmpa
        end if
        end do   ! icc
    end do   ! jcls
    call parampollu_tdx_entdet_diag01(   &
        2, lunaa,   &
        ifrom_where, ktau, k, kts, ktebnd, ktecen, ncls_use,   &
        ent_airamt_tot_sv1, ecls_aa, ecls_aaunasi, egrp_aa, egrp_aaunasi,   &
        det_airamt_tot_sv1, dcls_aa, dcls_aaunasi, dgrp_aa, dgrp_aaunasi,   &
        dcls_aalimit )


!   step 3
!   for up and dndraft detrainment, assign as much as possible of the det
!       as  "clear  up/dndraft" <--> "clear  quiescent"
!       and "cloudy up/dndraft" <--> "cloudy quiescent"
    do icc = 1, 2
        iccy = icc
        jclsy = jcls_qu
        jgrpy = jgrp_of_jcls(jclsy)

        ! tmpb = unassigned detrain from all up/dndraft
        tmpb = dgrp_aaunasi(icc,2) + dgrp_aaunasi(icc,3)
        ! tmpc = portion of tmpb that will be assigned in this step
        tmpc = min( tmpb, egrp_aaunasi(icc,1) )
        if (tmpc <= 0.0) cycle

        do jcls = 1, ncls_use
        if (jcls == jcls_qu) cycle
        if ( empty_old(icc ,jcls )  ) cycle
        if ( empty_new(iccy,jclsy)  ) cycle
        jgrp  = jgrp_of_jcls(jcls )

        ! tmpf is fraction of total-unassigned-draft detrainment due to this jcls
        tmpf = min( dcls_aaunasi(icc,jcls), tmpb ) / max( 1.0e-30_r8, tmpb )
        ! tmpa is portion of tmpc applied to this jcls
        tmpa = tmpf*tmpc
        tmpa = min( tmpa, dcls_aaunasi(icc ,jcls ), ecls_aaunasi(iccy,jclsy) )
        if (tmpa > 0.0) then
            dcls_aaunasi(icc ,jcls ) = dcls_aaunasi(icc ,jcls ) - tmpa
            ecls_aaunasi(iccy,jclsy) = ecls_aaunasi(iccy,jclsy) - tmpa
            dcls_aalimit(icc ,jcls ) = dcls_aalimit(icc ,jcls ) - tmpa
            dcls_aa(icc ,jcls ,iccy,jclsy) = dcls_aa(icc ,jcls ,iccy,jclsy) + tmpa
            ecls_aa(iccy,jclsy,icc ,jcls ) = ecls_aa(iccy,jclsy,icc ,jcls ) + tmpa

            dgrp_aaunasi(icc ,jgrp ) = dgrp_aaunasi(icc ,jgrp ) - tmpa
            egrp_aaunasi(iccy,jgrpy) = egrp_aaunasi(iccy,jgrpy) - tmpa
            dgrp_aa(icc ,jgrp ,iccy,jgrpy) = dgrp_aa(icc ,jgrp ,iccy,jgrpy) + tmpa
            egrp_aa(iccy,jgrpy,icc ,jgrp ) = egrp_aa(iccy,jgrpy,icc ,jgrp ) + tmpa
        end if
        end do   ! icc
    end do   ! jcls
    call parampollu_tdx_entdet_diag01(   &
        3, lunaa,   &
        ifrom_where, ktau, k, kts, ktebnd, ktecen, ncls_use,   &
        ent_airamt_tot_sv1, ecls_aa, ecls_aaunasi, egrp_aa, egrp_aaunasi,   &
        det_airamt_tot_sv1, dcls_aa, dcls_aaunasi, dgrp_aa, dgrp_aaunasi,   &
        dcls_aalimit )


!   step 4
!   for up and dndraft detrainment, assign any remaining detrainment to
!       quiescent based on the clear/cloudy quiescent areas

    ! tmpvecgg(1) = fraction of quiescent class that is clear (using new areas)
    tmpvecgg(1) = ardz_cen_new(k,1,jcls_qu)/ardz_cen_new(k,0,jcls_qu)
    tmpvecgg(1) = max( 0.0_r8, min( 1.0_r8, tmpvecgg(1) ) )
    ! tmpvecgg(2) = fraction of quiescent class that is cloudy (using new areas)
    tmpvecgg(2) = 1.0_r8 - tmpvecgg(1)
    tmpvecgg(2) = max( 0.0_r8, min( 1.0_r8, tmpvecgg(2) ) )

    ! tmpmatbb(0,0) = unassigned detrain from all up/dndraft
    ! tmpmatbb(1,0) = portion of tmpmatbb(0,0) from clear draft to all quiescent
    tmpmatbb(1,0) = sum( dgrp_aaunasi(1,2:3) )
    ! tmpmatbb(2,0) = portion of tmpmatbb(0,0) from cloudy draft to all quiescent
    tmpmatbb(2,0) = sum( dgrp_aaunasi(2,2:3) )
    tmpmatbb(0,0) = tmpmatbb(1,0) + tmpmatbb(2,0) 

    if (tmpmatbb(0,0) > 1.0e-30_r8) then

    ! tmpmatbb(0,1) = portion of tmpmatbb(0,0) from all draft to clear quiescent
    ! tmpmatbb(0,2) = portion of tmpmatbb(0,0) from all draft to cloudy quiescent
    tmpmatbb(0,1:2) = tmpmatbb(0,0)*tmpvecgg(1:2)

    ! this step can drive the ecls_aaunasi of a quiescent negative,
    !    and the negative entrainment gets converted to positive detrainment
    !    (from one quiescent subarea to the other)
    ! when doing area-change, check that this will not make 
    !    dcls_aaunasi exceed dcls_aalimit
    if (ifrom_where < 10) then
        tmpvecbb(1:2) = tmpmatbb(0,1:2)
        jclsy = jcls_qu
        do iccy = 2, 1, -1
        if (tmpvecbb(iccy) > ecls_aaunasi(iccy,jclsy)) then
            tmpd = dcls_aaunasi(iccy,jclsy)   &
                 + (tmpvecbb(iccy) - ecls_aaunasi(iccy,jclsy))
            if (tmpd > dcls_aalimit(iccy,jclsy)) then
            tmpvecbb(iccy) = tmpvecbb(iccy)    &
                           - (tmpd - dcls_aalimit(iccy,jclsy))
            tmpvecbb(iccy) = max( 0.0_r8, tmpvecbb(iccy) )
            tmpvecbb(3-iccy) = tmpmatbb(0,0) - tmpvecbb(iccy)
            end if
        end if
        end do
        tmpmatbb(0,1:2) = tmpvecbb(1:2)
    end if

    ! tmpmatbb(1,1) = portion of tmpmatbb(0,0) from clear draft to clear quiescent
    tmpmatbb(1,1) = min( tmpmatbb(0,1), tmpmatbb(1,0) )
    ! tmpmatbb(1,2) = portion of tmpmatbb(0,0) from clear draft to cloudy quiescent
    tmpmatbb(1,2) = max( 0.0_r8, (tmpmatbb(1,0) - tmpmatbb(1,1)) )

    ! tmpmatbb(2,2) = portion of tmpmatbb(0,0) from cloudy draft to cloudy quiescent
    tmpmatbb(2,2) = min( tmpmatbb(0,2), tmpmatbb(2,0) )
    ! tmpmatbb(2,1) = portion of tmpmatbb(0,0) from cloudy draft to clear quiescent
    tmpmatbb(2,1) = max( 0.0_r8, (tmpmatbb(2,0) - tmpmatbb(2,2)) )

    tmpmatff(1,2) = tmpmatbb(1,2) / max( 1.0e-37_r8, tmpmatbb(1,0) )
    tmpmatff(1,2) = max( 0.0_r8, min( 1.0_r8, tmpmatff(1,2) ) )
    tmpmatff(1,1) = 1.0_r8 - tmpmatff(1,2)
    tmpmatff(1,1) = max( 0.0_r8, min( 1.0_r8, tmpmatff(1,1) ) )

    tmpmatff(2,2) = tmpmatbb(2,2) / max( 1.0e-37_r8, tmpmatbb(2,0) )
    tmpmatff(2,2) = max( 0.0_r8, min( 1.0_r8, tmpmatff(2,2) ) )
    tmpmatff(2,1) = 1.0_r8 - tmpmatff(2,2)
    tmpmatff(2,1) = max( 0.0_r8, min( 1.0_r8, tmpmatff(2,1) ) )

! *** now need to apply these ***
    do jcls = 1, ncls_use
        if (jcls == jcls_qu) cycle   ! do jcls
        jgrp  = jgrp_of_jcls(jcls)
        jclsy = jcls_qu
        jgrpy  = jgrp_of_jcls(jclsy)
        do icc = 1, 2
        tmpc = dcls_aaunasi(icc,jcls)
        if (tmpc <= 0.0_r8) cycle   ! do icc

        do iccy = 1, 2
            if ( empty_old(icc,jcls) ) cycle   ! do iccy
            if ( empty_new(iccy,jclsy) ) cycle   ! do iccy

            tmpa = tmpmatff(icc,iccy) * tmpc
            if (tmpa <= 0.0_r8) cycle   ! do iccy

            dcls_aaunasi(icc ,jcls ) = dcls_aaunasi(icc ,jcls ) - tmpa
            ecls_aaunasi(iccy,jclsy) = ecls_aaunasi(iccy,jclsy) - tmpa
            dcls_aalimit(icc ,jcls ) = dcls_aalimit(icc ,jcls ) - tmpa
            dcls_aa(icc ,jcls ,iccy,jclsy) = dcls_aa(icc ,jcls ,iccy,jclsy) + tmpa
            ecls_aa(iccy,jclsy,icc ,jcls ) = ecls_aa(iccy,jclsy,icc ,jcls ) + tmpa

            dgrp_aaunasi(icc ,jgrp ) = dgrp_aaunasi(icc ,jgrp ) - tmpa
            egrp_aaunasi(iccy,jgrpy) = egrp_aaunasi(iccy,jgrpy) - tmpa
            dgrp_aa(icc ,jgrp ,iccy,jgrpy) = dgrp_aa(icc ,jgrp ,iccy,jgrpy) + tmpa
            egrp_aa(iccy,jgrpy,icc ,jgrp ) = egrp_aa(iccy,jgrpy,icc ,jgrp ) + tmpa

            ! if unassigned entrainment from quiescent goes negative,
            ! convert this to positive unassigned detrainment
            if (ecls_aaunasi(iccy,jclsy) < 0.0_r8) then
            dcls_aaunasi(iccy,jclsy) = dcls_aaunasi(iccy,jclsy) - ecls_aaunasi(iccy,jclsy)
            ecls_aaunasi(iccy,jclsy) = 0.0
            end if
            if (egrp_aaunasi(iccy,jgrpy) < 0.0_r8) then
            dgrp_aaunasi(iccy,jgrpy) = dgrp_aaunasi(iccy,jgrpy) - egrp_aaunasi(iccy,jgrpy)
            egrp_aaunasi(iccy,jgrpy) = 0.0
            end if
        end do   ! iccy
        end do   ! icc
    end do   ! jcls

    end if   ! (tmpmatbb(0,0) > 1.0e-30_r8)
    call parampollu_tdx_entdet_diag01(   &
        4, lunaa,   &
        ifrom_where, ktau, k, kts, ktebnd, ktecen, ncls_use,   &
        ent_airamt_tot_sv1, ecls_aa, ecls_aaunasi, egrp_aa, egrp_aaunasi,   &
        det_airamt_tot_sv1, dcls_aa, dcls_aaunasi, dgrp_aa, dgrp_aaunasi,   &
        dcls_aalimit )


!   step 5
!   up and dndraft entrainment
!   do this in a much simpler manner
!       all up and dndraft entrainment comes from quiescent
!       contributions from clear and cloudy quiescent are proportional to
!          their fractional areas (tmpvecgg(1) & tmpvecgg(2))
    ! tmpvecgg(1) = fraction of quiescent class that is clear (using old areas)
    tmpvecgg(1) = ardz_cen_old(k,1,jcls_qu)/ardz_cen_old(k,0,jcls_qu)
    tmpvecgg(1) = max( 0.0_r8, min( 1.0_r8, tmpvecgg(1) ) )
    ! tmpvecgg(2) = fraction of quiescent class that is cloudy (using old areas)
    tmpvecgg(2) = 1.0_r8 - tmpvecgg(1)
    tmpvecgg(2) = max( 0.0_r8, min( 1.0_r8, tmpvecgg(2) ) )

    jclsy = jcls_qu
    jgrpy = jgrp_of_jcls(jclsy)

    ! when doing area-change, check that this will not make 
    !    dcls_aalimit negative for either quiescent subarea
    if (ifrom_where < 10) then
        ! total unassigned entrainment to up/dndrafts
        tmpa = sum( egrp_aaunasi(1:2,2:3) )
        ! amount of detrainment that will come from quiescent iccy=1,2
        tmpvecbb(1:2) = tmpa*tmpvecgg(1:2)
        jclsy = jcls_qu
        do iccy = 2, 1, -1
        if (tmpvecbb(iccy) > dcls_aalimit(iccy,jclsy)) then
            tmpvecbb(iccy) = dcls_aalimit(iccy,jclsy)
            tmpvecbb(3-iccy) = tmpa - tmpvecbb(iccy)
        end if
        end do
        tmpvecgg(2) = tmpvecbb(2)/max( 1.0e-37_r8, tmpa )
        tmpvecgg(2) = max( 0.0_r8, min( 1.0_r8, tmpvecgg(2) ) )
        tmpvecgg(1) = 1.0_r8 - tmpvecgg(2)
        tmpvecgg(1) = max( 0.0_r8, min( 1.0_r8, tmpvecgg(1) ) )
    end if

    do jcls = 1, ncls_use
    do icc = 1, 2
        iccy = 0
        if (jcls == jcls_qu) cycle
        if ( empty_new(icc ,jcls )  ) cycle
        jgrp  = jgrp_of_jcls(jcls )

        ! tmpa is unassigned-draft entrainment due to this icc,jcls
        tmpa = ecls_aaunasi(icc,jcls)
        if (tmpa > 0.0) then
        do iccy = 1, 2
            if ( empty_old(iccy,jclsy)  ) cycle
            if (tmpvecgg(iccy) <= 0.0_r8) cycle
            ! tmpb is portion of tmpa coming from iccy,jclsy
            tmpb = tmpa*tmpvecgg(iccy)

            ecls_aaunasi(icc ,jcls ) = ecls_aaunasi(icc ,jcls ) - tmpb
            dcls_aaunasi(iccy,jclsy) = dcls_aaunasi(iccy,jclsy) - tmpb
            dcls_aalimit(iccy,jclsy) = dcls_aalimit(iccy,jclsy) - tmpb
            ecls_aa(icc ,jcls ,iccy,jclsy) = ecls_aa(icc ,jcls ,iccy,jclsy) + tmpb
            dcls_aa(iccy,jclsy,icc ,jcls ) = dcls_aa(iccy,jclsy,icc ,jcls ) + tmpb

            egrp_aaunasi(icc ,jgrp ) = egrp_aaunasi(icc ,jgrp ) - tmpb
            dgrp_aaunasi(iccy,jgrpy) = dgrp_aaunasi(iccy,jgrpy) - tmpb
            egrp_aa(icc ,jgrp ,iccy,jgrpy) = egrp_aa(icc ,jgrp ,iccy,jgrpy) + tmpb
            dgrp_aa(iccy,jgrpy,icc ,jgrp ) = dgrp_aa(iccy,jgrpy,icc ,jgrp ) + tmpb

            ! if unassigned detrainment from quiescent goes negative,
            ! convert this to positive unassigned entrainment
            if (dcls_aaunasi(iccy,jclsy) < 0.0_r8) then
            ecls_aaunasi(iccy,jclsy) = ecls_aaunasi(iccy,jclsy) - dcls_aaunasi(iccy,jclsy)
            dcls_aaunasi(iccy,jclsy) = 0.0
            end if
            if (dgrp_aaunasi(iccy,jgrpy) < 0.0_r8) then
            egrp_aaunasi(iccy,jgrpy) = egrp_aaunasi(iccy,jgrpy) - dgrp_aaunasi(iccy,jgrpy)
            dgrp_aaunasi(iccy,jgrpy) = 0.0
            end if
        end do   ! iccy
        end if   ! (tmpa > 0.0)
    end do   ! icc
    end do   ! jcls
    call parampollu_tdx_entdet_diag01(   &
        5, lunaa,   &
        ifrom_where, ktau, k, kts, ktebnd, ktecen, ncls_use,   &
        ent_airamt_tot_sv1, ecls_aa, ecls_aaunasi, egrp_aa, egrp_aaunasi,   &
        det_airamt_tot_sv1, dcls_aa, dcls_aaunasi, dgrp_aa, dgrp_aaunasi,   &
        dcls_aalimit )


!   step 6
!   quiescent clear <--> quiescent cloudy exchanges
!   if clear  is detraining and cloudy is entraining, then assign as much as 
!       possible of the det/ent as "clear quiescent" --> "cloudy quiescent"
!   if cloudy is detraining and clear  is entraining, then assign as much as 
!       possible of the det/ent as "cloudy quiescent" --> "clear quiescent"
    do jcls = 1, ncls_use
        if (jcls /= jcls_qu) cycle
        jgrp  = jgrp_of_jcls(jcls)
        jclsy = jcls
        jgrpy  = jgrp_of_jcls(jclsy)
        do icc = 1, 2
        iccy = 3 - icc
        if ( empty_old(icc ,jcls )  ) cycle
        if ( empty_new(iccy,jclsy)  ) cycle
        tmpa = min( dcls_aaunasi(icc,jcls), ecls_aaunasi(iccy,jcls) )
        if (tmpa > 0.0) then
            dcls_aaunasi(icc ,jcls ) = dcls_aaunasi(icc ,jcls ) - tmpa
            ecls_aaunasi(iccy,jclsy) = ecls_aaunasi(iccy,jclsy) - tmpa
            dcls_aa(icc ,jcls ,iccy,jclsy) = dcls_aa(icc ,jcls ,iccy,jclsy) + tmpa
            ecls_aa(iccy,jclsy,icc ,jcls ) = ecls_aa(iccy,jclsy,icc ,jcls ) + tmpa

            dgrp_aaunasi(icc ,jgrp ) = dgrp_aaunasi(icc ,jgrp ) - tmpa
            egrp_aaunasi(iccy,jgrpy) = egrp_aaunasi(iccy,jgrpy) - tmpa
            dgrp_aa(icc ,jgrp ,iccy,jgrpy) = dgrp_aa(icc ,jgrp ,iccy,jgrpy) + tmpa
            egrp_aa(iccy,jgrpy,icc ,jgrp ) = egrp_aa(iccy,jgrpy,icc ,jgrp ) + tmpa
        end if
        end do   ! icc
    end do   ! jcls
    call parampollu_tdx_entdet_diag01(   &
        6, lunaa,   &
        ifrom_where, ktau, k, kts, ktebnd, ktecen, ncls_use,   &
        ent_airamt_tot_sv1, ecls_aa, ecls_aaunasi, egrp_aa, egrp_aaunasi,   &
        det_airamt_tot_sv1, dcls_aa, dcls_aaunasi, dgrp_aa, dgrp_aaunasi,   &
        dcls_aalimit )



!   load the current-k ent/det values for each class into ent/det_airamt
    ent_airamt(:,:,:,:,k) = ecls_aa(:,:,:,:)
    det_airamt(:,:,:,:,k) = dcls_aa(:,:,:,:)

    ecls_aaunasi_sv2(:,:,k) = ecls_aaunasi(:,:)
    egrp_aaunasi_sv2(:,:,k) = egrp_aaunasi(:,:)
    dcls_aaunasi_sv2(:,:,k) = dcls_aaunasi(:,:)
    dgrp_aaunasi_sv2(:,:,k) = dgrp_aaunasi(:,:)


!   calc largest unassigned ent/det
    m = 1
    if (ifrom_where == 10) m = 2
    if (ifrom_where ==  2) m = 3
    do jcls = 1, ncls_use
    do icc = 1, 2
        if (abs(ecls_aaunasi(icc,jcls)) > abs(ecls_aaunasi_worst(m))) then
        ecls_aaunasi_worst(m) = ecls_aaunasi(icc,jcls)
        ecls_aaunasi_worst_i(m) = icc
        ecls_aaunasi_worst_j(m) = jcls
        ecls_aaunasi_worst_k(m) = k
        ecls_aaunasi_worst_ktau(m) = ktau
        end if
        if (abs(dcls_aaunasi(icc,jcls)) > abs(dcls_aaunasi_worst(m))) then
        dcls_aaunasi_worst(m) = dcls_aaunasi(icc,jcls)
        dcls_aaunasi_worst_i(m) = icc
        dcls_aaunasi_worst_j(m) = jcls
        dcls_aaunasi_worst_k(m) = k
        dcls_aaunasi_worst_ktau(m) = ktau
        end if
        jgrp = jcls
        if (jgrp > 3) cycle
        if (abs(egrp_aaunasi(icc,jgrp)) > abs(egrp_aaunasi_worst(m))) then
        egrp_aaunasi_worst(m) = egrp_aaunasi(icc,jgrp)
        egrp_aaunasi_worst_i(m) = icc
        egrp_aaunasi_worst_j(m) = jgrp
        egrp_aaunasi_worst_k(m) = k
        egrp_aaunasi_worst_ktau(m) = ktau
        end if
        if (abs(dgrp_aaunasi(icc,jgrp)) > abs(dgrp_aaunasi_worst(m))) then
        dgrp_aaunasi_worst(m) = dgrp_aaunasi(icc,jgrp)
        dgrp_aaunasi_worst_i(m) = icc
        dgrp_aaunasi_worst_j(m) = jgrp
        dgrp_aaunasi_worst_k(m) = k
        dgrp_aaunasi_worst_ktau(m) = ktau
        end if
    end do
    end do



    end do entdet_main_kloop_aa


!   now calc ent/det_airamt_tot
    do k = kts, ktecen

    do jcls = 1, ncls_use
    do icc = 1, 2
        tmpa = 0.0
        tmpb = 0.0
        if (k < ktebnd) then
        do jclsy = 1, ncls_use
        do iccy = 1, 2
            tmpa = tmpa + ent_airamt( icc,jcls, iccy,jclsy, k)
            tmpb = tmpb + det_airamt( icc,jcls, iccy,jclsy, k)
        end do
        end do
        end if
        ent_airamt_tot(icc,jcls,k) = tmpa
        det_airamt_tot(icc,jcls,k) = tmpb
    end do   ! icc
    end do   ! jcls

    end do   ! k


!   diagnostic output
    if (lunaa > 0) then
    do k = kts, ktecen

    write(lunaa,'(/a,3i5)') 'bb parampollu_tdx_entdet_sub1 - ktau, ifrom_where, k', ktau, ifrom_where, k

    write(lunaa,'(a)') 'ent_airamt_tot simple/full'
    write(lunaa,'(1p,10e11.3)') ent_airamt_tot_sv1(1:2,1:ncls_use,k)
    write(lunaa,'(1p,10e11.3)') ent_airamt_tot(    1:2,1:ncls_use,k)
    do jcls = 1, ncls_use
    write(lunaa,'(a,i3,a,i3)') 'ent_airamt simple/full for icc,jcls=  1', jcls, '   and  2', jcls
    write(lunaa,'(1p,6e11.3,4x,6e11.3)') (ent_airamt_sv1(icc,jcls,1:2,1:ncls_use,k), icc=1,2)
    write(lunaa,'(1p,6e11.3,4x,6e11.3)') (ent_airamt(    icc,jcls,1:2,1:ncls_use,k), icc=1,2)
    end do

    write(lunaa,'(a)') 'det_airamt_tot simple/full'
    write(lunaa,'(1p,10e11.3)') det_airamt_tot_sv1(1:2,1:ncls_use,k)
    write(lunaa,'(1p,10e11.3)') det_airamt_tot(    1:2,1:ncls_use,k)
    do jcls = 1, ncls_use
    write(lunaa,'(a,i3,a,i3)') 'det_airamt simple/full for icc,jcls=  1', jcls, '   and  2', jcls
    write(lunaa,'(1p,6e11.3,4x,6e11.3)') (det_airamt_sv1(icc,jcls,1:2,1:ncls_use,k), icc=1,2)
    write(lunaa,'(1p,6e11.3,4x,6e11.3)') (det_airamt(    icc,jcls,1:2,1:ncls_use,k), icc=1,2)
    end do

    write(lunaa,'(a)') 'final ecls_aaunasi & egrp_aaunasi // final dcls_aaunasi & dgrp_aaunasi'
    write(lunaa,'(1p,6e11.3,4x,6e11.3)') ecls_aaunasi_sv2(1:2,1:ncls_use,k), egrp_aaunasi_sv2(1:2,1:3,k)
    write(lunaa,'(1p,6e11.3,4x,6e11.3)') dcls_aaunasi_sv2(1:2,1:ncls_use,k), dgrp_aaunasi_sv2(1:2,1:3,k)

    end do   ! k = kts, kte
    end if   ! (lunaa > 0)


    lunbb = -1
    if ((parampollu_opt == 2223) .and. (ifrom_where ==  2)) lunbb = ldiagaa_ecpp(123)
    if ((parampollu_opt == 2220) .and. (ifrom_where == 10)) lunbb = ldiagaa_ecpp(123)
    lunbb = ldiagaa_ecpp(123)
    if (idiagaa_ecpp(123) <= 0) lunbb = -1

    if (lunbb > 0) then
    write(lunbb,'(/a,3i5)') 'parampollu_tdx_entdet_sub1 - ktau, ifrom_where', ktau, ifrom_where

    do m = 1, 3
    write(lunbb,'(a,i3)') 'm =', m
    write(lunbb,'(a,2(3x,3i3,i5,1p,e11.3))')   &
        'ecls_aaunasi_worst i/j/k/ktau/val  &  dcls',   &
        ecls_aaunasi_worst_i(m), ecls_aaunasi_worst_j(m), ecls_aaunasi_worst_k(m),   &
        ecls_aaunasi_worst_ktau(m), ecls_aaunasi_worst(m),   &
        dcls_aaunasi_worst_i(m), dcls_aaunasi_worst_j(m), dcls_aaunasi_worst_k(m),   &
        dcls_aaunasi_worst_ktau(m), dcls_aaunasi_worst(m)
    write(lunbb,'(a,2(3x,3i3,i5,1p,e11.3))')   &
        'egrp_aaunasi_worst i/j/k/ktau/val  &  dgrp',   &
        egrp_aaunasi_worst_i(m), egrp_aaunasi_worst_j(m), egrp_aaunasi_worst_k(m),   &
        egrp_aaunasi_worst_ktau(m), egrp_aaunasi_worst(m),   &
        dgrp_aaunasi_worst_i(m), dgrp_aaunasi_worst_j(m), dgrp_aaunasi_worst_k(m),   &
        dgrp_aaunasi_worst_ktau(m), dgrp_aaunasi_worst(m)
    end do

    end if   ! (lunbb > 0)


!   restore saved values
!   ent_airamt(:,:,:,:,:) = ent_airamt_sv1(:,:,:,:,:)
!   det_airamt(:,:,:,:,:) = det_airamt_sv1(:,:,:,:,:)
!   ent_airamt_tot(:,:,:) = ent_airamt_tot_sv1(:,:,:)
!   det_airamt_tot(:,:,:) = det_airamt_tot_sv1(:,:,:)


    return
    end subroutine parampollu_tdx_entdet_sub1



!-----------------------------------------------------------------------
    subroutine parampollu_tdx_entdet_diag01(   &
        istep, lun,   &
        ifrom_where, ktau, k, kts, ktebnd, ktecen, ncls_use,   &
        ent_airamt_tot_sv1, ecls_aa, ecls_aaunasi, egrp_aa, egrp_aaunasi,   &
        det_airamt_tot_sv1, dcls_aa, dcls_aaunasi, dgrp_aa, dgrp_aaunasi,   &
        dcls_aalimit )

    use module_data_ecpp1

    integer :: istep, lun, ifrom_where, ktau, k, kts, ktebnd, ktecen, ncls_use
    real(r8), dimension( 1:2, 1:maxcls_ecpp, kts:ktecen  ) :: &
        ent_airamt_tot_sv1, det_airamt_tot_sv1
    real(r8), dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp ) ::   &
        ecls_aa, dcls_aa
    real(r8), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        ecls_aaunasi, dcls_aaunasi, dcls_aalimit
    real(r8), dimension( 1:2, 1:3, 1:2, 1:3 ) ::   &
        egrp_aa, dgrp_aa
    real(r8), dimension( 1:2, 1:3 ) ::   &
        egrp_aaunasi, dgrp_aaunasi

    integer :: icc, jcls

    if (lun <= 0) return

    write(lun,'(/a,i1,a,3i5)') 'aa', istep, ' parampollu_tdx_entdet_sub1 - ktau, ifrom_where, k', ktau, ifrom_where, k

    write(lun,'(/i3,a)') istep, '=istep - ent_airamt_tot_sv1'
    write(lun,'(1p,10e16.8)') ent_airamt_tot_sv1(1:2,1:ncls_use,k)
    write(lun,'(i3,a)') istep, '=istep - ecls_aaunasi after'
    write(lun,'(1p,10e16.8)') ecls_aaunasi(1:2,1:ncls_use)
    write(lun,'(i3,a)') istep, '=istep - egrp_aaunasi after'
    write(lun,'(1p,10e16.8)') egrp_aaunasi(1:2,1:3)
    do jcls = 1, ncls_use
    write(lun,'(i3,a,i3,a,i3)') istep, '=istep - ecls_aa after for icc,jcls=  1', jcls, '   and  2', jcls
    write(lun,'(1p,6e16.8)') (ecls_aa(icc,jcls,1:2,1:ncls_use), icc=1,2)
    if (jcls > 3) cycle
    write(lun,'(i3,a,i3,a,i3)') istep, '=istep - egrp_aa after for icc,jcls=  1', jcls, '   and  2', jcls
    write(lun,'(1p,6e16.8)') (egrp_aa(icc,jcls,1:2,1:3), icc=1,2)
    end do

    write(lun,'(/i3,a)') istep, '=istep - det_airamt_tot_sv1'
    write(lun,'(1p,10e16.8)') det_airamt_tot_sv1(1:2,1:ncls_use,k)
    write(lun,'(i3,a)') istep, '=istep - dcls_aalimit after'
    write(lun,'(1p,10e16.8)') dcls_aalimit(1:2,1:ncls_use)
    write(lun,'(i3,a)') istep, '=istep - dcls_aaunasi after'
    write(lun,'(1p,10e16.8)') dcls_aaunasi(1:2,1:ncls_use)
    write(lun,'(i3,a)') istep, '=istep - dgrp_aaunasi after'
    write(lun,'(1p,10e16.8)') dgrp_aaunasi(1:2,1:3)
    do jcls = 1, ncls_use
    write(lun,'(i3,a,i3,a,i3)') istep, '=istep - dcls_aa after for icc,jcls=  1', jcls, '   and  2', jcls
    write(lun,'(1p,6e16.8)') (dcls_aa(icc,jcls,1:2,1:ncls_use), icc=1,2)
    if (jcls > 3) cycle
    write(lun,'(i3,a,i3,a,i3)') istep, '=istep - dgrp_aa after for icc,jcls=  1', jcls, '   and  2', jcls
    write(lun,'(1p,6e16.8)') (dgrp_aa(icc,jcls,1:2,1:3), icc=1,2)
    end do

    return
    end subroutine parampollu_tdx_entdet_diag01

!-----------------------------------------------------------------------
    subroutine set_of_aerosol_stuff(is_aerosol, &
        iphase_of_aerosol, isize_of_aerosol, itype_of_aerosol,   &
        inmw_of_aerosol, laicwpair_of_aerosol )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! sets following arrays 
!
! is_aerosol   : logical variable, whether it is an aeroosl speices or not
!
! iphase_of_aerosol(l) = 0              for non-aerosol species
!                      = ai/cw/..._phase    for aerosol species
! isize_of_aerosol(l)  = 0              for non-aerosol species
!                      = size/bin index     for aerosol species
! itype_of_aerosol(l)  = 0              for non-aerosol species
!                      = type index         for aerosol species
! inmw_of_aerosol(l)   = 0              for non-aerosol species
!                      = 1/2/3          for aerosol number/mass/water species
! laicwpair_of_aerosol(l) = -999888777  for non-aerosol species
!                      = species index of corresponding ai/cw species
!
!-----------------------------------------------------------------------

!   use module_configure, only:  chem_dname_table

    use module_data_ecpp1, only:  num_chem_ecpp, param_first_ecpp

    use module_data_mosaic_asect, only:  ai_phase, cw_phase,   &
        massptr_aer,   &
        ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer

!   arguments
    integer, intent(out), dimension( 1:num_chem_ecpp ) ::   &
        iphase_of_aerosol, isize_of_aerosol, itype_of_aerosol,   &
        inmw_of_aerosol, laicwpair_of_aerosol
        logical, intent(out) :: is_aerosol(1:num_chem_ecpp)

!   local variables
    integer :: j, j2, l, l2, ll, m, n
    integer, save :: ientry = 0
    character(len=16) :: tmpname

        is_aerosol (:) = .false.
    iphase_of_aerosol(:) = 0
    isize_of_aerosol(:) = 0
    itype_of_aerosol(:) = 0
    laicwpair_of_aerosol(:) = -999888777
    inmw_of_aerosol(:) = 0

    do j = 1, nphase_aer
    do n = 1, ntype_aer
    do m = 1, nsize_aer(n)
    do ll = 0, ncomp_aer(n)

        l = -999888777
        if (ll == 0) then
        l = numptr_aer(m,n,j)
        else if (ll <= ncomp_aer(n)) then
        l = massptr_aer(ll,m,n,j)
        end if
        if ((l >= param_first_ecpp) .and. (l <= num_chem_ecpp)) then
                is_aerosol(l) = .true.
        iphase_of_aerosol(l) = j
        isize_of_aerosol(l) = m
        itype_of_aerosol(l) = n
        if (ll == 0) then
            inmw_of_aerosol(l) = 1
        else if (ll <= ncomp_aer(n)) then
            inmw_of_aerosol(l) = 2
        else
            inmw_of_aerosol(l) = 3
        end if
        end if

        if ( (nphase_aer >= 2) .and.   &
             (ai_phase > 0) .and. (cw_phase > 0) ) then
        if (j == ai_phase) then
            j2 = cw_phase
        else if (j == cw_phase) then
            j2 = ai_phase
        else
            cycle
        end if
        end if
        if (ll == 0) then
        l2 = numptr_aer(m,n,j2)
        else if (ll <= ncomp_aer(n)) then
        l2 = massptr_aer(ll,m,n,j2)
        else
        cycle
        end if
        if ((l  >= param_first_ecpp) .and. (l  <= num_chem_ecpp) .and.   &
            (l2 >= param_first_ecpp) .and. (l2 <= num_chem_ecpp))   &
        laicwpair_of_aerosol(l) = l2

    end do
    end do
    end do
    end do

    if (ientry == 0) then
    do l = param_first_ecpp, num_chem_ecpp
!       tmpname = chem_dname_table(1,l)
!       write(*,'(2a,6i5)') 'iphase, isize, itype, inmw, l, laicw_pairptr   ', tmpname,   &
!==Guangxing Lin
!           write(*,'(2a,7i5)') 'iphase, isize, itype, inmw, l, laicw_pairptr   ',  &
           write(*,'(a,l2,7i5)') 'iphase, isize, itype, inmw, l, laicw_pairptr   ',     &
        is_aerosol(l), iphase_of_aerosol(l), isize_of_aerosol(l), itype_of_aerosol(l),   &
        inmw_of_aerosol(l), l, max(-999,laicwpair_of_aerosol(l))
    end do
    end if
    ientry = 1

    return
    end subroutine set_of_aerosol_stuff

!-----------------------------------------------------------------------
    subroutine parampollu_tdx_startup(                        &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        rhocen_bar, dzcen,                                &
        chem_bar, chem_cls,                               &
        ncls_ecpp,                                        &
        acen_tbeg,                                        &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use,                                         &
        chem_sub_beg,                                     &
        acen_tbeg_use, ardz_cen_tbeg, rhodz_cen,          &
        activate_onoff_use,                               &
        iphase_of_aerosol, laicwpair_of_aerosol           )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_startup does some "startup" calculations
!
!    re-initializes the acen_tbeg to all-quiescent and the 
!       chem_cls to chem_bar at the re-init time (if this is turned on)
!
!    calculates chem_sub from chem_cls (which involves some assumptions
!       for the interstial and activated aerosols)
!
!-----------------------------------------------------------------------

!   use module_state_descption, only:  &
!   p_num_a01, p_num_cw01, p_oin_a01, p_oin_cw01, &
!   p_num_a03, p_num_cw03, p_oin_a03, p_oin_cw03
!        use module_data_ecpp1, only:  &
!        p_num_a01, p_num_cw01, p_oin_a01, p_oin_cw01, &
!        p_num_a03, p_num_cw03, p_oin_a03, p_oin_cw03

    use module_data_radm2, only:  epsilc

    use module_data_mosaic_asect, only: ai_phase, cw_phase

    use module_data_ecpp1

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message   

!   arguments
    integer, intent(in) ::                  &
        ktau, ktau_pp,              &
        it, jt, kts, ktebnd, ktecen
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!   chem_driver and routines under it do calculations
!   over these spatial indices.

    integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

    real(r8), intent(in) :: dtstep, dtstep_pp
!   dtstep - main model time step (s)
!   dtstep_pp - time step (s) for "parameterized pollutants" calculations

    real(r8), intent(in), dimension( kts:ktecen ) ::   &
        rhocen_bar, dzcen
!   rhocen_bar, rhobnd_bar - dry air density (kg/m^3) at layer centers and boundaries
!   dzcen - layer thicknesses (m)
!
    real(r8), intent(in), dimension( kts:ktecen, 1:num_chem_ecpp ) :: &
        chem_bar
!   chem_bar - mixing ratios of trace gase (ppm) and aerosol species
!   (ug/kg for mass species, #/kg for number species)

    real(r8), intent(inout), dimension( kts:ktecen, 1:maxcls_ecpp, 1:num_chem_ecpp ) :: &
        chem_cls
 
    integer, intent(in) :: ncls_ecpp
!   ncls_ecpp - number of ecpp transport classes in the grid column
    
    real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
        acen_tbeg

    integer, intent(in)    :: ncls_use

    real(r8), intent(inout),   &
        dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
        chem_sub_beg

    real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) :: &
        acen_tbeg_use, ardz_cen_tbeg

    real(r8), intent(inout), dimension( kts:ktecen ) :: rhodz_cen

    integer, intent(in) :: activate_onoff_use

    integer, intent(in), dimension( 1:num_chem_ecpp ) ::   &
        iphase_of_aerosol, laicwpair_of_aerosol


!   local variables
    integer :: icc, itmpa, jcls, jclsbb
    integer :: k, l, la, laa, lbb, lc
    integer :: lun161, lun162, lun164
    integer :: p1st

    real(r8) :: tmpa, tmpb, tmpqbarold
    real(r8), dimension( 0:2 ) :: tmp_acen
    real(r8), dimension( 1:num_chem_ecpp ) :: tmp_chem_cls
    real(r8), dimension( 1:2, 1:num_chem_ecpp ) :: tmp_chem_sub



    p1st = param_first_ecpp
    lun161 = -1
    if (idiagaa_ecpp(161) > 0) lun161 = ldiagaa_ecpp(161)
    lun162 = -1
    if (idiagaa_ecpp(162) > 0) lun162 = ldiagaa_ecpp(162)
    lun164 = -1
    if (idiagaa_ecpp(164) > 0) lun164 = ldiagaa_ecpp(164)

!   do sums of fractional areas over clear/cloudy and classes
    do k = kts, ktecen
        do jcls = 1, ncls_use
        acen_tbeg(k,0,jcls) = sum( acen_tbeg(k,1:2,jcls) )
        end do
        do icc = 0, 2
        tmpa = 0.0
        do jclsbb = 2, ncls_use+1
        ! sum order is [2,3,...,ncls,1] instead of [1,2,...,ncls]
            jcls = mod(jclsbb-1,ncls_use) + 1
            tmpa = tmpa + acen_tbeg(k,icc,jcls)
        end do
        acen_tbeg(k,icc,0) = tmpa
        end do
    end do


!
!   with hybrid-time-dependent drafts, always do reinit calcs
!
!   set all chem_cls = chem_bar for all species and levels
    chem_cls(:,:,:) = 0.0
    do l = 1, num_chem_ecpp
    do jcls = 1, ncls_use
    do k = kts, ktecen
        chem_cls(k,jcls,l) = chem_bar(k,l)
    end do
    end do
    end do

!   set up/dndraft areas to zero
!   set quiescent areas to overall clear/cloudy fractions
    do k = kts, ktecen
        tmpa = acen_tbeg(k,1,0)   ! this is total clear area (all classes)
        tmpa = max( 0.0_r8, min( 1.0_r8, tmpa ) )

!   force 100%/0%/70%/30% clear when iflag_ecpp_test_fixed_fcloud = 2/3/4/5
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

        acen_tbeg(k,:,:) = 0.0
        acen_tbeg(k,0,jcls_qu) = 1.0_r8
        acen_tbeg(k,1,jcls_qu) = tmpa
        acen_tbeg(k,2,jcls_qu) = 1.0_r8-tmpa
        acen_tbeg(k,0:2,0) = acen_tbeg(k,0:2,jcls_qu)
    end do


!
!   update the chem_cls values based on "host-code" changes to chem_bar
!   when iflag_ecpp_startup_host_chemtend > 0
    if (iflag_ecpp_startup_host_chemtend > 0) then
        do l = p1st, num_chem_ecpp
        do k = kts, ktecen
        tmpa = 0.0
        tmpb = 0.0
        do jcls = 1, ncls_use
            tmpa = tmpa + acen_tbeg(k,0,jcls)*chem_cls(k,jcls,l)
            tmpb = tmpb + acen_tbeg(k,0,jcls)
        end do
        tmpqbarold = tmpa/max(tmpb,0.99_r8)
        if (tmpqbarold < 1.01_r8*max(epsilc,1.0e-20)) then
            chem_cls(k,1:ncls_use,l) = chem_bar(k,l)
        else if (chem_bar(k,l) > tmpqbarold) then
            chem_cls(k,1:ncls_use,l) = chem_cls(k,1:ncls_use,l) + (chem_bar(k,l)-tmpqbarold)
        else
            chem_cls(k,1:ncls_use,l) = chem_cls(k,1:ncls_use,l) * (chem_bar(k,l)/tmpqbarold)
        end if
        end do
        end do
    end if


!   do chem_sub_beg <-- chem_cls and acen_tbeg_use <-- acen_tbeg
!   TODO - for aerosols, special treatment for "a" and "cw" in clear/cloudy sub-classes
    acen_tbeg_use(:,:,:) = acen_tbeg(:,:,:)
    chem_sub_beg(:,:,:,:) = 0.0

    do k = kts, ktecen
        do jcls = 0, ncls_use
        ardz_cen_tbeg(k,0:2,jcls) = acen_tbeg_use(k,0:2,jcls)*rhodz_cen(k)
        end do
    end do

    do jcls = 1, ncls_use
    do k = kts, ktecen
        do l = p1st, num_chem_ecpp
        chem_sub_beg(k,1:2,jcls,l) = chem_cls(k,jcls,l)
        end do
    end do
    end do

!   for aerosols, special treatment for "a" and "cw" in clear/cloudy sub-classes
    if ((activate_onoff_use > 0) .and. (iflag_ecpp_startup_acw_partition > 0)) then

acwxx1_jcls_loop: &
    do jcls = 1, ncls_use
acwxx1_k_loop: &
    do k = kts, ktecen

        ! clear subarea ~= 0 --> all cloudy
        ! no special treatment in this case
        if (acen_tbeg_use(k,1,jcls) < afrac_cut_0p5) cycle acwxx1_k_loop

        ! cloudy subarea ~= 0 and clear subarea > 0 
        ! resuspend any cloudborne material
        if (acen_tbeg_use(k,2,jcls) < afrac_cut_0p5) then
        do la = p1st, num_chem_ecpp
            if (iphase_of_aerosol(la) /= ai_phase) cycle
            lc = laicwpair_of_aerosol(la)
            if (lc < p1st) cycle
            if (iphase_of_aerosol(lc) /= cw_phase) cycle

            tmpa = chem_cls(k,jcls,la) + chem_cls(k,jcls,lc)
            chem_sub_beg(k,1:2,jcls,la) = tmpa
            chem_sub_beg(k,1:2,jcls,lc) = 0.0
            chem_cls(k,jcls,la) = tmpa
            chem_cls(k,jcls,lc) = 0.0
        end do ! la
        cycle acwxx1_k_loop
        end if

        ! at this point, clear and cloudy subareas > 0 
        tmp_acen(0:2) = acen_tbeg_use(k,0:2,jcls)
        tmp_chem_cls(p1st:num_chem_ecpp) = chem_cls(k,jcls,p1st:num_chem_ecpp)
        tmp_chem_sub(1:2,p1st:num_chem_ecpp) = chem_sub_beg(k,1:2,jcls,p1st:num_chem_ecpp)

        if (lun164 > 0) &
        write(lun164,'(/a,8i5)') 'aa ktau,jcls,k         ', ktau,jcls,k
        call parampollu_tdx_partition_acw( &
        tmp_acen, tmp_chem_cls, tmp_chem_sub, &
        ktau, it, jt, k, jcls, lun164 )

        chem_sub_beg(k,1:2,jcls,p1st:num_chem_ecpp) = tmp_chem_sub(1:2,p1st:num_chem_ecpp)

    end do acwxx1_k_loop
    end do acwxx1_jcls_loop

    end if ! ((activate_onoff_use > 0) .and. (iflag_ecpp_startup_acw_partition > 0))

    if ((lun161 > 0) .and. (kts > -1)) then
!   la = p_num_a03 ; lc = p_num_cw03
!   write(lun161,'(/a,4i6)') 'startup - ktau, l_num_ac03', ktau, la, lc, laicwpair_of_aerosol(la)
!   la = p_oin_a03 ; lc = p_oin_cw03
!   if (lun162 > 0) &
!      write(lun162,'(/a,4i6)') 'startup - ktau, l_oin_ac03', ktau, la, lc, laicwpair_of_aerosol(la)
!   do k = min(10,ktecen), kts, -1

!       write(lun161,'(i2,2(1x,2l1),2(2x,  2x,2(2x,2f11.8)))') k, &
!       (( (acen_tbeg_use(k,icc,jcls)>afrac_cut_0p5), icc=1,2 ), jcls=1,2 ), &
!       (( acen_tbeg_use(k,icc,jcls), icc=1,2 ), jcls=1,2 )

!       la = p_num_a01 ; lc = p_num_cw01 ; tmpa = 1.0e-9
!       la = p_num_a03 ; lc = p_num_cw03 ; tmpa = 1.0e-6
!       write(lun161,'(i2,  1x,a5,   2(3x,f6.3,2(1x,3f6.3)))') k, 'num_3', &
!       ( tmpa*chem_bar(k,l), &
!       ( tmpa*chem_cls(k,jcls,l), tmpa*chem_sub_beg(k,1:2,jcls,l), jcls=1,2 ), &
!       l=la,lc,lc-la )
!       la = p_oin_a01 ; lc = p_oin_cw01 ; tmpa = 1.0
!       la = p_oin_a03 ; lc = p_oin_cw03 ; tmpa = 1.0
!       write(lun161,'(i2,  1x,a5,   2(3x,f6.3,2(1x,3f6.3)))') k, 'oin_3', &
!       ( tmpa*chem_bar(k,l), &
!       ( tmpa*chem_cls(k,jcls,l), tmpa*chem_sub_beg(k,1:2,jcls,l), jcls=1,2 ), &
!       l=la,lc,lc-la )

!   end do
    end if ! ((lun161 > 0) .and. (kts > -1))


    return
    end subroutine parampollu_tdx_startup


!-----------------------------------------------------------------------
    subroutine parampollu_tdx_partition_acw( &
        acen, chem_cls, chem_sub, &
        ktau, i, j, k, jcls, lun164 )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_partition_acw paritions interstitial ("a") and
!    activate/cloudborne ("cw") aerosol species to the clear and cloudy
!    fractions of a grid cell (or grid cell transport-class)
!
!-----------------------------------------------------------------------

    use module_data_mosaic_asect, only:  ai_phase, cw_phase, &
        ncomp_aer, nsize_aer, ntype_aer, &
        massptr_aer, numptr_aer,   & !waterptr_aer, &
        dens_aer, volumlo_sect, volumhi_sect

    use module_data_ecpp1

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message,   &
                                 parampollu_1clm_set_opts

!   arguments
    integer, intent(in) ::                  &
        ktau, i, j, k, jcls, lun164
!   ktau - time step number
!   [i, k, j] - spatial (x,z,y) indices for grid cell

    real(r8), intent(in), dimension( 0:2 ) :: acen

    real(r8), intent(in), dimension( 1:num_chem_ecpp ) :: chem_cls

    real(r8), intent(inout), dimension( 1:2, 1:num_chem_ecpp ) :: chem_sub


!   local variables
    integer :: iphase, isize, itmpa, itype
    integer :: la, lc, ll

    real(r8) :: fx, fy
    real(r8) :: q_a_x, q_a_y, q_a_bar, &
                q_c_x, q_c_y, q_c_bar, &
                q_ac_x, q_ac_y, q_ac_bar, &
                qn_a_x, qn_a_y, qn_a_bar, &
                qn_a_x_sv, qn_a_y_sv, &
                qv_a_x, qv_a_y, qv_a_bar
    real(r8) :: tmpa

    character(len=120) :: msg



    if (min(acen(1),acen(2)) < afrac_cut_0p5) then
        write(msg,'(a,i10,3i5,1p,2e12.4)') &
        '*** parampollu_tdx_partition_acw - bad acen(1:2)', &
        ktau, i, j, k, acen(1:2)
        call ecpp_message( lunout, msg )
        call ecpp_error_fatal( lunout, msg )
        return
    end if
    fy = acen(2)/(acen(1)+acen(2))
    fx = 1.0_r8 - fy

! main loops over aerosol types and sizes
    do itype = 1, ntype_aer
    do isize = 1, nsize_aer(itype)

! first partition number and dry-mass species
!    in a manner that attempts to get the "a+cw" mixing ratios
!    in clear and cloudy subareas to be equal the 
!    cell/class average (clear+cloudy) "a+cw" mixing ratios
        qv_a_x = 0.0 ; qv_a_y = 0.0
        do ll = 0, ncomp_aer(itype)
        if (ll == 0) then
            la = numptr_aer(isize,itype,ai_phase)
            lc = numptr_aer(isize,itype,cw_phase)
        else
            la = massptr_aer(ll,isize,itype,ai_phase)
            lc = massptr_aer(ll,isize,itype,cw_phase)
        end if

! nomenclature for q_...
!    a = interstitial; c = cloudborne; ac = a+c
!    x = in clear subarea; y = in cloudy subarea; 
!        bar = average over both subareas
!
! following always hold
!    q_ac_any == q_a_any + q_c_any
!    q_any_bar == q_any_x*fx + q_any_y*fy
!
        q_a_bar = max( 0.0_r8, chem_cls(la) )
        q_c_bar = max( 0.0_r8, chem_cls(lc) )
        q_ac_bar = q_a_bar + q_c_bar
        q_c_y = q_c_bar/fy
        q_c_x = 0.0
        q_a_y = max( 0.0_r8, (q_ac_bar - q_c_y) )
        q_a_x = max( 0.0_r8, (q_a_bar - q_a_y*fy)/fx )

!       if ((k <= 5) .and. (isize == 1) .and. (ll == 3)) then
        if ((k <= 5) .and. (isize == 3) .and. (ll==3 .or. ll==0)) then
          if (lun164 > 0) then
            write(lun164,'(/a,8i5)') 'bb ktau,jcls,k,isize,ll', ktau,jcls,k,isize,ll
            write(lun164,'(a,1p,8e12.4)') 'acen1/2, fx/y', acen(1:2), fx, fy
            write(lun164,'(a,1p,8e12.4)') 'chem_cls     ', chem_cls(la), chem_cls(lc)
            write(lun164,'(a,1p,8e12.4)') 'chem_sub old ', chem_sub(1:2,la), chem_sub(1:2,lc)
          end if
        end if
        chem_sub(1,la) = q_a_x
        chem_sub(2,la) = q_a_y
        chem_sub(1,lc) = q_c_x
        chem_sub(2,lc) = q_c_y
!       if ((k <= 5) .and. (isize == 1) .and. (ll == 3)) then
        if ((k <= 5) .and. (isize == 3) .and. (ll==3 .or. ll==0)) then
          if (lun164 > 0) &
            write(lun164,'(a,1p,8e12.4)') 'chem_sub new ', chem_sub(1:2,la), chem_sub(1:2,lc)
        end if

        if (ll == 0) then
            qn_a_x = q_a_x
            qn_a_y = q_a_y
        else
            qv_a_x = qv_a_x + q_a_x/dens_aer(ll,itype)
            qv_a_y = qv_a_y + q_a_y/dens_aer(ll,itype)
        end if
        end do
        qv_a_x = qv_a_x*1.0e-6   ! because mass mixratios are ug/kg,
        qv_a_y = qv_a_y*1.0e-6   ! and want volume mixratio in cm3-aerosol/kg

! now check that the partitioning has not produced an out-of-bounds size
!    (size = mean 1-particle volume) for interstitial in clear or cloudy subareas
! if this has occurred, then partition the number differently
        qv_a_bar = qv_a_x*fx + qv_a_y*fy
        qn_a_bar = qn_a_x*fx + qn_a_y*fy
        qn_a_x_sv = qn_a_x ; qn_a_y_sv = qn_a_y
        if ( (qv_a_bar <= 1.0e-30) .or. &
         (qv_a_bar <= qn_a_bar*volumlo_sect(isize,itype)) .or. &
         (qv_a_bar >= qn_a_bar*volumhi_sect(isize,itype)) ) then
            ! neglible dry volume, or size already out-of-bounds
        tmpa = max(qv_a_bar,1.0e-35_r8)
        qn_a_x = qn_a_bar * ( max(qv_a_x,0.5e-35_r8) / tmpa )
        qn_a_y = qn_a_bar * ( max(qv_a_y,0.5e-35_r8) / tmpa )
        if (qv_a_bar <= 1.0e-30) then
            itmpa = 1
        else if (qv_a_bar <= qn_a_bar*volumlo_sect(isize,itype)) then
            itmpa = 2
        else
            itmpa = 3
        end if

        else if (qv_a_x <= qn_a_x*volumlo_sect(isize,itype)) then
            ! size to small in clear subarea
        qn_a_x = qv_a_x/volumlo_sect(isize,itype)
        qn_a_y = max( 0.0_r8, (qn_a_bar - qn_a_x*fx)/fy )
        itmpa = 4
        else if (qv_a_y <= qn_a_y*volumlo_sect(isize,itype)) then
            ! size to small in cloudy subarea
        qn_a_y = qv_a_y/volumlo_sect(isize,itype)
        qn_a_x = max( 0.0_r8, (qn_a_bar - qn_a_y*fy)/fx )
        itmpa = 5

        else if (qv_a_x >= qn_a_x*volumhi_sect(isize,itype)) then
            ! size to large in clear subarea
        qn_a_x = qv_a_x/volumhi_sect(isize,itype)
        qn_a_y = max( 0.0_r8, (qn_a_bar - qn_a_x*fx)/fy )
        itmpa = 6
        else if (qv_a_y >= qn_a_y*volumhi_sect(isize,itype)) then
            ! size to large in cloudy subarea
        qn_a_y = qv_a_y/volumhi_sect(isize,itype)
        qn_a_x = max( 0.0_r8, (qn_a_bar - qn_a_y*fy)/fx )
        itmpa = 7
        else
        itmpa = 0
        end if
        la = numptr_aer(isize,itype,ai_phase)
        chem_sub(1,la) = qn_a_x
        chem_sub(2,la) = qn_a_y
        if ((k <= 5) .and. (isize == 3)) then
        if ((itmpa==5) .and. (qv_a_y>0.0)) itmpa=8
        if ((itmpa==5) .and. (qn_a_y>0.0)) itmpa=9
        if (lun164 > 0) then
            write(lun164,'(/i1,a,1p,8e12.4)') itmpa, ' final num_a', chem_sub(1:2,la)
            write(lun164,'(  13x,1p,8e12.4)') qn_a_x_sv, qn_a_y_sv, qn_a_bar, qv_a_x, qv_a_y
        end if
        end if

! aerosol water - do this for now, but it should be improved
! comment out now, need to check with Dick Easter. +++mhwang
!
!       la = waterptr_aer(isize,itype)
!       tmpa = max(qv_a_bar,1.0e-35)
!       chem_sub(1,la) = ( max(qv_a_x,0.5e-35) / tmpa ) * chem_cls(la)
!       chem_sub(2,la) = ( max(qv_a_y,0.5e-35) / tmpa ) * chem_cls(la)

    end do
    end do



    return
    end subroutine parampollu_tdx_partition_acw

!-----------------------------------------------------------------------
    subroutine parampollu_tdx_cleanup(                        &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        chem_bar, chem_cls,                               &
        ncls_ecpp,                                        &
        acen_tfin_ecpp,                                   &
        it,      jt,      kts,ktebnd,ktecen,              &
        ncls_use,                                         &
        chem_sub_beg, chem_sub_new,                       &
        del_chem_clm_cldchem, del_chem_clm_wetscav,       &
                del_cldchem3d,  del_rename3d,                     &
                del_wetdep3d, del_wetresu3d,      &
                del_activate3d, del_conv3d,                       &
        acen_tbeg_use, acen_tfin_use, rhodz_cen,          &
        activate_onoff_use,                               &
        iphase_of_aerosol, isize_of_aerosol,              &
        itype_of_aerosol, inmw_of_aerosol,                &
        laicwpair_of_aerosol                              )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_cleanup does some final "cleanup" calculations
!
!    calculates final chem_cls and chem_bar from the final chem_sub
!
!    calculates beginning and final column-average mixing ratios
!   and checks for mass conservation
!
!-----------------------------------------------------------------------

    use module_data_mosaic_asect, only:  ai_phase, cw_phase, &
                nsize_aer, massptr_aer, numptr_aer

    use module_data_radm2, only:  epsilc

    use module_data_ecpp1

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message

!   arguments
    integer, intent(in) ::                  &
        ktau, ktau_pp,              &
        it, jt, kts, ktebnd, ktecen
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!   chem_driver and routines under it do calculations
!   over these spatial indices.

    integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

    real(r8), intent(in) :: dtstep, dtstep_pp
!   dtstep - main model time step (s)
!   dtstep_pp - time step (s) for "parameterized pollutants" calculations

    real(r8), intent(inout), dimension( kts:ktecen, 1:num_chem_ecpp ) :: &
        chem_bar
!   chem_bar - mixing ratios of trace gase (ppm) and aerosol species
!   (ug/kg for mass species, #/kg for number species)

    real(r8), intent(inout), dimension( kts:ktecen, 1:maxcls_ecpp, 1:num_chem_ecpp ) :: &
        chem_cls
 
    integer, intent(in) :: ncls_ecpp
!   ncls_ecpp - number of ecpp transport classes in the grid column
    

    real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
        acen_tfin_ecpp

    integer, intent(in) :: ncls_use

    real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
        chem_sub_beg, chem_sub_new

    real(r8), intent(inout), dimension( 1:num_chem_ecpp ) ::   &
        del_chem_clm_cldchem, del_chem_clm_wetscav

        real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:num_chem_ecpp ) ::   &
                del_cldchem3d, del_rename3d, del_wetdep3d, del_wetresu3d

        real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
                del_activate3d   

        real(r8), intent(out), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
                del_conv3d

    real(r8), intent(in), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) :: &
        acen_tbeg_use, acen_tfin_use

    real(r8), intent(in), dimension( kts:ktecen ) :: rhodz_cen

    integer, intent(in) :: activate_onoff_use

    integer, intent(in), dimension( 1:num_chem_ecpp ) ::   &
        iphase_of_aerosol, isize_of_aerosol, itype_of_aerosol,   &
        inmw_of_aerosol, laicwpair_of_aerosol



!   local variables
    integer :: ia, ib, icc
    integer :: jcls, jclsbb
    integer :: k
    integer :: l, la, laa, lbb, lc, lewa, lewc, lun119, lun121
    integer :: laicwpair_flagaa
    integer, save :: ktaueww = 0

    real(r8) :: air_clmmass
    real(r8) :: chem_cutoff_aa
    real(r8) :: tmpa, tmpb, tmpe, tmpew, tmpx, tmpy, tmpz
    real(r8) :: tmpa_clmavg(1:6), tmpw_clmavg(1:6)
    real(r8) :: tmpveca( kts:ktecen ), tmpvecb( kts:ktecen )
    real(r8) :: tmpvece(1:6)
    real(r8), save :: tmpeww = 0.0

    real(r8), dimension( 1:6, 1:num_chem_ecpp ) :: chem_clmavg
    real(r8), dimension( kts:ktecen, 1:num_chem_ecpp ) :: chem_bar_beg



    lun121 = -1
    if (idiagaa_ecpp(121) > 0) lun121 = ldiagaa_ecpp(121)

        del_conv3d = 0.0
!   calculate initial clmmass  and del_conv3d
    air_clmmass = sum( rhodz_cen(kts:ktecen) )
    do l = param_first_ecpp, num_chem_ecpp
       tmpveca(:) = 0.0 ; tmpvecb(:) = 0.0
       do jcls = 1, ncls_use
       do icc = 1, 2
       do k = kts, ktecen
          tmpveca(k) = tmpveca(k) + acen_tbeg_use(k,icc,jcls)*chem_cls(    k,    jcls,l)
          tmpvecb(k) = tmpvecb(k) + acen_tbeg_use(k,icc,jcls)*chem_sub_beg(k,icc,jcls,l)
       end do
       end do
       end do
       chem_clmavg(1,l) = sum( rhodz_cen(kts:ktecen)*chem_bar(kts:ktecen,l) )
       chem_clmavg(2,l) = sum( rhodz_cen(kts:ktecen)*tmpveca(kts:ktecen) )
       chem_clmavg(3,l) = sum( rhodz_cen(kts:ktecen)*tmpvecb(kts:ktecen) )
    end do
    if ((ktau < 0) .and. (lun121 > 0)) then
        l = 17
        icc = 1
!       write(lun121,*) 'ktau, l, ncls_use', ktau, l, ncls_use
!       write(lun121,*) 'k, old chem_bar, old chem_cls, chem_sub_beg, acen_tbeg_use'
!       do k = ktecen, kts, -1
!       write(lun121,'(i3,1p,e12.5,3(3x,3e12.5))') k, chem_bar(k,l),   &
!           chem_cls(k,1:3,l), chem_sub_beg(k,icc,1:3,l), acen_tbeg_use(k,icc,1:3)
!       end do
    end if
!   if (ktau > 1) stop


!   do acen_tfin_ecpp <-- acen_tfin_use
    acen_tfin_ecpp(:,:,:) = acen_tfin_use(:,:,:)


! compute new chem_cls (class-avg mix ratios) and chem_bar (grid-avg mix ratios)
    chem_bar_beg(:,:) = chem_bar(:,:)
    do l = param_first_ecpp, num_chem_ecpp
    do k = kts, ktecen

        tmpa = 0.0 ; tmpb = 0.0
        do jcls = 1, ncls_use
        do icc = 1, 2
            tmpa = tmpa + acen_tfin_use(k,icc,jcls)* &
                             max(0.0_r8,chem_sub_new(k,icc,jcls,l))
            tmpb = tmpb + acen_tfin_use(k,icc,jcls)

                    del_conv3d(k,icc,jcls,l) = (acen_tfin_use(k,icc,jcls)*max(0.0_r8, chem_sub_new(k,icc,jcls,l))   &
                                               - acen_tbeg_use(k,icc,jcls)*chem_sub_beg(k,icc,jcls,l))              &
                                               - del_activate3d(k,icc,jcls,l)                                       &
                                               - del_cldchem3d(k,icc,jcls,1,l)-del_cldchem3d(k,icc,jcls,2,l)        &
                                               - del_rename3d(k,icc,jcls,1,l)-del_rename3d(k,icc,jcls,2,l)          &
                                               - del_wetdep3d(k,icc,jcls,1,l)-del_wetdep3d(k,icc,jcls,2,l)          &
                                               - del_wetresu3d(k,icc,jcls,1,l)-del_wetresu3d(k,icc,jcls,2,l)
        end do
        end do
!       chem_bar(k,l) = max(0.0_r8,tmpa)/tmpb
            chem_bar(k,l) = tmpa   ! chem_bar is used to calcualte q tendency at the MMF model, 
                                   ! so keep it consistent with del_conv3d

        do jcls = 1, ncls_use
        tmpa = 0.0 ; tmpb = 0.0
        do icc = 1, 2
            tmpa = tmpa + acen_tfin_use(k,icc,jcls)* &
                             max(0.0_r8,chem_sub_new(k,icc,jcls,l))
            tmpb = tmpb + acen_tfin_use(k,icc,jcls)
        end do
        if (tmpb >= afrac_cut_0p5) then
            chem_cls(k,jcls,l) = max(0.0_r8,tmpa)/tmpb
        else
            chem_cls(k,jcls,l) = chem_bar(k,l)
        end if
        end do

    end do
    end do


!   calculate final clmmass
    do l = param_first_ecpp, num_chem_ecpp
       tmpveca(:) = 0.0 ; tmpvecb(:) = 0.0
       do jcls = 1, ncls_use
       do icc = 1, 2
       do k = kts, ktecen
          tmpveca(k) = tmpveca(k) + acen_tfin_use(k,icc,jcls)*chem_cls(    k,    jcls,l)
          tmpvecb(k) = tmpvecb(k) + acen_tfin_use(k,icc,jcls)*chem_sub_new(k,icc,jcls,l)
       end do
       end do
       end do
       chem_clmavg(4,l) = sum( rhodz_cen(kts:ktecen)*chem_bar(kts:ktecen,l) )
       chem_clmavg(5,l) = sum( rhodz_cen(kts:ktecen)*tmpveca(kts:ktecen) )
       chem_clmavg(6,l) = sum( rhodz_cen(kts:ktecen)*tmpvecb(kts:ktecen) )
       chem_clmavg(1:6,l) = chem_clmavg(1:6,l)/air_clmmass
    end do
    if ((ktau < 0) .and. (lun121 > 0)) then
        l = 17
        icc = 1
!       write(lun121,*) 'ktau, l, ncls_use', ktau, l, ncls_use
!       write(lun121,*) 'k, new chem_bar, new chem_cls, chem_sub_new, acen_tfin_use'
        do k = ktecen, kts, -1
!       write(lun121,'(i3,1p,e12.5,3(3x,3e12.5))') k, chem_bar(k,l),   &
!           chem_cls(k,1:3,l), chem_sub_new(k,icc,1:3,l), acen_tfin_use(k,icc,1:3)
        end do
    end if
!   if (ktau > 5) stop
    if ((ktau < 5) .and. (lun121 > 0)) then
        l = 9
!       write(lun121,'(/a,3i5)') 'ktau, l, ncls_use ', ktau, l, ncls_use
!       write(lun121,'(a)') 'k, ((chem_sub_beg(k,icc,jcls,l), chem_sub_new(k,icc,jcls,l), icc=1,2), jcls=1,...) '
        do k = ktecen, kts, -1
!       write(lun121,'(i3,1p,6(2x,2e10.3))') k,   &
!           ((chem_sub_beg(k,icc,jcls,l), chem_sub_new(k,icc,jcls,l), icc=1,2), jcls=1,ncls_use)
        end do
    end if


!   diagnostic output to unit 121
    if (lun121 > 0) then

!   write(lun121,'(/a,2i6)') 'parampollu_1clm clmmass check - ktau, ktau_pp =',   &
!       ktau, ktau_pp
    lewa = 0
    lewc = 0
    tmpew = 0.0
    chem_cutoff_aa = 3.0*epsilc
    laicwpair_flagaa = 0
    if ( (activate_onoff_use > 0) .and.   &
         (activate_onoff_use /=100) ) laicwpair_flagaa = 2
    do la = param_first_ecpp, num_chem_ecpp
        l = -999888777
        lc = 0
        if (laicwpair_flagaa == 2) then
        if (iphase_of_aerosol(la) == ai_phase) then
            lc = laicwpair_of_aerosol(la)
        else if (iphase_of_aerosol(la) == cw_phase) then
            cycle
        end if
        end if
        if ((lc < param_first_ecpp) .or. (lc > num_chem_ecpp)) lc = 0

        ! these are the 3 initial and 3 final values of column-average mixing ratio
        ! for the current species (or species la-lc pair)
        tmpa_clmavg(1:6) = chem_clmavg(1:6,la)
        if (lc > 0) tmpa_clmavg(1:6) = tmpa_clmavg(1:6) + chem_clmavg(1:6,lc)

        ! for the 3 final values, subtract off the change from cldchem and wetscav
        tmpa = del_chem_clm_cldchem(la) + del_chem_clm_wetscav(la)
        if (lc > 0) tmpa = tmpa + del_chem_clm_cldchem(lc) + del_chem_clm_wetscav(lc)
        tmpa = tmpa/air_clmmass
        tmpa_clmavg(4:6) = tmpa_clmavg(4:6) - tmpa

        do ia = 1, 6
        ib = mod(ia,6) + 1
        tmpa = tmpa_clmavg(ia)
        tmpb = tmpa_clmavg(ib)
        tmpvece(ia) = abs( tmpa-tmpb )   &
                    / max( abs(tmpa), abs(tmpb), 1.0e-30_r8 )
        end do
        tmpx = maxval( tmpa_clmavg(1:6) )
        tmpy = minval( tmpa_clmavg(1:6) )
        tmpz = max( abs(tmpx), abs(tmpy), 1.0e-30_r8 )
        ! ignore species with max,min( clmavg mixratios ) < chem_cutoff_aa
        if (tmpz >= chem_cutoff_aa) then
        tmpe = abs( tmpx-tmpy ) / tmpz
        else
        tmpe = 0.0
        end if
        if (tmpe > tmpew) then
        tmpew = tmpe
        lewa = la
        lewc = lc
        tmpw_clmavg(:) = tmpa_clmavg(:)
        end if

        if (tmpe > 1.0e-12_r8 ) then
        write(lun121,'(a,2i3,1p,2(3x,6e10.2))') 'la/c=', la, lc,   &
            tmpa_clmavg(1:6), tmpvece(1:6)

                write(0,'(a,2i3,1p,2(3x,6e10.2))') 'mass convervation error in ecpp, la/c=', la, lc,   &
                    tmpa_clmavg(1:6), tmpvece(1:6)
                call endrun('mass convervation error in ecpp_cleanup')
        end if
    end do
    if (tmpew > tmpeww) then
        tmpeww = tmpew
        ktaueww = ktau
    end if
    if (lewa > 0) then
        write(lun121,'(a,2i3,1p,e10.2,10x,2i6,e10.2)') 'worst clmmass error - la/c=',   &
        lewa, lewc, tmpew, ktau, ktaueww, tmpeww
        write(lun121,'(a,1p,6e14.6)') 'chem_clmavg(1:6,l)', tmpw_clmavg(1:6)
    end if

    end if   ! (lun121 > 0)


!   diagnostic output to unit 119
    lun119 = -1
    if (idiagaa_ecpp(119) > 0) lun119 = ldiagaa_ecpp(119)
    if (lun119 > 0) then
    write(lun119,'(/a,2i5)') 'parampollu_1clm - pt2 ktau, ktau_pp =', ktau, ktau_pp
!   do laa = param_first_ecpp, num_chem_ecpp, 3
!       lbb = min( laa+2, num_chem_ecpp )
!   do laa = param_first_ecpp, num_chem_ecpp, 4
!       lbb = min( laa+3, num_chem_ecpp )
    do laa = 9, 9
        lbb = min( laa+3, num_chem_ecpp )
        write(lun119,'(/a,4i5)') 'ktau, ktau_pp, laa, lbb =', ktau, ktau_pp, laa, lbb
        write(lun119,'(a)') ' k,  chem_bar_beg, chem_bar'
        do k = ktecen, kts, -1
!       write(lun119,'(i2,4(2x,2f9.5))') k,   &
        write(lun119,'(i2,4(2x,1p,2e10.2))') k,   &
            (chem_bar_beg(k,l), chem_bar(k,l), l=laa, lbb)
        end do
!       write(lun119,'(i2,4(2x,2f9.5))') -1,  &
        write(lun119,'(i2,4(2x,1p,2e10.2))') -1,   &
        (chem_clmavg(2,l), chem_clmavg(5,l), l=laa, lbb)
!       write(lun119,'(i2,1p,4e20.5))') -2,  &
        write(lun119,'(i2,4(2x,1p,e20.2))') -2,   &
        ( (chem_clmavg(2,l)-chem_clmavg(5,l)), l=laa, lbb)
    end do
    end if   ! (lun119 > 0)


    return
    end subroutine parampollu_tdx_cleanup



!-----------------------------------------------------------------------
    subroutine parampollu_check_adjust_inputs(                &
        ipass_check_adjust_inputs,                        &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        idiagaa_ecpp, ldiagaa_ecpp,                       &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, zbnd, wbnd_bar,                       &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        kdraft_bot_ecpp, kdraft_top_ecpp,                 &
        mtype_updnenv_ecpp,                               &
        mfbnd_ecpp, abnd_tavg_ecpp,                       &
        acen_tavg_ecpp, acen_tfin_ecpp, acen_prec_ecpp,   &
        wbnd_bar_use,                                     &
        ncls_use,                                         &
        kdraft_bot_use, kdraft_top_use,                   &
        mtype_updnenv_use,                                &
        mfbnd_use, mfbnd_quiescn_up, mfbnd_quiescn_dn,    &
        abnd_tavg_use,                                    &
        acen_tavg_use, acen_tfin_use, acen_prec_use,      &
        rhodz_cen,                                        &
        it,      jt,      kts,ktebnd,ktecen               )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_check_adjust_inputs does checking and adjustment
!    of several of the ecpp arrays
!
!    fractional areas less than afrac_cut are set to zero
!    up and downdraft mass fluxes less than ... are set to zero
!    remaining fractional areas are adjusted so that the sum is 1.0
!
!    all mass fluxes are set to zero at/above k_max_wnonzero
!    up and downdraft mass fluxes and areas are set to zero at/above k_max_updndraft
!    cloud fractional areas are set to zero at/above k_max_clouds
!    
! the checks and adjustment are designed to eliminate "problems" in
!    the input/incoming arrays that might cause the rest of the
!    parampollu code to fail
!
!-----------------------------------------------------------------------

    use module_data_ecpp1

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message

!   arguments
    integer, intent(in) ::                  &
        ipass_check_adjust_inputs,      &
        ktau, ktau_pp,              &
        it, jt, kts, ktebnd, ktecen
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations

!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!   chem_driver and routines under it do calculations
!   over these spatial indices.

    integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

    real(r8), intent(in) :: dtstep, dtstep_pp
!   dtstep - main model time step (s)
!   dtstep_pp - time step (s) for "parameterized pollutants" calculations

    real(r8), intent(in), dimension( kts:ktecen ) ::   &
        tcen_bar, pcen_bar, rhocen_bar, dzcen
    real(r8), intent(in), dimension( kts:ktebnd ) ::   &
        rhobnd_bar, zbnd, wbnd_bar

    real(r8), intent(inout), dimension( kts:ktebnd ) ::   &
        wbnd_bar_use

    real(r8), intent(inout), dimension( kts:ktecen, 1:num_chem_ecpp ) :: &
        chem_bar

    integer, intent(in) :: ncls_ecpp
    integer, intent(inout) :: ncls_use
    
    integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        kdraft_bot_ecpp, kdraft_top_ecpp,   &
        mtype_updnenv_ecpp
    integer, intent(inout), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        kdraft_bot_use, kdraft_top_use,   &
        mtype_updnenv_use
    
    real(r8), intent(in), dimension( kts:ktebnd, 0:2, 0:maxcls_ecpp ) ::   &
        mfbnd_ecpp, abnd_tavg_ecpp 
        real(r8), intent(inout), dimension( kts:ktebnd, 0:2, 0:maxcls_ecpp ) ::   &
                mfbnd_use, abnd_tavg_use
    real(r8), intent(in), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
        acen_tavg_ecpp, acen_tfin_ecpp, acen_prec_ecpp
        real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
        acen_tavg_use,  acen_tfin_use,  acen_prec_use
    real(r8), intent(inout), dimension( kts:ktebnd, 0:2, 0:2 ) ::   &
        mfbnd_quiescn_up, mfbnd_quiescn_dn
    real(r8), intent(inout), dimension( kts:ktecen ) :: rhodz_cen


!   local variables
    integer :: k_max_updndraft
    integer :: k_max_clouds
    integer :: k_max_wnonzero

    integer :: i, icc, itmpa, itmpb
    integer :: ido_downdr_area_zeroout, ido_updndr_area_adjust, ipass_2_changes
    integer :: ispecial_check_acen_tfin
    integer :: ja, jb
    integer :: jcls, jclsbb 
    integer :: jclsicc, jclsicc_noc, jclsicc_cld
    integer :: k, ka, kb, ktmpa, ktmpb
    integer :: lun63, lun141, lun155
    integer :: ncls_noc, ncls_cld
    integer :: nchanges(10)
    integer :: kdraft_bot_tmp(1:2,1:maxcls_ecpp), kdraft_top_tmp(1:2,1:maxcls_ecpp)
    integer :: mtype_updnenv_tmp(1:2,1:maxcls_ecpp)

    real(r8) :: ardz_cut   ! sub-class fractional areas below this value are set to zero
    real(r8) :: arw_draft_cut   ! mass fluxes below this value are set to zero
    real(r8) :: a_sum_toleraa = 1.0e-5_r8   ! tolerance for abs(sum(axxx) - 1.0)
    real(r8) :: afrac_noc, afrac_cld
    real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf, tmpq, tmpu
    real(r8) :: tmp_afrac
    real(r8) :: tmp_mfa, tmp_mfb
    real(r8) :: tmp_tola, tmp_tolb
    real(r8) :: tmpvecaa(0:ktebnd), tmpvecbb(0:ktebnd), tmpvecdd(0:ktebnd)
    real(r8) :: tmp0202aa(0:2,0:2)
    real(r8) :: updndr_area_adjust

    character(len=100) :: msg
    character(len=10) :: area_name10(1:3) = &
                         (/ 'abnd_tavg ', 'acen_tavg ', 'acen_tfin ' /)


    lun63 = -1
    if (idiagaa_ecpp(63) > 0) lun63 = ldiagaa_ecpp(63)
    lun141 = -1
    if (idiagaa_ecpp(141) > 0) lun141 = ldiagaa_ecpp(141)
    lun155 = -1
    if (idiagaa_ecpp(155) > 0) lun155 = ldiagaa_ecpp(155)

    if ((ipass_check_adjust_inputs /= 1) .and.   &
        (ipass_check_adjust_inputs /= 2)) return


!   force w = 0 at kbnd >= k_max_wnonzero
!   (note - doing k_max_wnonzero = ktebnd-1 would probably be ok)
    k_max_wnonzero = ktebnd-1

!   force up/dn draft mf & afrac = 0 at kbnd,kcen >= k_max_updndraft
!   (note - currently set k_max_updndraft & _kclouds to almost top of domain)
    k_max_updndraft = ktebnd-1

!   force cloud fraction = 0 at kbnd,kcen >= k_max_clouds
    k_max_clouds = ktebnd-1

    nchanges(:) = 0


!-----------------------------------------------------
!   when ipass_check_adjust_inputs == 2, 
!   skip to he beginning of the special stuff for ipass_check_adjust_inputs == 2
    if (ipass_check_adjust_inputs == 2) goto 20000
!-----------------------------------------------------


!
!   copy from "_ecpp" arrays to "_use" arrays
!
    ncls_use = ncls_ecpp

    kdraft_bot_use(:,:) = kdraft_bot_ecpp(:,:)
    kdraft_top_use(:,:) = kdraft_top_ecpp(:,:)

    mtype_updnenv_use(:,:) = mtype_updnenv_ecpp(:,:)

    wbnd_bar_use(:) = wbnd_bar(:)

    mfbnd_use(:,:,:) = mfbnd_ecpp(:,:,:)
    abnd_tavg_use(:,:,:) = max( abnd_tavg_ecpp(:,:,:), 0.0_r8 )
    acen_tavg_use(:,:,:) = max( acen_tavg_ecpp(:,:,:), 0.0_r8 )
    acen_tfin_use(:,:,:) = max( acen_tfin_ecpp(:,:,:), 0.0_r8 )
!   acen_tavg_use(kte,:,:) = 0.0
!   acen_tfin_use(kte,:,:) = 0.0

!   calc rhodz_cen
    rhodz_cen(kts:ktecen) = rhocen_bar(kts:ktecen)*dzcen(kts:ktecen)


!   check that 
!   the mtype_updnenv_use are valid
!   there is exactly one of each quiescent transport class (cloudy, clear)
    jclsicc_noc = -1
    jclsicc_cld = -1
    ncls_noc = 0
    ncls_cld = 0
    msg = ' '

    do jcls = 1, ncls_use
    do icc = 1, 2
        jclsicc = jcls*10 + icc
        if ((mtype_updnenv_use(icc,jcls) == mtype_quiescn_ecpp) .and.   &
            (icc == 1)) then
        jclsicc_noc = jclsicc
        ncls_noc = ncls_noc + 1
        end if
        if ((mtype_updnenv_use(icc,jcls) == mtype_quiescn_ecpp) .and.   &
            (icc == 2)) then
        jclsicc_cld = jclsicc
        ncls_cld = ncls_cld + 1
        end if

        if ( ((jcls == jcls_qu) .and. &
              (mtype_updnenv_use(icc,jcls) /= mtype_quiescn_ecpp)) .or. &
             ((jcls /= jcls_qu) .and. &
              (mtype_updnenv_use(icc,jcls) /= mtype_updraft_ecpp) .and. &
              (mtype_updnenv_use(icc,jcls) /= mtype_dndraft_ecpp)) ) then
        write( msg, '(a,5(1x,i5))' ) &
        '*** parampollu_check_adjust_inputs - bad mtype_updnenv', &
        it, jt, jcls, icc, mtype_updnenv_use(icc,jcls)
        call ecpp_message( lunout, msg )
        end if
    end do
    end do

    if ((jclsicc_noc <= 0) .or. (ncls_noc > 1)) then
        write(msg,'(a,2(1x,i5))')   &
        '*** parampollu_check_adjust_inputs - bad jclsicc_noc, ncls_noc =',   &
        jclsicc_noc, ncls_noc
        call ecpp_message( lunout, msg )
    end if
    if ((jclsicc_cld <= 0) .or. (ncls_cld > 1)) then
        write(msg,'(a,2(1x,i5))')   &
        '*** parampollu_check_adjust_inputs - bad jclsicc_cld, ncls_cld =',   &
        jclsicc_cld, ncls_cld
        call ecpp_message( lunout, msg )
    end if
    if (msg /= ' ') call ecpp_error_fatal( lunout, msg )


    if ((ktau==4) .and. (lun155 > 0)) then
    write(lun155,'(/a,3i5)') 'aaa', ktau, ipass_check_adjust_inputs
    write(lun155,'(3(i5,i3,1pe16.8))') ((jcls,icc,acen_tavg_use(26,icc,jcls),icc=0,2),jcls=0,3)
    end if
!   *** this is for testing 
!       when iflag_ecpp_test_fixed_fcloud == 2/3/4/5, 
!           set clear  fractions to 1.0/0.0/0.7/0.3
!           set cloudy fractions to 0.0/1.0/0.3/0.7
!
!   *** also set k_max_clouds=kte+1 so that it has no effect
!
    if ((iflag_ecpp_test_fixed_fcloud >= 2) .and. &
        (iflag_ecpp_test_fixed_fcloud <= 5)) then
        k_max_clouds = ktebnd+1

        if      (iflag_ecpp_test_fixed_fcloud == 2) then
        tmpvecaa(1) = 1.0_r8
        else if (iflag_ecpp_test_fixed_fcloud == 3) then
        tmpvecaa(1) = 0.0_r8
        else if (iflag_ecpp_test_fixed_fcloud == 4) then
        tmpvecaa(1) = 0.7_r8
        else
        tmpvecaa(1) = 0.3_r8
        end if
        tmpvecaa(2) = 1.0_r8 - tmpvecaa(1)

        do k = kts, ktebnd
        do jcls = 1, ncls_use
        tmpa = sum( mfbnd_use(k,1:2,jcls) )
        mfbnd_use(k,1:2,jcls) = tmpa*tmpvecaa(1:2)

        tmpa = sum( abnd_tavg_use(k,1:2,jcls) )
        abnd_tavg_use(k,1:2,jcls) = tmpa*tmpvecaa(1:2)

        if (k > ktecen) cycle

        tmpa = sum( acen_tavg_use(k,1:2,jcls) )
        acen_tavg_use(k,1:2,jcls) = tmpa*tmpvecaa(1:2)

        tmpa = sum( acen_tfin_use(k,1:2,jcls) )
        acen_tfin_use(k,1:2,jcls) = tmpa*tmpvecaa(1:2)
        end do   ! jcls
        end do   ! k
    end if   ! ((iflag_ecpp_test_fixed_fcloud >= 2) .and. (iflag_ecpp_test_fixed_fcloud <= 5))


!   check that fractional areas sum to 1.0 (within small tolerance)
!   then normalize to exactly 1.0
!   also check and total quiescent areas are each >= a_quiescn_minaa
    do k = kts, ktebnd
        do jcls = 1, ncls_use
            abnd_tavg_use(k,0,jcls) = sum( abnd_tavg_use(k,1:2,jcls) )
            if (k > ktecen) cycle
            acen_tavg_use(k,0,jcls) = sum( acen_tavg_use(k,1:2,jcls) )
            acen_tfin_use(k,0,jcls) = sum( acen_tfin_use(k,1:2,jcls) )
        end do
        do icc = 0, 2
            abnd_tavg_use(k,icc,0) = sum( abnd_tavg_use(k,icc,1:ncls_use) )
            if (k > ktecen) cycle
            acen_tavg_use(k,icc,0) = sum( acen_tavg_use(k,icc,1:ncls_use) )
            acen_tfin_use(k,icc,0) = sum( acen_tfin_use(k,icc,1:ncls_use) )
        end do

        do i = 1, 3
            if ((i >= 2) .and. (k > ktecen)) cycle
            if (i == 1) then
                tmpa = abnd_tavg_use(k,0,0)
            else if (i == 2) then
                tmpa = acen_tavg_use(k,0,0)
            else
                tmpa = acen_tfin_use(k,0,0)
            end if
            if (abs(tmpa-1.0_r8) < a_sum_toleraa) cycle
            write(msg,'(2a,i5,1pe15.7)') &
                '*** parampollu_check_adjust_inputs - bad ', &
                area_name10(i), k, tmpa  
            call ecpp_message( lunout, msg )
            call ecpp_error_fatal( lunout, msg )
        end do

        tmpa = abnd_tavg_use(k,0,0)
        abnd_tavg_use(k,0:2,0:ncls_use) = abnd_tavg_use(k,0:2,0:ncls_use)/tmpa
        if (k <= ktecen) then
            tmpa = acen_tavg_use(k,0,0)
            acen_tavg_use(k,0:2,0:ncls_use) = acen_tavg_use(k,0:2,0:ncls_use)/tmpa
            tmpa = acen_tfin_use(k,0,0)
            acen_tfin_use(k,0:2,0:ncls_use) = acen_tfin_use(k,0:2,0:ncls_use)/tmpa
        end if

        do i = 1, 3
            if ((i >= 2) .and. (k > ktecen)) cycle
            jcls = jcls_qu
            if (i == 1) then
                tmpa = abnd_tavg_use(k,0,jcls)
            else if (i == 2) then
                tmpa = acen_tavg_use(k,0,jcls)
            else
                tmpa = acen_tfin_use(k,0,jcls)
            end if
            msg = ' '
            if (tmpa < a_quiescn_minaa) then
                write(msg,'(2a,i5,1p,2e10.2)') &
                '*** parampollu_check_adjust_inputs - a_quiescent(v1) too small ', &
                area_name10(i), k, tmpa, a_quiescn_minaa
                call ecpp_message( lunout, msg )
                call ecpp_error_fatal( lunout, msg )
            end if
        end do

    end do


!   eliminate cloudy subareas when k >= k_max_clouds
    do k = kts, ktebnd
        if (k < k_max_clouds) cycle
        mfbnd_use(    k,1,0:ncls_use) = mfbnd_use(    k,1,0:ncls_use) &
                                      + mfbnd_use(    k,2,0:ncls_use)
        mfbnd_use(    k,2,0:ncls_use) = 0.0
        abnd_tavg_use(k,1,0:ncls_use) = abnd_tavg_use(k,1,0:ncls_use) &
                                      + abnd_tavg_use(k,2,0:ncls_use)
        abnd_tavg_use(k,2,0:ncls_use) = 0.0
        if (k > ktecen) cycle
        acen_tavg_use(k,1,0:ncls_use) = acen_tavg_use(k,1,0:ncls_use) &
                                      + acen_tavg_use(k,2,0:ncls_use)
        acen_tavg_use(k,2,0:ncls_use) = 0.0
        acen_tfin_use(k,1,0:ncls_use) = acen_tfin_use(k,1,0:ncls_use) &
                                      + acen_tfin_use(k,2,0:ncls_use)
        acen_tfin_use(k,2,0:ncls_use) = 0.0
    end do


!   at k = kts and k >= k_max_wnonzero
!   set mfbnd and wbnd_bar = 0
!       set areas = 0 for drafts (at kts set abnd=0 but allow acen>0)
    do k = kts, ktebnd
        if ((k > kts) .and. (k < k_max_wnonzero)) cycle

        mfbnd_use(k,:,:) = 0.0
        wbnd_bar_use(k) = 0.0

        do jcls = 1, ncls_use
        if (jcls == jcls_qu) then
            abnd_tavg_use(k,0:2,jcls) = abnd_tavg_use(k,0:2,0)
            if ((k == kts) .or. (k > ktecen)) cycle
            acen_tavg_use(k,0:2,jcls) = acen_tavg_use(k,0:2,0)
            acen_tfin_use(k,0:2,jcls) = acen_tfin_use(k,0:2,0)
        else
            abnd_tavg_use(k,0:2,jcls) = 0.0
            if ((k == kts) .or. (k > ktecen)) cycle
            acen_tavg_use(k,0:2,jcls) = 0.0
            acen_tfin_use(k,0:2,jcls) = 0.0
        end if
        end do
    end do


!   at k >= k_max_updndraft
!   set mfbnd = 0 and areas = 0 for drafts
!       set mfbnd = abnd*wbnd_bar*rhobnd_bar for quiescents
    do k = kts, ktebnd
        if ((k < k_max_updndraft) .or. (k >= k_max_wnonzero)) cycle

        do jcls = 1, ncls_use
        if (jcls == jcls_qu) then
            abnd_tavg_use(k,0:2,jcls) = abnd_tavg_use(k,0:2,0)
            mfbnd_use(k,1:2,jcls) = &
            abnd_tavg_use(k,1:2,jcls)*wbnd_bar_use(k)*rhobnd_bar(k)
            if (k > ktecen) cycle
            acen_tavg_use(k,0:2,jcls) = acen_tavg_use(k,0:2,0)
            acen_tfin_use(k,0:2,jcls) = acen_tfin_use(k,0:2,0)
        else
            abnd_tavg_use(k,0:2,jcls) = 0.0
            mfbnd_use(k,0:2,jcls) = 0.0
            if (k > ktecen) cycle
            acen_tavg_use(k,0:2,jcls) = 0.0
            acen_tfin_use(k,0:2,jcls) = 0.0
        end if
        end do
    end do


    if ((ktau==4) .and. (lun155 > 0)) then
    write(lun155,'(/a,3i5)') 'bbb', ktau, ipass_check_adjust_inputs
    write(lun155,'(3(i5,i3,1pe16.8))') ((jcls,icc,acen_tavg_use(26,icc,jcls),icc=0,2),jcls=0,3)
    end if
!
!   check updraft/dndraft
!
    do 3590 jcls = 1, ncls_use
    if (jcls == jcls_qu) goto 3590

    do 3490 icc = 1, 2
    jclsicc = jcls*10 + icc

!   check   kts        <= kdraft_bot <= ktecen
!   and     kdraft_bot <  kdraft_top <= ktecen
    if ( (kdraft_bot_use(icc,jcls) < kts)   .or.   &
         (kdraft_bot_use(icc,jcls) > ktecen) .or.   &
         (kdraft_top_use(icc,jcls) <= kdraft_bot_use(icc,jcls)) .or.   &
         (kdraft_top_use(icc,jcls) > ktecen) ) then
        msg = '*** parampollu_check_adjust_inputs - ' //   &
        'bad up/dndraft kdraft_bot/_top'
        call ecpp_message( lunout, msg )
        write( msg, '(a,4(1x,i5))' ) 'it, jt, jclsicc, mtype_updnenv =',   &
        it, jt, jclsicc, icc, mtype_updnenv_use(icc,jcls)
        call ecpp_message( lunout, msg )
        write( msg, '(a,2(1x,i5),2(1x,i10))' )   &
        'kts, ktebnd, kdraft_bot, kdraft_top =',   &
        kts, ktebnd, kdraft_bot_use(icc,jcls), kdraft_top_use(icc,jcls)
        call ecpp_message( lunout, msg )
        call ecpp_error_fatal( lunout, msg )
    end if

!   check/adjust mbfnd_use and abnd_tavg_use
!      if either is below the cut-off value, set both to zero
!      also set both to zero outside of [kdraft_bot_use, kdraft_top_use]
!   set the kdraft_bot/top_use
!
!   note that kdraft_bot/top define bottom/top for layer centers
!   for layer boundaries, the up/dndraft mfbnd and abnd are zero
!      at the bottom of kdraft_bot and at the top of kdraft_top
!
    ktmpa = -999888777 ; ktmpb = -999888778
    do k = kts, ktebnd
        arw_draft_cut = aw_draft_cut*rhobnd_bar(k)

        tmp_mfa = mfbnd_use(k,icc,jcls)
        tmp_mfb = tmp_mfa
        tmpa = abnd_tavg_use(k,icc,jcls)
        tmpb = tmpa

        if ( (k <= kdraft_bot_use(icc,jcls)) .or.   &
             (k >  kdraft_top_use(icc,jcls)) .or.   &
             (k == kts) ) then
        tmp_mfb = 0.0_r8
        else
        if (mtype_updnenv_use(icc,jcls) == mtype_updraft_ecpp) then
            if ( tmp_mfa < arw_draft_cut) tmp_mfb = 0.0_r8
        else
            if (-tmp_mfa < arw_draft_cut) tmp_mfb = 0.0_r8
        end if
        if (abnd_tavg_use(k,icc,jcls) < afrac_cut) tmp_mfb = 0.0_r8
        end if

        if (tmp_mfb /= 0.0_r8) then
        tmpb = max( tmpb, afrac_cut )
        else
        tmpb = 0.0_r8
        end if

        mfbnd_use(k,icc,jcls) = tmp_mfb
        abnd_tavg_use(k,icc,jcls) = tmpb
        if (tmp_mfb /= 0.0_r8) then
        if (ktmpa <= 0) ktmpa = k-1
        ktmpb = k
        end if

!   set change counts
!   increment/decrement abnd of quiescent class if up/dndraft abnd has changed
        if (tmp_mfb /= tmp_mfa) then
        nchanges(1) = nchanges(1) + 1
        end if
        if (tmpb /= tmpa) then
        nchanges(2) = nchanges(2) + 1
        abnd_tavg_use(k,icc,jcls_qu) = abnd_tavg_use(k,icc,jcls_qu) &
                                   + (tmpa-tmpb)
        end if

    end do

    kdraft_bot_use(icc,jcls) = ktmpa
    kdraft_top_use(icc,jcls) = ktmpb

!   check/adjust acen_tavg_use
!       set acen_tavg to zero outside of kdraft_bot:kdraft_top
!       set acen_tavg to zero if abnd_tavg=0 at both layer boundaries (14-apr-2009)
    do k = kts, ktecen
        tmpa = acen_tavg_use(k,icc,jcls)
        tmpb = tmpa

        if ( (k < kdraft_bot_use(icc,jcls)) .or.   &
             (k > kdraft_top_use(icc,jcls)) ) then
        tmpb = 0.0_r8
        else
        tmpe = 0.5*( abnd_tavg_use(k,  icc,jcls) + &
                     abnd_tavg_use(k+1,icc,jcls) )
        if (tmpe > 0.0) then
            tmpb = max( afrac_cut, tmpe )
        else
            tmpb = 0.0
        end if
        end if

        if (tmpb /= tmpa) then
        nchanges(3) = nchanges(3) + 1
        acen_tavg_use(k,icc,jcls_qu) = &
            acen_tavg_use(k,icc,jcls_qu) + (tmpa-tmpb)
        end if

        acen_tavg_use(k,icc,jcls) = tmpb
    end do

!   check/adjust acen_tfin_use
!       set acen_tfin to zero if it is < afrac_cut or if k >= k_max_updndraft
!       set acen_tfin to zero if acen_tavg=0 (14-apr-2009)
!   for case of parampollu_opt == 2220, but iflag_ecpp_test_fixed_fcloud /= 2,3,4,5
!       do not allow acen_tfin=0 if acen_tavg>0
!   (14-apr-2009 -- do similar for all parampollu_opt)
    ispecial_check_acen_tfin = 0
    if (parampollu_opt == 2220) then
        ispecial_check_acen_tfin = 1
        if ((iflag_ecpp_test_fixed_fcloud >= 2) .and. &
            (iflag_ecpp_test_fixed_fcloud <= 5)) ispecial_check_acen_tfin = 0
    end if
    if (ispecial_check_acen_tfin <= 0) then
        ispecial_check_acen_tfin = 2
        if ((iflag_ecpp_test_fixed_fcloud >= 2) .and. &
            (iflag_ecpp_test_fixed_fcloud <= 5)) ispecial_check_acen_tfin = 0
    end if

    do k = kts, ktecen
        tmpa = acen_tfin_use(k,icc,jcls)
        tmpb = tmpa

        if ((tmpa < afrac_cut) .or. &
            (k >= k_max_updndraft)) then
        tmpb = 0.0_r8
        end if
        if (acen_tavg_use(k,icc,jcls) <= 0.0) then
        tmpb = 0.0_r8
        end if

        if (ispecial_check_acen_tfin > 0) then
        if (tmpb < afrac_cut) then
            if (acen_tavg_use(k,icc,jcls) >= afrac_cut) then
            if (ispecial_check_acen_tfin == 2) then
                tmpb = max( 0.5_r8*acen_tavg_use(k,icc,jcls), afrac_cut )
            else
                tmpb = acen_tavg_use(k,icc,jcls)
            end if
            end if
        end if
        end if

        if (tmpb /= tmpa) then
        nchanges(4) = nchanges(4) + 1
        acen_tfin_use(k,icc,jcls_qu) = &
            acen_tfin_use(k,icc,jcls_qu) + (tmpa-tmpb)
        end if

        acen_tfin_use(k,icc,jcls) = tmpb
    end do

!   for empty sub-class (mfbnd/abnd/acen=0 at all levels), 
!   set kdraft_bot/top_use to ktecen
    if ((kdraft_bot_use(icc,jcls) < -999888000) .and.   &
        (kdraft_top_use(icc,jcls) < -999888000)) then
        kdraft_bot_use(icc,jcls) = ktecen
        kdraft_top_use(icc,jcls) = ktecen
    end if

3490    continue

!   sum clear and cloudy mfbnd_use
    do k = kts, ktebnd
        mfbnd_use(k,0,jcls) = sum( mfbnd_use(k,1:2,jcls) ) 
    end do

3590    continue


!
!   check/adjust quiescent transport-class
!

    if ((ktau==4) .and. (lun155 > 0)) then
    write(lun155,'(/a,3i5)') 'ccc', ktau, ipass_check_adjust_inputs
    write(lun155,'(3(i5,i3,1pe16.8))') ((jcls,icc,acen_tavg_use(26,icc,jcls),icc=0,2),jcls=0,3)
    end if
!   first set to zero any areas that are < afrac_cut
    do k = kts, ktebnd
        do i = 1, 3
            jcls = jcls_qu
            if ((i >= 2) .and. (k > ktecen)) cycle

            if (i == 1) then
                tmpvecaa(0:2) = abnd_tavg_use(k,0:2,jcls)
            else if (i == 2) then
                tmpvecaa(0:2) = acen_tavg_use(k,0:2,jcls)
            else
                tmpvecaa(0:2) = acen_tfin_use(k,0:2,jcls)
            end if

            tmpvecbb(0:2) = tmpvecaa(0:2)
            tmpvecbb(0) = tmpvecbb(1) + tmpvecbb(2)
            do icc = 1, 2
                if (tmpvecbb(icc) < afrac_cut) then
                ! if (tmpvecbb(icc) < 0.0_r8) then  ! whannah - this didn't work
                    tmpvecbb(3-icc) = tmpvecbb(0)
                    tmpvecbb(icc) = 0.0_r8
                end if
            end do

            ! for case of parampollu_opt == 2220, but iflag_ecpp_test_fixed_fcloud /= 2,3,4,5
            ! do not allow acen_tfin=0 if acen_tavg>0
            if ((i == 3) .and. (ispecial_check_acen_tfin > 0)) then
                do icc = 1, 2
                if (tmpvecbb(icc) < afrac_cut) then
                    if (acen_tavg_use(k,icc,jcls) >= afrac_cut) then
                        tmpvecbb(icc) = acen_tavg_use(k,icc,jcls)
                        tmpvecbb(3-icc) = tmpvecbb(0) - tmpvecbb(icc)
                    end if
                end if
                end do
            end if

    ! whannah 
    if (tmpvecbb(1) < 0.0_r8) then  
        tmpvecbb(2) = tmpvecbb(0)
        tmpvecbb(1) = 0.0_r8
    end if
    if (tmpvecbb(2) < 0.0_r8) then  
        tmpvecbb(1) = tmpvecbb(0)
        tmpvecbb(2) = 0.0_r8
    end if
    ! whannah 

            if ((tmpvecbb(1) < 0.0_r8) .or. &
                (tmpvecbb(2) < 0.0_r8) .or. &
                (tmpvecbb(0) < a_quiescn_minbb)) then
                ! at this point, the total (adjusted) quiescent area is too small
                write(msg,'(a,1p,3e12.4)') &
                '    tmpvecaa(0:2) = v1 quiescent areas =', tmpvecaa(0:2)
                call ecpp_message( lunout, msg )
                write(msg,'(a,1p,3e12.4)') &
                '    tmpvecbb(0:2) = v2 quiescent areas =', tmpvecbb(0:2)
                call ecpp_message( lunout, msg )

                write(msg,'(2a,2i5)') &
                '*** parampollu_check_adjust_inputs - a_quiescent(v2) too small ', &
                area_name10(i), k, i
                call ecpp_message( lunout, msg )
                call ecpp_error_fatal( lunout, msg )
            end if

            if (i == 1) then
                abnd_tavg_use(k,0:2,jcls) = tmpvecbb(0:2)
            else if (i == 2) then
                acen_tavg_use(k,0:2,jcls) = tmpvecbb(0:2)
            else
                acen_tfin_use(k,0:2,jcls) = tmpvecbb(0:2)
            end if
        end do   ! i = 1, 3
    end do   ! k = kts, ktebnd


!   recalc summed area fractions
    do k = kts, ktebnd
        do jcls = 1, ncls_use
        abnd_tavg_use(k,0,jcls) = sum( abnd_tavg_use(k,1:2,jcls) )
        if (k > ktecen) cycle
        acen_tavg_use(k,0,jcls) = sum( acen_tavg_use(k,1:2,jcls) )
        acen_tfin_use(k,0,jcls) = sum( acen_tfin_use(k,1:2,jcls) )
        end do
        do icc = 0, 2
        abnd_tavg_use(k,icc,0) = sum( abnd_tavg_use(k,icc,1:ncls_use) )
        if (k > ktecen) cycle
        acen_tavg_use(k,icc,0) = sum( acen_tavg_use(k,icc,1:ncls_use) )
        acen_tfin_use(k,icc,0) = sum( acen_tfin_use(k,icc,1:ncls_use) )
        end do
    end do   ! k = kts, ktebnd


!   calc kdraft_bot_use & kdraft_top_use
    jcls = jcls_qu
    do icc = 1, 2
        ktmpa = -999888777 ; ktmpb = -999888778
        do k = kts, ktecen
        if (acen_tavg_use(k,icc,jcls) > 0.0_r8) then
            if (ktmpa <= 0) ktmpa = k
            ktmpb = k
        end if
        end do
        kdraft_bot_use(icc,jcls) = ktmpa
        kdraft_top_use(icc,jcls) = ktmpb
    end do

!   normally allow cloudy quiescent to be empty
!   if iflag_ecpp_test_fixed_fcloud=3 (special testing), allow clear quiescent to be empty
    icc = 2
    if (iflag_ecpp_test_fixed_fcloud == 3) icc = 1
    if ((kdraft_bot_use(icc,jcls) < -999888000) .and.   &
        (kdraft_top_use(icc,jcls) < -999888000)) then
         kdraft_bot_use(icc,jcls) = ktecen
         kdraft_top_use(icc,jcls) = ktecen
    end if

!   check for validity of kdraft_bot_use & kdraft_top_use
    ka = min( kdraft_bot_use(1,jcls), kdraft_bot_use(2,jcls) )
    kb = max( kdraft_top_use(1,jcls), kdraft_top_use(2,jcls) )
    do icc = 1, 2
        if ( (kdraft_bot_use(icc,jcls) <  kts) .or.   &
         (kdraft_bot_use(icc,jcls) > ktecen) .or.   &
         (kdraft_bot_use(icc,jcls) >  kdraft_top_use(icc,jcls)) .or.   &
         (kdraft_top_use(icc,jcls) > ktecen) .or.   &
         (ka /= kts)                           .or.   &
         (kb /= ktecen)                        ) then
        jclsicc = jcls*10 + icc
        msg = '*** parampollu_check_adjust_inputs - ' //   &
            'bad quiescent transport-class kdraft_bot/top_use'
        call ecpp_message( lunout, msg )
        write( msg, '(a,4(1x,i5))' ) 'it, jt, jclsicc, mtype_updnenv =',   &
            it, jt, jclsicc, mtype_updnenv_use(icc,jcls)
        call ecpp_message( lunout, msg )
        write( msg, '(a,2(1x,i5),2(1x,i10))' )   &
            'kts, ktebnd, kdraft_bot, kdraft_top =',   &
            kts, ktebnd, kdraft_bot_use(icc,jcls), kdraft_top_use(icc,jcls)
        call ecpp_message( lunout, msg )
        call ecpp_error_fatal( lunout, msg )
        end if
    end do


!-----------------------------------------------------
!   here ipass_check_adjust_inputs == 1
!   skip over the special stuff for ipass_check_adjust_inputs == 2
!-----------------------------------------------------
    if ((ktau==4) .and. (lun155 > 0)) then
    write(lun155,'(/a,3i5)') 'ddd', ktau, ipass_check_adjust_inputs
    write(lun155,'(3(i5,i3,1pe16.8))') ((jcls,icc,acen_tavg_use(26,icc,jcls),icc=0,2),jcls=0,3)
    end if
    goto 30000


!-----------------------------------------------------
!   special stuff for ipass_check_adjust_inputs == 2
!-----------------------------------------------------
20000   continue
    ipass_2_changes = 0


!   for testing only -- reduce up/dndraft areas
!   *** NOTE / TODO - in the "new" code, this may not work correctly
    ido_updndr_area_adjust = 0
    if (ido_updndr_area_adjust > 0) then
    ipass_2_changes = ipass_2_changes + 1

    updndr_area_adjust = 1.0_r8
    tmpb = 1.0_r8 - updndr_area_adjust
    do k = kts, ktebnd
        do icc = 0, 2
        do jcls = 1, ncls_use
            if (jcls == jcls_qu) cycle

            abnd_tavg_use(k,icc,jcls_qu) = abnd_tavg_use(k,icc,jcls_qu)   &
                                         + abnd_tavg_use(k,icc,jcls   )*tmpb
            abnd_tavg_use(k,icc,jcls   ) = abnd_tavg_use(k,icc,jcls   )*updndr_area_adjust

            if (k > ktecen) cycle

            acen_tavg_use(k,icc,jcls_qu) = acen_tavg_use(k,icc,jcls_qu)   &
                                         + acen_tavg_use(k,icc,jcls   )*tmpb
            acen_tavg_use(k,icc,jcls   ) = acen_tavg_use(k,icc,jcls   )*updndr_area_adjust

            acen_tfin_use(k,icc,jcls_qu) = acen_tfin_use(k,icc,jcls_qu)   &
                                         + acen_tfin_use(k,icc,jcls   )*tmpb
            acen_tfin_use(k,icc,jcls   ) = acen_tfin_use(k,icc,jcls   )*updndr_area_adjust

        end do
        end do
    end do
    end if   ! (ido_updndr_area_adjust > 0)


!   for testing only -- zero out downdraft
!   *** NOTE / TODO - in the "new" code, this may not work correctly
    ido_downdr_area_zeroout = 0
    if (ido_downdr_area_zeroout > 0) then
    ipass_2_changes = ipass_2_changes + 1

    do k = kts, ktebnd
        do icc = 0, 2
        do jcls = 1, ncls_use
            if (jcls == jcls_qu) cycle
            if (mtype_updnenv_use(icc,jcls) /= mtype_dndraft_ecpp) cycle

            abnd_tavg_use(k,icc,jcls_qu) = abnd_tavg_use(k,icc,jcls_qu)   &
                                         + abnd_tavg_use(k,icc,jcls   )
            abnd_tavg_use(k,icc,jcls   ) = 0.0

            mfbnd_use(    k,icc,jcls_qu) = mfbnd_use(    k,icc,jcls_qu)   &
                                         + mfbnd_use(    k,icc,jcls   )
            mfbnd_use(    k,icc,jcls   ) = 0.0

            if (k > ktecen) cycle

            acen_tavg_use(k,icc,jcls_qu) = acen_tavg_use(k,icc,jcls_qu)   &
                                         + acen_tavg_use(k,icc,jcls   )
            acen_tavg_use(k,icc,jcls   ) = 0.0

            acen_tfin_use(k,icc,jcls_qu) = acen_tfin_use(k,icc,jcls_qu)   &
                                         + acen_tfin_use(k,icc,jcls   )
            acen_tfin_use(k,icc,jcls   ) = 0.0
        end do
        end do
    end do
    end if   ! (ido_downdr_area_zeroout > 0)


!   if (ipass_2_changes == 0) return


!-----------------------------------------------------
!   common stuff for ipass_check_adjust_inputs == 1,2
!-----------------------------------------------------
30000   continue
!
!   check/adjust quiescent abnd_tavg_use (and mfbnd_use)
!
!   before 15-jul-2008 code
!      code above may have set afrac_bnd=0 in some transport-class
!      now adjust afrac_bnd in quiescent transport-class so that 
!          all-transport-class-sum = 1.0
!
!   on/after 15-jul-2008 code
!      the post-processor does not correctly identify the clear versus
!         cloudy parts of the quiescent abnd_tavg
!         (it calcs an average qcloud for 2 layers adjacent to the boundary,
!         and if qcloud in either layer exceeds cutoff, then the average
!         will too (almost always), so this is biased)
!      so instead, set these based on the clear/cloud quiescent acen_tavg_use
!      also, apportion the quiescent mfbnd_use similarly
!
    mfbnd_quiescn_up(:,:,:) = 0.0
    mfbnd_quiescn_dn(:,:,:) = 0.0

    jcls = jcls_qu
    do k = kts, ktecen
!   first calc tmpvecdd(k) = fraction of layer-k quiescent-class that is clear
        ardz_cut = afrac_cut*rhodz_cen(k)*0.3_r8
        if ((acen_tavg_use(k,1,jcls) >= ardz_cut) .and.   &
        (acen_tavg_use(k,2,jcls) >= ardz_cut)) then
        ! clear and cloudy both > 0
        tmpvecdd(k) = acen_tavg_use(k,1,jcls)/acen_tavg_use(k,0,jcls)
        tmpvecdd(k) = max( 0.0_r8, min( 1.0_r8, tmpvecdd(k) ) )
        else if (acen_tavg_use(k,2,jcls) >= ardz_cut) then
        ! only cloudy > 0
        tmpvecdd(k) = 0.0
        else
        ! only clear > 0
        tmpvecdd(k) = 1.0_r8
        end if
    end do


    do k = kts+1, ktecen
!   calc (total quiescent "w-prime" mass flux) = - (sum of up/dndraft mass fluxes)
        tmp_mfa = 0.0
        do jcls = 1, ncls_use
        if (jcls == jcls_qu) cycle
        mfbnd_use(k,0,jcls) = sum( mfbnd_use(k,1:2,jcls) ) 
        tmp_mfa = tmp_mfa + mfbnd_use(k,0,jcls)
        end do
        jcls = jcls_qu
        mfbnd_use(k,0,jcls) = -tmp_mfa

!   partition total quiescent mass flux to clear/cloudy using the
!   quiescent clear/cloud amounts in the "upwind" layer
        if (mfbnd_use(k,0,jcls) < 0.0) then
        tmpvecaa(1) = tmpvecdd(k)     !   upwind is layer above
        tmpvecbb(1) = tmpvecdd(k-1)   ! downwind is layer below
        else
        tmpvecaa(1) = tmpvecdd(k-1)   !   upwind is layer below
        tmpvecbb(1) = tmpvecdd(k)     ! downwind is layer above
        end if
        tmpvecaa(2) = 1.0_r8 - tmpvecaa(1)
        tmpvecbb(2) = 1.0_r8 - tmpvecbb(1)

        mfbnd_use(k,1:2,jcls) = mfbnd_use(k,0,jcls)*tmpvecaa(1:2)
!   same for abnd
        abnd_tavg_use(k,1:2,jcls) = abnd_tavg_use(k,0,jcls)*tmpvecaa(1:2)

!   do other sums
        do icc = 0, 2
        mfbnd_use(    k,icc,0) = sum( mfbnd_use(    k,icc,1:ncls_use) )
        abnd_tavg_use(k,icc,0) = sum( abnd_tavg_use(k,icc,1:ncls_use) )
        end do


!   now calculate more detailed up and down fluxes
!      mfbnd_quiescn_up(k,jccfrom,jcctooo) is mbbnd from (k,jccfrom) to (k+1,jcctooo)
!         with jccfrom=0/1/2=both/clear/cloudy; and jcctooo=0/1/2=similar
!
!      the clear-->both   and cloudy-->both   are already determined
!      the clear-->clear  and cloudy-->cloudy are calculated maximum overlap
!         of cloudy and clear regions
!      the clear-->cloudy and cloudy-->clear  are simply what is left
!
!      tmpvecaa holds clear/cloudy fractions of the   upwind layer
!      tmpvecbb holds clear/cloudy fractions of the downwind layer
        jcls = jcls_qu
        tmp0202aa(0:2,0) = mfbnd_use(k,0:2,jcls) 
        tmp0202aa(0:2,1:2) = 0.0
        do ja = 1, 2
        jb = 3-ja
        tmpa = 0.0
        if (tmpvecaa(ja) > 0.0) &
            tmpa = min(tmpvecbb(ja),tmpvecaa(ja))/tmpvecaa(ja)
        tmpa = max( 0.0_r8, min( 1.0_r8, tmpa ) )
        tmp0202aa(ja,ja) = tmp0202aa(ja,0)*tmpa
        tmp0202aa(ja,jb) = tmp0202aa(ja,0)*(1.0_r8-tmpa)
        end do
        do jb = 1, 2
        tmp0202aa(0,jb) = sum( tmp0202aa(1:2,jb) )
        end do
        if (mfbnd_use(k,0,jcls) < 0.0) then
        mfbnd_quiescn_dn(k,0:2,0:2) = tmp0202aa(0:2,0:2) 
        else if (mfbnd_use(k,0,jcls) > 0.0) then
        mfbnd_quiescn_up(k,0:2,0:2) = tmp0202aa(0:2,0:2) 
        end if

!   if ((ipass_check_adjust_inputs == 2) .and. (lun141 > 0)) then
!   if (k == kts+1) write( 141, '(/a,2i5)' )   &
!       'mfbnd_quiescn at ktau, ipass =', ktau, ipass_check_adjust_inputs
!   write( 141, '(i3,1p,2e11.3,2(2x,4e11.3))' ) k, mfbnd_use(k,1:2,jcls),   &
!       mfbnd_quiescn_up(k,1:2,1:2), mfbnd_quiescn_dn(k,1:2,1:2) 
!   end if

    end do   ! k = kts+1, ktecen


!   for "empty" drafts, reset the kbot & ktop, and also the mtype_updnenv_use
!
!   *** currently the reset of mtype_updnenv_use is deactivated
    kdraft_bot_tmp(:,:) = kdraft_bot_use(:,:)
    kdraft_top_tmp(:,:) = kdraft_top_use(:,:)
    mtype_updnenv_tmp(:,:) = mtype_updnenv_use(:,:)
    if (lun63 > 0) write(lun63,'(a/2a)') &
        'parampollu_check_adjust_inputs transport-class summary', &
        ' jcls  mcc, mf/af nonzero, mtype_tmp/use, ', &
        'kbase/top_inp, kbase/top_tmp, kbase/top_use'
    do jcls = 1, ncls_use
    do icc = 1, 2
        itmpa = 0
        itmpb = 0
        do k = kts, ktebnd
        if (mfbnd_use(k,icc,jcls)    /= 0.0) itmpa = itmpa + 1
        if (abnd_tavg_use(k,icc,jcls) /= 0.0) itmpb = itmpb + 1
        end do
        if (itmpa+itmpb <= 0) then
            kdraft_bot_use(icc,jcls) = ktecen
            kdraft_top_use(icc,jcls) = ktecen
!       if (mtype_updnenv_use(icc,jcls) == mtype_updraft_ecpp) then
!           mtype_updnenv_use(icc,jcls) = mtype_upempty_ecpp
!       else if (mtype_updnenv_use(icc,jcls) == mtype_dndraft_ecpp) then
!           mtype_updnenv_use(icc,jcls) = mtype_dnempty_ecpp
!       else
!           mtype_updnenv_use(icc,jcls) = mtype_quempty_ecpp
!       end if
        end if
        if (lun63 > 0) write(lun63,'(2i5,5(5x,2i5))')   &
        jcls, icc, itmpa, itmpb,   &
        mtype_updnenv_tmp(icc,jcls), mtype_updnenv_use(icc,jcls),   &
        kdraft_bot_ecpp(icc,jcls), kdraft_top_ecpp(icc,jcls),   &
        kdraft_bot_tmp(icc,jcls), kdraft_top_tmp(icc,jcls),   &
        kdraft_bot_use(icc,jcls), kdraft_top_use(icc,jcls)
    end do
    end do


!   now adjust area with precipitation
    acen_prec_use(:,:,:) = 0.0
    do jcls = 1, ncls_use
    do icc = 1, 2
    do k = kts, ktecen
        if (acen_tavg_use(k,icc,jcls) < afrac_cut) cycle
        if (acen_prec_ecpp(k,icc,jcls) < afrac_cut) cycle

        tmpa = acen_prec_ecpp(k,icc,jcls)   ! portion of sub-area with precip
        tmpb = acen_tavg_use(k,icc,jcls) - tmpa   ! portion of sub-area without precip
        if (tmpb < afrac_cut) tmpa = acen_tavg_use(k,icc,jcls)
        acen_prec_use(k,icc,jcls) = tmpa
    end do
    end do
    end do


!   final recalc summed area fractions
    do k = kts, ktebnd
        do jcls = 1, ncls_use
        abnd_tavg_use(k,0,jcls) = sum( abnd_tavg_use(k,1:2,jcls) )
        if (k > ktecen) cycle
        acen_tavg_use(k,0,jcls) = sum( acen_tavg_use(k,1:2,jcls) )
        acen_tfin_use(k,0,jcls) = sum( acen_tfin_use(k,1:2,jcls) )
        acen_prec_use(k,0,jcls) = sum( acen_prec_use(k,1:2,jcls) )
        end do
        do icc = 0, 2
        abnd_tavg_use(k,icc,0) = sum( abnd_tavg_use(k,icc,1:ncls_use) )
        if (k > ktecen) cycle
        acen_tavg_use(k,icc,0) = sum( acen_tavg_use(k,icc,1:ncls_use) )
        acen_tfin_use(k,icc,0) = sum( acen_tfin_use(k,icc,1:ncls_use) )
        acen_prec_use(k,icc,0) = sum( acen_prec_use(k,icc,1:ncls_use) )
        end do
    end do


    if (lun63 > 0) then
    write(lun63,'(a,i2)') 'parampollu_check_adjust_inputs -- ipass =',   &
        ipass_check_adjust_inputs
    do k = 1, 10
        write(lun63,'(a,i2,a,i10)') '    nchanges(', k, ') =', nchanges(k)
    end do
    end if ! (lun63 > 0)


    if ((ktau==4) .and. (lun155 > 0)) then
    write(lun155,'(/a,3i5)') 'eee', ktau, ipass_check_adjust_inputs
    write(lun155,'(3(i5,i3,1pe16.8))') ((jcls,icc,acen_tavg_use(26,icc,jcls),icc=0,2),jcls=0,3)
    end if

    return
    end subroutine parampollu_check_adjust_inputs



!-----------------------------------------------------------------------
    subroutine parampollu_1clm_dumpaa(                        &
        ktau, dtstep, ktau_pp, dtstep_pp,             &
        tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
        rhobnd_bar, zbnd, wbnd_bar,                       &
        chem_bar,                                         &
        ncls_ecpp,                                        &
        kdraft_bot_ecpp, kdraft_top_ecpp,                 &
        mtype_updnenv_ecpp,                               &
        mfbnd, abnd_tavg,                                 &
        acen_tavg, acen_tbeg, acen_tfin,                  &
        it,      jt,      kts,ktebnd,ktecen,              &
        lun                                               )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_1clm_dumpaa does a diagnostic print of
!    numerous ecpp arrays
!
!-----------------------------------------------------------------------

    use module_data_ecpp1

    use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message

!   arguments
    integer, intent(in) ::                  &
        ktau, ktau_pp,              &
        it, jt, kts, ktebnd, ktecen,    &
        lun
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations

!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!   chem_driver and routines under it do calculations
!   over these spatial indices.

    real(r8), intent(in) :: dtstep, dtstep_pp
!   dtstep - main model time step (s)
!   dtstep_pp - time step (s) for "parameterized pollutants" calculations

    real(r8), intent(in), dimension( kts:ktecen ) ::   &
        tcen_bar, pcen_bar, rhocen_bar, dzcen
    real(r8), intent(in), dimension( kts:ktebnd ) ::   &
        rhobnd_bar, zbnd, wbnd_bar

    real(r8), intent(in), dimension( kts:ktecen, 1:num_chem_ecpp ) :: &
        chem_bar

    integer, intent(in) :: ncls_ecpp
    
    integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
        kdraft_bot_ecpp, kdraft_top_ecpp,   &
        mtype_updnenv_ecpp
    
    real(r8), intent(in), dimension( kts:ktebnd, 0:2, 0:maxcls_ecpp ) ::   &
        mfbnd, abnd_tavg
    real(r8), intent(in), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) ::   &
        acen_tavg, acen_tbeg, acen_tfin
    
    character(len=8), dimension( kts:ktebnd ) :: dumchar8
    

!   local variables
    integer :: iclrcld
    integer :: itmp_mtype_clrcldy(1:2)
    integer :: jcls, jclsaa, jclsbb
    integer :: k, l

    real(r8) :: duma
    real(r8), dimension( kts:ktebnd ) :: dumarr1, dumarr2, dumarr3, dumarr4, dumarr5


!
!   output with same format as ppboxmakeinp01
!
9400    format( a )
9410    format( 5i15 )
9415    format( a, i10 )
9416    format( a, 5i10 )
!9420   format( 5(1pe15.7) )
9420    format( 5(1pe12.4) )

    if (lun <= 0) return


    itmp_mtype_clrcldy(1) = mtype_nocloud_ecpp
    itmp_mtype_clrcldy(2) = mtype_iscloud_ecpp

!   write(lun,9400) 'output from ppboxmakeinp01'
    write(lun,9400)
    write(lun,9400)
    write(lun,9416) 'output from ppboxmakeinp01 - ktau, ktau_pp',   &
        ktau, ktau_pp

    write(lun,9400) 'kts, kte, ncls_ecpp_clm'
    write(lun,9410) kts, ktebnd, ncls_ecpp

    write(lun,9400) 'num_moist, num_chem'
    write(lun,9410) num_moist_ecpp, num_chem_ecpp

    write(lun,9400) 'rho,z,w bnd'
    do k = kts, ktebnd
        write(lun,9420) rhobnd_bar(k),   &
        zbnd(k), wbnd_bar(k)
    end do

    write(lun,9400) 'p,t,rho cen'
    do k = kts, ktecen
        write(lun,9420) pcen_bar(k), tcen_bar(k), rhocen_bar(k)
    end do

    do l = 1, num_moist_ecpp
        write(lun,9415) 'moist', l
        write(lun,9420) (0.0, k=kts,ktecen)
    end do

    do l = 1, num_chem_ecpp
        write(lun,9415) 'chem ', l
        write(lun,9420) (chem_bar(k,l), k=kts,ktecen)
    end do

    do jcls = 1, ncls_ecpp
    do iclrcld = 1, 2
        write(lun,9416) 'jcls, iclrcld // mtype a,b,c; kdraft a,b', jcls, iclrcld
        write(lun,9410)   &
        mtype_updnenv_ecpp(iclrcld,jcls),   &
        itmp_mtype_clrcldy(iclrcld), mtype_noprecip_ecpp,   &
        kdraft_bot_ecpp(iclrcld,jcls), kdraft_top_ecpp(iclrcld,jcls)

        write(lun,9416) 'afrac', jcls, iclrcld
        write(lun,9420) (abnd_tavg(k,iclrcld,jcls), k=kts,ktebnd)

        write(lun,9416) 'mf', jcls, iclrcld
        write(lun,9420) (mfbnd(k,iclrcld,jcls), k=kts,ktebnd)
    end do
    end do
     

    write(lun,'(/a)') 'baraa'
    write(lun,'(a)') ' k     z(km)   p(mb)     rho    t(C) w(cm/s)'
    do k = ktebnd, kts, -1
        if (k < ktebnd) then
          duma = zbnd(k) + 0.5*dzcen(k)
          write(lun,'(i2,2x,f8.3,f8.1,f8.4,f8.1,  8x)')   &
        k, duma*1.0e-3, pcen_bar(k)*1.0e-2, rhocen_bar(k), tcen_bar(k)-273.16
        end if
        duma = k-0.5
        write(lun,'( f4.1,f8.3,  8x,f8.4,  8x,f8.2)')   &
        duma, zbnd(k)*1.0e-3, rhobnd_bar(k), wbnd_bar(k)*1.0e2
    end do
    write(lun,'(a)') ' k     z(km)   p(mb)     rho    t(C) w(cm/s)'

    write(lun,'(/a)') 'draftaa'
    do jcls = 1, ncls_ecpp
    do iclrcld = 1, 2
        write(lun,'(/a,7i5)') 'draftbb - ktau_pp, jcls, iclrcld, updn, clrcldy, top, bot =',   &
        ktau_pp, jcls, iclrcld,   &
        mtype_updnenv_ecpp(iclrcld,jcls), itmp_mtype_clrcldy(iclrcld),   &
        kdraft_bot_ecpp(iclrcld,jcls), kdraft_top_ecpp(iclrcld,jcls) 

        write(lun,'(a)') 'afrac'
        do k = kts, ktebnd
        duma = abnd_tavg(k,iclrcld,jcls)
        if (duma == 0.0) then
            dumchar8(k) = '  0.    '
        else if (abs(duma) >= 5.0e-5) then
            write(dumchar8(k),'(f8.4)') duma
        else
           write(dumchar8(k),'(1p,e8.0)') duma
        end if
        end do
        write(lun,'(15a)') (dumchar8(k), k=kts,ktebnd)

        do k = kts, ktebnd
        duma = max( 1.0e-10_r8, abnd_tavg(k,iclrcld,jcls) )
        dumarr1(k) = mfbnd(k,iclrcld,jcls)/(rhobnd_bar(k)*duma)
        end do
        write(lun,'(a)') 'w'
        write(lun,'(15f8.4)') (dumarr1(k), k=kts,ktebnd)

        write(lun,'(a)') 'mfbnd'
        write(lun,'(1p,10e12.5)') (mfbnd(k,iclrcld,jcls), k=kts,ktebnd)

        write(lun,'(a)') 'abnd_tavg'
        write(lun,'(1p,10e12.5)') (abnd_tavg(k,iclrcld,jcls), k=kts,ktebnd)
!       write(lun,'(1p,15e8.1 )') (abnd_tavg(k,iclrcld,jcls), k=kts,ktebnd)

        write(lun,'(a)') 'acen_tavg'
        write(lun,'(1p,10e12.5)') (acen_tavg(k,iclrcld,jcls), k=kts,ktecen)

        write(lun,'(a)') 'acen_tbeg'
        write(lun,'(1p,10e12.5)') (acen_tbeg(k,iclrcld,jcls), k=kts,ktecen)

        write(lun,'(a)') 'acen_tfin'
        write(lun,'(1p,10e12.5)') (acen_tfin(k,iclrcld,jcls), k=kts,ktecen)

    end do
    end do

    do k = kts, ktebnd
        dumarr1(k) = 0.0
        dumarr2(k) = 0.0
        dumarr3(k) = 0.0
        dumarr4(k) = 0.0
        dumarr5(k) = 0.0
        do jcls = 1, ncls_ecpp
        do iclrcld = 1, 2
        dumarr1(k) = dumarr1(k) + mfbnd(k,iclrcld,jcls)
        dumarr2(k) = dumarr2(k) + abnd_tavg(k,iclrcld,jcls)
        if (k > ktecen) cycle
        dumarr3(k) = dumarr3(k) + acen_tavg(k,iclrcld,jcls)
        dumarr4(k) = dumarr4(k) + acen_tbeg(k,iclrcld,jcls)
        dumarr5(k) = dumarr5(k) + acen_tfin(k,iclrcld,jcls)
        end do
        end do
        duma = max( 1.0e-10_r8, dumarr2(k) )
        dumarr1(k) = dumarr1(k)/(rhobnd_bar(k)*duma)
    end do
    write(lun,'(/a,4i5)') 'draftbb - ktau_pp, all subs            =',   &
        ktau_pp
    write(lun,'(a)') 'wbar'
    write(lun,'(12f10.5)') (wbnd_bar(k), k=kts,ktebnd)
    write(lun,'(a)') '(mfbnd summed over all subs)/rhobnd'
    write(lun,'(12f10.5)') (dumarr1(k), k=kts,ktebnd)
    write(lun,'(a)') '(abnd_tavg-1) summed over all subs'
    write(lun,'(1p,12e10.2)') ((dumarr2(k)-1.0), k=kts,ktebnd)
    write(lun,'(a)') '(acen_tavg-1) summed over all subs'
    write(lun,'(1p,12e10.2)') ((dumarr3(k)-1.0), k=kts,ktecen)
    write(lun,'(a)') '(acen_tbeg-1) summed over all subs'
    write(lun,'(1p,12e10.2)') ((dumarr4(k)-1.0), k=kts,ktecen)
    write(lun,'(a)') '(acen_tfin-1) summed over all subs'
    write(lun,'(1p,12e10.2)') ((dumarr5(k)-1.0), k=kts,ktecen)


    do jclsaa = 1, ncls_ecpp, 3
    jclsbb = min( jclsaa+2, ncls_ecpp )
    write(lun,'(/a,3i5)') 'draftcc - ktau_pp, jclsaa, jclsbb',   &
        ktau_pp, jclsaa, jclsbb
    write(lun,'(a)')   &
        'k,  acen_tavg(k,1:2,jclsaa:jclsbb),  mfbnd(k+1,1:2,jclsaa:jclsbb)'
    do k = ktecen, kts, -1
       write(lun,'(i3,2x,3(1x,2f8.5),2x,1p,3(1x,2e10.2))') k,   &
        acen_tavg(k,1:2,jclsaa:jclsbb),   &
        mfbnd(k,1:2,jclsaa:jclsbb) 
    end do
    end do



    return
    end subroutine parampollu_1clm_dumpaa



!-----------------------------------------------------------------------
    end module module_ecpp_td2clm
