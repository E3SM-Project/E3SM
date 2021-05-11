module ecpp_modal_wetscav

!-----------------------------------------------------------------
!  Module interface for cloud chemistry used in the ECPP treatment 
! in the MMF model
!  Adopted the similar one used in the ECPP 
!  for the WRF-chem model written by Dick Easter
!
!  Minghuai Wang, 2009-11
!------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use perf_mod
!==Guangxing Lin
!   use abortutils,   only: endrun
   use cam_abortutils,   only: endrun
!==Guangxing Lin

   implicit none

   public parampollu_tdx_wetscav_2

contains

!-----------------------------------------------------------------------
	subroutine parampollu_tdx_wetscav_2(                       &
		ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
		itstep_hybrid,                                     &
		idiagaa_ecpp, ldiagaa_ecpp, idiagbb_wetscav,       &
		tcen_bar, pcen_bar, rhocen_bar, dzcen,             &
!		rhobnd_bar, zbnd, wbnd_bar,                        &   not needed ?
!		chem_bar,                                          &   not needed ?
!		ncls_ecpp,                                         &
		it,      jt,      kts,ktebnd,ktecen,               &
		ncls_use,                                          &
!		kdraft_bot_use, kdraft_top_use,                    &   not needed ?
!		mtype_updnenv_use,                                 &   not needed ?
		chem_sub_new,                                      &
		del_chem_clm_wetscav,                              &
		del_wetscav3d, del_wetresu3d,                      &
!		ardz_cen_old, ardz_cen_new,                        &   not needed ?
		rhodz_cen,             &
		acen_tavg_use, acen_prec_use,                      &
		rh_sub2, qcloud_sub2, qlsink_sub2,                 &
		precr_sub2, precs_sub2,                            &
!		chem_bar_iccfactor,                                &   not needed ?
		activate_onoff_use,                                &
		iphase_of_aerosol, isize_of_aerosol,               &
		itype_of_aerosol, inmw_of_aerosol,                 &
		laicwpair_of_aerosol                               )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_wetscav_2 does wet scavenging of aerosols only
!    for one main-integ time sub-step
!
! incoming chem_sub_new holds current sub-class mixing ratios
! outgoing chem_sub_new holds updated sub-class mixing ratios
!
!-----------------------------------------------------------------------

!	use module_state_description, only:  p_qv, p_qc

!	use module_data_radm2, only:  epsilc

!	use module_data_mosaic_asect, only:  ai_phase, cw_phase,   &
!		massptr_aer, maxd_asize, maxd_atype,   &
!		ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer,   &
!		waterptr_aer
	use module_data_mosaic_asect, only:  &
		ai_phase, cw_phase,   &
		massptr_aer, maxd_asize, maxd_atype,   &
		ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer

	use module_data_ecpp1

!	use module_ecpp_hoststuff, only:  config_flags_ecpp

!	use module_mosaic_wetscav, only:  wetscav_cbmz_mosaic

!	use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message,   &
!	                             parampollu_1clm_set_opts

	implicit none

!   arguments
	integer, intent(in) ::                  &
		ktau, ktau_pp, itstep_sub,   &
		it, jt, kts, ktebnd, ktecen
	integer, intent(in) :: itstep_hybrid
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!	chem_driver and routines under it do calculations
!	over these spatial indices.

	integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199), &
		idiagbb_wetscav

	real(r8), intent(in) :: dtstep, dtstep_sub
!   dtstep - main model time step (s)
!   dtstep_sub - sub time step (s) currently used in ecpp main-integ routine

	real(r8), intent(in), dimension( kts:ktecen ) ::   &
		tcen_bar, pcen_bar, rhocen_bar, dzcen
!	real(r8), intent(in), dimension( kts:ktebnd ) ::   &
!		rhobnd_bar, wbnd_bar, zbnd
!   tcen_bar - temperature (K) at layer centers
!   rhocen_bar, rhobnd_bar - dry air density (kg/m^3) at layer centers and boundaries
!   pcen_bar - air pressure (Pa) at layer centers
!   wbnd_bar - vertical velocity (m/s) at layer boundaries
!   zbnd - elevation (m) at layer boundaries
!   dzcen - layer thicknesses (m)

!	real(r8), intent(in), dimension( kts:ktecen, 1:num_chem_ecpp ) :: &
!		chem_bar
!   chem_bar - mixing ratios of trace gase (ppm) and aerosol species
!	(ug/kg for mass species, #/kg for number species)

	integer, intent(in) :: ncls_use
!	integer, intent(in) :: ncls_ecpp

!	integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
!		kdraft_bot_use, kdraft_top_use,   &
!		mtype_updnenv_use

	real(r8), intent(inout), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
		chem_sub_new

	real(r8), intent(inout), dimension( 1:num_chem_ecpp ) :: &
		del_chem_clm_wetscav
! del_chem_clm_wetscav(l) = &
! sum( rhodz_cen(kts:ktecen) * ( del_wetscav3d(kts:ktecen,1:2,1:ncls_use,1:2,l) &
!                              + del_wetresu3d(kts:ktecen,1:2,1:ncls_use,1:2,l) ) )

	real(r8), intent(inout), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:num_chem_ecpp ) ::   &
		del_wetscav3d, del_wetresu3d
! del_wetscav3d = acen * (change to chem_sub due to uptake by precip)
!    the change for the current time sub-step is added to this array, so the array holds
!       the cummulative change over multiple time steps
!    this is always negative (or zero), and units are (kg/m^2)
! del_wetresu3d = acen * (change to chem_sub due to resuspension from precip evaporation)
!    this is always positive (or zero), and units are (kg/m^2)
!
! units for del_wetscav/resu3d will be (kg/m^2) or (#/m^2) in cam, 
!    where all tracer mixing ratios are (kg/kgair)
! in wrfchem, units are (ug/m^2) and (#/m^2) for aerosol mass and number
!    for gases, they are (mg/m^2) AFTER one applies a molecular weight ratio
! the important thing is that their sum is always equal to the column burden change,
!    where column burden = sum_over_k[ (mixing ratio)*(air density, kg/m^3)*(dz, m) ]

	real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) :: &
		acen_tavg_use, acen_prec_use
!		ardz_cen_old, ardz_cen_new, 

	real(r8), intent(inout), dimension( kts:ktecen ) :: rhodz_cen

	real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		rh_sub2, qcloud_sub2, qlsink_sub2, precr_sub2, precs_sub2

!	real(r8), intent(in), dimension( 1:2, num_chem_ecpp ) :: chem_bar_iccfactor

	integer, intent(in) :: activate_onoff_use

	integer, intent(in), dimension( 1:num_chem_ecpp ) ::   &
		iphase_of_aerosol, isize_of_aerosol, itype_of_aerosol,   &
		inmw_of_aerosol, laicwpair_of_aerosol



!   local variables
	integer, parameter :: nwdt = 1

	integer :: icc, icc_g, icc_l, iphase, ipp, ipp_l, ipp_g
	integer :: jcls, jcls_g, jcls_l
	integer :: k, kk, km1, kp1
	integer :: l, ll, lun142
	integer :: lgas_scav(1:num_chem_ecpp)
	integer :: m, mwdt
	integer :: n
	integer, parameter :: maxgas_scav = 4
	integer :: ngas_scav
	integer :: p1st
        integer :: inwdt

	logical :: skip_aer_resu, skip_gas_scav
	logical, dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		is_active, is_precp, is_ptgain, is_ptloss, is_rain
! is_active = .true. if sub-subarea has acen > afrac_cut_0p5
! is_precp  = .true. if sub-subarea has prtb > prsmall
! is_ptgain = .true. prtb increases from k+1 to k for the sub-subarea
! is_ptloss = .true. prtb decreases from k+1 to k for the sub-subarea
	logical, dimension( 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		ltmp_aa3d

	real(r8) :: delprtb_gtot, delprtb_ltot, delprtb_xtot
	real(r8) :: dt_scav
	real(r8) :: flxdt, flxdt_kp1
	real(r8) :: qgcx, qgcx_bgn
	real(r8) :: frac_scav
	real(r8) :: prsmall
	real(r8) :: rate_scav
	real(r8) :: scavcoef
	real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf, tmpq
	real(r8) :: tmpx, tmpx2, tmpy, tmpy2
	real(r8) :: tmpa1, tmpa2, tmpb1, tmpb2, tmpq1, tmpq2
	real(r8) :: tmp_ardzcen, tmpvol

	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:num_chem_ecpp ) ::   &
		chem_tmpa, chem_tmpb
	real(r8), dimension( 1:num_chem_ecpp ) :: curdel_chem_clm_wetscav
	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:max_wetdiagtype, 1:num_chem_ecpp ) ::   &
		delchem_wetscav, delchem_wetresu 
! delchem_wetscav = [ change to chem from wet scavenging over dt_scav ] ] * acen_tmp * rhodz_cen
!                   so units are (kg/m^2)
! delchem_wetresu = similar, but change from resuspension (due to precip evap)
	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:max_wetdiagtype, 1:num_chem_ecpp ) ::   &
		chem_prflxdt,  chem_prflxdt_xfer
! chem_prflxdt = [ downwards flux of precip-borne-tracers (kg/m^2/s) for subarea
!                  if it were spread over the entire host-code grid cell area ] * dt_scav
!                so units are (kg/m^2)
! chem_prflxdt_xfer = net transfer of chem_prflxdt into subarea from other subareas

	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) :: acen_tmp
! acen_tmp = fractional at layer centers for all 2 X 3 X 2 sub-subareas
	real(r8), dimension( kts:ktebnd, 1:2, 1:maxcls_ecpp, 1:2 ) :: prra, prsa, prta, prtb
! prta = total (liquid + solid) precip rate (kg/m^2/s) within the subarea
! prra, prsa =  liquid,  solid  precip rate (kg/m^2/s) within the subarea
! prtb = prta*acen_tmp = subarea precip rate 
!                           if it were spread over the entire host-code grid cell area
	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) :: &
		delprtb, delprtb_g, delprtb_l
! depprtb   = change in prtb from k+1 to k (kg/m^2/s)
! depprtb_g = increase in prtb from k+1 to k
! depprtb_l = abs( decrease in prtb from k+1 to k )
	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) :: frac_evap_prtb
	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, &
	                                 1:2, 1:maxcls_ecpp, 1:2 ) :: frac_xfer_prtb
	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, &
	                                 1:2, 1:maxcls_ecpp, 1:2 ) :: fxaa_evap_prtb
! frac_evap_prtb = fraction of precip (and precip-borne aerosols) entering the 
!                  top of a subarea that is evaporated/resuspended
! frac_xfer_prtb = fraction of precip (and precip-borne aerosols) entering the 
!                  top of a subarea that is transferred to another subarea
!                  (the first  set of icc,jcls,ipp indices are the "xfer from" subarea)
!                  (the second set of icc,jcls,ipp indices are the "xfer  to " subarea)

	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:maxd_asize, 1:maxd_atype ) ::   &
		scavcoef_num, scavcoef_vol
!	scavcoef_vol = below-cloud scavenging coeficient for volume (1/mm)
!	scavcoef_num = below-cloud scavenging coeficient for number (1/mm)
!       when precip rate = xxx kg/m2/s == xxx mm/s, the scavenging rate (1/s) = scavcoef*xxx 

	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, maxgas_scav ) ::   &
		gasscav_aa, gasscav_bb

	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		gasscav_cc


        call t_startf('ecpp_wetscav_init') 
!	write(*,'(a)') 'wetscav_2 doing part 1 stuff'

	lun142 = -1
	if (idiagaa_ecpp(142) > 0) lun142 = ldiagaa_ecpp(142)
	if (idiagbb_wetscav <= 0) lun142 = -1

	p1st = param_first_ecpp
	dt_scav = dtstep_sub

	mwdt = 1

	skip_gas_scav = .false.   ! flag for gas scavenging on/off
	if (wetscav_onoff_ecpp < 400) skip_gas_scav = .true.
	skip_aer_resu = .false.   ! flag for aerosol resuspension on/off
	if (wetscav_onoff_ecpp == 310) skip_aer_resu = .true.
	if (wetscav_onoff_ecpp == 410) skip_aer_resu = .true.

! load chem_tmpa array
        chem_tmpa = 0.0_r8
	do l = p1st, num_chem_ecpp
	do jcls = 1, ncls_use
	do icc = 1, 2
	do k = kts, ktecen
	    chem_tmpa(k,icc,jcls,1:2,l) = chem_sub_new(k,icc,jcls,l)
	end do
	end do
	end do
	end do
	chem_tmpb(:,:,:,:,:) = chem_tmpa(:,:,:,:,:)

	curdel_chem_clm_wetscav(:) = 0.0_r8
	delchem_wetscav(:,:,:,:,:,:) = 0.0_r8
	delchem_wetresu(:,:,:,:,:,:) = 0.0_r8
	chem_prflxdt(:,:,:,:,:,:) = 0.0_r8
	chem_prflxdt_xfer(:,:,:,:,:,:) = 0.0_r8


! precip rates -- 1.0 kgwtr/m^2/s = 1.0e-3 m3wtr/m^2/s = 1.0e-3 m/s
!     7.06e-5 kg/m^2/s = 7.06e-8  m/s = 0.01 inch/h
!     1.00e-7 kg/m^2/s = 1.00e-10 m/s = (0.01 inch/h) * 0.0014 is a very small precip rate!
	prsmall = 1.0e-7

! load precip rates for each icc,jcls,ipp subarea
	prta(:,:,:,:) = 0.0_r8
	prtb(:,:,:,:) = 0.0_r8
	prra(:,:,:,:) = 0.0_r8
	prsa(:,:,:,:) = 0.0_r8
	acen_tmp(:,:,:,:) = 0.0_r8

	is_active(:,:,:,:) = .false.
	is_precp(:,:,:,:) = .false.
	is_ptgain(:,:,:,:) = .false.
	is_ptloss(:,:,:,:) = .false.
	is_rain(:,:,:,:) = .false.

	do jcls = 1, ncls_use
	do icc = 1, 2
	do k = kts, ktecen
	    tmpa = max( 0.0_r8, acen_tavg_use(k,icc,jcls) )
	    tmpb = max( 0.0_r8, acen_prec_use(k,icc,jcls) )
	    tmpb = min( tmpa, tmpb )

	    if (tmpa <= afrac_cut_0p5) then   ! both ipp=1&2 have near-zero area
		continue
	    else if (tmpb <= afrac_cut_0p5) then   ! ipp=2 has near-zero area
		is_active(k,icc,jcls,1) = .true.
		acen_tmp(k,icc,jcls,1) = tmpa
		prta(k,icc,jcls,1) = precr_sub2(k,icc,jcls,1) + precs_sub2(k,icc,jcls,1)
		prtb(k,icc,jcls,1) = prta(k,icc,jcls,1)*acen_tmp(k,icc,jcls,1)
	    else if (tmpa-tmpb <= afrac_cut_0p5) then   ! ipp=1 has near-zero area
		is_active(k,icc,jcls,2) = .true.
		acen_tmp(k,icc,jcls,2) = tmpb
		prta(k,icc,jcls,2) = precr_sub2(k,icc,jcls,2) + precs_sub2(k,icc,jcls,2)
		prtb(k,icc,jcls,2) = prta(k,icc,jcls,2)*acen_tmp(k,icc,jcls,2)
	    else   ! both ipp=1&2 have areas > threshold
		is_active(k,icc,jcls,1) = .true.
		acen_tmp(k,icc,jcls,1) = tmpa-tmpb
		prta(k,icc,jcls,1) = precr_sub2(k,icc,jcls,1) + precs_sub2(k,icc,jcls,1)
		prtb(k,icc,jcls,1) = prta(k,icc,jcls,1)*acen_tmp(k,icc,jcls,1)
		is_active(k,icc,jcls,2) = .true.
		acen_tmp(k,icc,jcls,2) = tmpb
		prta(k,icc,jcls,2) = precr_sub2(k,icc,jcls,2) + precs_sub2(k,icc,jcls,2)
		prtb(k,icc,jcls,2) = prta(k,icc,jcls,2)*acen_tmp(k,icc,jcls,2)
	    end if

	    do ipp = 1, 2
		if ( is_active(k,icc,jcls,ipp) ) then
		    prtb(k,icc,jcls,ipp) = prta(k,icc,jcls,ipp)*acen_tmp(k,icc,jcls,ipp)
		    if (prtb(k,icc,jcls,ipp) > prsmall) then
			is_precp(k,icc,jcls,ipp) = .true.
			prsa(k,icc,jcls,ipp) = precs_sub2(k,icc,jcls,ipp)
			if (precr_sub2(k,icc,jcls,ipp)*acen_tmp(k,icc,jcls,ipp) > prsmall) then
			    prra(k,icc,jcls,ipp) = precr_sub2(k,icc,jcls,ipp)
			    is_rain(k,icc,jcls,ipp) = .true.
			end if
		    else
			prta(k,icc,jcls,ipp) = 0.0
			prtb(k,icc,jcls,ipp) = 0.0
		    end if
		end if
	    end do
	end do
	end do
	end do
        call t_stopf('ecpp_wetscav_init')


!
! calculate the fractions of precip (and precip-borne aerosols) 
! entering the top of a subarea that are either 
!	> evaporated/resuspended or 
!	> transferred to another subarea
!
        call t_startf('ecpp_wetscav_precip_evap')
	call wetscav_2_precip_evap_xfer(                           &
		ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
		idiagaa_ecpp, ldiagaa_ecpp, idiagbb_wetscav,       &
		it,      jt,      kts,ktebnd,ktecen,               &
		ncls_use,                                          &
		is_active, is_precp, is_ptgain, is_ptloss,         &
		acen_tmp, prtb, frac_evap_prtb, frac_xfer_prtb,    &
		fxaa_evap_prtb                                     )
        call t_stopf('ecpp_wetscav_precip_evap')


!
! calculate below-cloud scavenging coeficients for interstitial aerosols
!
        call t_startf('ecpp_wetscav_bcscav')
	call wetscav_2_bcscavcoef(                                 &
		ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
		idiagaa_ecpp, ldiagaa_ecpp, idiagbb_wetscav,       &
		tcen_bar, pcen_bar, rhocen_bar,                    &
		it,      jt,      kts,ktebnd,ktecen,               &
		ncls_use,                                          &
		rh_sub2,                                           &
		is_active, is_precp,                               &
		chem_tmpa, scavcoef_num, scavcoef_vol              )
        call t_stopf('ecpp_wetscav_bcscav')


!
! calculate stuff for below-cloud gas scavenging
!
        call t_startf('ecpp_wetscav_gascav')
	call wetscav_2_gasscav(                                    &
		ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
		dt_scav,                                           &
		idiagaa_ecpp, ldiagaa_ecpp, idiagbb_wetscav,       &
		tcen_bar, pcen_bar, rhocen_bar, dzcen,             &
		it,      jt,      kts,ktebnd,ktecen,               &
		ncls_use,                                          &
		is_active, is_precp, is_rain,                      &
		maxgas_scav, ngas_scav, lgas_scav,                 &
		acen_tmp, prra,                                    &
		qcloud_sub2, qlsink_sub2,                          &
		gasscav_aa, gasscav_bb, gasscav_cc                 )
        call t_stopf('ecpp_wetscav_gascav')


!
!
! now calculate
!    in-cloud & below-cloud aerosol wet removal
!    below-cloud resuspension from evaporating precip
!
!
        call t_startf('ecpp_wetscav_main')
wetscav_main_kloop_aa: &
	do k = ktecen, kts, -1


! set precip-borne_flux to that of layer above
	if (k < ktecen) then
	    chem_prflxdt(k,:,:,:,:,:) = chem_prflxdt(k+1,:,:,:,:,:)
	end if
	if (wetscav_onoff_ecpp < 200) cycle wetscav_main_kloop_aa


!
! do transfer of precip-borne tracers between subareas
! and resuspension from evaporation
!
	if (k < ktecen) then
! loop over the "losing" subareas
	    do jcls_l = 1, ncls_use
	    do ipp_l = 1, 2
	    do icc_l = 1, 2
		if ( .not. is_ptloss(k,icc_l,jcls_l,ipp_l) ) cycle

! loop over the "gaining" subareas,
! transferring chem_prflxdt from losing to gaining subarea
		do jcls_g = 1, ncls_use
		do ipp_g = 1, 2
		do icc_g = 1, 2
		    if ( .not. is_ptgain(k,icc_g,jcls_g,ipp_g) ) cycle
		    tmpa = frac_xfer_prtb(k,icc_l,jcls_l,ipp_l, icc_g,jcls_g,ipp_g)
		    if (tmpa <= 0.0_r8) cycle
		    do l = p1st, num_chem_ecpp
			if ( skip_gas_scav .and. (inmw_of_aerosol(l) <= 0)) cycle
			tmpb = chem_prflxdt(k+1,icc_l,jcls_l,ipp_l,mwdt,l)*tmpa
			chem_prflxdt(k  ,icc_g,jcls_g,ipp_g,mwdt,l) = &
			chem_prflxdt(k  ,icc_g,jcls_g,ipp_g,mwdt,l) + tmpb
			chem_prflxdt(k  ,icc_l,jcls_l,ipp_l,mwdt,l) = &
			chem_prflxdt(k  ,icc_l,jcls_l,ipp_l,mwdt,l) - tmpb

			chem_prflxdt_xfer(k  ,icc_g,jcls_g,ipp_g,mwdt,l) = &
			chem_prflxdt_xfer(k  ,icc_g,jcls_g,ipp_g,mwdt,l) + tmpb
			chem_prflxdt_xfer(k  ,icc_l,jcls_l,ipp_l,mwdt,l) = &
			chem_prflxdt_xfer(k  ,icc_l,jcls_l,ipp_l,mwdt,l) - tmpb
		    end do
		end do ! icc_g
		end do ! ipp_g
		end do ! jcls_g

! do resuspension from evaporation here
		tmpa = frac_evap_prtb(k,icc_l,jcls_l,ipp_l)
		if (tmpa <=  0.0_r8) cycle

		tmp_ardzcen = acen_tmp(k,icc_l,jcls_l,ipp_l)*rhodz_cen(k)
		do l = p1st, num_chem_ecpp
		    if ( skip_gas_scav .and. (inmw_of_aerosol(l) <= 0)) cycle
		    if ( skip_aer_resu .and. (inmw_of_aerosol(l) >  0)) cycle
		    tmpd = chem_prflxdt(k+1,icc_l,jcls_l,ipp_l,mwdt,l)*tmpa
		    delchem_wetresu(k,icc_l,jcls_l,ipp_l,mwdt,l) = &
		    delchem_wetresu(k,icc_l,jcls_l,ipp_l,mwdt,l) + tmpd
		    chem_prflxdt(k,icc_l,jcls_l,ipp_l,mwdt,l) = &
		    chem_prflxdt(k,icc_l,jcls_l,ipp_l,mwdt,l) - tmpd

		    if ( is_active(k,icc_l,jcls_l,ipp_l) ) then
! normally resuspend into (k,icc_l,jcls_l,ipp_l)
			chem_tmpb(k,icc_l,jcls_l,ipp_l,l) = &
			chem_tmpb(k,icc_l,jcls_l,ipp_l,l) + tmpd/tmp_ardzcen
		    else
! if (k,icc_l,jcls_l,ipp_l) is not active (acen_tmp ~= 0), then resuspend
!    uniformly across all active subareas
! (tmpd/rhodz_cen(k)) is the delta(chem) spread over the entire grid area
			do jcls_g = 1, ncls_use
			do ipp_g = 1, 2
			do icc_g = 1, 2
			    tmpf = fxaa_evap_prtb(k,icc_l,jcls_l,ipp_l, icc_g,jcls_g,ipp_g)
			    if (tmpf <= afrac_cut_0p5) cycle
			    chem_tmpb(k,icc_g,jcls_g,ipp_g,l) = &
			    chem_tmpb(k,icc_g,jcls_g,ipp_g,l) + tmpd/(tmpf*rhodz_cen(k))
			end do ! icc_g
			end do ! ipp_g
			end do ! jcls_g
		    end if
		end do ! l

	    end do ! icc_l
	    end do ! ipp_l
	    end do ! jcls_l
	end if ! (k < kte_cen)


!
! do additional resuspension for gases
!    currently gases are only in rain (none in solid precip), 
!        and the previous resuspension involves total precip
!    if rain ~= zero in a subarea, then resuspend any rainborne gases
!
	if ((k < ktecen) .and. ( .not. skip_gas_scav )) then
	    do jcls_l = 1, ncls_use
	    do ipp_l = 1, 2
	    do icc_l = 1, 2
		if ( is_rain(k,icc_l,jcls_l,ipp_l) ) cycle

		tmp_ardzcen = acen_tmp(k,icc_l,jcls_l,ipp_l)*rhodz_cen(k)
		if ( .not. is_active(k,icc_l,jcls_l,ipp_l) ) then
		    tmpf = 0.0
		    ltmp_aa3d(:,:,:) = .false.
		    do jcls_g = 1, ncls_use
		    do ipp_g = 1, 2
		    do icc_g = 1, 2
			if ( .not. is_active(k,icc_g,jcls_g,ipp_g) ) cycle
			if ((jcls_g == jcls_l) .and. &
			    (ipp_g == ipp_l) .and. (icc_g == icc_l)) cycle
			tmpf = tmpf + acen_tmp(k,icc_g,jcls_g,ipp_g)
			ltmp_aa3d(icc_g,jcls_g,ipp_g) = .true.
		    end do ! icc_g
		    end do ! ipp_g
		    end do ! jcls_g
		end if

		do ll = 1, ngas_scav
		    l = lgas_scav(ll)
		    if ((l < p1st) .or. (l > num_chem_ecpp)) cycle
		    tmpd = chem_prflxdt(k,icc_l,jcls_l,ipp_l,mwdt,l)
		    if (tmpd <= 0.0_r8) cycle

		    delchem_wetresu(k,icc_l,jcls_l,ipp_l,mwdt,l) = &
		    delchem_wetresu(k,icc_l,jcls_l,ipp_l,mwdt,l) + tmpd
		    chem_prflxdt(k,icc_l,jcls_l,ipp_l,mwdt,l) = 0.0_r8

		    if ( is_active(k,icc_l,jcls_l,ipp_l) ) then
! resuspend into (k,icc_l,jcls_l,ipp_l)
			chem_tmpb(k,icc_l,jcls_l,ipp_l,l) = &
			chem_tmpb(k,icc_l,jcls_l,ipp_l,l) + tmpd/tmp_ardzcen
		    else
! (k,icc_l,jcls_l,ipp_l) is not active, so resuspend across all active subareas
			do jcls_g = 1, ncls_use
			do ipp_g = 1, 2
			do icc_g = 1, 2
			    if ( .not. ltmp_aa3d(icc_g,jcls_g,ipp_g) ) cycle
			    chem_tmpb(k,icc_g,jcls_g,ipp_g,l) = &
			    chem_tmpb(k,icc_g,jcls_g,ipp_g,l) + tmpd/rhodz_cen(k)
			end do ! icc_g
			end do ! ipp_g
			end do ! jcls_g
		    end if
		end do ! ll

	    end do ! icc_l
	    end do ! ipp_l
	    end do ! jcls_l
	end if ! ((k < ktecen) .and. ( .not. skip_gas_scav ))


!
! calc in-cloud scavenging of activated aerosols
!
	do jcls = 1, ncls_use
	do ipp = 1, 2
	do icc = 1, 2
!	    cycle   ! *** skip for testing
	    if ( .not. is_active(k,icc,jcls,ipp) ) cycle
	    if ( .not. is_precp( k,icc,jcls,ipp) ) cycle

	    frac_scav = max( 0.0_r8, min( 1.0_r8, qlsink_sub2(k,icc,jcls,ipp)*dt_scav ) )
	    tmp_ardzcen = acen_tmp(k,icc,jcls,ipp)*rhodz_cen(k)

	    iphase = cw_phase
	    do n = 1, ntype_aer
	    do m = 1, nsize_aer(n)
	    do ll = 0, ncomp_aer(n)
		if (ll == 0) then
		    l = numptr_aer(m,n,iphase)
		else
		    l = massptr_aer(ll,m,n,iphase)
		end if
		if ((l < p1st) .or. (l > num_chem_ecpp)) cycle

		tmpa = frac_scav*chem_tmpb(k,icc,jcls,ipp,l)
		chem_tmpb(k,icc,jcls,ipp,l) = chem_tmpb(k,icc,jcls,ipp,l) - tmpa

		tmpb = tmpa*tmp_ardzcen
		delchem_wetscav(k,icc,jcls,ipp,mwdt,l) = &
		delchem_wetscav(k,icc,jcls,ipp,mwdt,l) - tmpb
		chem_prflxdt(k,icc,jcls,ipp,mwdt,l) = &
		chem_prflxdt(k,icc,jcls,ipp,mwdt,l) + tmpb
	    end do ! ll
	    end do ! m
	    end do ! n
	end do ! icc
	end do ! ipp
	end do ! jcls


!
! calc below-cloud scavenging of interstitial aerosols
!
	do jcls = 1, ncls_use
	do ipp = 1, 2
	do icc = 1, 2
!	    cycle   ! *** skip for testing
	    if ( .not. is_active(k,icc,jcls,ipp) ) cycle
	    if ( .not. is_precp( k,icc,jcls,ipp) ) cycle

	    tmp_ardzcen = acen_tmp(k,icc,jcls,ipp)*rhodz_cen(k)

	    iphase = ai_phase
	    do n = 1, ntype_aer
	    do m = 1, nsize_aer(n)
	    do ll = 0, ncomp_aer(n)
		if (ll == 0) then
		    l = numptr_aer(m,n,iphase)
		    scavcoef = scavcoef_num(k,icc,jcls,ipp,m,n)
		else
		    l = massptr_aer(ll,m,n,iphase)
		    scavcoef = scavcoef_vol(k,icc,jcls,ipp,m,n)
		end if
		if ((l < p1st) .or. (l > num_chem_ecpp)) cycle
!		scavcoef = 0.01_r8  ! use simple constant value
!		scavcoef = 0.0_r8   ! turn off below-cloud scav

		rate_scav = prta(k,icc,jcls,ipp)*scavcoef
		frac_scav = 1.0_r8 - exp( -rate_scav*dt_scav )
		frac_scav = max( 0.0_r8, min( 1.0_r8, frac_scav ) )

		tmpa = frac_scav*chem_tmpb(k,icc,jcls,ipp,l)
		chem_tmpb(k,icc,jcls,ipp,l) = chem_tmpb(k,icc,jcls,ipp,l) - tmpa

		tmpb = tmpa*tmp_ardzcen
		delchem_wetscav(k,icc,jcls,ipp,mwdt,l) = &
		delchem_wetscav(k,icc,jcls,ipp,mwdt,l) - tmpb
		chem_prflxdt(k,icc,jcls,ipp,mwdt,l) = &
		chem_prflxdt(k,icc,jcls,ipp,mwdt,l) + tmpb
	    end do ! ll
	    end do ! m
	    end do ! n
	end do ! icc
	end do ! ipp
	end do ! jcls


!
! calc gas scavenging
!
	if ( .not. skip_gas_scav ) then
	    do jcls = 1, ncls_use
	    do ipp = 1, 2
	    do icc = 1, 2
!	    cycle   ! *** skip for testing
		if ( .not. is_rain(k,icc,jcls,ipp) ) cycle
		tmp_ardzcen = acen_tmp(k,icc,jcls,ipp)*rhodz_cen(k)

		do ll = 1, ngas_scav
		    l = lgas_scav(ll)
		    if ((l < p1st) .or. (l > num_chem_ecpp)) cycle

		    flxdt_kp1 = chem_prflxdt(k,icc,jcls,ipp,mwdt,l)
		    qgcx_bgn = chem_tmpb(k,icc,jcls,ipp,l)
		    tmpa = gasscav_aa(k,icc,jcls,ipp,ll)
		    tmpb = gasscav_bb(k,icc,jcls,ipp,ll)
		    tmpc = gasscav_cc(k,icc,jcls,ipp)
		    tmpe = tmpb + tmpc + tmpa*tmpc

! this is the solution to the 2 final equations in subr wetscav_2_gasscav
		    flxdt = flxdt_kp1*((1.0_r8 + tmpa)*tmpc/tmpe) + qgcx_bgn*(tmpa/tmpe)
		    qgcx = qgcx_bgn*((1.0_r8 + tmpa*(tmpb/tmpe))/(1.0_r8 + tmpa)) &
		         + flxdt_kp1*(tmpc*(tmpb/tmpe))

		    chem_tmpb(k,icc,jcls,ipp,l) = qgcx
		    chem_prflxdt(k,icc,jcls,ipp,mwdt,l) = flxdt
		    tmpf = (qgcx - qgcx_bgn)*tmp_ardzcen
		    if (tmpf > 0.0_r8) then
			delchem_wetresu(k,icc,jcls,ipp,mwdt,l) = &
			delchem_wetresu(k,icc,jcls,ipp,mwdt,l) + tmpf
		    else
			delchem_wetscav(k,icc,jcls,ipp,mwdt,l) = &
			delchem_wetscav(k,icc,jcls,ipp,mwdt,l) + tmpf
		    end if
		end do ! ll

	    end do ! icc
	    end do ! ipp
	    end do ! jcls
	end if ! ( .not. skip_gas_scav )


	end do wetscav_main_kloop_aa
        call t_stopf('ecpp_wetscav_main')


        call t_startf('ecpp_wetscav_endcopy')
!
! load new chem mixratios into chem_sub_new (only if wetscav_onoff_ecpp >= 300)
! calc overall changes to column burdens    (only if wetscav_onoff_ecpp >= 200) 
!
	if (wetscav_onoff_ecpp >= 200) then

	do l = p1st, num_chem_ecpp
	    if ( skip_gas_scav .and. (inmw_of_aerosol(l) <= 0)) cycle
	    tmpx = 0.0_r8 ; tmpx2 = 0.0_r8
	    do k = kts, ktecen
		tmpy = 0.0_r8 ; tmpy2 = 0.0_r8
		do jcls = 1, ncls_use
		do icc = 1, 2
		    tmpb = 0.0_r8
		    tmpc = 0.0_r8
		    do ipp = 1, 2
		        tmpa = acen_tmp(k,icc,jcls,ipp)
		        if ( is_active(k,icc,jcls,ipp) ) then
			    tmpy = tmpy + tmpa*(chem_tmpb(k,icc,jcls,ipp,l) &
			                      - chem_tmpa(k,icc,jcls,ipp,l))
			    tmpb = tmpb + tmpa*chem_tmpb(k,icc,jcls,ipp,l)
			    tmpc = tmpc + tmpa
		        end if
                        tmpd = 0.0_r8
                        do inwdt=1, nwdt
			   tmpd = tmpd + delchem_wetscav(k,icc,jcls,ipp,inwdt,l) /rhodz_cen(k)
                        end do
			del_wetscav3d(k,icc,jcls,ipp,l) = del_wetscav3d(k,icc,jcls,ipp,l) + tmpd
                        tmpe = 0.0_r8
                        do inwdt=1, nwdt
                           tmpe = tmpe + delchem_wetresu(k,icc,jcls,ipp,inwdt,l) /rhodz_cen(k)
                        end do
			del_wetresu3d(k,icc,jcls,ipp,l) = del_wetresu3d(k,icc,jcls,ipp,l) + tmpe
			tmpy2 = tmpy2 + tmpd + tmpe
		    end do ! ipp
		    if ((acen_tavg_use(k,icc,jcls) > afrac_cut_0p5) .and. &
		        (tmpc > 0.0_r8) .and. (wetscav_onoff_ecpp >= 300)) then
			chem_sub_new(k,icc,jcls,l) = max( 0.0_r8, tmpb )/tmpc
		    end if
		end do ! icc
		end do ! jcls
		tmpx = tmpx + tmpy*rhodz_cen(k)
		tmpx2 = tmpx2 + tmpy2*rhodz_cen(k)
	    end do ! k
	    curdel_chem_clm_wetscav(l) = tmpx
	    ! *** increment del_chem_clm_wetscav with tmpx2 (new way) 
            !     instead of tmpx (old way)
	    del_chem_clm_wetscav(l) = del_chem_clm_wetscav(l) + tmpx2
	end do ! l

	end if ! (wetscav_onoff_ecpp >= 200)
        call t_stopf('ecpp_wetscav_endcopy')
   
        call t_startf('ecpp_wetscav_enddiag')
!
! diagnostic checks on the new arrays to see that they are "making sense"
!
	if (lun142 > 0) then

	do l = p1st, num_chem_ecpp

	write(lun142,'(//a,i5)') 'diags for species l =', l

	if (lun142 == -999888777) then   ! *** skip for testing

	write(lun142,'(a,i5)') 'chem_tmpa for icc=ipp=2 & grid-avg;  chem_tmpb ...;  b-a ...'
	icc = 2 ; ipp = 2
	do k = ktecen, kts, -1
	    write(lun142,'(i3,1p,3(2x,4e10.3))') k, &
		   ( chem_tmpa(k,icc,jcls,ipp,l), jcls=1,ncls_use ), &
		sum( chem_tmpa(k,1:2,1:ncls_use,1:2,l)* &
		      acen_tmp(k,1:2,1:ncls_use,1:2) ), &
		   ( chem_tmpb(k,icc,jcls,ipp,l), jcls=1,ncls_use ), &
		sum( chem_tmpb(k,1:2,1:ncls_use,1:2,l)* &
		      acen_tmp(k,1:2,1:ncls_use,1:2) ), &
		   ( (chem_tmpb(k,icc,jcls,ipp,l) - &
		      chem_tmpa(k,icc,jcls,ipp,l)), jcls=1,ncls_use ), &
		sum( ( chem_tmpb(k,1:2,1:ncls_use,1:2,l) - &
		       chem_tmpa(k,1:2,1:ncls_use,1:2,l) )* &
		        acen_tmp(k,1:2,1:ncls_use,1:2) )
	end do

	write(lun142,'(/a,i5)') &
	'delchem_wetscav for icc=ipp=2 & grid-avg;  delchem_wetresu ...;  chem_prflxdt_xfer ...'
	icc = 2 ; ipp = 2
	do k = ktecen, kts, -1
	    write(lun142,'(i3,1p,3(2x,4e10.3))') k, &
		   ( delchem_wetscav(  k,icc,jcls,ipp,mwdt,l), jcls=1,ncls_use ), &
		sum( delchem_wetscav(  k,1:2,1:ncls_use,1:2,1:nwdt,l) ), &
		   ( delchem_wetresu(  k,icc,jcls,ipp,mwdt,l), jcls=1,ncls_use ), &
		sum( delchem_wetresu(  k,1:2,1:ncls_use,1:2,1:nwdt,l) ), &
		   ( chem_prflxdt_xfer(k,icc,jcls,ipp,mwdt,l), jcls=1,ncls_use ), &
		sum( chem_prflxdt_xfer(k,1:2,1:ncls_use,1:2,1:nwdt,l) )
	end do

	write(lun142,'(/a,i5)') &
	'chem_prflxdt for icc=ipp=2 & grid-avg;  conserve check stuff'
	icc = 2 ; ipp = 2
	do k = ktecen, kts, -1
	    kp1 = min(k+1,ktecen) ; tmpa = kp1 - k
	    write(lun142,'(i3,1p,3(2x,4e10.3))') k, &
		   ( chem_prflxdt(   k,icc,jcls,ipp,mwdt,l), jcls=1,ncls_use ), &
		sum( chem_prflxdt(   k,1:2,1:ncls_use,1:2,1:nwdt,l) ), &
		   ( chem_prflxdt(   kp1,icc,jcls,ipp,mwdt,l)*tmpa &
		   - chem_prflxdt(     k,icc,jcls,ipp,mwdt,l) &
		   - delchem_wetscav(  k,icc,jcls,ipp,mwdt,l) &
		   - delchem_wetresu(  k,icc,jcls,ipp,mwdt,l) &
		   + chem_prflxdt_xfer(k,icc,jcls,ipp,mwdt,l), jcls=1,ncls_use ), &
		sum( chem_prflxdt(   kp1,1:2,1:ncls_use,1:2,1:nwdt,l)*tmpa &
		   - chem_prflxdt(     k,1:2,1:ncls_use,1:2,1:nwdt,l) &
		   - delchem_wetscav(  k,1:2,1:ncls_use,1:2,1:nwdt,l) &
		   - delchem_wetresu(  k,1:2,1:ncls_use,1:2,1:nwdt,l) &
		   + chem_prflxdt_xfer(k,1:2,1:ncls_use,1:2,1:nwdt,l) )
	end do

	end if ! (lun142 == -999888777) 

	write(lun142,'(/2a,i5)') &
	'sum( delchem_wetscav ),  sum( delchem_wetresu ),  sum( both ),', &
	'  curdel_chem_clm_wetscav,  (4)-(5)/max(...)'
	tmpa = sum( delchem_wetscav(  kts:ktecen,1:2,1:ncls_use,1:2,1:nwdt,l) )
	tmpb = sum( delchem_wetresu(  kts:ktecen,1:2,1:ncls_use,1:2,1:nwdt,l) )
	tmpc = tmpa + tmpb
	tmpd = curdel_chem_clm_wetscav(l)
	tmpe = (tmpc - tmpd)/max( abs(tmpc), abs(tmpd), 1.0e-38_r8 )
	write(lun142,'(1p,3(2x,2e11.3))') &
	    tmpa, tmpb, tmpc, tmpd, tmpe
!	if (l == 2) write(lun142,'(3a)') 'qakee - ktau, it_hyb, it_sub, l', &
!	'sum( delchem_wetscav ),  sum( delchem_wetresu ),  sum( both ),', &
!	'  curdel_chem_clm_wetscav,  (4)-(5)/max(...)'
!	if (l >= 39) write(lun142,'(a,4i4,1p,3(2x,2e11.3))') &
!	    'qakee', ktau, itstep_hybrid, itstep_sub, l, &
!	    tmpa, tmpb, tmpc, tmpd, tmpe

	write(lun142,'(/2a,i5)') &
	'sum( del_wetscav3d ),  sum( del_wetresu3d ),  sum( both ),', &
	'  del_chem_clm_wetscav,  (4)-(5)/max(...)'
	tmpa = 0.0_r8 ; tmpb = 0.0_r8
	do k = kts, ktecen
	    tmpa = tmpa + sum( del_wetscav3d(k,1:2,1:ncls_use,1:2,l) ) * rhodz_cen(k)
	    tmpb = tmpb + sum( del_wetresu3d(k,1:2,1:ncls_use,1:2,l) ) * rhodz_cen(k)
	end do
	tmpc = tmpa + tmpb
	tmpd = del_chem_clm_wetscav(l)
	tmpe = (tmpc - tmpd)/max( abs(tmpc), abs(tmpd), 1.0e-38_r8 )
	write(lun142,'(1p,3(2x,2e11.3))') &
	    tmpa, tmpb, tmpc, tmpd, tmpe
!	if (l == 2) write(lun142,'(3a)') 'qakff - ktau, it_hyb, it_sub, l', &
!	'sum( del_wetscav3d ),  sum( del_wetresu3d ),  sum( both ),', &
!	'  del_chem_clm_wetscav,  (4)-(5)/max(...)'
!	if (l >= 39) write(lun142,'(a,4i4,1p,3(2x,2e11.3))') &
!	    'qakff', ktau, itstep_hybrid, itstep_sub, l, &
!	    tmpa, tmpb, tmpc, tmpd, tmpe

	end do ! l

	write(lun142,'(//a,i5)') 'qlsink*dt_scav for icc=ipp=2;  qcloud ...;  ardzcen ...'
	icc = 2 ; ipp = 2
	do k = ktecen, kts, -1
	    write(lun142,'(i3,1p,4(2x,3e10.3))') k, &
		( qlsink_sub2(k,icc,jcls,ipp)*dt_scav, jcls=1,ncls_use ), &
		( qcloud_sub2(k,icc,jcls,ipp), jcls=1,ncls_use ), &
		( acen_tmp(k,icc,jcls,ipp)*rhodz_cen(k), jcls=1,ncls_use )
	end do

	write(lun142,'(//a,i5)') 'prta for icc=ipp=2;  prtb ...;  delprtb ...'
	icc = 2 ; ipp = 2
	do k = ktecen, kts, -1
	    write(lun142,'(i3,1p,4(2x,3e10.3))') k, &
		( prta(k,icc,jcls,ipp), jcls=1,ncls_use ), &
		( prtb(k,icc,jcls,ipp), jcls=1,ncls_use ), &
		( prtb(k,icc,jcls,ipp)-prtb(k+1,icc,jcls,ipp), jcls=1,ncls_use )
	end do

	end if ! (lun142 > 0)

        call t_stopf('ecpp_wetscav_enddiag')


!	write(*,'(a)') 'wetscav_2 DONE'
	return
	end subroutine parampollu_tdx_wetscav_2



!-----------------------------------------------------------------------
	subroutine wetscav_2_gasscav(                              &
		ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
		dt_scav,                                           &
		idiagaa_ecpp, ldiagaa_ecpp, idiagbb_wetscav,       &
		tcen_bar, pcen_bar, rhocen_bar, dzcen,             &
		it,      jt,      kts,ktebnd,ktecen,               &
		ncls_use,                                          &
		is_active, is_precp, is_rain,                      &
		maxgas_scav, ngas_scav, lgas_scav,                 &
		acen_tmp, prra,                                    &
		qcloud_sub2, qlsink_sub2,                          &
		gasscav_aa, gasscav_bb, gasscav_cc                 )


!-----------------------------------------------------------------------
! DESCRIPTION
!
! wetscav_2_gasscav does pre-calculations for in-cloud and below-cloud 
!    of gases (h2o2, so2, and nh3) by rain
! the results are applied in subr parampollu_tdx_wetscav_2
!
! main assumptions
!    reversible scavenging of gases
!    prescribed pH for rainwater and cloudwater
!    no aqueous phase reactions are treated here
!-----------------------------------------------------------------------

!	use module_state_description, only:  p_qv, p_qc

!	use module_data_radm2, only:  epsilc

	use module_data_mosaic_asect, only:  &
		ai_phase, dens_aer, hygro_aer,   &
		massptr_aer, maxd_asize, maxd_atype,   &
		ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer,   &
		dcen_sect, dhi_sect, dlo_sect, sigmag_aer,   &
		volumhi_sect, volumlo_sect

	use module_data_ecpp1

	use constituents, only:  cnst_get_ind

	use module_ecpp_util, only:  ecpp_error_fatal

	implicit none

!   arguments
!   ( for definitions see subr parampollu_tdx_wetscav_2 )
	integer, intent(in) ::                  &
		ktau, ktau_pp, itstep_sub,   &
		it, jt, kts, ktebnd, ktecen

	integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199), &
		idiagbb_wetscav

	real(r8), intent(in) :: dtstep, dtstep_sub, dt_scav

	real(r8), intent(in), dimension( kts:ktecen ) ::   &
		tcen_bar, pcen_bar, rhocen_bar, dzcen

	integer, intent(in) :: ncls_use

	logical, intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		is_active, is_precp, is_rain

	integer, intent(in) :: maxgas_scav
	integer, intent(out)  :: ngas_scav
	integer, intent(out), dimension( 1:maxgas_scav ) ::   &
		lgas_scav

	real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		acen_tmp, qcloud_sub2, qlsink_sub2

	real(r8), intent(in),  dimension( kts:ktebnd, 1:2, 1:maxcls_ecpp, 1:2 ) :: &
		prra

	real(r8), intent(out), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, maxgas_scav ) ::   &
		gasscav_aa, gasscav_bb

	real(r8), intent(out), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		gasscav_cc



!   local variables
	integer :: icc, ipp
	integer :: itmpa
	integer :: jcls
	integer :: k
	integer :: ll, lun143
	integer :: m
	integer :: n
	integer :: p1st

!	real(r8), parameter :: piover6 = 3.14159265358979323846_r8/6.0_r8
	real(r8), parameter :: qcldwtr_cutoff = 1.0e-6
	real(r8), parameter :: tmp8over9 = 8.0_r8/9.0_r8

	real(r8) :: frac_c, frac_g
	real(r8) :: hen1c(maxgas_scav), hen1r(maxgas_scav)
	real(r8) :: hen2c(maxgas_scav), hen2r(maxgas_scav)
	real(r8) :: heffcx(maxgas_scav), heffrx(maxgas_scav)
	real(r8) :: hionc, hionr
	real(r8) :: kxf_cr, kxf_gcr, kxf_gr(maxgas_scav)
	real(r8) :: qcwtr, qrwtr
	real(r8) :: scavrate_hno3
	real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf, tmpg
	real(r8) :: vfallr


!   pointers for gases that are scavenged
	p1st = param_first_ecpp

	ngas_scav = 4
	if (ngas_scav > maxgas_scav) then
	    write(*,*) 'subr wetscav_2_gasscav -- ngas_scav > maxgas_scav', &
		ngas_scav > maxgas_scav
	    call ecpp_error_fatal( lunout, &
		'subr wetscav_2_gasscav -- ngas_scav > maxgas_scav' )
	end if

	lgas_scav(:) = -1

	call cnst_get_ind(                 'so2', itmpa, .false. )
	if (itmpa <= 0) call cnst_get_ind( 'SO2', itmpa, .false. )
	if ((itmpa >= p1st) .and. (itmpa <= num_chem_ecpp)) lgas_scav(1) = itmpa

	call cnst_get_ind(                 'h2o2', itmpa, .false. )
	if (itmpa <= 0) call cnst_get_ind( 'H2O2', itmpa, .false. )
	if ((itmpa >= p1st) .and. (itmpa <= num_chem_ecpp)) lgas_scav(2) = itmpa

	call cnst_get_ind(                 'nh3', itmpa, .false. )
	if (itmpa <= 0) call cnst_get_ind( 'NH3', itmpa, .false. )
	if ((itmpa >= p1st) .and. (itmpa <= num_chem_ecpp)) lgas_scav(3) = itmpa

	call cnst_get_ind(                 'h2so4', itmpa, .false. )
	if (itmpa <= 0) call cnst_get_ind( 'H2SO4', itmpa, .false. )
	if ((itmpa >= p1st) .and. (itmpa <= num_chem_ecpp)) lgas_scav(4) = itmpa

!	write(*,'(a,10i5)') 'wetscav_2_gasscav - ngas_scav', ngas_scav
!	write(*,'(a,10i5)') 'wetscav_2_gasscav - lgas_scav', lgas_scav(1:maxgas_scav)
!	if (ngas_scav /= -13579) stop


!
!   treatment of gas scavenging (by rain)
!
!   primary assumptions are
!      gases are reversibly scavenging in rain (e.g., transfer from gas to rain
!         and transfer from rain to gas are both treated)
!      rainborne gases are treated a locally steady-state, but vary with height
!      cloudborne gases in equilibrium with the "interstitial gases" 
!         and are collected by rain
!      pH for the cloud and rainwater are prescribed
!      aqueous chemical reaction in rain are not treated
!
!   define
!      qrx = mixing ratio of rainborne  species x (kg-x/kg-air)
!      qcx = mixing ratio of cloudborne species x (kg-x/kg-air)
!      qgx = mixing ratio of gaseous    species x (kg-x/kg-air)
!      qgcx = qgx + qcx
!
!      the above are defined for each vertical layer and each ecpp subarea
!      (in wrf-chem, they are units are actually mg-x/kg-air after a molecular weight
!       ratio is applied, but the equations work anyway)
!
!   basic equations:
!
!      d[acen*rho*qgcx]/dt = acen*rho*[ -kxf_gr*(qgx - qrx/heffrx) - kxf_ct*qcx ]
!
!      d[acen*rho*vfallr*qrx]/dz = acen*rho*[ -kxf_gr*(qgx - qrx/heffrx) - kxf_ct*qcx ]
!
!      qcx = heffcx*qgx
!
!   where
!      acen = fractional area of subarea
!      rho = air density (kg-air/m^3)
!      vfallr = rain fall velocity (m/s, and positive)
!      kxf_gr = mass transfer coefficient for gas <--> rain (1/s)
!             a power-law curve fit to Schwarz and Levine (19xx) is used
!      kxf_ct = rate of collection of cloudwater by rainwater (1/s) == qlsink
!      heffrx, heffcx = gaseous-rainborne and gaseous-cloudborne equilibirum partitioning 
!             coefficients (i.e., modified effective henry law constants) with units of
!             [(mol-x/kg-h2o)/(mol-x/kg-air)] == [(kg-x/kg-h2o)/(kg-x/kg-air)]
!
!   define
!      frac_c = heffcx/(1 + heffgx)  so qcx = frac_c*qgcx
!      frac_g = 1 - frac_c           so qgx = frac_g*qgcx
!      kxf_gcr = frac_g*kxf_gr + frac_c*kxf_cr
!
!   then
!
!      d[acen*rho*qgcx]/dt = acen*rho*[ -kxf_gcr*qgcx + kxf_gr*qrx/heffrx) ]
!
!      d[acen*rho*vfallr*qrx]/dz = acen*rho*[ -kxf_gcr*qgcx + kxf_gr*qrx/heffrx) ]
!
!   define
!      dt = time step ( = ecpp sub time step )
!      flxdt = acen*rho*vfallr*qrx*dt = chem_prflxdt of subr parampollu_tdx_wetscav_2
!
!   then
!
!      d[acen*rho*qgcx]/dt = -[acen*rho*kxf_gcr]*qgcx + [kxf_gr/(heffrx*vfallr*dt)]*flxdt
!
!      d[flxdt]/dz = -[acen*rho*kxf_gcr*dt]*qgcx + [kxf_gr/(heffrx*vfallr)]*flxdt
!
!   now define
!      dt = time step (s)
!      dz = thickness of layer k (m)
!      qgcx     = qgcx in layer k at end of time step
!      qgcx_bgn = qgcx in layer k at beginning of time step
!      flxdt     = flxdt in layer k   at end of time step
!      flxdt_kp1 = flxdt in layer k+1 at end of time step
!
!   and use the following finite differencing which is implicit in time
!
!      (acen*rho)*(qgcx - qgcx_o)/dt = -[acen*rho*kxf_gcr]*qgcx + [kxf_gr/(heffrx*vfallr*dt)]*flxdt
!   which yields
!      qgcx*[1 + kxf_gcr*dt] + flxdt*[-kxf_gr/(heffrx*vfallr*acen*rho)] = qgcx_bgn
!    
!      (flxdt+kp1 - flxdt)/dz = -[acen*rho*kxf_gcr*dt]*qgcx + [kxf_gr/(heffrx*vfallr)]*flxdt
!   which yields
!      qgcx*[-kxf_gcr*dt] + flxdt*[1/(dz*acen*rho) + kxf_gr/(heffrx*vfallr*acen*rho)] = flxdt_kp1*[1/(dz*acen*rho)]
!
!   define
!      aa = kxf_gcr*dt
!      bb = kxf_gr/(heffrx*vfallr*acen*rho)
!      cc = 1/(dz*acen*rho)
!
!   then
!      qgcx*[1 + aa] + flxdt*[-bb] = qgcx_bgn
!      qgcx*[-aa] + flxdt*[cc + bb] = flxdt_kp1*[cc]
!
!   these 2 equations are solved in the gas-scavenging section of subr parampollu_tdx_wetscav_2,
!   starting at ktecen (where flxdt_kp1 = )
!   the purpose of this routine (subr wetscav_2_gasscav) is to provide the aa, bb, and cc
!


	lun143 = -1
	if (idiagaa_ecpp(143) > 0) lun143 = ldiagaa_ecpp(143)
	if (idiagbb_wetscav /= 1) lun143 = -1

!   hionr, hionc = prescribed hydrogen ion concentrations (mol/liter-h2o)
!   for rainwater and cloudwater
	hionr = 10.0_r8**(-5.0_r8)
	hionc = 10.0_r8**(-4.5_r8)

!   calculate information needed for the gas scavenging equations
main_kloop_aa: &
	do k = kts, ktecen

	do ipp = 1, 2
	do jcls = 1, ncls_use
	do icc = 1, 2

	if ( .not. is_rain(k,icc,jcls,ipp) ) cycle
	if (lun143 > 0) write(lun143,'(/a,5i5)') 'wetscav_2_gasscav', &
		ktau, k, icc, jcls, ipp
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') 'aaaa stuff  ', &
	    tcen_bar(k), pcen_bar(k), rhocen_bar(k), dzcen(k), dt_scav


!   calculate rain fallspeed and rainwater mixing ratio using Kessler (1969)
	tmpa = prra(k,icc,jcls,ipp)   ! rain precip rate (kg/m^2/s)
	tmpb = sqrt( 1.22_r8/rhocen_bar(k) )   ! density factor for fallspeed
!   tmpc = first guess rain water conc (kg/m^3) from Kessler (1969)
	tmpc = (tmpa/(12.11_r8*tmpb))**tmp8over9   
!   vfallr = rain mean fallspeed (m/s) from its definition, but forced to >= 1 m/s
	vfallr = max( 1.0_r8, (tmpa/tmpc) )   
!   qrwtr = rain water mixing ratio (kg/kgair) from its definition
	qrwtr = tmpa/(vfallr*rhocen_bar(k))   
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') 'rain stuff  ', &
	    prra(k,icc,jcls,ipp), acen_tmp(k,icc,jcls,ipp), &
		tmpa, tmpb, tmpc, vfallr, qrwtr

!   qcwtr = cloud water mixing ratio (kg/kgair) from its definition
	qcwtr = qcloud_sub2(k,icc,jcls,ipp)
	if (qcwtr > qcldwtr_cutoff) then
	    kxf_cr = max( 0.0_r8, qlsink_sub2(k,icc,jcls,ipp) )
	else
	    qcwtr = 0.0_r8
	    kxf_cr = 0.0_r8
	end if


!   gas-liquid partitioning coefficients
!
!   hen1 = effective henry law constant at prescribed ph
!   [(mol-x/liter-h2o)/atm] = [(mol-x/kg-h2o)/atm]
	hen1r(:) = 0.0_r8
	hen1c(:) = 0.0_r8
	tmpa = (1.0_r8/tcen_bar(k)) - (1.0_r8/298.16_r8)
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') '0000 hen1   ', &
		tcen_bar(k), tmpa, qrwtr, qcwtr
!   so2
	tmpb = 1.23_r8*exp(3150.0_r8*tmpa)     ! henry law constant
	tmpc = 1.3e-2_r8*exp(1960.0_r8*tmpa)   ! 1st dissociation constant
	hen1r(1) = tmpb*(1.0_r8 + tmpc/hionr)  ! effective henry
	hen1c(1) = tmpb*(1.0_r8 + tmpc/hionc)  ! effective henry
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') 'so2  hen1   ', &
		tmpb, tmpc, hen1r(1), hen1c(1)
!   h2o2
	tmpb = 7.45e4_r8*exp(7300.0_r8*tmpa)   ! henry law constant
	hen1r(2) = tmpb
	hen1c(2) = tmpb
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') 'h2o2 hen1   ', &
		tmpb, 0.0, hen1r(2), hen1c(2)
!+++mhwang
! set hen1r and hen1d of so2 to be the same as H2O2, which is what used 
! in the conventional NCAR CAM. 
! Minghuai Wang (Minghuai.Wang@pnl.gov), 2010-02
!        hen1r(1) = hen1r(2)
!        hen1c(1) = hen1c(2)
!---mhwang

!   nh3
	tmpb = 6.21e1_r8*exp(4110.0_r8*tmpa)     ! henry law constant
	tmpc = 1.7e-5_r8*exp(-450.0_r8*tmpa)     ! 1st dissociation constant
	tmpd = 1.0e-14_r8*exp(-6710.0_r8*tmpa)   ! water dissociation constant
	hen1r(3) = tmpb*(1.0_r8 + (tmpc/tmpd)*hionr)  ! effective henry
	hen1c(3) = tmpb*(1.0_r8 + (tmpc/tmpd)*hionc)  ! effective henry
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') 'nh3  hen1   ', &
		tmpb, tmpc, hen1r(3), hen1c(3)
!   h2so4 (values are from CAPRAM website)
	tmpb = 8.7e11_r8                       ! henry law constant
	tmpc = 1.0e3_r8                        ! 1st dissociation constant
	hen1r(4) = tmpb*(1.0_r8 + tmpc/hionr)  ! effective henry
	hen1c(4) = tmpb*(1.0_r8 + tmpc/hionc)  ! effective henry
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') 'h2so4  hen1 ', &
		tmpb, tmpc, hen1r(4), hen1c(4)

!   hen2 = like hen1 but units = [(mol-x/kg-h2o)/(mol-x/kg-air)]
!   ax atm of x = ax*p0 Pa of x = ax*p0/pair (mol-x/mol-air) 
!                               = ax*p0/(pair*0.029) (mol-x/kg-air) 
	tmpa = (pcen_bar(k)/1.01325e5)*0.028966
	hen2r(1:ngas_scav) = hen1r(1:ngas_scav)*tmpa
	hen2c(1:ngas_scav) = hen1c(1:ngas_scav)*tmpa

!   heffrx,cx units = [(mol-x/kg-air)/(mol-x/kg-air)] and includes 
!   rainwater,cloudwater mixing ratio factor
	heffrx(1:ngas_scav) = hen2r(1:ngas_scav)*qrwtr
	heffcx(1:ngas_scav) = hen2c(1:ngas_scav)*qcwtr
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') 'heffrx,cx   ', &
		heffrx(1:4), heffcx(1:4)


!   gas-rain mass transfer rates
!
!   scavrate_hno3 = rain scavenging rate for hno3 (1/s)
!   this is power law fit to levine and schwartz (1982, atmos environ) 
!      results, with temperature and pressure adjustments 
	tmpa = prra(k,icc,jcls,ipp)*3600.0   ! precip rate in mm/hr = kg/m^2/hr
        scavrate_hno3 = 6.262e-5*(tmpa**0.7366)   &
                * ((tcen_bar(k)/298.0)**1.12)   &
                * ((1.01325e5/pcen_bar(k))**.75)
!   for other gases, multiply hno3 rate by ratio of gas diffusivities
        kxf_gr(1) = scavrate_hno3*1.08   ! so2
        kxf_gr(2) = scavrate_hno3*1.38   ! h2o2
        kxf_gr(3) = scavrate_hno3*1.59   ! nh3
        kxf_gr(4) = scavrate_hno3*0.80   ! h2so4
	if (lun143 > 0) write(lun143,'(a,1p,8e11.3)') 'kxf_gr,cr   ', &
		kxf_gr(1:4), kxf_cr


!   aa, bb, and cc coefficients of the 2 final equations
	tmpa = acen_tmp(k,icc,jcls,ipp)*rhocen_bar(k)
!	cc = 1/(dz*acen*rho)
	gasscav_cc(k,icc,jcls,ipp) = 1.0_r8/(dzcen(k)*tmpa)

	do ll = 1, ngas_scav
	    frac_c = heffcx(ll)/(1.0_r8 + heffcx(ll))
	    frac_g = 1.0_r8 - frac_c
	    kxf_gcr = frac_g*kxf_gr(ll) + frac_c*kxf_cr
!	    aa = kxf_gcr*dt
	    gasscav_aa(k,icc,jcls,ipp,ll) = kxf_gcr*dt_scav
                                                                
!	    bb = kxf_gr/(heffrx*vfallr*acen*rho)
	    gasscav_bb(k,icc,jcls,ipp,ll) = kxf_gr(ll)/(heffrx(ll)*vfallr*tmpa)
!	    setting gasscav_bb=0 (heffrx = infinity) gives irreversible scavenging
!	    gasscav_bb(k,icc,jcls,ipp,ll) = 0.0   

	    if (lun143 > 0) write(lun143,'(a,i1,1p,8e11.3)') 'aa/bb/cc   ', &
		ll, gasscav_aa(k,icc,jcls,ipp,ll), gasscav_bb(k,icc,jcls,ipp,ll), &
		gasscav_cc(k,icc,jcls,ipp), frac_g, frac_c, kxf_gcr
	end do ! l



	end do ! icc
	end do ! jcls
	end do ! ipp

	end do main_kloop_aa


	return
	end subroutine wetscav_2_gasscav



!-----------------------------------------------------------------------
	subroutine wetscav_2_bcscavcoef(                           &
		ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
		idiagaa_ecpp, ldiagaa_ecpp, idiagbb_wetscav,       &
		tcen_bar, pcen_bar, rhocen_bar,                    &
		it,      jt,      kts,ktebnd,ktecen,               &
		ncls_use,                                          &
		rh_sub2,                                           &
		is_active, is_precp,                               &
		chem_tmpa, scavcoef_num, scavcoef_vol              )


!-----------------------------------------------------------------------
! DESCRIPTION
!
! wetscav_2_bcscavcoef calculates below-cloud scavenging coefficents
!    similar to subr modal_aero_bcscavcoef_get
!
!-----------------------------------------------------------------------

!	use module_state_description, only:  p_qv, p_qc

!	use module_data_radm2, only:  epsilc

!	use module_data_mosaic_asect, only:  ai_phase, cw_phase,   &
!		massptr_aer, maxd_asize, maxd_atype,   &
!		ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer,   &
!		waterptr_aer
	use module_data_mosaic_asect, only:  &
		ai_phase, dens_aer, hygro_aer,   &
		massptr_aer, maxd_asize, maxd_atype,   &
		ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer,   &
		dcen_sect, dhi_sect, dlo_sect, sigmag_aer,   &
		volumhi_sect, volumlo_sect

	use modal_aero_wateruptake, only:  modal_aero_kohler

	 use aero_model, only:  &
		calc_1_impact_rate,   &
		dlndg_nimptblgrow, nimptblgrow_mind, nimptblgrow_maxd,   &
		scavimptblnum, scavimptblvol

	 use module_data_ecpp1
!==Guangxing Lin

!	use module_ecpp_hoststuff, only:  config_flags_ecpp

!	use module_mosaic_wetscav, only:  wetscav_cbmz_mosaic

!	use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message,   &
!	                             parampollu_1clm_set_opts

	implicit none

!   arguments
!   ( for definitions see subr parampollu_tdx_wetscav_2 )
	integer, intent(in) ::                  &
		ktau, ktau_pp, itstep_sub,   &
		it, jt, kts, ktebnd, ktecen

	integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199), &
		idiagbb_wetscav

	real(r8), intent(in) :: dtstep, dtstep_sub

	real(r8), intent(in), dimension( kts:ktecen ) ::   &
		tcen_bar, pcen_bar, rhocen_bar

	integer, intent(in) :: ncls_use

	real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		rh_sub2

	logical, intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		is_active, is_precp

	real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:num_chem_ecpp ) ::   &
		chem_tmpa

	real(r8), intent(inout),   &
		dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:maxd_asize, 1:maxd_atype ) ::   &
		scavcoef_num, scavcoef_vol
!	scavcoef_vol = below-cloud scavenging coeficient for volume (1/mm)
!	scavcoef_num = below-cloud scavenging coeficient for number (1/mm)
!       when precip rate = xxx kg/m2/s == xxx mm/s, the scavenging rate (1/s) = scavcoef*xxx 


!   local variables
	integer :: icc, ipp
	integer :: jcls, jgrow
	integer :: k
	integer :: l, ll, lun142
	integer :: m
	integer :: n
	integer :: p1st

	real(r8) :: dgratio
	real(r8) :: dry_dens, dry_diam, dry_mass, dry_volu
	real(r8) :: dry_mass_cut, dry_volu_cut
	real(r8) :: fact_leng, fact_mass
	real(r8), parameter :: onethird = 1.0_r8/3.0_r8
	real(r8), parameter :: piover6 = 3.14159265358979323846_r8/6.0_r8
	real(r8) :: scavimpnum, scavimpvol
	real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf, tmpg
	real(r8) :: tmpflo, tmpfhi
	real(r8) :: tmp_hygro, tmp_num, tmp_rdry, tmp_rwet, tmp_rh
	real(r8) :: watr_mass, wet_dens, wet_diam, wet_volu
	real(r8) :: xgrow
        real(r8) :: rdry_in_mak(1), hygro_mak(1), s_mak(1), rwet_out_mak(1)

! NOTE ON UNITS
!
! hostcode		wrfchem		cam
! mass mixing ratios	ug/kg		kg/kg
! dry/wet_mass		g/kgair		kg/kgair
! dens_aer		g/cm^3		kg/m^3
! dry/wet_volu		cm^3/kgair	m^3/kgair
! volumlo/hi_sect	cm^3		m^3
! dcen_sect		cm		m
! dry/wet_diam		cm		m
!
	    if ( hostcode_is_wrfchem ) then
		fact_mass = 1.0e-6_r8   ! ug/kgair --> g/kgair
		fact_leng = 1.0e-2_r8   ! cm --> m
		dry_mass_cut = 1.0e-26_r8   ! g/kgair = 1.0e-20 ug/kgair
		dry_volu_cut = 1.0e-26_r8   ! cm^3/kgair
	    else
		fact_mass = 1.0_r8      ! kg/kgair, unchanged
		fact_leng = 1.0_r8      ! m, unchanged
		dry_mass_cut = 1.0e-29_r8   ! kg/kgair = 1.0e-20 ug/kgair
		dry_volu_cut = 1.0e-32_r8   ! m^3/kgair
	    end if

!
! calc below-cloud scavenging coefficients of interstitial aerosols
!
	scavcoef_num(:,:,:,:,:,:) = 0.0_r8
	scavcoef_vol(:,:,:,:,:,:) = 0.0_r8

	do k = kts, ktecen
	do jcls = 1, ncls_use
	do ipp = 1, 2
icc_loop: &
	do icc = 1, 2
	    if ( .not. is_active(k,icc,jcls,ipp) ) cycle
!	    if ( .not. is_precp( k,icc,jcls,ipp) ) cycle

	    lun142 = 0
!	    if ((ktau == 1) .and. (k == 5)) lun142 = 142
	    if (k == 5) lun142 = 142
	    if (idiagbb_wetscav <= 0) lun142 = -1


! calc below-cloud scavenging coefficients for each aerosol mode
	    do n = 1, ntype_aer
	    do m = 1, nsize_aer(n)

!               calc dry mass and dry volume mixing ratios
		dry_volu = 0.0
		dry_mass = 0.0
		tmp_hygro = 0.0
		do l = 1, ncomp_aer(n)
		    tmpa = chem_tmpa(k,icc,jcls,ipp,massptr_aer(l,m,n,ai_phase))
		    dry_mass = dry_mass + tmpa
		    dry_volu = dry_volu + tmpa/dens_aer(l,n)
		    tmp_hygro = tmp_hygro + (tmpa/dens_aer(l,n))*hygro_aer(l,n)
		end do
		dry_mass = dry_mass*fact_mass   ! g/kgair    OR kg/kgair
		dry_volu = dry_volu*fact_mass   ! cm^3/kgair OR m^3/kgair

!		if negligible aerosol is present at this size and type, cycle
		if ((dry_mass < dry_mass_cut) .or. (dry_volu < dry_volu_cut)) then
		    ! BUT FIRST set dgn_dry/wet and chem_sub( ... water ... ) to default values
		    cycle
		end if

!               calc volume-mean dry diameter
		tmp_num = chem_tmpa(k,icc,jcls,ipp,numptr_aer(m,n,ai_phase))
		if (dry_volu <= tmp_num*volumlo_sect(m,n)) then
		    dry_diam = dlo_sect(m,n)
		else if (dry_volu >= tmp_num*volumhi_sect(m,n)) then
		    dry_diam = dhi_sect(m,n)
		else
		    dry_diam = (dry_volu/(tmp_num*piover6))**onethird
		end if

!               calc volume-mean wet diameter
		tmp_hygro = tmp_hygro*fact_mass/dry_volu
		tmp_rh = max( 0.0_r8, min( 0.99_r8, rh_sub2(k,icc,jcls,ipp) ) )
		tmp_rdry = dry_diam*0.5_r8*fact_leng   ! cm OR m --> m
		tmp_rwet = tmp_rdry
 
                rdry_in_mak(1) = tmp_rdry
                hygro_mak(1) = tmp_hygro
                s_mak(1) = tmp_rh
                rwet_out_mak(1) = tmp_rwet
!             call modal_aero_kohler( tmp_rdry, tmp_hygro, tmp_rh, tmp_rwet, 1, 1 )
!                call modal_aero_kohler( rdry_in_mak, hygro_mak, s_mak, rwet_out_mak, 1, 1)
!==Guangxing Lin            
                call modal_aero_kohler( rdry_in_mak, hygro_mak, s_mak, rwet_out_mak, 1)
!==Guangxing Lin            
                tmp_rwet = rwet_out_mak(1)

		wet_diam = tmp_rwet*2.0_r8/fact_leng   ! m --> cm OR m
		wet_diam = min( wet_diam, dry_diam*100.0_r8, 50.0e-6_r8/fact_leng )
		wet_diam = max( wet_diam, dry_diam )

!		wet_diam = dry_diam   ! force water == 0 (for testing)

		wet_volu = dry_volu * (wet_diam/dry_diam)**3   ! cm^3/kgair
		watr_mass = max( 0.0_r8, (wet_volu-dry_volu) )   ! g/kgair, as rho_water = 1.0 g/cm^3
! *** eventually should store this in some array that can be used by cam3
!     for now, leave it alone
!		chem_tmpa(k,icc,jcls,ipp,waterptr_aer(m,n)) = watr_mass/fact_mass

		wet_dens = (dry_mass + watr_mass)/wet_volu
		dry_dens = dry_mass/dry_volu

!   compute impaction scavenging removal amount for volume
!		interpolate table values using log of (actual-wet-size)/(base-dry-size)

!   in the bcscavcoef_get routine, dgratio = dgnum_wet/dgnum_amode
!   BUT dgnum_wet/dgnum_amode = (b*dgnum_wet)/(b*dgnum_amode) = dvolmean_wet/dcen_sect
!       where b = exp( 1.5 * (log(sigmag)**2) )
!		dgratio = ((wet_volu/dry_volu)**onethird) * (dry_diam/dcen_sect(m,n))
		dgratio = wet_diam/dcen_sect(m,n)
                                                                                                                                            
		if ((dgratio .ge. 0.99) .and. (dgratio .le. 1.01)) then
		    scavimpvol = scavimptblvol(0,m)
		    scavimpnum = scavimptblnum(0,m)
		else
		    xgrow = log( dgratio ) / dlndg_nimptblgrow
		    jgrow = int( xgrow )
		    if (xgrow .lt. 0.) jgrow = jgrow - 1
		    if (jgrow .lt. nimptblgrow_mind) then
			jgrow = nimptblgrow_mind
			xgrow = jgrow
		    else
			jgrow = min( jgrow, nimptblgrow_maxd-1 )
		    end if
                                                                                                                                            
		    tmpfhi = xgrow - jgrow
		    tmpfhi = max( 0.0_r8, min( 1.0_r8, tmpfhi ) )
		    tmpflo = 1.0_r8 - tmpfhi
		    scavimpvol = tmpflo*scavimptblvol(jgrow,m) +   &
		                 tmpfhi*scavimptblvol(jgrow+1,m)
		    scavimpnum = tmpflo*scavimptblnum(jgrow,m) +   &
		                 tmpfhi*scavimptblnum(jgrow+1,m)
		end if
                                                                                                                                            
		!impaction scavenging removal amount for volume
		scavcoef_vol(k,icc,jcls,ipp,m,n) = exp( scavimpvol )
		!impaction scavenging removal amount to number
		scavcoef_num(k,icc,jcls,ipp,m,n) = exp( scavimpnum )
                                                                                                                                            
! test diagnostics
		if (lun142 > 0) then
		    write(lun142,'(/a,8i4)') 'wetscav_2_bcscavcoef diags', &
			ktau, k, jcls, ipp, icc, n, m
		    tmpb = sigmag_aer(m,n)
		    tmpg = log( sigmag_aer(m,n) )
		    tmpg = exp( 1.5_r8*tmpg*tmpg )
		    tmpa = dcen_sect(m,n)*dgratio/tmpg
		    tmpc = dens_aer(1,n) ! bcscavcoef_init uses this
		    if ( .not. hostcode_is_wrfchem ) then
			tmpa = tmpa*1.0e2_r8    ! m --> cm
			tmpc = tmpc*1.0e-3_r8   ! kg/m^3 --> g/cm^3
		    end if
		    tmpd = 273.16_r8     ! bcscavcoef_init uses this
		    tmpe = 0.75e6_r8     ! bcscavcoef_init uses this
!		    call calc_1_impact_rate(             &
!			dg0, sigmag, rhoaero, temp, press, &
!			scavratenum, scavratevol, lunerr )
		    call calc_1_impact_rate(             &
			tmpa, tmpb, tmpc, tmpd, tmpe, &
			tmpf, tmpg, lun142 )
		    write(lun142,'(1p,8e11.3)') dgratio, &
			tmpa, tmpb, tmpc, tmpd, tmpe
		    write(lun142,'(1p,8e11.3)') &
			scavcoef_num(k,icc,jcls,ipp,m,n), tmpf, &
			scavcoef_vol(k,icc,jcls,ipp,m,n), tmpg 
		    write(lun142,'(1p,8e11.3)') &
			dry_mass, dry_volu, wet_volu, dry_diam, wet_diam, tmp_rh, &
			chem_tmpa(k,icc,jcls,ipp,numptr_aer(m,n,ai_phase))
		end if

	    end do ! m
	    end do ! n

	end do icc_loop ! icc
	end do ! ipp
	end do ! jcls
	end do ! k


!	write(*,'(a)') 'wetscav_2_bcscavcoef DONE'
	return
	end subroutine wetscav_2_bcscavcoef



!-----------------------------------------------------------------------
	subroutine wetscav_2_precip_evap_xfer(                     &
		ktau, dtstep, ktau_pp, itstep_sub, dtstep_sub, &
		idiagaa_ecpp, ldiagaa_ecpp, idiagbb_wetscav,       &
		it,      jt,      kts,ktebnd,ktecen,               &
		ncls_use,                                          &
		is_active, is_precp, is_ptgain, is_ptloss,         &
		acen_tmp, prtb, frac_evap_prtb, frac_xfer_prtb,    &
		fxaa_evap_prtb                                     )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! wetscav_2_precip_evap_xfer calculates the fractions of precip 
! (and precip-borne aerosols) entering the top of a subarea that are either 
!	> evaporated/resuspended or 
!	> transferred to another subarea
!
!-----------------------------------------------------------------------

!	use module_state_description, only:  p_qv, p_qc

!	use module_data_radm2, only:  epsilc

!	use module_data_mosaic_asect, only:  ai_phase, cw_phase,   &
!		massptr_aer, maxd_asize, maxd_atype,   &
!		ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer,   &
!		waterptr_aer

	use module_data_ecpp1

	implicit none

!   subr arguments
	integer, intent(in) ::                  &
		ktau, ktau_pp, itstep_sub,   &
		it, jt, kts, ktebnd, ktecen
	integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199), &
		idiagbb_wetscav
	integer, intent(in) :: ncls_use

	real(r8), intent(in) :: dtstep, dtstep_sub

	logical,  intent(in),  dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		is_active, is_precp
	logical,  intent(out), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) ::   &
		is_ptgain, is_ptloss

	real(r8), intent(in),  dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) :: acen_tmp
	real(r8), intent(in),  dimension( kts:ktebnd, 1:2, 1:maxcls_ecpp, 1:2 ) :: prtb

	real(r8), intent(out), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) :: frac_evap_prtb
	real(r8), intent(out), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, &
	                                              1:2, 1:maxcls_ecpp, 1:2 ) :: frac_xfer_prtb
	real(r8), intent(out), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, &
	                                              1:2, 1:maxcls_ecpp, 1:2 ) :: fxaa_evap_prtb
! frac_evap_prtb = fraction of precip (and precip-borne aerosols) entering the 
!                  top of a subarea that is evaporated/resuspended
! frac_xfer_prtb = fraction of precip (and precip-borne aerosols) entering the 
!                  top of a subarea that is transferred to another subarea
!                  (the first  set of icc,jcls,ipp indices are the "xfer from" subarea)
!                  (the second set of icc,jcls,ipp indices are the "xfer  to " subarea)

!   local variables
	integer :: icc, icc_g, icc_l, iphase, ipp, ipp_l, ipp_g
	integer :: jcls, jcls_g, jcls_l
	integer :: k, km1
	integer :: lun141
	integer :: m

	real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf, tmpg, tmph
	real(r8) :: tmpvecb(100), tmpvece(100), tmpvecf(100)

	real(r8), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2 ) :: &
		delprtb, delprtb_g, delprtb_l
	real(r8), dimension( kts:ktecen ) :: &
		delprtb_gtot, delprtb_ltot, delprtb_xtot, &
		frac_evap, frac_xferg, frac_xferl


	lun141 = -1
	if (idiagaa_ecpp(141) > 0) lun141 = ldiagaa_ecpp(141)
	if (idiagbb_wetscav <= 0) lun141 = -1

	is_ptloss(kts:ktecen,1:2,1:ncls_use,1:2) = .false.
	is_ptgain(kts:ktecen,1:2,1:ncls_use,1:2) = .false.

	frac_evap_prtb(kts:ktecen,1:2,1:ncls_use,1:2) = 0.0
	frac_xfer_prtb(kts:ktecen,1:2,1:ncls_use,1:2, 1:2,1:ncls_use,1:2) = 0.0
	fxaa_evap_prtb(kts:ktecen,1:2,1:ncls_use,1:2, 1:2,1:ncls_use,1:2) = 0.0

	delprtb_gtot(:) = 0.0 ; delprtb_ltot(:) = 0.0 ; delprtb_xtot(:) = 0.0
	frac_evap(:) = 0.0 ; frac_xferg(:) = 0.0 ; frac_xferl(:) = 0.0

main_kloop_aa: &
	do k = ktecen, kts, -1

!
! calculate the fractions of precip (and precip-borne aerosols) 
!    entering the top of a subarea that are either 
!	> evaporated/resuspended or 
!	> transferred to another subarea
!
! this is a bit tricky because we do not have evaporation information,
!    and a decrease in precip from k+1 to k for one subarea 
!    can be due to that precip being classified as in another subarea
!
! approach here is to calculate precip loss and gains (from k+1 to k)
!    for each subarea, then try to balance them out
! any "unbalanced" loss is treated as true evaporation
!

	do ipp = 1, 2
	do jcls = 1, ncls_use
	do icc = 1, 2
! delprtb   = change in subarea precip from k+1 to k
! delprtb_g = gain   in subarea precip from k+1 to k
! delprtb_l = loss   in subarea precip from k+1 to k, but sign is positive
	    delprtb(k,icc,jcls,ipp) = prtb(k,  icc,jcls,ipp) &
	                            - prtb(k+1,icc,jcls,ipp)
	    delprtb_g(k,icc,jcls,ipp) = max( 0.0_r8,  delprtb(k,icc,jcls,ipp) )
	    delprtb_l(k,icc,jcls,ipp) = max( 0.0_r8, -delprtb(k,icc,jcls,ipp) )
	    if (delprtb_g(k,icc,jcls,ipp) > 0.0_r8) is_ptgain(k,icc,jcls,ipp) = .true.
	    if (delprtb_l(k,icc,jcls,ipp) > 0.0_r8) is_ptloss(k,icc,jcls,ipp) = .true.
	end do
	end do
	end do

! delprtb_gtot = sum of delprtb_g over all subareas ; similar for depltrb_ltot
	delprtb_gtot(k) = sum( delprtb_g(k,1:2,1:ncls_use,1:2) )
	delprtb_ltot(k) = sum( delprtb_l(k,1:2,1:ncls_use,1:2) )
! delprtb_xtot = is amount of precip loss that can be balance by precip gain
	delprtb_xtot(k) = min( delprtb_gtot(k), delprtb_ltot(k) )

	if (delprtb_gtot(k) > 0.0_r8) then
	    frac_xferg(k) = delprtb_xtot(k) / delprtb_gtot(k)
	    frac_xferg(k) = max( 0.0_r8, min( 1.0_r8, frac_xferg(k) ) )
	end if

	if (delprtb_ltot(k) <= 0.0_r8) cycle main_kloop_aa   ! bypass next steps if no loss

	frac_xferl(k) = delprtb_xtot(k) / delprtb_ltot(k)
	frac_xferl(k) = max( 0.0_r8, min( 1.0_r8, frac_xferl(k) ) )
	frac_evap(k) = 1.0_r8 - frac_xferl(k)


! do calcs associated with balancing of precip loss and gain
! current approach is that there is no preferred pairing of 
!    "losing" and "gaining" subareas
! one might want to pair the clear and cloud subareas of a
!    transport class first --  something to think about in the future
! *** this code is incomplete ***
!
! loop over the "losing" subareas
	do jcls_l = 1, ncls_use
	do ipp_l = 1, 2
	do icc_l = 1, 2
	    if ( .not. is_ptloss(k,icc_l,jcls_l,ipp_l) ) cycle
	    tmpa = delprtb_l(k,icc_l,jcls_l,ipp_l)/prtb(k+1,icc_l,jcls_l,ipp_l)
	    frac_evap_prtb(k,icc_l,jcls_l,ipp_l) = frac_evap(k)*tmpa

! loop over the "gaining" subareas
	    if (frac_xferl(k) <= 1.0e-7_r8) cycle
	    do jcls_g = 1, ncls_use
	    do ipp_g = 1, 2
	    do icc_g = 1, 2
		if ( .not. is_ptgain(k,icc_g,jcls_g,ipp_g) ) cycle
		tmpb = delprtb_g(k,icc_g,jcls_g,ipp_g)/delprtb_gtot(k)
	        frac_xfer_prtb(k,icc_l,jcls_l,ipp_l, icc_g,jcls_g,ipp_g) = &
		    frac_xferl(k)*tmpa*tmpb
	    end do ! icc_g
	    end do ! ipp_g
	    end do ! jcls_g

! if a subarea exists ( is_active ) and has precip>0 at k+1,
!    but does not exist at k, then the evaporated/resuspended material
!    from the losing subarea must go to other subareas
! the fxaa_evap_prtb are used for this
	    if ( .not. is_active(k,icc_l,jcls_l,ipp_l) ) then
		tmpf = 0.0_r8
		do jcls_g = 1, ncls_use
		do ipp_g = 1, 2
		do icc_g = 1, 2
		    if ( .not. is_active(k,icc_g,jcls_g,ipp_g) ) cycle
		    if ((jcls_g == jcls_l) .and. &
		        (ipp_g == ipp_l) .and. (icc_g == icc_l)) cycle
		    tmpf = tmpf + acen_tmp(k,icc_g,jcls_g,ipp_g)
                    fxaa_evap_prtb(k,icc_l,jcls_l,ipp_l, icc_g,jcls_g,ipp_g) = 1.0_r8
		end do ! icc_g
		end do ! ipp_g
		end do ! jcls_g
                fxaa_evap_prtb(k,icc_l,jcls_l,ipp_l, 1:2,1:ncls_use,1:2) = &
                fxaa_evap_prtb(k,icc_l,jcls_l,ipp_l, 1:2,1:ncls_use,1:2)*tmpf
	    end if

	end do ! icc_l
	end do ! ipp_l
	end do ! jcls_l


	end do main_kloop_aa


!
! diagnostics for testing
!
! first set shows main arrays that can be inspected visually 
	if (lun141 > 0) then

	tmph = 3600.0
	do k = ktecen, kts, -1

	write(lun141,'(a,i3)') 'k =', k

	tmpa = delprtb_ltot(k) - delprtb_xtot(k)
	tmpb = sum( frac_evap_prtb(k,1:2,1:ncls_use,1:2)*prtb(k+1,1:2,1:ncls_use,1:2) )
	write(lun141,'(a,2f9.5,2x,3f9.5,2x,2f9.5)') 'frac_xferl/evap, delg/l/xtot=', frac_xferl(k), frac_evap(k), &
		3600.0*delprtb_gtot(k), 3600.0*delprtb_ltot(k), 3600.0*delprtb_xtot(k), &
		3600.0*tmpa, 3600.0*tmpb

	write(lun141,'(a,3(2x,4f9.5))') 'acen       =', (((acen_tmp(k,icc,jcls,ipp), icc=1,2), ipp=1,2), jcls=1,ncls_use)
	write(lun141,'(a,3(2x,4f9.5))') 'prtb       =', (((prtb(k,icc,jcls,ipp)*tmph, icc=1,2), ipp=1,2), jcls=1,ncls_use)
	write(lun141,'(a,3(2x,4f9.5))') 'delprtb    =', (((delprtb(k,icc,jcls,ipp)*tmph, icc=1,2), ipp=1,2), jcls=1,ncls_use)
	write(lun141,'(a,3(2x,4f9.5))') 'delprtb_g  =', (((delprtb_g(k,icc,jcls,ipp)*tmph, icc=1,2), ipp=1,2), jcls=1,ncls_use)
	write(lun141,'(a,3(2x,4f9.5))') 'delprtb_l  =', (((delprtb_l(k,icc,jcls,ipp)*tmph, icc=1,2), ipp=1,2), jcls=1,ncls_use)
	icc_l = 2 ; ipp_l = 2 ; icc_g = 2 ; ipp_g = 2
	write(lun141,'(a,3(2x,4f9.5))') 'frac_ev/xf =', ( frac_evap_prtb(k,icc_l,jcls_l,ipp_l), &
		( frac_xfer_prtb(k,icc_l,jcls_l,ipp_l, icc_g,jcls_g,ipp_g), jcls_g=1,ncls_use), jcls_l=1,ncls_use) 

	end do


! second set does "conservation checks" 
! is prtb(k) equal to [prtb(k+1) + gains - losses] ?
	do k = ktecen, kts, -1

	write(lun141,'(a,i3)') 'k =', k

! here check sum( prtb ) over all subareas
	tmpa = sum( prtb(k+1,1:2,1:ncls_use,1:2) )
	tmpb = sum( prtb(k  ,1:2,1:ncls_use,1:2) )
	tmpc = tmpa + delprtb_gtot(k) - delprtb_ltot(k)
	tmpd = tmpa + delprtb_gtot(k)*(1.0_r8 - frac_xferg(k)) - delprtb_ltot(k)*(1.0_r8 - frac_xferl(k))
	tmpe = (tmpb-tmpc)*tmph                       ! absolute error in mm/h
	tmpf = (tmpb-tmpd)*tmph
	tmpe = (tmpb-tmpc)/max(tmpa,tmpb,1.0e-30_r8)  ! relative error
	tmpf = (tmpb-tmpd)/max(tmpa,tmpb,1.0e-30_r8)
	write(lun141,'(a,1p,2e10.2)') 'relerr1/2 =', tmpe, tmpf

! here check prtb for each subarea
	m = 0
	do jcls = 1, ncls_use
	do ipp = 1, 2
	do icc = 1, 2
	    tmpa = prtb(k+1,icc,jcls,ipp)
	    tmpb = prtb(k  ,icc,jcls,ipp)
	    tmpc = tmpa + delprtb_g(k,icc,jcls,ipp) - delprtb_l(k,icc,jcls,ipp)
	    if ( is_ptgain(k,icc,jcls,ipp) ) then
		tmpd = tmpa + delprtb_g(k,icc,jcls,ipp)*(1.0_r8 - frac_xferg(k)) &
		+ sum( frac_xfer_prtb(k,1:2,1:ncls_use,1:2,icc,jcls,ipp)*prtb(k+1,1:2,1:ncls_use,1:2) )
	    else if ( is_ptloss(k,icc,jcls,ipp) ) then
		tmpd = tmpa - prtb(k+1,icc,jcls,ipp)*( frac_evap_prtb(k,icc,jcls,ipp) &
			+ sum( frac_xfer_prtb(k,icc,jcls,ipp,1:2,1:ncls_use,1:2) ) )
	    else
		tmpd = tmpb
	    end if
	    tmpe = (tmpb-tmpc)*tmph                       ! absolute error in mm/h
	    tmpf = (tmpb-tmpd)*tmph
	    tmpe = (tmpb-tmpc)/max(tmpa,tmpb,1.0e-30_r8)  ! relative error
	    tmpf = (tmpb-tmpd)/max(tmpa,tmpb,1.0e-30_r8)
	    m = m + 1
	    tmpvece(m) = tmpe
	    tmpvecf(m) = tmpf
	    tmpvecb(m) = tmpb*tmph
	end do
	end do
	end do
	write(lun141,'(a,1p,3(2x,4e10.2))') 'tmpvecb   =', tmpvecb(1:m)
	write(lun141,'(a,1p,3(2x,4e10.2))') 'tmpvece   =', tmpvece(1:m)
	write(lun141,'(a,1p,3(2x,4e10.2))') 'tmpvecf   =', tmpvecf(1:m)
	 
	end do ! k = ktecen, kts, -1

	end if ! (lun141 > 0)


	end subroutine wetscav_2_precip_evap_xfer


end module ecpp_modal_wetscav

