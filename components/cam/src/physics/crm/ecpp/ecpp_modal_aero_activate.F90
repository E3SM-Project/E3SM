module ecpp_modal_aero_activate

!-----------------------------------------------------------------
!  Module interface of aerosol activaiton used in the ECPP treatment 
! in the MMF model
!  Adopted from ndrop.F90 and from the similar one used in the ECPP 
!  for the WRF-chem model written by Dick Easter
!
!  Minghuai Wang, 2009-11
!------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
!==Guangxing Lin    
   !use abortutils,   only: endrun
   use cam_abortutils,   only: endrun
!==Guangxing Lin    
 
   implicit none

   public parampollu_tdx_activate1
   public parampollu_tdx_activate_intface

contains

!-----------------------------------------------------------------------
	subroutine parampollu_tdx_activate1(                      &
		ktau, dtstep, ktau_pp, dtstep_pp,             &
		idiagaa_ecpp, ldiagaa_ecpp,                       &
		tcen_bar, pcen_bar, rhocen_bar, dzcen,            &
		rhobnd_bar, wbnd_bar,                             &
		ncls_ecpp,                                        &
		it,      jt,      kts,ktebnd,ktecen,              &
		ncls_use, ifrom_where, activate_onoff_use,        &
		kdraft_bot_use, kdraft_top_use,                   &
		mtype_updnenv_use,                                &
		chem_sub_old,                                     &
		mfbnd_use,                                        &
		ar_bnd_tavg,                                      &
		ent_airamt,                                       &
		ido_actres_horz, fmact_horz, fnact_horz,          &
		fmact_vert, fnact_vert, mfbnd_quiescn_up          )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_activate1 calculates number and mass activation
!    fractions associated with vertical and horizontal transfer
!    between subclasses
!
!-----------------------------------------------------------------------

	use module_data_mosaic_asect, only:  maxd_asize, maxd_atype,   &
		nsize_aer, ntype_aer

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
!	chem_driver and routines under it do calculations
!	over these spatial indices.

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
!
	integer, intent(in) :: ncls_ecpp
!   ncls_ecpp - number of ecpp transport classes in the grid column
	integer, intent(in)    :: ifrom_where
!	1,2 - from area_change; 10 - from main_integ
	integer, intent(in)    :: activate_onoff_use
!	 1-99 - calc real fmact,fnact
!	  200 - set fmact = fmact_testa, ...
!	other - set fmact,fnact = 0.0
!       ALSO, ido_actres_horz is set correctly when activate_onoff_use > 0
!          but is set to zero when activate_onoff_use <= 0

	integer, intent(in)    :: ncls_use

	integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
		kdraft_bot_use, kdraft_top_use,   &
		mtype_updnenv_use

	real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
		chem_sub_old

	real(r8), intent(inout), dimension( kts:ktebnd, 0:2, 0:maxcls_ecpp ) :: &
		ar_bnd_tavg, mfbnd_use

	real(r8), intent(in), dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
		ent_airamt

	integer, intent(out), dimension( 1:2, 1:maxcls_ecpp, 1:2, 1:maxcls_ecpp ) ::   &
		ido_actres_horz
!   ido_actres_horz(iccaa,jclsaa,iccbb,jclsbb) is associated with air moving 
!       into sub-class (iccaa,jclsaa) from sub-class (iccbb,jclsbb)
!   ido_actres_horz = +1 or +2 if activation, -1 if resuspension, 0 otherwise
!   note that its values are independent of k (i.e., they only depend on the source and 
!       destination sub-classes)
!   the fnact and fmact do depend on k

	real(r8), intent(out), dimension( 1:maxd_asize, 1:maxd_atype, 1:maxcls_ecpp,   &
		                          1:2, 1:maxcls_ecpp, kts:ktecen ) ::   &
		fmact_horz, fnact_horz
!   fmact_horz(m,n,jclsaa,iccbb,jclsbb,k) and fnact(...) are associated with air moving 
!       into sub-class (icc=2,jclsaa,k) from sub-class (iccbb,jclsbb,k)

	real(r8), optional, intent(out), dimension( 1:maxd_asize, 1:maxd_atype, kts:ktecen ) ::   &
		fmact_vert, fnact_vert
!   fnact_vert(m,n,k) and fmact(...) are associated with (quiescent, clear, layer k-1) air moving 
!      into (quiescent, cloudy, layer k)

	real(r8), optional, intent(in), dimension( kts:ktebnd, 0:2, 0:2 ) ::   &
		mfbnd_quiescn_up


!   local variables
	integer :: icc, iccb, iccy, ido_actres_tmp, ihorzvert, itmpa
	integer :: jcls, jclsy, jj
	integer :: k, l
	integer :: m, n

	real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpt
	real(r8) :: wbar_tmp, wmix_tmp

	real(r8), dimension( 1:maxd_asize, 1:maxd_atype ) ::   &
		fnact_tmp, fmact_tmp
	real(r8), dimension( 1:maxd_asize, 1:maxd_atype, 2 ) ::   &
		fnact_testa, fmact_testa


!   initialize fnact/fmact to zero
	ido_actres_horz(:,:,:,:) = 0
	fmact_horz(:,:,:,:,:,:) = 0.0
	fnact_horz(:,:,:,:,:,:) = 0.0
	if ( present(fmact_vert) ) fmact_vert(:,:,:) = 0.0
	if ( present(fnact_vert) ) fnact_vert(:,:,:) = 0.0

	if (activate_onoff_use <= 0) return


!   temporary values for testing purposes
	fmact_testa(:,:,:) = 0.0
	fnact_testa(:,:,:) = 0.0

        fmact_testa(1,1:3,1) = (/ 0.50, 0.90, 0.95 /) ! updraft
        fnact_testa(1,1:3,1) = (/ 0.40, 0.80, 0.90 /)
        fmact_testa(1,1:3,2) = (/ 0.30, 0.80, 0.90 /) ! quiescent
        fnact_testa(1,1:3,2) = (/ 0.20, 0.60, 0.80 /)

!
!   horizontal transfer
!

!   first set ido_actres_horz
!   note again:  ido_actres_horz(icc,jcls,iccy,jclsy) is from iccy,jclsy to icc,jcls
	do jclsy = 1, ncls_use
	do iccy = 1, 2
	do jcls = 1, ncls_use
	do icc = 1, 2

	if (icc == 1) then
	    if (iccy == 1) then
	    ! clear --> clear -- do nothing (no activation or resuspension)
		cycle
	    else
	    ! cloudy --> clear -- do resuspension
		ido_actres_horz(icc,jcls,iccy,jclsy) = -1
	    end if

	else
	    if (iccy == 1) then
	    ! clear --> cloudy -- do activation for into updrafts & quiescent
	    !                     do nothing for into downdrafts
		if (mtype_updnenv_use(icc,jcls) /= mtype_dndraft_ecpp)   &
		    ido_actres_horz(icc,jcls,iccy,jclsy) = 1
	    else
	    ! cloudy --> cloudy -- do (re)activation for into updrafts
	    !                      do nothing for into downdrafts & quiescent
	    !	if (mtype_updnenv_use(icc,jcls) == mtype_updraft_ecpp)   &
	    !		    ido_actres_horz(icc,jcls,iccy,jclsy) = 2
	    end if
	end if

	end do   ! icc
	end do   ! jcls
	end do   ! iccy
	end do   ! jclsy



!   next calc activation fractions
horz_k_loop:   &
	do k = kts, ktecen

horz_jcls_loop:   &
	do jcls = 1, ncls_use
	icc = 2

horz_jclsy_loop:   &
	do jclsy = 1, ncls_use

horz_iccy_loop:   &
	do iccy = 1, 2

	if (ent_airamt(icc,jcls,iccy,jclsy,k) <= 0.0) cycle horz_iccy_loop

	if (jcls == jcls_qu) then
!   quiescent class 
!      it can entrain from quiescent, updraft, dndraft
!      do activation for entrain from clear-any
	    if (iccy == 2) cycle horz_iccy_loop   ! only activate clear --> cloudy
	
	else if (mtype_updnenv_use(icc,jcls) == mtype_dndraft_ecpp) then
!   downdraft class
!      it can entrain from quiescent, dndraft
!      do activation for none of these
	    cycle horz_iccy_loop

	else
!   updraft class
!      it can entrain from quiescent, updraft
!      do activation for entrain from any-quiescent and clear-updraft
	    if (jclsy == jcls_qu) then
		continue
	    else if ( (iccy == 1) .and.   &
	              (mtype_updnenv_use(iccy,jclsy) ==   &
	               mtype_updraft_ecpp) ) then
		continue
	    else
		cycle horz_iccy_loop
	    end if
	end if

	if (activate_onoff_use == 200) then   ! use the fmnact_tst values
	    jj = 1
	    if (jcls == jcls_qu) jj = 2
	    fmact_horz(:,:,jcls,iccy,jclsy,k) = fmact_testa(:,:,jj)
	    fnact_horz(:,:,jcls,iccy,jclsy,k) = fnact_testa(:,:,jj)
	end if

	if (activate_onoff_use < 100) then   ! calculate "real" values
!	    stop '*** parampollu_tdx_activate1 - cannot do activate_onoff_use < 100'

	    tmpa = 0.5*(mfbnd_use(k,icc,jcls)+mfbnd_use(k+1,icc,jcls))
	    tmpb = 0.5*(ar_bnd_tavg(k,icc,jcls)+ar_bnd_tavg(k+1,icc,jcls))
	    if (tmpb > 0.0) then
		if (abs(tmpa) > abs(tmpb)*w_draft_max) then
		    wbar_tmp = w_draft_max
		else
	            wbar_tmp = tmpa/tmpb
		end if
	    else
		wbar_tmp = 0.0
	    end if
	    wbar_tmp = wbar_tmp + 0.5*(wbnd_bar(k)+wbnd_bar(k+1))
	    wmix_tmp = 0.0
	    if (max(wbar_tmp,wmix_tmp) <= 0.0) cycle horz_iccy_loop

	    ido_actres_tmp = ido_actres_horz(icc,jcls,iccy,jclsy)
	    ihorzvert = 1

	    call parampollu_tdx_activate_intface(                 &
		ktau, ktau_pp,                                &
		idiagaa_ecpp, ldiagaa_ecpp,                       &
		ncls_ecpp, ncls_use,                              &
		it,      jt,      kts,ktebnd,ktecen,              &
		k, iccy, jclsy, jcls,                             &
		activate_onoff_use, ido_actres_tmp,               &
		ihorzvert, ifrom_where,                           &
		chem_sub_old,                                     &
		tcen_bar(k), rhocen_bar(k),                       &
		wbar_tmp, wmix_tmp,                               &
		fmact_testa, fnact_testa,                         &
		fmact_tmp, fnact_tmp                              )

	    fmact_horz(:,:,jcls,iccy,jclsy,k) = fmact_tmp(:,:)
	    fnact_horz(:,:,jcls,iccy,jclsy,k) = fnact_tmp(:,:)
	end if

	end do horz_iccy_loop
	end do horz_jclsy_loop
	end do horz_jcls_loop
	end do horz_k_loop

!	write(*,'(a,i4,1p,4e10.2)') 'tdx_activate1 horz min/max', ifrom_where,   &
!	    minval(fmact_horz(:,:,:,:,:,:)), maxval(fmact_horz(:,:,:,:,:,:)),   &
!	    minval(fnact_horz(:,:,:,:,:,:)), maxval(fnact_horz(:,:,:,:,:,:))


!
! vertical transfer 
!    in up/dndrafts, vertical transport is clear<-->clear or cloudy<-->cloudy
!       so no activation
!    in quiescent, can have clear<-->cloudy
!       do activation for clear(k-1)-->cloud(k)
!
	if ( present(fmact_vert) .and. present(fnact_vert) ) then

vert_k_loop:   &
	do k = kts, ktecen
	if (k == kts) cycle vert_k_loop

	jcls = jcls_qu
	icc = 2
	jclsy = jcls_qu
	iccy = 1

!	mfbnd_quiescn_up(k,iccy,icc) is upwards mass flux from iccy to icc 
!		at bottom of layer k
	if (mfbnd_quiescn_up(k,iccy,icc) <= 0.0) cycle vert_k_loop

	if (activate_onoff_use == 200) then   ! use the fmnact_tst values
	    jj = 2
	    fmact_vert(:,:,k) = fmact_testa(:,:,jj)
	    fnact_vert(:,:,k) = fnact_testa(:,:,jj)
	end if

	if (activate_onoff_use < 100) then   ! calculate "real" values
!	    stop '*** parampollu_tdx_activate1 - cannot do activate_onoff_use < 100'

	    tmpa = mfbnd_use(k,iccy,jclsy)
	    tmpb = ar_bnd_tavg(k,iccy,jclsy)
	    if (tmpb > 0.0) then
		if (abs(tmpa) > abs(tmpb)*w_draft_max) then
		    wbar_tmp = w_draft_max
		else
	            wbar_tmp = tmpa/tmpb
		end if
	    else
		wbar_tmp = 0.0
	    end if
	    wbar_tmp = wbar_tmp + wbnd_bar(k)
	    wmix_tmp = 0.0
	    if (max(wbar_tmp,wmix_tmp) <= 0.0) cycle vert_k_loop

	    ido_actres_tmp = 1

	    tmpt = 0.5*( tcen_bar(k) + tcen_bar(max(k-1,kts)) )

	    ido_actres_tmp = 1
	    ihorzvert = 2

	    call parampollu_tdx_activate_intface(                 &
		ktau, ktau_pp,                                &
		idiagaa_ecpp, ldiagaa_ecpp,                       &
		ncls_ecpp, ncls_use,                              &
		it,      jt,      kts,ktebnd,ktecen,              &
		k-1, iccy, jclsy, jcls,                           &
		activate_onoff_use, ido_actres_tmp,               &
		ihorzvert, ifrom_where,                           &
		chem_sub_old,                                     &
		tmpt, rhobnd_bar(k),                              &
		wbar_tmp, wmix_tmp,                               &
		fmact_testa, fnact_testa,                         &
		fmact_tmp, fnact_tmp                              )

	    fmact_vert(:,:,k) = fmact_tmp(:,:)
	    fnact_vert(:,:,k) = fnact_tmp(:,:)
	end if

	end do vert_k_loop

!	write(*,'(a,i4,1p,4e10.2)') 'tdx_activate1 vert min/max', ifrom_where,   &
!	    minval(fmact_vert(:,:,:)), maxval(fmact_vert(:,:,:)),   &
!	    minval(fnact_vert(:,:,:)), maxval(fnact_vert(:,:,:))

	end if   ! ( present(fmact_vert) .and. present(fnact_vert) )



	return
	end subroutine parampollu_tdx_activate1



!-----------------------------------------------------------------------
	subroutine parampollu_tdx_activate_intface(               &
		ktau, ktau_pp,                                &
		idiagaa_ecpp, ldiagaa_ecpp,                       &
		ncls_ecpp, ncls_use,                              &
		i,       j,       kts,ktebnd,ktecen,              &
		k, iccy, jclsy, jcls,                             &
		activate_onoff_use, ido_actres,                   &
		ihorzvert, ifrom_where,                           &
		chem_sub_old,                                     &
		tempair_in, rhoair_in,                            &
		wbar_in, wmix_in,                                 &
		fmact_testa, fnact_testa,                         &
		fmact, fnact                                      )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_activate1 calculates number and mass activation
!    fractions associated with vertical and horizontal transfer
!    between subclasses
!
!-----------------------------------------------------------------------

	use module_data_mosaic_asect, only:    &
		maxd_acomp, maxd_asize, maxd_atype,   &
		ncomp_aer, nsize_aer, ntype_aer,   &
		nphase_aer, ai_phase, cw_phase,   &
		numptr_aer, massptr_aer, sigmag_aer

	use module_data_ecpp1

	use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message

        use ndrop,           only: activate_modal
        use constituents,    only: pcnst

!   arguments
	integer, intent(in) ::                  &
		ktau, ktau_pp,              &
		i, j, kts, ktebnd, ktecen
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!	chem_driver and routines under it do calculations
!	over these spatial indices.

	integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

	integer, intent(in) ::                  &
		k, iccy, jclsy, jcls

	real(r8), intent(in) :: tempair_in, rhoair_in, wbar_in, wmix_in
!   tempair - temperature (k)
!   rhoair - air density (kg/m3)

	integer, intent(in) :: ncls_ecpp
!   ncls_ecpp - number of ecpp transport classes in the grid column
	integer, intent(in)    :: ncls_use

	integer, intent(in)    :: activate_onoff_use
!	 1-99 - calc real fmact,fnact
!	  200 - set fmact = fmact_testa, ...
!	other - set fmact,fnact = 0.0
	integer, intent(in)    :: ido_actres, ihorzvert, ifrom_where

	real(r8), intent(in), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
		chem_sub_old

	real(r8), intent(in), dimension( 1:maxd_asize, 1:maxd_atype, 2 ) ::   &
		fnact_testa, fmact_testa

	real(r8), intent(out), dimension( 1:maxd_asize, 1:maxd_atype ) ::   &
		fmact, fnact


!   local variables
	integer :: iphase, jj, l, ll, lun, m, n
	integer, save :: ifrom_where_save, ktau_save
	data ifrom_where_save, ktau_save / -1, -1 /

	real(r8) :: factscale, flux_fullact
	real(r8) :: rhoair
	real(r8) :: sumhygro, sumvol
	real(r8) :: tempair, tmpc
	real(r8) :: wbar, wdiab, wmin, wmax, wmix, wmixmin

	real(r8) :: raercol( 1:1, 1:num_chem_ecpp )

        real(r8) :: raer (1:pcnst)    ! interstitial aerosols
        real(r8) :: qqcw (1:pcnst)    ! interstitial aerosols

	real(r8), dimension( 1:maxd_asize, 1:maxd_atype ) ::   &
	    fn, fs, fm, fluxn, fluxs, fluxm, hygro, &
	    maerosol_tot, maerosol_totcw, &
	    naerosol, naerosolcw, &
	    vaerosol, vaerosolcw, sigmag

	real(r8), dimension( 1:maxd_acomp, 1:maxd_asize, 1:maxd_atype ) ::   &
	    maerosol, maerosolcw



!   initialize fnact/fmact to zero
	fmact(:,:) = 0.0
	fnact(:,:) = 0.0

!   special testing cases
	if ((activate_onoff_use <= 0) .or. (activate_onoff_use >= 100)) then
	    return
	else if (activate_onoff_use == 81) then
	    return
	else if (activate_onoff_use == 82) then
	    jj = 1
	    if (jcls == jcls_qu) jj = 2
	    fmact(:,:) = fmact_testa(:,:,jj)
	    fnact(:,:) = fnact_testa(:,:,jj)
	    return
	end if

!
!   calc activation fractions
!
	tempair = tempair_in
	rhoair = rhoair_in
	wbar = wbar_in
	wmix = wmix_in

	wmixmin = 0.2
	! do single updraft, forced to wbar >= wmixmin
	wbar = max( wbar+wmix, wmixmin )
	wmix = 0.0

	wmin = 0.0
	wmax = 50.0
	wdiab = 0.0

!   load raercol (with units conversion) and calculate hygro
	raercol(:,:) = 0.0
        
        raer(1:pcnst) = chem_sub_old(k,iccy,jclsy,1:pcnst)
        qqcw(1:pcnst) = chem_sub_old(k,iccy,jclsy,pcnst+1:2*pcnst)

!   do loadaer calls
	do n=1,ntype_aer
	do m=1,nsize_aer(n)
	
          if(ido_actres ==2 ) then
            iphase = 3
          else
            iphase = 1
          end if
         call loadaer0D (raer, qqcw, n, rhoair, ai_phase,  &
               naerosol(m,n), vaerosol(m,n), hygro(m,n))
          sigmag(m, n) = sigmag_aer(m,n)
	enddo ! m
	enddo ! n

!   do activate call
        m = 1    ! for the CAM modal aeosol, nsize_aer is always 1.
        call activate_modal( wbar, wmix, wdiab, wmin, wmax, tempair, rhoair, &
            naerosol(m,:), ntype_aer, &
            vaerosol(m,:), hygro(m,:),         &
            fn(m,:), fm(m,:), fluxn(m,:), fluxm(m,:), flux_fullact )

!   load results
	fmact(:,:) = fm(:,:)
	fnact(:,:) = fn(:,:)

!   diagnostics
	lun = ldiagaa_ecpp(125)
	if ((idiagaa_ecpp(125) > 0) .and. (lun > 0)) then

	if ((ktau /= ktau_save) .or. (ifrom_where /= ifrom_where_save)) &
		write(lun,'(//a,4i8)') &
		'activate_intface - ktau, ifrom_where =', ktau, ifrom_where 
	ktau_save = ktau
	ifrom_where_save = ifrom_where

	write(lun,'(2i3,2x,2i2,2x,4i2, 1p,2x,3e8.1, 0p,3x,3f7.3, 2(3x,4f6.3))') &
		jcls, k, jclsy, iccy, ido_actres, ihorzvert, maxd_asize, maxd_atype, & 
		naerosol(1,1:3)*1.0e-6, wbar_in, wmix_in, wbar, fmact(1,1:3), fnact(1,1:3)
	write(lun,'(8x,a, 1p,2x,4e10.2)') '  vaerosol', vaerosol(1,1:3)
	write(lun,'(8x,a, 1p,2x,4e10.2)') '  hygro   ', hygro(1,1:3)
	write(lun,'(8x,a, 1p,2x,6e10.2)') '  t,rho', tempair, rhoair

	end if


	return
	end subroutine parampollu_tdx_activate_intface
!==========================================================================================================

!----------------------------------------------------------------------------------------------------------
      subroutine loadaer0D(raer,qqcw,m,cs, phase, &
                         naerosol, vaerosol,  hygro )
!-------------------------------------------------------------------------
!     This subroutine is adopted from loadaer in ndrop.F90. It is 2D in ndrop.F90, 
!     but it is 0D here (single point). So that we do not need to define arrays with
!     pcols, pver.
!     Minghuai Wang, 2009-11
!-------------------------------------------------------------------------
      use modal_aero_data

      implicit none

!      load aerosol number, volume concentrations, and bulk hygroscopicity

       real(r8), intent(in) :: raer(pcnst) ! aerosol mass, number mixing ratios
       real(r8), intent(in) :: qqcw(pcnst) ! cloud-borne aerosol mass, number mixing ratios
       integer, intent(in) ::  m          ! m=mode index
       real(r8), intent(in) :: cs  ! air density (kg/m3)
       integer, intent(in) :: phase ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum
       real(r8), intent(out) :: naerosol                ! interstitial number conc (/m3)
       real(r8), intent(out) :: vaerosol       ! interstitial+activated volume conc (m3/m3)
       real(r8), intent(out) :: hygro ! bulk hygroscopicity of mode

!      internal

       real(r8) vol ! aerosol volume mixing ratio
       integer i,lnum,lnumcw,l,ltype,lmass,lmasscw

          vaerosol=0._r8
          hygro=0._r8

          do l=1,nspec_amode(m)
             lmass=lmassptr_amode(l,m) ! interstitial
             lmasscw=lmassptrcw_amode(l,m) ! cloud-borne
             ltype=lspectype_amode(l,m)
             if(phase.eq.3)then
                vol=max(raer(lmass)+qqcw(lmasscw),0._r8)/specdens_amode(ltype)
             elseif(phase.eq.2)then
                vol=max(qqcw(lmasscw),0._r8)/specdens_amode(ltype)
             elseif(phase.eq.1)then
                vol=max(raer(lmass),0._r8)/specdens_amode(ltype)
             else
                write(6,*)'phase=',phase,' in loadaer'
                call endrun('phase error in loadaer')
             endif
             vaerosol=vaerosol+vol
             hygro=hygro+vol*spechygro(ltype)
          enddo
          if (vaerosol > 1.0e-30_r8) then   ! +++xl add 8/2/2007
             hygro=hygro/(vaerosol)
             vaerosol=vaerosol*cs
          else
             hygro=0.0_r8
             vaerosol=0.0_r8
          endif

          lnum=numptr_amode(m)
          lnumcw=numptrcw_amode(m)
!            aerosol number predicted
             if(phase.eq.3)then
                naerosol=(raer(lnum)+qqcw(lnumcw))*cs
             elseif(phase.eq.2)then
                naerosol=qqcw(lnumcw)*cs
             else
                naerosol=raer(lnum)*cs
             endif
!            adjust number so that dgnumlo < dgnum < dgnumhi
             naerosol = max( naerosol, vaerosol*voltonumbhi_amode(m) )
             naerosol = min( naerosol, vaerosol*voltonumblo_amode(m) )

       return
       end subroutine loadaer0D
!============================================================================================

end module ecpp_modal_aero_activate
