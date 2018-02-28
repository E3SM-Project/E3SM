! modal_aero_coag.F90


!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_coag --- modal aerosol coagulation
!
! !INTERFACE:
   module modal_aero_coag

! !USES:
   use shr_kind_mod,    only:  r8 => shr_kind_r8
   use cam_logfile,     only:  iulog
   use chem_mods,       only:  gas_pcnst
   use modal_aero_data, only:  maxd_aspectype

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public :: modal_aero_coag_sub, modal_aero_coag_init, getcoags_wrapper_f

! !PUBLIC DATA MEMBERS:
  integer, parameter :: pcnstxx = gas_pcnst

#if ( defined MODAL_AERO_9MODE || defined MODAL_AERO_7MODE || defined MODAL_AERO_4MODE || defined MODAL_AERO_4MODE_MOM )
  integer, parameter, public :: pair_option_acoag = 3
#elif ( defined MODAL_AERO_3MODE )
  integer, parameter, public :: pair_option_acoag = 1
#endif
! specifies pairs of modes for which coagulation is calculated
!   1 -- [aitken-->accum] 
!   2 -- [aitken-->accum], and [pcarbon-->accum]
!   3 -- [aitken-->accum], [pcarbon-->accum], 
!        and [aitken-->pcarbon--(aging)-->accum] 
! other -- do no coag

  integer, parameter, public :: maxpair_acoag = 10
  integer, parameter, public :: maxspec_acoag = maxd_aspectype

  integer, public :: npair_acoag
  integer, public :: modefrm_acoag(maxpair_acoag)
  integer, public :: modetoo_acoag(maxpair_acoag)
  integer, public :: modetooeff_acoag(maxpair_acoag)
  integer, public :: nspecfrm_acoag(maxpair_acoag)
  integer, public :: lspecfrm_acoag(maxspec_acoag,maxpair_acoag)
  integer, public :: lspectoo_acoag(maxspec_acoag,maxpair_acoag)

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
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_coag_sub --- ...
!
! !INTERFACE:
   subroutine modal_aero_coag_sub(                               &
                        lchnk,    ncol,     nstep,               &
                        loffset,  deltat_main,                   &
                        t,        pmid,     pdel,                &
                        q,                                       &
                        dgncur_a,           dgncur_awet,         &
                        wetdens_a                                )


!----------------------------------------------------------------------
!   Authors: R. Easter
!----------------------------------------------------------------------

! !USES:
   use mo_constants,     only: pi
   use modal_aero_data
   use modal_aero_gasaerexch, only:  n_so4_monolayers_pcage, &
                                     soa_equivso4_factor

   use cam_abortutils,   only: endrun
   use cam_history,      only: outfld, fieldname_len
   use chem_mods,        only: adv_mass
   use constituents,     only: pcnst, cnst_name
   use physconst,        only: gravit, mwdry, r_universal
   use ppgrid,           only: pcols, pver
   use spmd_utils,       only: iam, masterproc
   use ref_pres,         only: top_lev => clim_modal_aero_top_lev

   implicit none

! !PARAMETERS:
   integer, intent(in)  :: lchnk            ! chunk identifier
   integer, intent(in)  :: ncol             ! number of columns in chunk
   integer, intent(in)  :: nstep            ! model step
   integer, intent(in)  :: loffset          ! offset applied to modal aero "pointers"

   real(r8), intent(in) :: deltat_main      ! model timestep (s)

   real(r8), intent(in) :: t(pcols,pver)    ! temperature (K)
   real(r8), intent(in) :: pmid(pcols,pver) ! pressure at model levels (Pa)
   real(r8), intent(in) :: pdel(pcols,pver) ! pressure thickness of levels (Pa)

   real(r8), intent(inout) :: q(ncol,pver,pcnstxx) 
                                            ! tracer mixing ratio (TMR) array
                                            ! *** MUST BE mol/mol-air or #/mol-air
                                            ! *** NOTE ncol & pcnstxx dimensions
   real(r8), intent(in) :: dgncur_a(pcols,pver,ntot_amode)
                                 ! dry geo. mean dia. (m) of number distrib.
   real(r8), intent(in) :: dgncur_awet(pcols,pver,ntot_amode)
                                 ! wet geo. mean dia. (m) of number distrib.
   real(r8), intent(in) :: wetdens_a(pcols,pver,ntot_amode) 
                                 ! density of wet aerosol (kg/m3)

! !DESCRIPTION: 
!   computes changes due to coagulation involving
!	aitken mode (modeptr_aitken) with accumulation mode (modeptr_accum)
!   this version will 
!	compute changes to mass and number, but not to surface area
!	calculates coagulation rate coefficients using either
!	    new CMAQ V4.6 fast method
!	    older cmaq slow method (direct gauss-hermite quadrature)
!
! !REVISION HISTORY:
!   RCE 07.04.15:  Adapted from MIRAGE2 code and CMAQ V4.6 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
	integer :: i, iok, ipair, ip_aitacc, ip_aitpca, ip_pcaacc, iq
	integer :: idomode(ntot_amode), iselfcoagdone(ntot_amode)
	integer :: jfreqcoag
	integer :: k
	integer :: l, l1, l2, la, lmz, lsfrm, lstoo, lunout
	integer :: modefrm, modetoo, mait, macc, mpca
	integer ::  n, nfreqcoag


	integer, save :: nerr = 0       ! number of errors for entire run
	integer, save :: nerrmax = 9999 ! maximum number of errors before abort
	integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1

	logical, parameter :: fastcoag_flag = .true. ! selects coag rate-coef method

	real(r8) :: aircon
      	real(r8) :: deltat, deltatinv_main
      	real(r8) :: dr_so4_monolayers_pcage
	real(r8) :: dryvol_a(pcols,pver,ntot_amode)
      	real(r8) :: dumexp, dumloss, dumprod
	real(r8) :: dumsfc_frm_old, dumsfc_frm_new
	real(r8) :: dum_m2v
	real(r8) :: fac_m2v_aitage(maxd_aspectype), fac_m2v_pcarbon(maxd_aspectype)
	real(r8) :: fac_volsfc_pcarbon
	real(r8) :: lnsg_frm, lnsg_too
	real(r8) :: sg_frm, sg_too
	real(r8) :: tmpa, tmpb, tmpc, tmpf, tmpg, tmph, tmpn
	real(r8) :: tmp1, tmp2
	real(r8) :: tmp_qold
	real(r8) :: v2ncur_a_tmp
	real(r8) :: vol_core, vol_shell
	real(r8) :: wetdens_frm, wetdens_too, wetdgnum_frm, wetdgnum_too
	real(r8) :: xbetaij0, xbetaij2i, xbetaij2j, xbetaij3, &
                    xbetaii0, xbetaii2,  xbetajj0,  xbetajj2     
      	real(r8) :: xferamt, xferfracvol, xferfrac_pcage, xferfrac_max
	real(r8) :: xnumbconc(ntot_amode)
	real(r8) :: xnumbconcavg(ntot_amode), xnumbconcnew(ntot_amode)
	real(r8) :: ybetaij0(maxpair_acoag), ybetaij3(maxpair_acoag)
	real(r8) :: ybetaii0(maxpair_acoag), ybetajj0(maxpair_acoag)

        real(r8) :: dqdt(ncol,pver,pcnstxx)  ! TMR "dq/dt" array - NOTE dims
        logical  :: dotend(pcnst)            ! identifies the species that
                                             ! tendencies are computed for
	real(r8) :: qsrflx(pcols)

        character(len=fieldname_len)   :: tmpname
        character(len=fieldname_len+3) :: fieldname

! begin
!   check if any coagulation pairs exist
	if (npair_acoag <= 0) return

!--------------------------------------------------------------------------------
!!$   if (ldiag1 > 0) then
!!$   if (nstep <= 3) then
!!$   do i = 1, ncol
!!$   if (lonndx(i) /= 37) cycle
!!$   if (latndx(i) /= 23) cycle
!!$   if (nstep > 3)       cycle
!!$   write( *, '(/a,i7,i5,2(2x,2i5))' )   &
!!$         '*** modal_aero_coag_sub -- nstep, iam, lat, lon, pcols, ncol =',   &
!!$         nstep, iam, latndx(i), lonndx(i), pcols, ncol
!!$   end do
!!$   end if
!!$!  if (ncol /= -999888777) return
!!$   if (nstep > 3) call endrun( 'modal_aero_coag_sub -- nstep>3 testing halt' )
!!$   end if   ! (ldiag1 > 0)
!--------------------------------------------------------------------------------

	dotend(:) = .false.
	dqdt(1:ncol,:,:) = 0.0_r8

	lunout = iulog


!
!   determine if coagulation will be done on this time-step
!   currently coagulation is done every 3 hours
!
!	deltat = 3600.0*3.0
	deltat = deltat_main
	nfreqcoag = max( 1, nint( deltat/deltat_main ) )
	jfreqcoag = nfreqcoag/2
	xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

	if (nfreqcoag .gt. 1) then
	    if ( mod(nstep,nfreqcoag) .ne. jfreqcoag ) return
	end if

!
!   set idomode
!
	idomode(:) = 0
	do ipair = 1, npair_acoag
	    idomode(modefrm_acoag(ipair)) = 1
	    idomode(modetoo_acoag(ipair)) = 1
	end do

!
!   other init
!
	macc = modeptr_accum
	mait = modeptr_aitken
	mpca = modeptr_pcarbon

	fac_m2v_aitage(:) = 0.0_r8
	fac_m2v_pcarbon(:) = 0.0_r8
	if (pair_option_acoag == 3) then
!   following ipair definitions MUST BE CONSISTENT with
!   the coding in modal_aero_coag_init for pair_option_acoag == 3
	    ip_aitacc = 1
	    ip_pcaacc = 2
	    ip_aitpca = 3

	! use 1 mol (bi-)sulfate = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
	    dr_so4_monolayers_pcage = n_so4_monolayers_pcage * 4.76e-10_r8

	    ipair = ip_aitpca
	    do iq = 1, nspecfrm_acoag(ipair)
		lsfrm = lspecfrm_acoag(iq,ipair)
		if (lsfrm == lptr_so4_a_amode(mait)) then
		    fac_m2v_aitage(iq) = specmw_so4_amode / specdens_so4_amode
		else if (lsfrm == lptr_nh4_a_amode(mait)) then
		    fac_m2v_aitage(iq) = specmw_nh4_amode / specdens_nh4_amode
		else if (lsfrm == lptr_soa_a_amode(mait)) then
		    fac_m2v_aitage(iq) = soa_equivso4_factor*   &
                                        (specmw_soa_amode / specdens_soa_amode)
!   for soa, the soa_equivso4_factor converts the soa volume into an
!	so4(+nh4) volume that has same hygroscopicity contribution as soa
!   this allows aging calculations to be done in terms of the amount
!	of (equivalent) so4(+nh4) in the shell
!   (see modal_aero_gasaerexch)
		end if
	    end do
	    
	    do l = 1, nspec_amode(mpca)
		l2 = lspectype_amode(l,mpca)
!   fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!		[m3-AP/kmol-AP]    = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
		fac_m2v_pcarbon(l) = specmw_amode(l2) / specdens_amode(l2)
	    end do

	    fac_volsfc_pcarbon = exp( 2.5_r8*(alnsg_amode(mpca)**2) )
	else
	    ip_aitacc = -999888777
	    ip_pcaacc = -999888777
	    ip_aitpca = -999888777
	end if

!
!   loop over levels and columns to calc the coagulation
!
!   integrate coagulation changes over deltat = nfreqcoag*deltat_main
!   then compute tendencies as
!      dqdt = (q(t+deltat) - q(t))/deltat_main
!   because tendencies are applied (in physics_update) over deltat_main
!
	deltat = nfreqcoag*deltat_main
	deltatinv_main = 1.0_r8/(deltat_main*(1.0_r8 + 1.0e-15_r8))

main_k: do k = top_lev, pver
main_i: do i = 1, ncol

!   air molar density (kmol/m3)
	aircon = (pmid(i,k)/(r_universal*t(i,k)))

!   calculate number conc. (#/m3) for modes doing coagulation
	do n = 1, ntot_amode
	    if (idomode(n) .gt. 0) then
		xnumbconc(n) = q(i,k,numptr_amode(n)-loffset)*aircon
		xnumbconc(n) = max( 0.0_r8, xnumbconc(n) ) 
	    end if
	    iselfcoagdone(n) = 0
	end do

!
!   calculate coagulation rates for each pair
!
main_ipair1: do ipair = 1, npair_acoag

	modefrm = modefrm_acoag(ipair)
	modetoo = modetoo_acoag(ipair)

!
! compute coagulation rates using cmaq "fast" method
!    (based on E. Whitby's approximation approach)
! here subr. arguments are all in mks unit
!
        call getcoags_wrapper_f(                                       &
          t(i,k), pmid(i,k),                                           &
          dgncur_awet(i,k,modefrm),     dgncur_awet(i,k,modetoo),      &
          sigmag_amode(modefrm),        sigmag_amode(modetoo),         &
          alnsg_amode(modefrm),         alnsg_amode(modetoo),          &
          wetdens_a(i,k,modefrm),       wetdens_a(i,k,modetoo),        &
          xbetaij0, xbetaij2i, xbetaij2j, xbetaij3,                    &
          xbetaii0, xbetaii2,  xbetajj0,  xbetajj2                     )


!   test diagnostics begin --------------------------------------------
!!$ 	if (ldiag2 > 0) then
!!$ 	if (nstep <= 3) then
!!$ 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) then
!!$ 	if ((mod(k-1,5) == 0) .or. (k>=23)) then
!!$
!!$	wetdgnum_frm = dgncur_awet(i,k,modefrm)
!!$	wetdgnum_too = dgncur_awet(i,k,modetoo)
!!$	wetdens_frm  = wetdens_a(i,k,modefrm)
!!$	wetdens_too  = wetdens_a(i,k,modetoo)
!!$	sg_frm   = sigmag_amode(modefrm)
!!$	sg_too   = sigmag_amode(modetoo)
!!$	lnsg_frm = alnsg_amode(modefrm)
!!$	lnsg_too = alnsg_amode(modetoo)
!!$
!!$        call getcoags_wrapper_f(                                       &
!!$          t(i,k), pmid(i,k),                                           &
!!$          wetdgnum_frm,                   wetdgnum_too,                &
!!$          sg_frm,                         sg_too,                      &
!!$          lnsg_frm,                       lnsg_too,                    &
!!$          wetdens_frm,                    wetdens_too,                 &
!!$          xbetaij0, xbetaij2i, xbetaij2j, xbetaij3,                    &
!!$          xbetaii0, xbetaii2,  xbetajj0,  xbetajj2                     )
!!$
!!$
!!$ 	    write(lunout,9801)
!!$ 	    write(lunout,9810) 'nstep,lat,lon,k,ipair   ',   &
!!$ 		nstep, latndx(i), lonndx(i), k, ipair
!!$ 	    write(lunout,9820) 'tk, pmb, aircon, pdel   ',   &
!!$ 		t(i,k), pmid(i,k)*1.0e-2_r8, aircon, pdel(i,k)*1.0e-2_r8
!!$ 	    write(lunout,9820) 'wetdens-cgs, sg      f/t',   &
!!$ 		wetdens_frm*1.0e-3_r8, wetdens_too*1.0e-3_r8,   &
!!$ 		sg_frm, sg_too
!!$ 	    write(lunout,9820) 'dgnwet-um, dgndry-um f/t',   &
!!$ 		1.0e6_r8*wetdgnum_frm, 1.0e6_r8*wetdgnum_too,   &
!!$ 		1.0e6_r8*dgncur_a(i,k,modefrm), 1.0e6_r8*dgncur_a(i,k,modetoo)
!!$ 	    write(lunout,9820) 'xbeta ij0, ij3, ii0, jj0',   &
!!$ 		xbetaij0, xbetaij3, xbetaii0, xbetajj0
!!$ 	    write(lunout,9820) 'xbeta ij2i & j, ii2, jj2',   &
!!$ 		xbetaij2i, xbetaij2j, xbetaii2, xbetajj2
!!$ 	    write(lunout,9820) 'numbii, numbjj, deltat  ',   &
!!$ 		xnumbconc(modefrm), xnumbconc(modetoo), deltat
!!$ 	    write(lunout,9820) 'loss ij3, ii0, jj0      ',   &
!!$ 		(xbetaij3*xnumbconc(modetoo)*deltat),   &
!!$ 		(xbetaij0*xnumbconc(modetoo)*deltat+    &
!!$ 		 xbetaii0*xnumbconc(modefrm)*deltat),   &
!!$ 		(xbetajj0*xnumbconc(modetoo)*deltat)
!!$ 9801	format( / 72x, 'ACOAG' )
!!$ 9810	format( 'ACOAG ', a, 2i8, 3i7, 3(1pe15.6) )
!!$ 9820	format( 'ACOAG ', a, 4(1pe15.6) )
!!$ 9830	format( 'ACOAG ', a, i1, a, 4(1pe15.6) )
!!$ 	end if
!!$ 	end if
!!$ 	end if
!!$ 	end if   ! (ldiag2 > 0)
!   test diagnostics end ----------------------------------------------

	ybetaij0(ipair) = xbetaij0
	ybetaij3(ipair) = xbetaij3
	ybetaii0(ipair) = xbetaii0
	ybetajj0(ipair) = xbetajj0

	end do main_ipair1



	if ( (pair_option_acoag == 1) .or.   &
	     (pair_option_acoag == 2) ) then
!
!   calculate number and mass changes for pair_option_acoag == 1,2
!
main_ipair2: do ipair = 1, npair_acoag

	modefrm = modefrm_acoag(ipair)
	modetoo = modetoo_acoag(ipair)

!   calculate number changes
!   apply self-coagulation losses only once to a mode (when iselfcoagdone=0)
!	first calc change to "too" mode
!	next  calc change to "frm" mode, using average number conc of "too"
	if ( (mprognum_amode(modetoo) > 0) .and.   &
	     (iselfcoagdone(modetoo) <= 0) ) then
	    iselfcoagdone(modetoo) = 1
	    tmpn = xnumbconc(modetoo)
	    xnumbconcnew(modetoo) = tmpn/(1.0_r8 + deltat*ybetajj0(ipair)*tmpn)
	    xnumbconcavg(modetoo) = 0.5_r8*(xnumbconcnew(modetoo) + tmpn)
	    lstoo = numptr_amode(modetoo) - loffset
	    q(i,k,lstoo) = xnumbconcnew(modetoo)/aircon
	    dqdt(i,k,lstoo) = (xnumbconcnew(modetoo)-tmpn)*deltatinv_main/aircon
	end if

	if ( (mprognum_amode(modefrm) > 0) .and.   &
	     (iselfcoagdone(modefrm) <= 0) ) then
	    iselfcoagdone(modefrm) = 1
	    tmpn = xnumbconc(modefrm)
	    tmpa = deltat*ybetaij0(ipair)*xnumbconcavg(modetoo)
	    tmpb = deltat*ybetaii0(ipair)
	    tmpc = tmpa + tmpb*tmpn
	    if (abs(tmpc) < 0.01_r8) then
		xnumbconcnew(modefrm) = tmpn*exp(-tmpc)
	    else if (abs(tmpa) < 0.001_r8) then
		xnumbconcnew(modefrm) =   &
		    exp(-tmpa)*tmpn/(1.0_r8 + tmpb*tmpn)
	    else
		tmpf = tmpb*tmpn/tmpc
		tmpg = exp(-tmpa)
		tmph = tmpg*(1.0_r8 - tmpf)/(1.0_r8 - tmpg*tmpf)
		xnumbconcnew(modefrm) = tmpn*max( 0.0_r8, min( 1.0_r8, tmph ) )
	    end if
	    xnumbconcavg(modefrm) = 0.5_r8*(xnumbconcnew(modefrm) + tmpn)
	    lsfrm = numptr_amode(modefrm) - loffset
	    q(i,k,lsfrm) = xnumbconcnew(modefrm)/aircon
	    dqdt(i,k,lsfrm) = (xnumbconcnew(modefrm)-tmpn)*deltatinv_main/aircon
	end if

!   calculate mass changes
!     xbetaij3*xnumbconc(modetoo) = first order loss rate for modefrm volume
!     xferfracvol = fraction of modefrm volume transferred to modetoo over deltat
	dumloss = ybetaij3(ipair)*xnumbconcavg(modetoo)
	xferfracvol = 1.0_r8 - exp( -dumloss*deltat )
	xferfracvol = max( 0.0_r8, min( xferfrac_max, xferfracvol ) )

	do iq = 1, nspecfrm_acoag(ipair)
	    lsfrm = lspecfrm_acoag(iq,ipair) - loffset
	    lstoo = lspectoo_acoag(iq,ipair) - loffset
	    if (lsfrm > 0) then
		xferamt = q(i,k,lsfrm)*xferfracvol
		dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferamt*deltatinv_main
		q(i,k,lsfrm) = q(i,k,lsfrm) - xferamt
		if (lstoo > 0) then
		    dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferamt*deltatinv_main
		    q(i,k,lstoo) = q(i,k,lstoo) + xferamt
		end if
	    end if
	end do

	end do main_ipair2


	else if (pair_option_acoag == 3) then
!
!   calculate number and mass changes for pair_option_acoag == 3
!

!   calculate number changes to accum mode
	if (mprognum_amode(macc) > 0) then
	    tmpn = xnumbconc(macc)
	    xnumbconcnew(macc) = tmpn/(1.0_r8 + deltat*ybetajj0(ip_aitacc)*tmpn)
	    xnumbconcavg(macc) = 0.5_r8*(xnumbconcnew(macc) + tmpn)
	    lstoo = numptr_amode(macc) - loffset
	    q(i,k,lstoo) = xnumbconcnew(macc)/aircon
	    dqdt(i,k,lstoo) = (xnumbconcnew(macc)-tmpn)*deltatinv_main/aircon
	end if

!   calculate number changes to primary carbon mode
	modefrm = modeptr_pcarbon
	if (mprognum_amode(mpca) > 0) then
	    tmpn = xnumbconc(mpca)
	    tmpa = deltat*ybetaij0(ip_pcaacc)*xnumbconcavg(macc)
	    tmpb = deltat*ybetaii0(ip_pcaacc)
	    tmpc = tmpa + tmpb*tmpn
	    if (abs(tmpc) < 0.01_r8) then
		xnumbconcnew(mpca) = tmpn*exp(-tmpc)
	    else if (abs(tmpa) < 0.001_r8) then
		xnumbconcnew(mpca) =   &
		    exp(-tmpa)*tmpn/(1.0_r8 + tmpb*tmpn)
	    else
		tmpf = tmpb*tmpn/tmpc
		tmpg = exp(-tmpa)
		tmph = tmpg*(1.0_r8 - tmpf)/(1.0_r8 - tmpg*tmpf)
		xnumbconcnew(mpca) = tmpn*max( 0.0_r8, min( 1.0_r8, tmph ) )
	    end if
	    xnumbconcavg(mpca) = 0.5_r8*(xnumbconcnew(mpca) + tmpn)
	    lsfrm = numptr_amode(mpca) - loffset
	    q(i,k,lsfrm) = xnumbconcnew(mpca)/aircon
	    dqdt(i,k,lsfrm) = (xnumbconcnew(mpca)-tmpn)*deltatinv_main/aircon
	end if

!   calculate number changes to aitken mode
	if (mprognum_amode(mait) > 0) then
	    tmpn = xnumbconc(mait)
	    tmpa = deltat*( ybetaij0(ip_aitacc)*xnumbconcavg(macc)   &
	                  + ybetaij0(ip_aitpca)*xnumbconcavg(mpca) )
	    tmpb = deltat*ybetaii0(ip_aitacc)
	    tmpc = tmpa + tmpb*tmpn
	    if (abs(tmpc) < 0.01_r8) then
		xnumbconcnew(mait) = tmpn*exp(-tmpc)
	    else if (abs(tmpa) < 0.001_r8) then
		xnumbconcnew(mait) =   &
		    exp(-tmpa)*tmpn/(1.0_r8 + tmpb*tmpn)
	    else
		tmpf = tmpb*tmpn/tmpc
		tmpg = exp(-tmpa)
		tmph = tmpg*(1.0_r8 - tmpf)/(1.0_r8 - tmpg*tmpf)
		xnumbconcnew(mait) = tmpn*max( 0.0_r8, min( 1.0_r8, tmph ) )
	    end if
	    xnumbconcavg(mait) = 0.5_r8*(xnumbconcnew(mait) + tmpn)
	    lsfrm = numptr_amode(mait) - loffset
	    q(i,k,lsfrm) = xnumbconcnew(mait)/aircon
	    dqdt(i,k,lsfrm) = (xnumbconcnew(mait)-tmpn)*deltatinv_main/aircon
	end if


!   calculate mass changes from aitken-->accum direct coagulation and
!	aitken-->pcarbon-->accum coagulation/aging
!   also calc volume of shell material (so4 & nh4 from aitken-->pcarbon)
	dumloss = ybetaij3(ip_aitacc)*xnumbconcavg(macc)   &
	        + ybetaij3(ip_aitpca)*xnumbconcavg(mpca)
	tmpa = ybetaij3(ip_aitpca)*xnumbconcavg(mpca)/max( dumloss, 1.0e-37_r8 )
	xferfracvol = 1.0_r8 - exp( -dumloss*deltat )
	xferfracvol = max( 0.0_r8, min( xferfrac_max, xferfracvol ) )
	vol_shell = 0.0_r8

	ipair = ip_aitacc
	do iq = 1, nspecfrm_acoag(ipair)
	    lsfrm = lspecfrm_acoag(iq,ipair) - loffset
	    lstoo = lspectoo_acoag(iq,ipair) - loffset
	    if (lsfrm > 0) then
		xferamt = q(i,k,lsfrm)*xferfracvol
		dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferamt*deltatinv_main
		q(i,k,lsfrm) = q(i,k,lsfrm) - xferamt
		if (lstoo > 0) then
		    dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferamt*deltatinv_main
		    q(i,k,lstoo) = q(i,k,lstoo) + xferamt
		end if
		vol_shell = vol_shell + xferamt*tmpa*fac_m2v_aitage(iq)
	    end if
	end do


!   now calculate aging transfer fraction for pcarbon-->accum
!   this duplicates the code in modal_aero_gasaerexch
	vol_core = 0.0_r8
	do l = 1, nspec_amode(mpca)
	    vol_core = vol_core + &
		q(i,k,lmassptr_amode(l,mpca)-loffset)*fac_m2v_pcarbon(l)
	end do
	tmp1 = vol_shell*dgncur_a(i,k,mpca)*fac_volsfc_pcarbon
	tmp2 = 6.0_r8*dr_so4_monolayers_pcage*vol_core
	tmp2 = max( tmp2, 0.0_r8 )
	if (tmp1 >= tmp2) then
	    xferfrac_pcage = xferfrac_max
	else
	    xferfrac_pcage = min( tmp1/tmp2, xferfrac_max )
	end if


!   calculate mass changes from pcarbon-->accum by direct coagulation
!   and aging
	dumloss = ybetaij3(ip_pcaacc)*xnumbconcavg(macc)
	xferfracvol = 1.0_r8 - exp( -dumloss*deltat )
	xferfracvol = xferfracvol + xferfrac_pcage
	xferfracvol = max( 0.0_r8, min( xferfrac_max, xferfracvol ) )

	ipair = ip_pcaacc
	do iq = 1, nspecfrm_acoag(ipair)
	    lsfrm = lspecfrm_acoag(iq,ipair) - loffset
	    lstoo = lspectoo_acoag(iq,ipair) - loffset
	    if (lsfrm > 0) then
		xferamt = q(i,k,lsfrm)*xferfracvol
		dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferamt*deltatinv_main
		q(i,k,lsfrm) = q(i,k,lsfrm) - xferamt
		if (lstoo > 0) then
		    dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferamt*deltatinv_main
		    q(i,k,lstoo) = q(i,k,lstoo) + xferamt
		end if
	    end if
	end do

	lsfrm = numptr_amode(mpca) - loffset
	lstoo = numptr_amode(macc) - loffset
	if (lsfrm > 0) then
	    xferamt = q(i,k,lsfrm)*xferfrac_pcage
	    dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferamt*deltatinv_main
	    q(i,k,lsfrm) = q(i,k,lsfrm) - xferamt
	    if (lstoo > 0) then
		dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferamt*deltatinv_main
		q(i,k,lstoo) = q(i,k,lstoo) + xferamt
	    end if
	end if



	else   ! (pair_option_acoag /= 1,2,3) then

	write(lunout,*) '*** modal_aero_coag_sub error'
	write(lunout,*) '    cannot do _coag_sub error pair_option_acoag =', &
		pair_option_acoag
	call endrun( 'modal_aero_coag_sub error' )


	end if   ! (pair_option_acoag == ...)


!   test diagnostics begin --------------------------------------------
!!$ 	if (ldiag3 > 0) then
!!$ 	if (nstep <= 3) then
!!$ 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) then
!!$ 	if ((mod(k-1,5) == 0) .or. (k>=23)) then
!!$ 	   if (pair_option_acoag == 3) then
!!$ 		write(*,*)
!!$ 		write(lunout,9820) 'xnumbconcavg ait,acc,pca', &
!!$ 		    xnumbconcavg(mait), xnumbconcavg(macc), xnumbconcavg(mpca)
!!$ 		write(lunout,9820) 'vshell, core            ', &
!!$ 		    vol_shell, vol_core
!!$ 		write(lunout,9820) 'dr_mono, dgn            ', &
!!$ 		    dr_so4_monolayers_pcage, dgncur_a(i,k,mpca)
!!$ 		write(lunout,9820) 'tmp1, tmp2              ', tmp1, tmp2
!!$ 		write(lunout,9820) 'xferfrac_age            ', xferfrac_pcage
!!$ 	   end if
!!$
!!$ 	   do ipair = 1, npair_acoag
!!$ 	   modefrm = modefrm_acoag(ipair)
!!$ 	   modetoo = modetoo_acoag(ipair)
!!$ 	   if (npair_acoag > 1) then
!!$ 		write(lunout,*)
!!$ 		write(lunout,9810) 'ipair =   ', ipair
!!$ 	   end if
!!$
!!$ 	   do iq = 1, nspecfrm_acoag(ipair)
!!$ 	   lsfrm = lspecfrm_acoag(iq,ipair) - loffset
!!$ 	   lstoo = lspectoo_acoag(iq,ipair) - loffset
!!$ 	   if (lsfrm > 0) then
!!$ 	   tmp_qold = q(i,k,lsfrm) - dqdt(i,k,lsfrm)*deltat_main
!!$!	   write(lunout,9820) 'm1 frm dqdt/q0,dqdt,q0/1',   &
!!$ 	   write(lunout,9830) 'm', iq,   &
!!$ 	                        ' frm dqdt/q0,dqdt,q0/1',   &
!!$ 		dqdt(i,k,lsfrm)/tmp_qold, dqdt(i,k,lsfrm), tmp_qold, q(i,k,lsfrm)
!!$ 	   end if
!!$ 	   if (lstoo > 0) then
!!$ 	   tmp_qold = q(i,k,lstoo) - dqdt(i,k,lstoo)*deltat_main
!!$ 	   write(lunout,9830) 'm', iq,   &
!!$ 	                        ' too dqdt/q0,dqdt,q0/1',   &
!!$ 		dqdt(i,k,lstoo)/tmp_qold, dqdt(i,k,lstoo), tmp_qold, q(i,k,lstoo)
!!$ 	   end if
!!$ 	   end do   ! iq
!!$
!!$ 	   lsfrm = numptr_amode(modefrm) - loffset
!!$ 	   lstoo = numptr_amode(modetoo) - loffset
!!$ 	   if (lsfrm > 0) then
!!$ 	   tmp_qold = q(i,k,lsfrm) - dqdt(i,k,lsfrm)*deltat_main
!!$ 	   write(lunout,9820) 'n  frm dqdt/q0,dqdt,q0/1',   &
!!$ 		dqdt(i,k,lsfrm)/tmp_qold, dqdt(i,k,lsfrm), tmp_qold, q(i,k,lsfrm)
!!$ 	   end if
!!$ 	   if (lstoo > 0) then
!!$ 	   tmp_qold = q(i,k,lstoo) - dqdt(i,k,lstoo)*deltat_main
!!$ 	   write(lunout,9820) 'n  too dqdt/q0,dqdt,q0/1',   &
!!$ 		dqdt(i,k,lstoo)/tmp_qold, dqdt(i,k,lstoo), tmp_qold, q(i,k,lstoo)
!!$ 	   end if
!!$
!!$ 	   end do   ! ipair
!!$ 	end if
!!$ 	end if
!!$ 	end if
!!$ 	end if   ! (ldiag3 > 0)
!   test diagnostics end ----------------------------------------------



	end do main_i
	end do main_k


! set dotend's
	do ipair = 1, npair_acoag
	    modefrm = modefrm_acoag(ipair)
	    modetoo = modetoo_acoag(ipair)

	    do iq = 1, nspecfrm_acoag(ipair)
		lsfrm = lspecfrm_acoag(iq,ipair) - loffset
		lstoo = lspectoo_acoag(iq,ipair) - loffset
		if (lsfrm > 0) dotend(lsfrm) = .true.
		if (lstoo > 0) dotend(lstoo) = .true.
	    end do

	    if (mprognum_amode(modefrm) > 0) then
		lsfrm = numptr_amode(modefrm) - loffset
		if (lsfrm > 0) dotend(lsfrm) = .true.
	    end if
	    if (mprognum_amode(modetoo) > 0) then
		lstoo = numptr_amode(modetoo) - loffset
		if (lstoo > 0) dotend(lstoo) = .true.
	    end if

	end do


!   do history file column-tendency fields
	do l = loffset+1, pcnst
	    lmz = l - loffset
	    if ( .not. dotend(lmz) ) cycle

	    qsrflx(:) = 0.0_r8
	    do k = top_lev, pver
	    do i = 1, ncol
		qsrflx(i) = qsrflx(i) + dqdt(i,k,lmz)*pdel(i,k)
	    end do
	    end do
	    qsrflx(:) = qsrflx(:)*(adv_mass(lmz)/(gravit*mwdry))
	    fieldname = trim(cnst_name(l)) // '_sfcoag1'
	    call outfld( fieldname, qsrflx, pcols, lchnk )
!	    if (( masterproc ) .and. (nstep < 1)) &
!		write(*,'(2(a,2x),1p,e11.3)') &
!		'modal_aero_coag_sub outfld', fieldname, adv_mass(lmz)
	end do ! l = ...


	return


!EOC
	end subroutine modal_aero_coag_sub


!----------------------------------------------------------------------
!----------------------------------------------------------------------
	subroutine modal_aero_coag_init
!
!   computes pointers for species transfer during coagulation
!
	use modal_aero_data
	use modal_aero_gasaerexch, only:  &
		modefrm_pcage, nspecfrm_pcage, lspecfrm_pcage, lspectoo_pcage

	use cam_abortutils,      only: endrun
	use cam_history,     only: addfld, horiz_only, add_default, fieldname_len
	use constituents,    only: pcnst, cnst_name
	use spmd_utils,      only: masterproc
        use phys_control,    only: phys_getopts

	implicit none

!   local variables
	integer :: ipair, iq, iqfrm, iqfrm_aa, iqtoo, iqtoo_aa
	integer :: l, lsfrm, lstoo, lunout
	integer :: m, mfrm, mtoo, mtef
	integer :: nsamefrm, nsametoo, nspec

	character(len=fieldname_len)   :: tmpname
	character(len=fieldname_len+3) :: fieldname
	character(128)                 :: long_name
	character(8)                   :: unit

	logical :: dotend(pcnst)
        logical :: history_aerosol      ! Output the MAM aerosol tendencies
 
        !-----------------------------------------------------------------------     
    
        call phys_getopts( history_aerosol_out        = history_aerosol   )

	lunout = iulog
!
!   define "from mode" and "to mode" for each coagulation pairing
!	currently just a2-->a1 coagulation
!
	if (pair_option_acoag == 1) then
	    npair_acoag = 1
	    modefrm_acoag(1) = modeptr_aitken
	    modetoo_acoag(1) = modeptr_accum
	    modetooeff_acoag(1) = modeptr_accum
	else if (pair_option_acoag == 2) then
	    npair_acoag = 2
	    modefrm_acoag(1) = modeptr_aitken
	    modetoo_acoag(1) = modeptr_accum
	    modetooeff_acoag(1) = modeptr_accum
	    modefrm_acoag(2) = modeptr_pcarbon
	    modetoo_acoag(2) = modeptr_accum
	    modetooeff_acoag(2) = modeptr_accum
	else if (pair_option_acoag == 3) then
	    npair_acoag = 3
	    modefrm_acoag(1) = modeptr_aitken
	    modetoo_acoag(1) = modeptr_accum
	    modetooeff_acoag(1) = modeptr_accum
	    modefrm_acoag(2) = modeptr_pcarbon
	    modetoo_acoag(2) = modeptr_accum
	    modetooeff_acoag(2) = modeptr_accum
	    modefrm_acoag(3) = modeptr_aitken
	    modetoo_acoag(3) = modeptr_pcarbon
	    modetooeff_acoag(3) = modeptr_accum
	    if (modefrm_pcage <= 0) then
		write(iulog,*) '*** modal_aero_coag_init error'
		write(iulog,*) '    pair_option_acoag, modefrm_pcage mismatch'
		write(iulog,*) '    pair_option_acoag, modefrm_pcage =', &
		    pair_option_acoag, modefrm_pcage
		call endrun( 'modal_aero_coag_init error' )
	    end if
	else
	    npair_acoag = 0
	    return
	end if

!
!   define species involved in each coagulation pairing
!	(include aerosol water)
!
aa_ipair: do ipair = 1, npair_acoag

	mfrm = modefrm_acoag(ipair)
	mtoo = modetoo_acoag(ipair)
	mtef = modetooeff_acoag(ipair)
	if ( (mfrm < 1) .or. (mfrm > ntot_amode) .or.   &
	     (mtoo < 1) .or. (mtoo > ntot_amode) .or.   &
	     (mtef < 1) .or. (mtef > ntot_amode) ) then
	    write(iulog,*) '*** modal_aero_coag_init error'
	    write(iulog,*) '    ipair, ntot_amode =', ipair, ntot_amode
	    write(iulog,*) '    mfrm, mtoo, mtef  =', mfrm, mtoo, mtef
	    call endrun( 'modal_aero_coag_init error' )
	end if


	mtoo = mtef   ! effective modetoo
	nspec = 0
aa_iqfrm: do iqfrm = 1, nspec_amode(mfrm)
	    lsfrm = lmassptr_amode(iqfrm,mfrm)
	    if ((lsfrm .lt. 1) .or. (lsfrm .gt. pcnst)) cycle aa_iqfrm

! find "too" species having same lspectype_amode as the "frm" species
! several species in a mode may have the same lspectype_amode, so also
!    use the ordering as a criterion (e.g., 1st <--> 1st, 2nd <--> 2nd)
	    iqfrm_aa = 1
	    iqtoo_aa = 1
	    if (iqfrm .gt. nspec_amode(mfrm)) then
		iqfrm_aa = nspec_amode(mfrm) + 1
		iqtoo_aa = nspec_amode(mtoo) + 1
	    end if
	    nsamefrm = 0
	    do iq = iqfrm_aa, iqfrm
		if ( lspectype_amode(iq   ,mfrm) .eq.   &
      		     lspectype_amode(iqfrm,mfrm) ) then
		    nsamefrm = nsamefrm + 1
		end if
	    end do
	    nsametoo = 0
	    lstoo = 0
	    do iqtoo = iqtoo_aa, nspec_amode(mtoo)
		if ( lspectype_amode(iqtoo,mtoo) .eq.   &
      		     lspectype_amode(iqfrm,mfrm) ) then
		    nsametoo = nsametoo + 1
		    if (nsametoo .eq. nsamefrm) then
			lstoo = lmassptr_amode(iqtoo,mtoo)
			exit
		    end if
		end if
	    end do

	    nspec = nspec + 1
	    lspecfrm_acoag(nspec,ipair) = lsfrm
	    lspectoo_acoag(nspec,ipair) = lstoo
	end do aa_iqfrm

!       lsfrm = lwaterptr_amode(mfrm)
!       if ((lsfrm .ge. 1) .and. (lsfrm .le. pcnst)) then
!           lstoo = lwaterptr_amode(mtoo)
!           if ((lstoo .lt. 1) .or. (lstoo .gt. pcnst)) lstoo = 0
!           nspec = nspec + 1
!           lspecfrm_acoag(nspec,ipair) = lsfrm
!           lspectoo_acoag(nspec,ipair) = lstoo
!       end if

	nspecfrm_acoag(ipair) = nspec
	end do aa_ipair

!
!   output results
!
	if ( masterproc ) then

	write(lunout,9310)

	do ipair = 1, npair_acoag
	  mfrm = modefrm_acoag(ipair)
	  mtoo = modetoo_acoag(ipair)
	  mtef = modetooeff_acoag(ipair)
	  write(lunout,9320) ipair, mfrm, mtoo, mtef

	  do iq = 1, nspecfrm_acoag(ipair)
	    lsfrm = lspecfrm_acoag(iq,ipair)
	    lstoo = lspectoo_acoag(iq,ipair)
	    if (lstoo .gt. 0) then
		write(lunout,9330) lsfrm, cnst_name(lsfrm),   &
      			lstoo, cnst_name(lstoo)
	    else
		write(lunout,9340) lsfrm, cnst_name(lsfrm)
	    end if
	  end do

	end do ! ipair = ...
	write(lunout,*)

	end if ! ( masterproc ) 

9310	format( / 'subr. modal_aero_coag_init' )
9320	format( 'pair', i3, 5x, 'mode', i3, &
		' ---> mode', i3, '   eff', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )


!
!   create history file column-tendency fields
!
	dotend(:) = .false.
	do ipair = 1, npair_acoag
	  do iq = 1, nspecfrm_acoag(ipair)
	    l = lspecfrm_acoag(iq,ipair)
	    if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
	    l = lspectoo_acoag(iq,ipair)
	    if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
	  end do

	  m = modefrm_acoag(ipair)
	  if ((m > 0) .and. (m <= ntot_amode)) then
	    l = numptr_amode(m)
	    if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
	  end if
	  m = modetoo_acoag(ipair)
	  if ((m > 0) .and. (m <= ntot_amode)) then
	    l = numptr_amode(m)
	    if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
	  end if
	end do ! ipair = ...

	if (pair_option_acoag == 3) then
	   do iq = 1, nspecfrm_pcage
	      lsfrm = lspecfrm_pcage(iq)
	      lstoo = lspectoo_pcage(iq)
	      if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
	         dotend(lsfrm) = .true.
	         if ((lstoo > 0) .and. (lstoo <= pcnst)) then
	            dotend(lstoo) = .true.
	         end if
	      end if
	   end do
	end if

	do l = 1, pcnst
	    if ( .not. dotend(l) ) cycle
	    tmpname = cnst_name(l)
	    unit = 'kg/m2/s'
	    do m = 1, ntot_amode
	        if (l == numptr_amode(m)) unit = '#/m2/s'
	    end do
	    fieldname = trim(tmpname) // '_sfcoag1'
	    long_name = trim(tmpname) // ' modal_aero coagulation column tendency'
	    call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
	    endif
	    if ( masterproc ) write(iulog,'(3(a,2x))') &
		'modal_aero_coag_init addfld', fieldname, unit
	end do ! l = ...
	if ( masterproc ) write(iulog,'(a)') &
		'modal_aero_coag_init ALL DONE'


	return
	end subroutine modal_aero_coag_init

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine getcoags_wrapper_f(              &
          airtemp, airprs,                        &
          dgatk, dgacc,                           &
          sgatk, sgacc,                           &
          xxlsgat, xxlsgac,                       &
          pdensat, pdensac,                       &
          betaij0, betaij2i, betaij2j, betaij3,   &
          betaii0, betaii2, betajj0, betajj2      )
        use physconst, only: p0 => pstd, &
                             tmelt, &
                             boltz
!
! interface to subr. getcoags
!
! interface code adapted from subr. aeroproc of cmaq v4.6,
!     with some of the parameter values from module aero_info_ae4
!
      implicit none

! *** arguments

      real(r8), intent(in) :: airtemp  ! air temperature [ k ]
      real(r8), intent(in) :: airprs   ! air pressure in [ pa ]

      real(r8), intent(in) :: dgatk    ! aitken mode geometric mean diameter [m]
      real(r8), intent(in) :: dgacc    ! accumulation mode geometric mean diam [m]

      real(r8), intent(in) :: sgatk    ! aitken mode geometric standard deviation
      real(r8), intent(in) :: sgacc    ! accumulation mode geometric standard deviation

      real(r8), intent(in) :: xxlsgat  ! natural log of geometric standard
      real(r8), intent(in) :: xxlsgac  !  deviations

      real(r8), intent(in) :: pdensat  ! aitken mode particle density [ kg / m**3 ]
      real(r8), intent(in) :: pdensac  ! accumulation mode density [ kg / m**3 ]

      real(r8), intent(out) :: betaij0, betaij2i, betaij2j, betaij3,   &
                               betaii0, betaii2,  betajj0,  betajj2


! *** local parameters
      real(r8) :: t0  ! standard surface temperature (15 deg C) [ k ]
      real(r8), parameter :: two3 = 2.0_r8/3.0_r8

! *** local variables
      real(r8) amu            ! atmospheric dynamic viscosity [ kg/m s ]
      real(r8) sqrt_temp      ! square root of ambient temperature
      real(r8) lamda          ! mean free path [ m ]

! *** intramodal coagulation rates [ m**3/s ] ( 0th & 2nd moments )
      real(r8)    batat( 2 )  ! aitken mode
      real(r8)    bacac( 2 )  ! accumulation mode
! *** intermodal coagulation rates [ m**3/s ] ( 0th & 2nd moments )
      real(r8)    batac( 2 )  ! aitken to accumulation
      real(r8)    bacat( 2 )  ! accumulation from aitken
! *** intermodal coagulation rate [ m**3/s ] ( 3rd moment )
      real(r8)    c3ij        ! aitken to accumulation
! *** 3rd moment intermodal transfer rate by coagulation
      real(r8)    c30atac     ! aitken to accumulation

! *** near continnuum regime (independent of mode)
      real(r8)    knc         ! knc = two3 * boltz *  airtemp / amu
! *** free molecular regime (depends upon modal density)
      real(r8)    kfmat       ! kfmat = sqrt(3.0*boltz*airtemp/pdensat)
      real(r8)    kfmac       ! kfmac = sqrt(3.0*boltz*airtemp/pdensac)
      real(r8)    kfmatac     ! kfmatac = sqrt( 6.0 * boltz * airtemp /
                              !                ( pdensat + pdensac ) )

      real(r8)    dumacc2, dumatk2, dumatk3

      t0 = tmelt + 15._r8

      sqrt_temp = sqrt( airtemp)

! *** calculate mean free path [ m ]:
!     6.6328e-8 is the sea level value given in table i.2.8
!     on page 10 of u.s. standard atmosphere 1962
      lamda = 6.6328e-8_r8 * p0 * airtemp  / ( t0 * airprs )

! *** calculate dynamic viscosity [ kg m**-1 s**-1 ]:
!     u.s. standard atmosphere 1962 page 14 expression
!     for dynamic viscosity is:
!     dynamic viscosity =  beta * t * sqrt(t) / ( t + s)
!     where beta = 1.458e-6 [ kg sec^-1 k**-0.5 ], s = 110.4 [ k ].
      amu = 1.458e-6_r8 * airtemp * sqrt_temp / ( airtemp + 110.4_r8 )

! *** coagulation
!     calculate coagulation coefficients using a method dictated by
!     the value of fastcoag_flag.  if true, the computationally-
!     efficient getcoags routine is used.  if false, the more intensive
!     gauss-hermite numerical quadrature method is used.  see section
!     2.1 of bhave et al. (2004) for further discussion.

! *** calculate term used in equation a6 of binkowski & shankar (1995)
      knc      = two3 * boltz *  airtemp / amu
! *** calculate terms used in equation a5 of binkowski & shankar (1995)
      kfmat    = sqrt( 3.0_r8 * boltz * airtemp / pdensat )
      kfmac    = sqrt( 3.0_r8 * boltz * airtemp / pdensac )
      kfmatac  = sqrt( 6.0_r8 * boltz * airtemp / ( pdensat + pdensac ) )

! *** transfer of number to accumulation mode from aitken mode is zero
      bacat(1) = 0.0_r8

! *** calculate intermodal and intramodal coagulation coefficients
!     for zeroth and second moments, and intermodal coagulation
!     coefficient for third moment
        call getcoags( lamda, kfmatac, kfmat, kfmac, knc,   &
                       dgatk,   dgacc,   sgatk,   sgacc,     &
                       xxlsgat,  xxlsgac,     &
                       batat(2), batat(1), bacac(2), bacac(1),   &
                       batac(2), bacat(2), batac(1), c3ij )

! convert from the "cmaq" coag rate parameters 
! to the "mirage2" parameters
        dumacc2 = ( (dgacc**2) * exp( 2.0_r8*xxlsgac*xxlsgac ) )
        dumatk2 = ( (dgatk**2) * exp( 2.0_r8*xxlsgat*xxlsgat ) )
        dumatk3 = ( (dgatk**3) * exp( 4.5_r8*xxlsgat*xxlsgat ) )

        betaii0  = max( 0.0_r8, batat(1) )
        betajj0  = max( 0.0_r8, bacac(1) )
        betaij0  = max( 0.0_r8, batac(1) )
        betaij3  = max( 0.0_r8, c3ij / dumatk3 )

        betajj2  = max( 0.0_r8, bacac(2) / dumacc2 )
        betaii2  = max( 0.0_r8, batat(2) / dumatk2 )
        betaij2i = max( 0.0_r8, batac(2) / dumatk2 )
        betaij2j = max( 0.0_r8, bacat(2) / dumatk2 )


      return
      end subroutine getcoags_wrapper_f



! //////////////////////////////////////////////////////////////////
!  subroutine getcoags calculates the coagulation rates using a new
!     approximate algorithm for the 2nd moment.  the 0th and 3rd moments
!     are done by analytic expressions from whitby et al. (1991).  the
!     correction factors are also similar to those from whitby et al.
!     (1991), but are derived from the gauss-hermite numerical
!     quadratures used by binkowski and roselle (2003).
!
!     called from aerostep as:
!     call getcoags( lamda, kfmatac, kfmat, kfmac, knc,
!                    dgat,dgac, sgatk, sgacc, xxlsgat,xxlsgac,
!                    batat(2), batat(1), bacac(2), bacac(1),
!                    batac(2), bacat(2), batac(1), c3ij )
!     where all input and outputs are real*8
!
!  revision history:
!   fsb 08/25/03 coded by dr. francis s. binkowksi
!
!   fsb 08/25/04 added in-line documentation
!
!   rce 04/15/2007 
!       code taken from cmaq v4.6 code; converted to f90;
!	added "intent" to subr arguments;
!       renamed "r4" & "r8" variables to "rx4" & "rx8";
!       changed "real*N" declarations to "real(rN)" (N = 4 or 8)
!
!  references:
!   1. whitby, e. r., p. h. mcmurry, u. shankar, and f. s. binkowski,
!   modal aerosol dynamics modeling, rep. 600/3-91/020, atmospheric
!   research and exposure assessment laboratory,
!   u.s. environmental protection agency, research triangle park, n.c.,
!   (ntis pb91-161729/as), 1991
!
!   2. binkowski, f.s. an u. shankar, the regional particulate matter
!   model 1. model decsription and preliminary results, journal of
!   geophysical research, 100, d12, pp 26,191-26,209,
!   december 20, 1995.
!
!   3. binkowski, f.s. and s.j. roselle, models-3 community
!      multiscale air quality (cmaq) model aerosol component 1:
!      model description.  j. geophys. res., vol 108, no d6, 4183
!      doi:10.1029/2001jd001409, 2003.


      subroutine getcoags( lamda, kfmatac, kfmat, kfmac, knc,   &
                           dgatk, dgacc, sgatk, sgacc, xxlsgat,xxlsgac,   &
                           qs11, qn11, qs22, qn22,   &
                           qs12, qs21, qn12, qv12 )

      implicit none

      real(r8), intent(in) ::  lamda     ! mean free path [ m ]

! *** coefficients for free molecular regime
      real(r8), intent(in) ::  kfmat     ! aitken mode
      real(r8), intent(in) ::  kfmac     ! accumulation mode
      real(r8), intent(in) ::  kfmatac   ! aitken to accumulation mode

      real(r8), intent(in) ::  knc   ! coefficient for near continnuum regime

! *** modal geometric mean diameters: [ m ]
      real(r8), intent(in) :: dgatk          ! aitken mode
      real(r8), intent(in) :: dgacc          ! accumulation mode

! *** modal geometric standard deviation
      real(r8), intent(in) :: sgatk          ! atken mode
      real(r8), intent(in) :: sgacc          ! accumulation mode

! *** natural log of modal geometric standard deviation
      real(r8), intent(in) :: xxlsgat         ! aitken mode
      real(r8), intent(in) :: xxlsgac         ! accumulation mode

! *** coagulation coefficients
      real(r8), intent(out) :: qs11, qn11, qs22, qn22,   &
                               qs12, qs21, qn12, qv12

      integer ibeta, n1, n2a, n2n ! indices for correction factors

      real(r8)  i1fm_at
      real(r8)  i1nc_at
      real(r8)  i1_at

      real(r8)  i1fm_ac
      real(r8)  i1nc_ac
      real(r8)  i1_ac

      real(r8)  i1fm
      real(r8)  i1nc
      real(r8)  i1

      real(r8) constii

      real(r8)    kngat, kngac
      real(r8)    one, two, half
       parameter( one = 1.0_r8, two = 2.0_r8, half = 0.5_r8 )
      real(r8)    a
!       parameter( a = 2.492_r8)
      parameter( a = 1.246_r8)
      real(r8)      two3rds
       parameter( two3rds = 2._r8 / 3._r8)

      real(r8)   sqrttwo  !  sqrt(two)
      real(r8)   dlgsqt2  !  1/ln( sqrt( 2 ) )


      real(r8)    esat01         ! aitken mode exp( log^2( sigmag )/8 )
      real(r8)    esac01         ! accumulation mode exp( log^2( sigmag )/8 )

      real(r8)    esat04
      real(r8)    esac04

      real(r8)    esat05
      real(r8)    esac05

      real(r8)    esat08
      real(r8)    esac08

      real(r8)    esat09
      real(r8)    esac09

      real(r8)    esat16
      real(r8)    esac16

      real(r8)    esat20
      real(r8)    esac20

      real(r8)    esat24
      real(r8)    esac24

      real(r8)    esat25
      real(r8)    esac25

      real(r8)    esat36
      real(r8)    esac36

      real(r8)    esat49

      real(r8)    esat64
      real(r8)    esac64

      real(r8)    esat100

      real(r8) dgat2, dgac2, dgat3, dgac3
      real(r8) sqdgat, sqdgac
      real(r8) sqdgat5, sqdgac5
      real(r8) sqdgat7
      real(r8) r, r2, r3, rx4, r5, r6, rx8
      real(r8) ri1, ri2, ri3, ri4
      real(r8) rat
      real(r8) coagfm0, coagnc0
      real(r8) coagfm3, coagnc3
      real(r8) coagfm_at, coagfm_ac
      real(r8) coagnc_at, coagnc_ac
      real(r8) coagatat0
      real(r8) coagacac0
      real(r8) coagatat2
      real(r8) coagacac2
      real(r8) coagatac0, coagatac3
      real(r8) coagatac2
      real(r8) coagacat2
      real(r8) xm2at, xm3at, xm2ac, xm3ac

! *** correction factors for coagulation rates
      real(r8), save :: bm0( 10 )          ! m0 intramodal fm - rpm values
      real(r8), save :: bm0ij( 10, 10, 10 ) ! m0 intermodal fm
      real(r8), save :: bm3i( 10, 10, 10 ) ! m3 intermodal fm- rpm values
      real(r8), save :: bm2ii(10) ! m2 intramodal fm
      real(r8), save :: bm2iitt(10) ! m2 intramodal total
      real(r8), save :: bm2ij(10,10,10) ! m2 intermodal fm i to j
      real(r8), save :: bm2ji(10,10,10) ! m2 total intermodal  j from i

! *** populate the arrays for the correction factors.

! rpm 0th moment correction factors for unimodal fm coagulation  rates
      data      bm0  /   &
            0.707106785165097_r8, 0.726148960080488_r8, 0.766430744110958_r8,   &
            0.814106389441342_r8, 0.861679526483207_r8, 0.903600509090092_r8,   &
            0.936578814219156_r8, 0.960098926735545_r8, 0.975646823342881_r8,   &
            0.985397173215326_r8   /


! fsb new fm correction factors for m0 intermodal coagulation

      data (bm0ij (  1,  1,ibeta), ibeta = 1,10) /   &
        0.628539_r8,  0.639610_r8,  0.664514_r8,  0.696278_r8,  0.731558_r8,   &
        0.768211_r8,  0.804480_r8,  0.838830_r8,  0.870024_r8,  0.897248_r8/
      data (bm0ij (  1,  2,ibeta), ibeta = 1,10) /   &
        0.639178_r8,  0.649966_r8,  0.674432_r8,  0.705794_r8,  0.740642_r8,   &
        0.776751_r8,  0.812323_r8,  0.845827_r8,  0.876076_r8,  0.902324_r8/
      data (bm0ij (  1,  3,ibeta), ibeta = 1,10) /   &
        0.663109_r8,  0.673464_r8,  0.697147_r8,  0.727637_r8,  0.761425_r8,   &
        0.796155_r8,  0.829978_r8,  0.861419_r8,  0.889424_r8,  0.913417_r8/
      data (bm0ij (  1,  4,ibeta), ibeta = 1,10) /   &
        0.693693_r8,  0.703654_r8,  0.726478_r8,  0.755786_r8,  0.787980_r8,   &
        0.820626_r8,  0.851898_r8,  0.880459_r8,  0.905465_r8,  0.926552_r8/
      data (bm0ij (  1,  5,ibeta), ibeta = 1,10) /   &
        0.727803_r8,  0.737349_r8,  0.759140_r8,  0.786870_r8,  0.816901_r8,   &
        0.846813_r8,  0.874906_r8,  0.900060_r8,  0.921679_r8,  0.939614_r8/
      data (bm0ij (  1,  6,ibeta), ibeta = 1,10) /   &
        0.763461_r8,  0.772483_r8,  0.792930_r8,  0.818599_r8,  0.845905_r8,   &
        0.872550_r8,  0.897051_r8,  0.918552_r8,  0.936701_r8,  0.951528_r8/
      data (bm0ij (  1,  7,ibeta), ibeta = 1,10) /   &
        0.799021_r8,  0.807365_r8,  0.826094_r8,  0.849230_r8,  0.873358_r8,   &
        0.896406_r8,  0.917161_r8,  0.935031_r8,  0.949868_r8,  0.961828_r8/
      data (bm0ij (  1,  8,ibeta), ibeta = 1,10) /   &
        0.833004_r8,  0.840514_r8,  0.857192_r8,  0.877446_r8,  0.898147_r8,   &
        0.917518_r8,  0.934627_r8,  0.949106_r8,  0.960958_r8,  0.970403_r8/
      data (bm0ij (  1,  9,ibeta), ibeta = 1,10) /   &
        0.864172_r8,  0.870734_r8,  0.885153_r8,  0.902373_r8,  0.919640_r8,   &
        0.935494_r8,  0.949257_r8,  0.960733_r8,  0.970016_r8,  0.977346_r8/
      data (bm0ij (  1, 10,ibeta), ibeta = 1,10) /   &
        0.891658_r8,  0.897227_r8,  0.909343_r8,  0.923588_r8,  0.937629_r8,   &
        0.950307_r8,  0.961151_r8,  0.970082_r8,  0.977236_r8,  0.982844_r8/
      data (bm0ij (  2,  1,ibeta), ibeta = 1,10) /   &
        0.658724_r8,  0.670587_r8,  0.697539_r8,  0.731890_r8,  0.769467_r8,   &
        0.807391_r8,  0.843410_r8,  0.875847_r8,  0.903700_r8,  0.926645_r8/
      data (bm0ij (  2,  2,ibeta), ibeta = 1,10) /   &
        0.667070_r8,  0.678820_r8,  0.705538_r8,  0.739591_r8,  0.776758_r8,   &
        0.814118_r8,  0.849415_r8,  0.881020_r8,  0.908006_r8,  0.930121_r8/
      data (bm0ij (  2,  3,ibeta), ibeta = 1,10) /   &
        0.686356_r8,  0.697839_r8,  0.723997_r8,  0.757285_r8,  0.793389_r8,   &
        0.829313_r8,  0.862835_r8,  0.892459_r8,  0.917432_r8,  0.937663_r8/
      data (bm0ij (  2,  4,ibeta), ibeta = 1,10) /   &
        0.711425_r8,  0.722572_r8,  0.747941_r8,  0.780055_r8,  0.814518_r8,   &
        0.848315_r8,  0.879335_r8,  0.906290_r8,  0.928658_r8,  0.946526_r8/
      data (bm0ij (  2,  5,ibeta), ibeta = 1,10) /   &
        0.739575_r8,  0.750307_r8,  0.774633_r8,  0.805138_r8,  0.837408_r8,   &
        0.868504_r8,  0.896517_r8,  0.920421_r8,  0.939932_r8,  0.955299_r8/
      data (bm0ij (  2,  6,ibeta), ibeta = 1,10) /   &
        0.769143_r8,  0.779346_r8,  0.802314_r8,  0.830752_r8,  0.860333_r8,   &
        0.888300_r8,  0.913014_r8,  0.933727_r8,  0.950370_r8,  0.963306_r8/
      data (bm0ij (  2,  7,ibeta), ibeta = 1,10) /   &
        0.798900_r8,  0.808431_r8,  0.829700_r8,  0.855653_r8,  0.882163_r8,   &
        0.906749_r8,  0.928075_r8,  0.945654_r8,  0.959579_r8,  0.970280_r8/
      data (bm0ij (  2,  8,ibeta), ibeta = 1,10) /   &
        0.827826_r8,  0.836542_r8,  0.855808_r8,  0.878954_r8,  0.902174_r8,   &
        0.923316_r8,  0.941345_r8,  0.955989_r8,  0.967450_r8,  0.976174_r8/
      data (bm0ij (  2,  9,ibeta), ibeta = 1,10) /   &
        0.855068_r8,  0.862856_r8,  0.879900_r8,  0.900068_r8,  0.919956_r8,   &
        0.937764_r8,  0.952725_r8,  0.964726_r8,  0.974027_r8,  0.981053_r8/
      data (bm0ij (  2, 10,ibeta), ibeta = 1,10) /   &
        0.879961_r8,  0.886755_r8,  0.901484_r8,  0.918665_r8,  0.935346_r8,   &
        0.950065_r8,  0.962277_r8,  0.971974_r8,  0.979432_r8,  0.985033_r8/
      data (bm0ij (  3,  1,ibeta), ibeta = 1,10) /   &
        0.724166_r8,  0.735474_r8,  0.761359_r8,  0.794045_r8,  0.828702_r8,   &
        0.862061_r8,  0.891995_r8,  0.917385_r8,  0.937959_r8,  0.954036_r8/
      data (bm0ij (  3,  2,ibeta), ibeta = 1,10) /   &
        0.730416_r8,  0.741780_r8,  0.767647_r8,  0.800116_r8,  0.834344_r8,   &
        0.867093_r8,  0.896302_r8,  0.920934_r8,  0.940790_r8,  0.956237_r8/
      data (bm0ij (  3,  3,ibeta), ibeta = 1,10) /   &
        0.745327_r8,  0.756664_r8,  0.782255_r8,  0.814026_r8,  0.847107_r8,   &
        0.878339_r8,  0.905820_r8,  0.928699_r8,  0.946931_r8,  0.960977_r8/
      data (bm0ij (  3,  4,ibeta), ibeta = 1,10) /   &
        0.765195_r8,  0.776312_r8,  0.801216_r8,  0.831758_r8,  0.863079_r8,   &
        0.892159_r8,  0.917319_r8,  0.937939_r8,  0.954145_r8,  0.966486_r8/
      data (bm0ij (  3,  5,ibeta), ibeta = 1,10) /   &
        0.787632_r8,  0.798347_r8,  0.822165_r8,  0.850985_r8,  0.880049_r8,   &
        0.906544_r8,  0.929062_r8,  0.947218_r8,  0.961288_r8,  0.971878_r8/
      data (bm0ij (  3,  6,ibeta), ibeta = 1,10) /   &
        0.811024_r8,  0.821179_r8,  0.843557_r8,  0.870247_r8,  0.896694_r8,   &
        0.920365_r8,  0.940131_r8,  0.955821_r8,  0.967820_r8,  0.976753_r8/
      data (bm0ij (  3,  7,ibeta), ibeta = 1,10) /   &
        0.834254_r8,  0.843709_r8,  0.864356_r8,  0.888619_r8,  0.912245_r8,   &
        0.933019_r8,  0.950084_r8,  0.963438_r8,  0.973530_r8,  0.980973_r8/
      data (bm0ij (  3,  8,ibeta), ibeta = 1,10) /   &
        0.856531_r8,  0.865176_r8,  0.883881_r8,  0.905544_r8,  0.926290_r8,   &
        0.944236_r8,  0.958762_r8,  0.969988_r8,  0.978386_r8,  0.984530_r8/
      data (bm0ij (  3,  9,ibeta), ibeta = 1,10) /   &
        0.877307_r8,  0.885070_r8,  0.901716_r8,  0.920729_r8,  0.938663_r8,   &
        0.953951_r8,  0.966169_r8,  0.975512_r8,  0.982442_r8,  0.987477_r8/
      data (bm0ij (  3, 10,ibeta), ibeta = 1,10) /   &
        0.896234_r8,  0.903082_r8,  0.917645_r8,  0.934069_r8,  0.949354_r8,   &
        0.962222_r8,  0.972396_r8,  0.980107_r8,  0.985788_r8,  0.989894_r8/
      data (bm0ij (  4,  1,ibeta), ibeta = 1,10) /   &
        0.799294_r8,  0.809144_r8,  0.831293_r8,  0.858395_r8,  0.885897_r8,   &
        0.911031_r8,  0.932406_r8,  0.949642_r8,  0.963001_r8,  0.973062_r8/
      data (bm0ij (  4,  2,ibeta), ibeta = 1,10) /   &
        0.804239_r8,  0.814102_r8,  0.836169_r8,  0.862984_r8,  0.890003_r8,   &
        0.914535_r8,  0.935274_r8,  0.951910_r8,  0.964748_r8,  0.974381_r8/
      data (bm0ij (  4,  3,ibeta), ibeta = 1,10) /   &
        0.815910_r8,  0.825708_r8,  0.847403_r8,  0.873389_r8,  0.899185_r8,   &
        0.922275_r8,  0.941543_r8,  0.956826_r8,  0.968507_r8,  0.977204_r8/
      data (bm0ij (  4,  4,ibeta), ibeta = 1,10) /   &
        0.831348_r8,  0.840892_r8,  0.861793_r8,  0.886428_r8,  0.910463_r8,   &
        0.931614_r8,  0.948993_r8,  0.962593_r8,  0.972872_r8,  0.980456_r8/
      data (bm0ij (  4,  5,ibeta), ibeta = 1,10) /   &
        0.848597_r8,  0.857693_r8,  0.877402_r8,  0.900265_r8,  0.922180_r8,   &
        0.941134_r8,  0.956464_r8,  0.968298_r8,  0.977143_r8,  0.983611_r8/
      data (bm0ij (  4,  6,ibeta), ibeta = 1,10) /   &
        0.866271_r8,  0.874764_r8,  0.892984_r8,  0.913796_r8,  0.933407_r8,   &
        0.950088_r8,  0.963380_r8,  0.973512_r8,  0.981006_r8,  0.986440_r8/
      data (bm0ij (  4,  7,ibeta), ibeta = 1,10) /   &
        0.883430_r8,  0.891216_r8,  0.907762_r8,  0.926388_r8,  0.943660_r8,   &
        0.958127_r8,  0.969499_r8,  0.978070_r8,  0.984351_r8,  0.988872_r8/
      data (bm0ij (  4,  8,ibeta), ibeta = 1,10) /   &
        0.899483_r8,  0.906505_r8,  0.921294_r8,  0.937719_r8,  0.952729_r8,   &
        0.965131_r8,  0.974762_r8,  0.981950_r8,  0.987175_r8,  0.990912_r8/
      data (bm0ij (  4,  9,ibeta), ibeta = 1,10) /   &
        0.914096_r8,  0.920337_r8,  0.933373_r8,  0.947677_r8,  0.960579_r8,   &
        0.971111_r8,  0.979206_r8,  0.985196_r8,  0.989520_r8,  0.992597_r8/
      data (bm0ij (  4, 10,ibeta), ibeta = 1,10) /   &
        0.927122_r8,  0.932597_r8,  0.943952_r8,  0.956277_r8,  0.967268_r8,   &
        0.976147_r8,  0.982912_r8,  0.987882_r8,  0.991450_r8,  0.993976_r8/
      data (bm0ij (  5,  1,ibeta), ibeta = 1,10) /   &
        0.865049_r8,  0.872851_r8,  0.889900_r8,  0.909907_r8,  0.929290_r8,   &
        0.946205_r8,  0.959991_r8,  0.970706_r8,  0.978764_r8,  0.984692_r8/
      data (bm0ij (  5,  2,ibeta), ibeta = 1,10) /   &
        0.868989_r8,  0.876713_r8,  0.893538_r8,  0.913173_r8,  0.932080_r8,   &
        0.948484_r8,  0.961785_r8,  0.972080_r8,  0.979796_r8,  0.985457_r8/
      data (bm0ij (  5,  3,ibeta), ibeta = 1,10) /   &
        0.878010_r8,  0.885524_r8,  0.901756_r8,  0.920464_r8,  0.938235_r8,   &
        0.953461_r8,  0.965672_r8,  0.975037_r8,  0.982005_r8,  0.987085_r8/
      data (bm0ij (  5,  4,ibeta), ibeta = 1,10) /   &
        0.889534_r8,  0.896698_r8,  0.912012_r8,  0.929395_r8,  0.945647_r8,   &
        0.959366_r8,  0.970227_r8,  0.978469_r8,  0.984547_r8,  0.988950_r8/
      data (bm0ij (  5,  5,ibeta), ibeta = 1,10) /   &
        0.902033_r8,  0.908713_r8,  0.922848_r8,  0.938648_r8,  0.953186_r8,   &
        0.965278_r8,  0.974729_r8,  0.981824_r8,  0.987013_r8,  0.990746_r8/
      data (bm0ij (  5,  6,ibeta), ibeta = 1,10) /   &
        0.914496_r8,  0.920599_r8,  0.933389_r8,  0.947485_r8,  0.960262_r8,   &
        0.970743_r8,  0.978839_r8,  0.984858_r8,  0.989225_r8,  0.992348_r8/
      data (bm0ij (  5,  7,ibeta), ibeta = 1,10) /   &
        0.926281_r8,  0.931761_r8,  0.943142_r8,  0.955526_r8,  0.966600_r8,   &
        0.975573_r8,  0.982431_r8,  0.987485_r8,  0.991128_r8,  0.993718_r8/
      data (bm0ij (  5,  8,ibeta), ibeta = 1,10) /   &
        0.937029_r8,  0.941877_r8,  0.951868_r8,  0.962615_r8,  0.972112_r8,   &
        0.979723_r8,  0.985488_r8,  0.989705_r8,  0.992725_r8,  0.994863_r8/
      data (bm0ij (  5,  9,ibeta), ibeta = 1,10) /   &
        0.946580_r8,  0.950819_r8,  0.959494_r8,  0.968732_r8,  0.976811_r8,   &
        0.983226_r8,  0.988047_r8,  0.991550_r8,  0.994047_r8,  0.995806_r8/
      data (bm0ij (  5, 10,ibeta), ibeta = 1,10) /   &
        0.954909_r8,  0.958581_r8,  0.966049_r8,  0.973933_r8,  0.980766_r8,   &
        0.986149_r8,  0.990166_r8,  0.993070_r8,  0.995130_r8,  0.996577_r8/
      data (bm0ij (  6,  1,ibeta), ibeta = 1,10) /   &
        0.914182_r8,  0.919824_r8,  0.931832_r8,  0.945387_r8,  0.957999_r8,   &
        0.968606_r8,  0.976982_r8,  0.983331_r8,  0.988013_r8,  0.991407_r8/
      data (bm0ij (  6,  2,ibeta), ibeta = 1,10) /   &
        0.917139_r8,  0.922665_r8,  0.934395_r8,  0.947580_r8,  0.959792_r8,   &
        0.970017_r8,  0.978062_r8,  0.984138_r8,  0.988609_r8,  0.991843_r8/
      data (bm0ij (  6,  3,ibeta), ibeta = 1,10) /   &
        0.923742_r8,  0.928990_r8,  0.940064_r8,  0.952396_r8,  0.963699_r8,   &
        0.973070_r8,  0.980381_r8,  0.985866_r8,  0.989878_r8,  0.992768_r8/
      data (bm0ij (  6,  4,ibeta), ibeta = 1,10) /   &
        0.931870_r8,  0.936743_r8,  0.946941_r8,  0.958162_r8,  0.968318_r8,   &
        0.976640_r8,  0.983069_r8,  0.987853_r8,  0.991330_r8,  0.993822_r8/
      data (bm0ij (  6,  5,ibeta), ibeta = 1,10) /   &
        0.940376_r8,  0.944807_r8,  0.954004_r8,  0.963999_r8,  0.972928_r8,   &
        0.980162_r8,  0.985695_r8,  0.989779_r8,  0.992729_r8,  0.994833_r8/
      data (bm0ij (  6,  6,ibeta), ibeta = 1,10) /   &
        0.948597_r8,  0.952555_r8,  0.960703_r8,  0.969454_r8,  0.977181_r8,   &
        0.983373_r8,  0.988067_r8,  0.991507_r8,  0.993977_r8,  0.995730_r8/
      data (bm0ij (  6,  7,ibeta), ibeta = 1,10) /   &
        0.956167_r8,  0.959648_r8,  0.966763_r8,  0.974326_r8,  0.980933_r8,   &
        0.986177_r8,  0.990121_r8,  0.992993_r8,  0.995045_r8,  0.996495_r8/
      data (bm0ij (  6,  8,ibeta), ibeta = 1,10) /   &
        0.962913_r8,  0.965937_r8,  0.972080_r8,  0.978552_r8,  0.984153_r8,   &
        0.988563_r8,  0.991857_r8,  0.994242_r8,  0.995938_r8,  0.997133_r8/
      data (bm0ij (  6,  9,ibeta), ibeta = 1,10) /   &
        0.968787_r8,  0.971391_r8,  0.976651_r8,  0.982148_r8,  0.986869_r8,   &
        0.990560_r8,  0.993301_r8,  0.995275_r8,  0.996675_r8,  0.997657_r8/
      data (bm0ij (  6, 10,ibeta), ibeta = 1,10) /   &
        0.973822_r8,  0.976047_r8,  0.980523_r8,  0.985170_r8,  0.989134_r8,   &
        0.992215_r8,  0.994491_r8,  0.996124_r8,  0.997277_r8,  0.998085_r8/
      data (bm0ij (  7,  1,ibeta), ibeta = 1,10) /   &
        0.947410_r8,  0.951207_r8,  0.959119_r8,  0.967781_r8,  0.975592_r8,   &
        0.981981_r8,  0.986915_r8,  0.990590_r8,  0.993266_r8,  0.995187_r8/
      data (bm0ij (  7,  2,ibeta), ibeta = 1,10) /   &
        0.949477_r8,  0.953161_r8,  0.960824_r8,  0.969187_r8,  0.976702_r8,   &
        0.982831_r8,  0.987550_r8,  0.991057_r8,  0.993606_r8,  0.995434_r8/
      data (bm0ij (  7,  3,ibeta), ibeta = 1,10) /   &
        0.954008_r8,  0.957438_r8,  0.964537_r8,  0.972232_r8,  0.979095_r8,   &
        0.984653_r8,  0.988907_r8,  0.992053_r8,  0.994330_r8,  0.995958_r8/
      data (bm0ij (  7,  4,ibeta), ibeta = 1,10) /   &
        0.959431_r8,  0.962539_r8,  0.968935_r8,  0.975808_r8,  0.981882_r8,   &
        0.986759_r8,  0.990466_r8,  0.993190_r8,  0.995153_r8,  0.996552_r8/
      data (bm0ij (  7,  5,ibeta), ibeta = 1,10) /   &
        0.964932_r8,  0.967693_r8,  0.973342_r8,  0.979355_r8,  0.984620_r8,   &
        0.988812_r8,  0.991974_r8,  0.994285_r8,  0.995943_r8,  0.997119_r8/
      data (bm0ij (  7,  6,ibeta), ibeta = 1,10) /   &
        0.970101_r8,  0.972517_r8,  0.977428_r8,  0.982612_r8,  0.987110_r8,   &
        0.990663_r8,  0.993326_r8,  0.995261_r8,  0.996644_r8,  0.997621_r8/
      data (bm0ij (  7,  7,ibeta), ibeta = 1,10) /   &
        0.974746_r8,  0.976834_r8,  0.981055_r8,  0.985475_r8,  0.989280_r8,   &
        0.992265_r8,  0.994488_r8,  0.996097_r8,  0.997241_r8,  0.998048_r8/
      data (bm0ij (  7,  8,ibeta), ibeta = 1,10) /   &
        0.978804_r8,  0.980591_r8,  0.984187_r8,  0.987927_r8,  0.991124_r8,   &
        0.993617_r8,  0.995464_r8,  0.996795_r8,  0.997739_r8,  0.998403_r8/
      data (bm0ij (  7,  9,ibeta), ibeta = 1,10) /   &
        0.982280_r8,  0.983799_r8,  0.986844_r8,  0.989991_r8,  0.992667_r8,   &
        0.994742_r8,  0.996273_r8,  0.997372_r8,  0.998149_r8,  0.998695_r8/
      data (bm0ij (  7, 10,ibeta), ibeta = 1,10) /   &
        0.985218_r8,  0.986503_r8,  0.989071_r8,  0.991711_r8,  0.993945_r8,   &
        0.995669_r8,  0.996937_r8,  0.997844_r8,  0.998484_r8,  0.998932_r8/
      data (bm0ij (  8,  1,ibeta), ibeta = 1,10) /   &
        0.968507_r8,  0.970935_r8,  0.975916_r8,  0.981248_r8,  0.985947_r8,   &
        0.989716_r8,  0.992580_r8,  0.994689_r8,  0.996210_r8,  0.997297_r8/
      data (bm0ij (  8,  2,ibeta), ibeta = 1,10) /   &
        0.969870_r8,  0.972210_r8,  0.977002_r8,  0.982119_r8,  0.986619_r8,   &
        0.990219_r8,  0.992951_r8,  0.994958_r8,  0.996405_r8,  0.997437_r8/
      data (bm0ij (  8,  3,ibeta), ibeta = 1,10) /   &
        0.972820_r8,  0.974963_r8,  0.979339_r8,  0.983988_r8,  0.988054_r8,   &
        0.991292_r8,  0.993738_r8,  0.995529_r8,  0.996817_r8,  0.997734_r8/
      data (bm0ij (  8,  4,ibeta), ibeta = 1,10) /   &
        0.976280_r8,  0.978186_r8,  0.982060_r8,  0.986151_r8,  0.989706_r8,   &
        0.992520_r8,  0.994636_r8,  0.996179_r8,  0.997284_r8,  0.998069_r8/
      data (bm0ij (  8,  5,ibeta), ibeta = 1,10) /   &
        0.979711_r8,  0.981372_r8,  0.984735_r8,  0.988263_r8,  0.991309_r8,   &
        0.993706_r8,  0.995499_r8,  0.996801_r8,  0.997730_r8,  0.998389_r8/
      data (bm0ij (  8,  6,ibeta), ibeta = 1,10) /   &
        0.982863_r8,  0.984292_r8,  0.987172_r8,  0.990174_r8,  0.992750_r8,   &
        0.994766_r8,  0.996266_r8,  0.997352_r8,  0.998125_r8,  0.998670_r8/
      data (bm0ij (  8,  7,ibeta), ibeta = 1,10) /   &
        0.985642_r8,  0.986858_r8,  0.989301_r8,  0.991834_r8,  0.993994_r8,   &
        0.995676_r8,  0.996923_r8,  0.997822_r8,  0.998460_r8,  0.998910_r8/
      data (bm0ij (  8,  8,ibeta), ibeta = 1,10) /   &
        0.988029_r8,  0.989058_r8,  0.991116_r8,  0.993240_r8,  0.995043_r8,   &
        0.996440_r8,  0.997472_r8,  0.998214_r8,  0.998739_r8,  0.999108_r8/
      data (bm0ij (  8,  9,ibeta), ibeta = 1,10) /   &
        0.990046_r8,  0.990912_r8,  0.992640_r8,  0.994415_r8,  0.995914_r8,   &
        0.997073_r8,  0.997925_r8,  0.998536_r8,  0.998968_r8,  0.999271_r8/
      data (bm0ij (  8, 10,ibeta), ibeta = 1,10) /   &
        0.991732_r8,  0.992459_r8,  0.993906_r8,  0.995386_r8,  0.996633_r8,   &
        0.997592_r8,  0.998296_r8,  0.998799_r8,  0.999154_r8,  0.999403_r8/
      data (bm0ij (  9,  1,ibeta), ibeta = 1,10) /   &
        0.981392_r8,  0.982893_r8,  0.985938_r8,  0.989146_r8,  0.991928_r8,   &
        0.994129_r8,  0.995783_r8,  0.996991_r8,  0.997857_r8,  0.998473_r8/
      data (bm0ij (  9,  2,ibeta), ibeta = 1,10) /   &
        0.982254_r8,  0.983693_r8,  0.986608_r8,  0.989673_r8,  0.992328_r8,   &
        0.994424_r8,  0.995998_r8,  0.997146_r8,  0.997969_r8,  0.998553_r8/
      data (bm0ij (  9,  3,ibeta), ibeta = 1,10) /   &
        0.984104_r8,  0.985407_r8,  0.988040_r8,  0.990798_r8,  0.993178_r8,   &
        0.995052_r8,  0.996454_r8,  0.997474_r8,  0.998204_r8,  0.998722_r8/
      data (bm0ij (  9,  4,ibeta), ibeta = 1,10) /   &
        0.986243_r8,  0.987386_r8,  0.989687_r8,  0.992087_r8,  0.994149_r8,   &
        0.995765_r8,  0.996971_r8,  0.997846_r8,  0.998470_r8,  0.998913_r8/
      data (bm0ij (  9,  5,ibeta), ibeta = 1,10) /   &
        0.988332_r8,  0.989313_r8,  0.991284_r8,  0.993332_r8,  0.995082_r8,   &
        0.996449_r8,  0.997465_r8,  0.998200_r8,  0.998723_r8,  0.999093_r8/
      data (bm0ij (  9,  6,ibeta), ibeta = 1,10) /   &
        0.990220_r8,  0.991053_r8,  0.992721_r8,  0.994445_r8,  0.995914_r8,   &
        0.997056_r8,  0.997902_r8,  0.998513_r8,  0.998947_r8,  0.999253_r8/
      data (bm0ij (  9,  7,ibeta), ibeta = 1,10) /   &
        0.991859_r8,  0.992561_r8,  0.993961_r8,  0.995403_r8,  0.996626_r8,   &
        0.997574_r8,  0.998274_r8,  0.998778_r8,  0.999136_r8,  0.999387_r8/
      data (bm0ij (  9,  8,ibeta), ibeta = 1,10) /   &
        0.993250_r8,  0.993837_r8,  0.995007_r8,  0.996208_r8,  0.997223_r8,   &
        0.998007_r8,  0.998584_r8,  0.998999_r8,  0.999293_r8,  0.999499_r8/
      data (bm0ij (  9,  9,ibeta), ibeta = 1,10) /   &
        0.994413_r8,  0.994903_r8,  0.995878_r8,  0.996876_r8,  0.997716_r8,   &
        0.998363_r8,  0.998839_r8,  0.999180_r8,  0.999421_r8,  0.999591_r8/
      data (bm0ij (  9, 10,ibeta), ibeta = 1,10) /   &
        0.995376_r8,  0.995785_r8,  0.996597_r8,  0.997425_r8,  0.998121_r8,   &
        0.998655_r8,  0.999048_r8,  0.999328_r8,  0.999526_r8,  0.999665_r8/
      data (bm0ij ( 10,  1,ibeta), ibeta = 1,10) /   &
        0.989082_r8,  0.989991_r8,  0.991819_r8,  0.993723_r8,  0.995357_r8,   &
        0.996637_r8,  0.997592_r8,  0.998286_r8,  0.998781_r8,  0.999132_r8/
      data (bm0ij ( 10,  2,ibeta), ibeta = 1,10) /   &
        0.989613_r8,  0.990480_r8,  0.992224_r8,  0.994039_r8,  0.995594_r8,   &
        0.996810_r8,  0.997717_r8,  0.998375_r8,  0.998845_r8,  0.999178_r8/
      data (bm0ij ( 10,  3,ibeta), ibeta = 1,10) /   &
        0.990744_r8,  0.991523_r8,  0.993086_r8,  0.994708_r8,  0.996094_r8,   &
        0.997176_r8,  0.997981_r8,  0.998564_r8,  0.998980_r8,  0.999274_r8/
      data (bm0ij ( 10,  4,ibeta), ibeta = 1,10) /   &
        0.992041_r8,  0.992716_r8,  0.994070_r8,  0.995470_r8,  0.996662_r8,   &
        0.997591_r8,  0.998280_r8,  0.998778_r8,  0.999133_r8,  0.999383_r8/
      data (bm0ij ( 10,  5,ibeta), ibeta = 1,10) /   &
        0.993292_r8,  0.993867_r8,  0.995015_r8,  0.996199_r8,  0.997205_r8,   &
        0.997985_r8,  0.998564_r8,  0.998981_r8,  0.999277_r8,  0.999487_r8/
      data (bm0ij ( 10,  6,ibeta), ibeta = 1,10) /   &
        0.994411_r8,  0.994894_r8,  0.995857_r8,  0.996847_r8,  0.997685_r8,   &
        0.998334_r8,  0.998814_r8,  0.999159_r8,  0.999404_r8,  0.999577_r8/
      data (bm0ij ( 10,  7,ibeta), ibeta = 1,10) /   &
        0.995373_r8,  0.995776_r8,  0.996577_r8,  0.997400_r8,  0.998094_r8,   &
        0.998630_r8,  0.999026_r8,  0.999310_r8,  0.999512_r8,  0.999654_r8/
      data (bm0ij ( 10,  8,ibeta), ibeta = 1,10) /   &
        0.996181_r8,  0.996516_r8,  0.997181_r8,  0.997861_r8,  0.998435_r8,   &
        0.998877_r8,  0.999202_r8,  0.999435_r8,  0.999601_r8,  0.999717_r8/
      data (bm0ij ( 10,  9,ibeta), ibeta = 1,10) /   &
        0.996851_r8,  0.997128_r8,  0.997680_r8,  0.998242_r8,  0.998715_r8,   &
        0.999079_r8,  0.999346_r8,  0.999538_r8,  0.999673_r8,  0.999769_r8/
      data (bm0ij ( 10, 10,ibeta), ibeta = 1,10) /   &
        0.997402_r8,  0.997632_r8,  0.998089_r8,  0.998554_r8,  0.998945_r8,   &
        0.999244_r8,  0.999464_r8,  0.999622_r8,  0.999733_r8,  0.999811_r8/


! rpm....   3rd moment nuclei mode corr. fac. for bimodal fm coag rate

       data (bm3i( 1, 1,ibeta ), ibeta=1,10)/   &
       0.70708_r8,0.71681_r8,0.73821_r8,0.76477_r8,0.79350_r8,0.82265_r8,0.85090_r8,0.87717_r8,   &
       0.90069_r8,0.92097_r8/
       data (bm3i( 1, 2,ibeta ), ibeta=1,10)/   &
       0.72172_r8,0.73022_r8,0.74927_r8,0.77324_r8,0.79936_r8,0.82601_r8,0.85199_r8,0.87637_r8,   &
       0.89843_r8,0.91774_r8/
       data (bm3i( 1, 3,ibeta ), ibeta=1,10)/   &
       0.78291_r8,0.78896_r8,0.80286_r8,0.82070_r8,0.84022_r8,0.85997_r8,0.87901_r8,0.89669_r8,   &
       0.91258_r8,0.92647_r8/
       data (bm3i( 1, 4,ibeta ), ibeta=1,10)/   &
       0.87760_r8,0.88147_r8,0.89025_r8,0.90127_r8,0.91291_r8,0.92420_r8,0.93452_r8,0.94355_r8,   &
       0.95113_r8,0.95726_r8/
       data (bm3i( 1, 5,ibeta ), ibeta=1,10)/   &
       0.94988_r8,0.95184_r8,0.95612_r8,0.96122_r8,0.96628_r8,0.97085_r8,0.97467_r8,0.97763_r8,   &
       0.97971_r8,0.98089_r8/
       data (bm3i( 1, 6,ibeta ), ibeta=1,10)/   &
       0.98318_r8,0.98393_r8,0.98551_r8,0.98728_r8,0.98889_r8,0.99014_r8,0.99095_r8,0.99124_r8,   &
       0.99100_r8,0.99020_r8/
       data (bm3i( 1, 7,ibeta ), ibeta=1,10)/   &
       0.99480_r8,0.99504_r8,0.99551_r8,0.99598_r8,0.99629_r8,0.99635_r8,0.99611_r8,0.99550_r8,   &
       0.99450_r8,0.99306_r8/
       data (bm3i( 1, 8,ibeta ), ibeta=1,10)/   &
       0.99842_r8,0.99848_r8,0.99858_r8,0.99861_r8,0.99850_r8,0.99819_r8,0.99762_r8,0.99674_r8,   &
       0.99550_r8,0.99388_r8/
       data (bm3i( 1, 9,ibeta ), ibeta=1,10)/   &
       0.99951_r8,0.99951_r8,0.99949_r8,0.99939_r8,0.99915_r8,0.99872_r8,0.99805_r8,0.99709_r8,   &
       0.99579_r8,0.99411_r8/
       data (bm3i( 1,10,ibeta ), ibeta=1,10)/   &
       0.99984_r8,0.99982_r8,0.99976_r8,0.99962_r8,0.99934_r8,0.99888_r8,0.99818_r8,0.99719_r8,   &
       0.99587_r8,0.99417_r8/
       data (bm3i( 2, 1,ibeta ), ibeta=1,10)/   &
       0.72957_r8,0.73993_r8,0.76303_r8,0.79178_r8,0.82245_r8,0.85270_r8,0.88085_r8,0.90578_r8,   &
       0.92691_r8,0.94415_r8/
       data (bm3i( 2, 2,ibeta ), ibeta=1,10)/   &
       0.72319_r8,0.73320_r8,0.75547_r8,0.78323_r8,0.81307_r8,0.84287_r8,0.87107_r8,0.89651_r8,   &
       0.91852_r8,0.93683_r8/
       data (bm3i( 2, 3,ibeta ), ibeta=1,10)/   &
       0.74413_r8,0.75205_r8,0.76998_r8,0.79269_r8,0.81746_r8,0.84258_r8,0.86685_r8,0.88938_r8,   &
       0.90953_r8,0.92695_r8/
       data (bm3i( 2, 4,ibeta ), ibeta=1,10)/   &
       0.82588_r8,0.83113_r8,0.84309_r8,0.85825_r8,0.87456_r8,0.89072_r8,0.90594_r8,0.91972_r8,   &
       0.93178_r8,0.94203_r8/
       data (bm3i( 2, 5,ibeta ), ibeta=1,10)/   &
       0.91886_r8,0.92179_r8,0.92831_r8,0.93624_r8,0.94434_r8,0.95192_r8,0.95856_r8,0.96409_r8,   &
       0.96845_r8,0.97164_r8/
       data (bm3i( 2, 6,ibeta ), ibeta=1,10)/   &
       0.97129_r8,0.97252_r8,0.97515_r8,0.97818_r8,0.98108_r8,0.98354_r8,0.98542_r8,0.98665_r8,   &
       0.98721_r8,0.98709_r8/
       data (bm3i( 2, 7,ibeta ), ibeta=1,10)/   &
       0.99104_r8,0.99145_r8,0.99230_r8,0.99320_r8,0.99394_r8,0.99439_r8,0.99448_r8,0.99416_r8,   &
       0.99340_r8,0.99217_r8/
       data (bm3i( 2, 8,ibeta ), ibeta=1,10)/   &
       0.99730_r8,0.99741_r8,0.99763_r8,0.99779_r8,0.99782_r8,0.99762_r8,0.99715_r8,0.99636_r8,   &
       0.99519_r8,0.99363_r8/
       data (bm3i( 2, 9,ibeta ), ibeta=1,10)/   &
       0.99917_r8,0.99919_r8,0.99921_r8,0.99915_r8,0.99895_r8,0.99856_r8,0.99792_r8,0.99698_r8,   &
       0.99570_r8,0.99404_r8/
       data (bm3i( 2,10,ibeta ), ibeta=1,10)/   &
       0.99973_r8,0.99973_r8,0.99968_r8,0.99955_r8,0.99928_r8,0.99883_r8,0.99814_r8,0.99716_r8,   &
       0.99584_r8,0.99415_r8/
       data (bm3i( 3, 1,ibeta ), ibeta=1,10)/   &
       0.78358_r8,0.79304_r8,0.81445_r8,0.84105_r8,0.86873_r8,0.89491_r8,0.91805_r8,0.93743_r8,   &
       0.95300_r8,0.96510_r8/
       data (bm3i( 3, 2,ibeta ), ibeta=1,10)/   &
       0.76412_r8,0.77404_r8,0.79635_r8,0.82404_r8,0.85312_r8,0.88101_r8,0.90610_r8,0.92751_r8,   &
       0.94500_r8,0.95879_r8/
       data (bm3i( 3, 3,ibeta ), ibeta=1,10)/   &
       0.74239_r8,0.75182_r8,0.77301_r8,0.79956_r8,0.82809_r8,0.85639_r8,0.88291_r8,0.90658_r8,   &
       0.92683_r8,0.94350_r8/
       data (bm3i( 3, 4,ibeta ), ibeta=1,10)/   &
       0.78072_r8,0.78758_r8,0.80317_r8,0.82293_r8,0.84437_r8,0.86589_r8,0.88643_r8,0.90526_r8,   &
       0.92194_r8,0.93625_r8/
       data (bm3i( 3, 5,ibeta ), ibeta=1,10)/   &
       0.87627_r8,0.88044_r8,0.88981_r8,0.90142_r8,0.91357_r8,0.92524_r8,0.93585_r8,0.94510_r8,   &
       0.95285_r8,0.95911_r8/
       data (bm3i( 3, 6,ibeta ), ibeta=1,10)/   &
       0.95176_r8,0.95371_r8,0.95796_r8,0.96297_r8,0.96792_r8,0.97233_r8,0.97599_r8,0.97880_r8,   &
       0.98072_r8,0.98178_r8/
       data (bm3i( 3, 7,ibeta ), ibeta=1,10)/   &
       0.98453_r8,0.98523_r8,0.98670_r8,0.98833_r8,0.98980_r8,0.99092_r8,0.99160_r8,0.99179_r8,   &
       0.99145_r8,0.99058_r8/
       data (bm3i( 3, 8,ibeta ), ibeta=1,10)/   &
       0.99534_r8,0.99555_r8,0.99597_r8,0.99637_r8,0.99662_r8,0.99663_r8,0.99633_r8,0.99569_r8,   &
       0.99465_r8,0.99318_r8/
       data (bm3i( 3, 9,ibeta ), ibeta=1,10)/   &
       0.99859_r8,0.99864_r8,0.99872_r8,0.99873_r8,0.99860_r8,0.99827_r8,0.99768_r8,0.99679_r8,   &
       0.99555_r8,0.99391_r8/
       data (bm3i( 3,10,ibeta ), ibeta=1,10)/   &
       0.99956_r8,0.99956_r8,0.99953_r8,0.99942_r8,0.99918_r8,0.99875_r8,0.99807_r8,0.99711_r8,   &
       0.99580_r8,0.99412_r8/
       data (bm3i( 4, 1,ibeta ), ibeta=1,10)/   &
       0.84432_r8,0.85223_r8,0.86990_r8,0.89131_r8,0.91280_r8,0.93223_r8,0.94861_r8,0.96172_r8,   &
       0.97185_r8,0.97945_r8/
       data (bm3i( 4, 2,ibeta ), ibeta=1,10)/   &
       0.82299_r8,0.83164_r8,0.85101_r8,0.87463_r8,0.89857_r8,0.92050_r8,0.93923_r8,0.95443_r8,   &
       0.96629_r8,0.97529_r8/
       data (bm3i( 4, 3,ibeta ), ibeta=1,10)/   &
       0.77870_r8,0.78840_r8,0.81011_r8,0.83690_r8,0.86477_r8,0.89124_r8,0.91476_r8,0.93460_r8,   &
       0.95063_r8,0.96316_r8/
       data (bm3i( 4, 4,ibeta ), ibeta=1,10)/   &
       0.76386_r8,0.77233_r8,0.79147_r8,0.81557_r8,0.84149_r8,0.86719_r8,0.89126_r8,0.91275_r8,   &
       0.93116_r8,0.94637_r8/
       data (bm3i( 4, 5,ibeta ), ibeta=1,10)/   &
       0.82927_r8,0.83488_r8,0.84756_r8,0.86346_r8,0.88040_r8,0.89704_r8,0.91257_r8,0.92649_r8,   &
       0.93857_r8,0.94874_r8/
       data (bm3i( 4, 6,ibeta ), ibeta=1,10)/   &
       0.92184_r8,0.92481_r8,0.93136_r8,0.93925_r8,0.94724_r8,0.95462_r8,0.96104_r8,0.96634_r8,   &
       0.97048_r8,0.97348_r8/
       data (bm3i( 4, 7,ibeta ), ibeta=1,10)/   &
       0.97341_r8,0.97457_r8,0.97706_r8,0.97991_r8,0.98260_r8,0.98485_r8,0.98654_r8,0.98760_r8,   &
       0.98801_r8,0.98777_r8/
       data (bm3i( 4, 8,ibeta ), ibeta=1,10)/   &
       0.99192_r8,0.99229_r8,0.99305_r8,0.99385_r8,0.99449_r8,0.99486_r8,0.99487_r8,0.99449_r8,   &
       0.99367_r8,0.99239_r8/
       data (bm3i( 4, 9,ibeta ), ibeta=1,10)/   &
       0.99758_r8,0.99768_r8,0.99787_r8,0.99800_r8,0.99799_r8,0.99777_r8,0.99727_r8,0.99645_r8,   &
       0.99527_r8,0.99369_r8/
       data (bm3i( 4,10,ibeta ), ibeta=1,10)/   &
       0.99926_r8,0.99928_r8,0.99928_r8,0.99921_r8,0.99900_r8,0.99860_r8,0.99795_r8,0.99701_r8,   &
       0.99572_r8,0.99405_r8/
       data (bm3i( 5, 1,ibeta ), ibeta=1,10)/   &
       0.89577_r8,0.90190_r8,0.91522_r8,0.93076_r8,0.94575_r8,0.95876_r8,0.96932_r8,0.97751_r8,   &
       0.98367_r8,0.98820_r8/
       data (bm3i( 5, 2,ibeta ), ibeta=1,10)/   &
       0.87860_r8,0.88547_r8,0.90052_r8,0.91828_r8,0.93557_r8,0.95075_r8,0.96319_r8,0.97292_r8,   &
       0.98028_r8,0.98572_r8/
       data (bm3i( 5, 3,ibeta ), ibeta=1,10)/   &
       0.83381_r8,0.84240_r8,0.86141_r8,0.88425_r8,0.90707_r8,0.92770_r8,0.94510_r8,0.95906_r8,   &
       0.96986_r8,0.97798_r8/
       data (bm3i( 5, 4,ibeta ), ibeta=1,10)/   &
       0.78530_r8,0.79463_r8,0.81550_r8,0.84127_r8,0.86813_r8,0.89367_r8,0.91642_r8,0.93566_r8,   &
       0.95125_r8,0.96347_r8/
       data (bm3i( 5, 5,ibeta ), ibeta=1,10)/   &
       0.79614_r8,0.80332_r8,0.81957_r8,0.84001_r8,0.86190_r8,0.88351_r8,0.90368_r8,0.92169_r8,   &
       0.93718_r8,0.95006_r8/
       data (bm3i( 5, 6,ibeta ), ibeta=1,10)/   &
       0.88192_r8,0.88617_r8,0.89565_r8,0.90728_r8,0.91931_r8,0.93076_r8,0.94107_r8,0.94997_r8,   &
       0.95739_r8,0.96333_r8/
       data (bm3i( 5, 7,ibeta ), ibeta=1,10)/   &
       0.95509_r8,0.95698_r8,0.96105_r8,0.96583_r8,0.97048_r8,0.97460_r8,0.97796_r8,0.98050_r8,   &
       0.98218_r8,0.98304_r8/
       data (bm3i( 5, 8,ibeta ), ibeta=1,10)/   &
       0.98596_r8,0.98660_r8,0.98794_r8,0.98943_r8,0.99074_r8,0.99172_r8,0.99227_r8,0.99235_r8,   &
       0.99192_r8,0.99096_r8/
       data (bm3i( 5, 9,ibeta ), ibeta=1,10)/   &
       0.99581_r8,0.99600_r8,0.99637_r8,0.99672_r8,0.99691_r8,0.99687_r8,0.99653_r8,0.99585_r8,   &
       0.99478_r8,0.99329_r8/
       data (bm3i( 5,10,ibeta ), ibeta=1,10)/   &
       0.99873_r8,0.99878_r8,0.99884_r8,0.99883_r8,0.99869_r8,0.99834_r8,0.99774_r8,0.99684_r8,   &
       0.99558_r8,0.99394_r8/
       data (bm3i( 6, 1,ibeta ), ibeta=1,10)/   &
       0.93335_r8,0.93777_r8,0.94711_r8,0.95764_r8,0.96741_r8,0.97562_r8,0.98210_r8,0.98701_r8,   &
       0.99064_r8,0.99327_r8/
       data (bm3i( 6, 2,ibeta ), ibeta=1,10)/   &
       0.92142_r8,0.92646_r8,0.93723_r8,0.94947_r8,0.96096_r8,0.97069_r8,0.97842_r8,0.98431_r8,   &
       0.98868_r8,0.99186_r8/
       data (bm3i( 6, 3,ibeta ), ibeta=1,10)/   &
       0.88678_r8,0.89351_r8,0.90810_r8,0.92508_r8,0.94138_r8,0.95549_r8,0.96693_r8,0.97578_r8,   &
       0.98243_r8,0.98731_r8/
       data (bm3i( 6, 4,ibeta ), ibeta=1,10)/   &
       0.83249_r8,0.84124_r8,0.86051_r8,0.88357_r8,0.90655_r8,0.92728_r8,0.94477_r8,0.95880_r8,   &
       0.96964_r8,0.97779_r8/
       data (bm3i( 6, 5,ibeta ), ibeta=1,10)/   &
       0.79593_r8,0.80444_r8,0.82355_r8,0.84725_r8,0.87211_r8,0.89593_r8,0.91735_r8,0.93566_r8,   &
       0.95066_r8,0.96255_r8/
       data (bm3i( 6, 6,ibeta ), ibeta=1,10)/   &
       0.84124_r8,0.84695_r8,0.85980_r8,0.87575_r8,0.89256_r8,0.90885_r8,0.92383_r8,0.93704_r8,   &
       0.94830_r8,0.95761_r8/
       data (bm3i( 6, 7,ibeta ), ibeta=1,10)/   &
       0.92721_r8,0.93011_r8,0.93647_r8,0.94406_r8,0.95166_r8,0.95862_r8,0.96460_r8,0.96949_r8,   &
       0.97326_r8,0.97595_r8/
       data (bm3i( 6, 8,ibeta ), ibeta=1,10)/   &
       0.97573_r8,0.97681_r8,0.97913_r8,0.98175_r8,0.98421_r8,0.98624_r8,0.98772_r8,0.98860_r8,   &
       0.98885_r8,0.98847_r8/
       data (bm3i( 6, 9,ibeta ), ibeta=1,10)/   &
       0.99271_r8,0.99304_r8,0.99373_r8,0.99444_r8,0.99499_r8,0.99528_r8,0.99522_r8,0.99477_r8,   &
       0.99390_r8,0.99258_r8/
       data (bm3i( 6,10,ibeta ), ibeta=1,10)/   &
       0.99782_r8,0.99791_r8,0.99807_r8,0.99817_r8,0.99813_r8,0.99788_r8,0.99737_r8,0.99653_r8,   &
       0.99533_r8,0.99374_r8/
       data (bm3i( 7, 1,ibeta ), ibeta=1,10)/   &
       0.95858_r8,0.96158_r8,0.96780_r8,0.97460_r8,0.98073_r8,0.98575_r8,0.98963_r8,0.99252_r8,   &
       0.99463_r8,0.99615_r8/
       data (bm3i( 7, 2,ibeta ), ibeta=1,10)/   &
       0.95091_r8,0.95438_r8,0.96163_r8,0.96962_r8,0.97688_r8,0.98286_r8,0.98751_r8,0.99099_r8,   &
       0.99353_r8,0.99536_r8/
       data (bm3i( 7, 3,ibeta ), ibeta=1,10)/   &
       0.92751_r8,0.93233_r8,0.94255_r8,0.95406_r8,0.96473_r8,0.97366_r8,0.98070_r8,0.98602_r8,   &
       0.98994_r8,0.99278_r8/
       data (bm3i( 7, 4,ibeta ), ibeta=1,10)/   &
       0.88371_r8,0.89075_r8,0.90595_r8,0.92351_r8,0.94028_r8,0.95474_r8,0.96642_r8,0.97544_r8,   &
       0.98220_r8,0.98715_r8/
       data (bm3i( 7, 5,ibeta ), ibeta=1,10)/   &
       0.82880_r8,0.83750_r8,0.85671_r8,0.87980_r8,0.90297_r8,0.92404_r8,0.94195_r8,0.95644_r8,   &
       0.96772_r8,0.97625_r8/
       data (bm3i( 7, 6,ibeta ), ibeta=1,10)/   &
       0.81933_r8,0.82655_r8,0.84279_r8,0.86295_r8,0.88412_r8,0.90449_r8,0.92295_r8,0.93890_r8,   &
       0.95215_r8,0.96281_r8/
       data (bm3i( 7, 7,ibeta ), ibeta=1,10)/   &
       0.89099_r8,0.89519_r8,0.90448_r8,0.91577_r8,0.92732_r8,0.93820_r8,0.94789_r8,0.95616_r8,   &
       0.96297_r8,0.96838_r8/
       data (bm3i( 7, 8,ibeta ), ibeta=1,10)/   &
       0.95886_r8,0.96064_r8,0.96448_r8,0.96894_r8,0.97324_r8,0.97701_r8,0.98004_r8,0.98228_r8,   &
       0.98371_r8,0.98435_r8/
       data (bm3i( 7, 9,ibeta ), ibeta=1,10)/   &
       0.98727_r8,0.98786_r8,0.98908_r8,0.99043_r8,0.99160_r8,0.99245_r8,0.99288_r8,0.99285_r8,   &
       0.99234_r8,0.99131_r8/
       data (bm3i( 7,10,ibeta ), ibeta=1,10)/   &
       0.99621_r8,0.99638_r8,0.99671_r8,0.99700_r8,0.99715_r8,0.99707_r8,0.99670_r8,0.99599_r8,   &
       0.99489_r8,0.99338_r8/
       data (bm3i( 8, 1,ibeta ), ibeta=1,10)/   &
       0.97470_r8,0.97666_r8,0.98064_r8,0.98491_r8,0.98867_r8,0.99169_r8,0.99399_r8,0.99569_r8,   &
       0.99691_r8,0.99779_r8/
       data (bm3i( 8, 2,ibeta ), ibeta=1,10)/   &
       0.96996_r8,0.97225_r8,0.97693_r8,0.98196_r8,0.98643_r8,0.99003_r8,0.99279_r8,0.99482_r8,   &
       0.99630_r8,0.99735_r8/
       data (bm3i( 8, 3,ibeta ), ibeta=1,10)/   &
       0.95523_r8,0.95848_r8,0.96522_r8,0.97260_r8,0.97925_r8,0.98468_r8,0.98888_r8,0.99200_r8,   &
       0.99427_r8,0.99590_r8/
       data (bm3i( 8, 4,ibeta ), ibeta=1,10)/   &
       0.92524_r8,0.93030_r8,0.94098_r8,0.95294_r8,0.96397_r8,0.97317_r8,0.98038_r8,0.98582_r8,   &
       0.98981_r8,0.99270_r8/
       data (bm3i( 8, 5,ibeta ), ibeta=1,10)/   &
       0.87576_r8,0.88323_r8,0.89935_r8,0.91799_r8,0.93583_r8,0.95126_r8,0.96377_r8,0.97345_r8,   &
       0.98072_r8,0.98606_r8/
       data (bm3i( 8, 6,ibeta ), ibeta=1,10)/   &
       0.83078_r8,0.83894_r8,0.85705_r8,0.87899_r8,0.90126_r8,0.92179_r8,0.93950_r8,0.95404_r8,   &
       0.96551_r8,0.97430_r8/
       data (bm3i( 8, 7,ibeta ), ibeta=1,10)/   &
       0.85727_r8,0.86294_r8,0.87558_r8,0.89111_r8,0.90723_r8,0.92260_r8,0.93645_r8,0.94841_r8,   &
       0.95838_r8,0.96643_r8/
       data (bm3i( 8, 8,ibeta ), ibeta=1,10)/   &
       0.93337_r8,0.93615_r8,0.94220_r8,0.94937_r8,0.95647_r8,0.96292_r8,0.96840_r8,0.97283_r8,   &
       0.97619_r8,0.97854_r8/
       data (bm3i( 8, 9,ibeta ), ibeta=1,10)/   &
       0.97790_r8,0.97891_r8,0.98105_r8,0.98346_r8,0.98569_r8,0.98751_r8,0.98879_r8,0.98950_r8,   &
       0.98961_r8,0.98912_r8/
       data (bm3i( 8,10,ibeta ), ibeta=1,10)/   &
       0.99337_r8,0.99367_r8,0.99430_r8,0.99493_r8,0.99541_r8,0.99562_r8,0.99551_r8,0.99501_r8,   &
       0.99410_r8,0.99274_r8/
       data (bm3i( 9, 1,ibeta ), ibeta=1,10)/   &
       0.98470_r8,0.98594_r8,0.98844_r8,0.99106_r8,0.99334_r8,0.99514_r8,0.99650_r8,0.99749_r8,   &
       0.99821_r8,0.99872_r8/
       data (bm3i( 9, 2,ibeta ), ibeta=1,10)/   &
       0.98184_r8,0.98330_r8,0.98624_r8,0.98934_r8,0.99205_r8,0.99420_r8,0.99582_r8,0.99701_r8,   &
       0.99787_r8,0.99848_r8/
       data (bm3i( 9, 3,ibeta ), ibeta=1,10)/   &
       0.97288_r8,0.97498_r8,0.97927_r8,0.98385_r8,0.98789_r8,0.99113_r8,0.99360_r8,0.99541_r8,   &
       0.99673_r8,0.99766_r8/
       data (bm3i( 9, 4,ibeta ), ibeta=1,10)/   &
       0.95403_r8,0.95741_r8,0.96440_r8,0.97202_r8,0.97887_r8,0.98444_r8,0.98872_r8,0.99190_r8,   &
       0.99421_r8,0.99586_r8/
       data (bm3i( 9, 5,ibeta ), ibeta=1,10)/   &
       0.91845_r8,0.92399_r8,0.93567_r8,0.94873_r8,0.96076_r8,0.97079_r8,0.97865_r8,0.98457_r8,   &
       0.98892_r8,0.99206_r8/
       data (bm3i( 9, 6,ibeta ), ibeta=1,10)/   &
       0.86762_r8,0.87533_r8,0.89202_r8,0.91148_r8,0.93027_r8,0.94669_r8,0.96013_r8,0.97062_r8,   &
       0.97855_r8,0.98441_r8/
       data (bm3i( 9, 7,ibeta ), ibeta=1,10)/   &
       0.84550_r8,0.85253_r8,0.86816_r8,0.88721_r8,0.90671_r8,0.92490_r8,0.94083_r8,0.95413_r8,   &
       0.96481_r8,0.97314_r8/
       data (bm3i( 9, 8,ibeta ), ibeta=1,10)/   &
       0.90138_r8,0.90544_r8,0.91437_r8,0.92513_r8,0.93602_r8,0.94615_r8,0.95506_r8,0.96258_r8,   &
       0.96868_r8,0.97347_r8/
       data (bm3i( 9, 9,ibeta ), ibeta=1,10)/   &
       0.96248_r8,0.96415_r8,0.96773_r8,0.97187_r8,0.97583_r8,0.97925_r8,0.98198_r8,0.98394_r8,   &
       0.98514_r8,0.98559_r8/
       data (bm3i( 9,10,ibeta ), ibeta=1,10)/   &
       0.98837_r8,0.98892_r8,0.99005_r8,0.99127_r8,0.99232_r8,0.99306_r8,0.99339_r8,0.99328_r8,   &
       0.99269_r8,0.99161_r8/
       data (bm3i(10, 1,ibeta ), ibeta=1,10)/   &
       0.99080_r8,0.99158_r8,0.99311_r8,0.99471_r8,0.99607_r8,0.99715_r8,0.99795_r8,0.99853_r8,   &
       0.99895_r8,0.99925_r8/
       data (bm3i(10, 2,ibeta ), ibeta=1,10)/   &
       0.98910_r8,0.99001_r8,0.99182_r8,0.99371_r8,0.99533_r8,0.99661_r8,0.99757_r8,0.99826_r8,   &
       0.99876_r8,0.99912_r8/
       data (bm3i(10, 3,ibeta ), ibeta=1,10)/   &
       0.98374_r8,0.98506_r8,0.98772_r8,0.99051_r8,0.99294_r8,0.99486_r8,0.99630_r8,0.99736_r8,   &
       0.99812_r8,0.99866_r8/
       data (bm3i(10, 4,ibeta ), ibeta=1,10)/   &
       0.97238_r8,0.97453_r8,0.97892_r8,0.98361_r8,0.98773_r8,0.99104_r8,0.99354_r8,0.99538_r8,   &
       0.99671_r8,0.99765_r8/
       data (bm3i(10, 5,ibeta ), ibeta=1,10)/   &
       0.94961_r8,0.95333_r8,0.96103_r8,0.96941_r8,0.97693_r8,0.98303_r8,0.98772_r8,0.99119_r8,   &
       0.99371_r8,0.99551_r8/
       data (bm3i(10, 6,ibeta ), ibeta=1,10)/   &
       0.90943_r8,0.91550_r8,0.92834_r8,0.94275_r8,0.95608_r8,0.96723_r8,0.97600_r8,0.98263_r8,   &
       0.98751_r8,0.99103_r8/
       data (bm3i(10, 7,ibeta ), ibeta=1,10)/   &
       0.86454_r8,0.87200_r8,0.88829_r8,0.90749_r8,0.92630_r8,0.94300_r8,0.95687_r8,0.96785_r8,   &
       0.97626_r8,0.98254_r8/
       data (bm3i(10, 8,ibeta ), ibeta=1,10)/   &
       0.87498_r8,0.88048_r8,0.89264_r8,0.90737_r8,0.92240_r8,0.93642_r8,0.94877_r8,0.95917_r8,   &
       0.96762_r8,0.97429_r8/
       data (bm3i(10, 9,ibeta ), ibeta=1,10)/   &
       0.93946_r8,0.94209_r8,0.94781_r8,0.95452_r8,0.96111_r8,0.96704_r8,0.97203_r8,0.97602_r8,   &
       0.97900_r8,0.98106_r8/
       data (bm3i(10,10,ibeta ), ibeta=1,10)/   &
       0.97977_r8,0.98071_r8,0.98270_r8,0.98492_r8,0.98695_r8,0.98858_r8,0.98970_r8,0.99027_r8,   &
       0.99026_r8,0.98968_r8/

! fsb fm correction for intramodal m2 coagulation
       data bm2ii /   &
        0.707107_r8,  0.720583_r8,  0.745310_r8,  0.748056_r8,  0.696935_r8,   &
        0.604164_r8,  0.504622_r8,  0.416559_r8,  0.343394_r8,  0.283641_r8/

! *** total correction for intramodal m2 coagulation

      data bm2iitt /   &
        1.000000_r8,  0.907452_r8,  0.680931_r8,  0.409815_r8,  0.196425_r8,   &
        0.078814_r8,  0.028473_r8,  0.009800_r8,  0.003322_r8,  0.001129_r8/


! fsb fm correction for m2 i to j coagulation

      data (bm2ij (  1,  1,ibeta), ibeta = 1,10) /   &
        0.707107_r8,  0.716828_r8,  0.738240_r8,  0.764827_r8,  0.793610_r8,   &
        0.822843_r8,  0.851217_r8,  0.877670_r8,  0.901404_r8,  0.921944_r8/
      data (bm2ij (  1,  2,ibeta), ibeta = 1,10) /   &
        0.719180_r8,  0.727975_r8,  0.747638_r8,  0.772334_r8,  0.799234_r8,   &
        0.826666_r8,  0.853406_r8,  0.878482_r8,  0.901162_r8,  0.920987_r8/
      data (bm2ij (  1,  3,ibeta), ibeta = 1,10) /   &
        0.760947_r8,  0.767874_r8,  0.783692_r8,  0.803890_r8,  0.826015_r8,   &
        0.848562_r8,  0.870498_r8,  0.891088_r8,  0.909823_r8,  0.926400_r8/
      data (bm2ij (  1,  4,ibeta), ibeta = 1,10) /   &
        0.830926_r8,  0.836034_r8,  0.847708_r8,  0.862528_r8,  0.878521_r8,   &
        0.894467_r8,  0.909615_r8,  0.923520_r8,  0.935959_r8,  0.946858_r8/
      data (bm2ij (  1,  5,ibeta), ibeta = 1,10) /   &
        0.903643_r8,  0.907035_r8,  0.914641_r8,  0.924017_r8,  0.933795_r8,   &
        0.943194_r8,  0.951806_r8,  0.959449_r8,  0.966087_r8,  0.971761_r8/
      data (bm2ij (  1,  6,ibeta), ibeta = 1,10) /   &
        0.954216_r8,  0.956094_r8,  0.960211_r8,  0.965123_r8,  0.970068_r8,   &
        0.974666_r8,  0.978750_r8,  0.982277_r8,  0.985268_r8,  0.987775_r8/
      data (bm2ij (  1,  7,ibeta), ibeta = 1,10) /   &
        0.980546_r8,  0.981433_r8,  0.983343_r8,  0.985568_r8,  0.987751_r8,   &
        0.989735_r8,  0.991461_r8,  0.992926_r8,  0.994150_r8,  0.995164_r8/
      data (bm2ij (  1,  8,ibeta), ibeta = 1,10) /   &
        0.992142_r8,  0.992524_r8,  0.993338_r8,  0.994272_r8,  0.995174_r8,   &
        0.995981_r8,  0.996675_r8,  0.997257_r8,  0.997740_r8,  0.998137_r8/
      data (bm2ij (  1,  9,ibeta), ibeta = 1,10) /   &
        0.996868_r8,  0.997026_r8,  0.997361_r8,  0.997742_r8,  0.998106_r8,   &
        0.998430_r8,  0.998705_r8,  0.998935_r8,  0.999125_r8,  0.999280_r8/
      data (bm2ij (  1, 10,ibeta), ibeta = 1,10) /   &
        0.998737_r8,  0.998802_r8,  0.998939_r8,  0.999094_r8,  0.999241_r8,   &
        0.999371_r8,  0.999481_r8,  0.999573_r8,  0.999648_r8,  0.999709_r8/
      data (bm2ij (  2,  1,ibeta), ibeta = 1,10) /   &
        0.729600_r8,  0.739948_r8,  0.763059_r8,  0.791817_r8,  0.822510_r8,   &
        0.852795_r8,  0.881000_r8,  0.905999_r8,  0.927206_r8,  0.944532_r8/
      data (bm2ij (  2,  2,ibeta), ibeta = 1,10) /   &
        0.727025_r8,  0.737116_r8,  0.759615_r8,  0.787657_r8,  0.817740_r8,   &
        0.847656_r8,  0.875801_r8,  0.901038_r8,  0.922715_r8,  0.940643_r8/
      data (bm2ij (  2,  3,ibeta), ibeta = 1,10) /   &
        0.738035_r8,  0.746779_r8,  0.766484_r8,  0.791340_r8,  0.818324_r8,   &
        0.845546_r8,  0.871629_r8,  0.895554_r8,  0.916649_r8,  0.934597_r8/
      data (bm2ij (  2,  4,ibeta), ibeta = 1,10) /   &
        0.784185_r8,  0.790883_r8,  0.806132_r8,  0.825501_r8,  0.846545_r8,   &
        0.867745_r8,  0.888085_r8,  0.906881_r8,  0.923705_r8,  0.938349_r8/
      data (bm2ij (  2,  5,ibeta), ibeta = 1,10) /   &
        0.857879_r8,  0.862591_r8,  0.873238_r8,  0.886539_r8,  0.900645_r8,   &
        0.914463_r8,  0.927360_r8,  0.939004_r8,  0.949261_r8,  0.958125_r8/
      data (bm2ij (  2,  6,ibeta), ibeta = 1,10) /   &
        0.925441_r8,  0.928304_r8,  0.934645_r8,  0.942324_r8,  0.950181_r8,   &
        0.957600_r8,  0.964285_r8,  0.970133_r8,  0.975147_r8,  0.979388_r8/
      data (bm2ij (  2,  7,ibeta), ibeta = 1,10) /   &
        0.966728_r8,  0.968176_r8,  0.971323_r8,  0.975027_r8,  0.978705_r8,   &
        0.982080_r8,  0.985044_r8,  0.987578_r8,  0.989710_r8,  0.991485_r8/
      data (bm2ij (  2,  8,ibeta), ibeta = 1,10) /   &
        0.986335_r8,  0.986980_r8,  0.988362_r8,  0.989958_r8,  0.991511_r8,   &
        0.992912_r8,  0.994122_r8,  0.995143_r8,  0.995992_r8,  0.996693_r8/
      data (bm2ij (  2,  9,ibeta), ibeta = 1,10) /   &
        0.994547_r8,  0.994817_r8,  0.995391_r8,  0.996046_r8,  0.996677_r8,   &
        0.997238_r8,  0.997719_r8,  0.998122_r8,  0.998454_r8,  0.998727_r8/
      data (bm2ij (  2, 10,ibeta), ibeta = 1,10) /   &
        0.997817_r8,  0.997928_r8,  0.998163_r8,  0.998429_r8,  0.998683_r8,   &
        0.998908_r8,  0.999099_r8,  0.999258_r8,  0.999389_r8,  0.999497_r8/
      data (bm2ij (  3,  1,ibeta), ibeta = 1,10) /   &
        0.783612_r8,  0.793055_r8,  0.814468_r8,  0.841073_r8,  0.868769_r8,   &
        0.894963_r8,  0.918118_r8,  0.937527_r8,  0.953121_r8,  0.965244_r8/
      data (bm2ij (  3,  2,ibeta), ibeta = 1,10) /   &
        0.772083_r8,  0.781870_r8,  0.803911_r8,  0.831238_r8,  0.859802_r8,   &
        0.887036_r8,  0.911349_r8,  0.931941_r8,  0.948649_r8,  0.961751_r8/
      data (bm2ij (  3,  3,ibeta), ibeta = 1,10) /   &
        0.755766_r8,  0.765509_r8,  0.787380_r8,  0.814630_r8,  0.843526_r8,   &
        0.871670_r8,  0.897443_r8,  0.919870_r8,  0.938557_r8,  0.953576_r8/
      data (bm2ij (  3,  4,ibeta), ibeta = 1,10) /   &
        0.763816_r8,  0.772145_r8,  0.790997_r8,  0.814784_r8,  0.840434_r8,   &
        0.865978_r8,  0.890034_r8,  0.911671_r8,  0.930366_r8,  0.945963_r8/
      data (bm2ij (  3,  5,ibeta), ibeta = 1,10) /   &
        0.813597_r8,  0.819809_r8,  0.833889_r8,  0.851618_r8,  0.870640_r8,   &
        0.889514_r8,  0.907326_r8,  0.923510_r8,  0.937768_r8,  0.950003_r8/
      data (bm2ij (  3,  6,ibeta), ibeta = 1,10) /   &
        0.886317_r8,  0.890437_r8,  0.899643_r8,  0.910955_r8,  0.922730_r8,   &
        0.934048_r8,  0.944422_r8,  0.953632_r8,  0.961624_r8,  0.968444_r8/
      data (bm2ij (  3,  7,ibeta), ibeta = 1,10) /   &
        0.944565_r8,  0.946855_r8,  0.951872_r8,  0.957854_r8,  0.963873_r8,   &
        0.969468_r8,  0.974438_r8,  0.978731_r8,  0.982372_r8,  0.985424_r8/
      data (bm2ij (  3,  8,ibeta), ibeta = 1,10) /   &
        0.976358_r8,  0.977435_r8,  0.979759_r8,  0.982467_r8,  0.985125_r8,   &
        0.987540_r8,  0.989642_r8,  0.991425_r8,  0.992916_r8,  0.994150_r8/
      data (bm2ij (  3,  9,ibeta), ibeta = 1,10) /   &
        0.990471_r8,  0.990932_r8,  0.991917_r8,  0.993048_r8,  0.994142_r8,   &
        0.995121_r8,  0.995964_r8,  0.996671_r8,  0.997258_r8,  0.997740_r8/
      data (bm2ij (  3, 10,ibeta), ibeta = 1,10) /   &
        0.996199_r8,  0.996389_r8,  0.996794_r8,  0.997254_r8,  0.997694_r8,   &
        0.998086_r8,  0.998420_r8,  0.998699_r8,  0.998929_r8,  0.999117_r8/
      data (bm2ij (  4,  1,ibeta), ibeta = 1,10) /   &
        0.844355_r8,  0.852251_r8,  0.869914_r8,  0.891330_r8,  0.912823_r8,   &
        0.932259_r8,  0.948642_r8,  0.961767_r8,  0.971897_r8,  0.979510_r8/
      data (bm2ij (  4,  2,ibeta), ibeta = 1,10) /   &
        0.831550_r8,  0.839954_r8,  0.858754_r8,  0.881583_r8,  0.904592_r8,   &
        0.925533_r8,  0.943309_r8,  0.957647_r8,  0.968779_r8,  0.977185_r8/
      data (bm2ij (  4,  3,ibeta), ibeta = 1,10) /   &
        0.803981_r8,  0.813288_r8,  0.834060_r8,  0.859400_r8,  0.885285_r8,   &
        0.909286_r8,  0.930084_r8,  0.947193_r8,  0.960714_r8,  0.971078_r8/
      data (bm2ij (  4,  4,ibeta), ibeta = 1,10) /   &
        0.781787_r8,  0.791080_r8,  0.811931_r8,  0.837749_r8,  0.864768_r8,   &
        0.890603_r8,  0.913761_r8,  0.933477_r8,  0.949567_r8,  0.962261_r8/
      data (bm2ij (  4,  5,ibeta), ibeta = 1,10) /   &
        0.791591_r8,  0.799355_r8,  0.816916_r8,  0.838961_r8,  0.862492_r8,   &
        0.885595_r8,  0.907003_r8,  0.925942_r8,  0.942052_r8,  0.955310_r8/
      data (bm2ij (  4,  6,ibeta), ibeta = 1,10) /   &
        0.844933_r8,  0.850499_r8,  0.863022_r8,  0.878593_r8,  0.895038_r8,   &
        0.911072_r8,  0.925939_r8,  0.939227_r8,  0.950765_r8,  0.960550_r8/
      data (bm2ij (  4,  7,ibeta), ibeta = 1,10) /   &
        0.912591_r8,  0.916022_r8,  0.923607_r8,  0.932777_r8,  0.942151_r8,   &
        0.951001_r8,  0.958976_r8,  0.965950_r8,  0.971924_r8,  0.976965_r8/
      data (bm2ij (  4,  8,ibeta), ibeta = 1,10) /   &
        0.959859_r8,  0.961617_r8,  0.965433_r8,  0.969924_r8,  0.974382_r8,   &
        0.978472_r8,  0.982063_r8,  0.985134_r8,  0.987716_r8,  0.989865_r8/
      data (bm2ij (  4,  9,ibeta), ibeta = 1,10) /   &
        0.983377_r8,  0.984162_r8,  0.985844_r8,  0.987788_r8,  0.989681_r8,   &
        0.991386_r8,  0.992860_r8,  0.994104_r8,  0.995139_r8,  0.995991_r8/
      data (bm2ij (  4, 10,ibeta), ibeta = 1,10) /   &
        0.993343_r8,  0.993672_r8,  0.994370_r8,  0.995169_r8,  0.995937_r8,   &
        0.996622_r8,  0.997209_r8,  0.997700_r8,  0.998106_r8,  0.998439_r8/
      data (bm2ij (  5,  1,ibeta), ibeta = 1,10) /   &
        0.895806_r8,  0.901918_r8,  0.915233_r8,  0.930783_r8,  0.945768_r8,   &
        0.958781_r8,  0.969347_r8,  0.977540_r8,  0.983697_r8,  0.988225_r8/
      data (bm2ij (  5,  2,ibeta), ibeta = 1,10) /   &
        0.885634_r8,  0.892221_r8,  0.906629_r8,  0.923540_r8,  0.939918_r8,   &
        0.954213_r8,  0.965873_r8,  0.974951_r8,  0.981794_r8,  0.986840_r8/
      data (bm2ij (  5,  3,ibeta), ibeta = 1,10) /   &
        0.860120_r8,  0.867858_r8,  0.884865_r8,  0.904996_r8,  0.924724_r8,   &
        0.942177_r8,  0.956602_r8,  0.967966_r8,  0.976616_r8,  0.983043_r8/
      data (bm2ij (  5,  4,ibeta), ibeta = 1,10) /   &
        0.827462_r8,  0.836317_r8,  0.855885_r8,  0.879377_r8,  0.902897_r8,   &
        0.924232_r8,  0.942318_r8,  0.956900_r8,  0.968222_r8,  0.976774_r8/
      data (bm2ij (  5,  5,ibeta), ibeta = 1,10) /   &
        0.805527_r8,  0.814279_r8,  0.833853_r8,  0.857892_r8,  0.882726_r8,   &
        0.906095_r8,  0.926690_r8,  0.943938_r8,  0.957808_r8,  0.968615_r8/
      data (bm2ij (  5,  6,ibeta), ibeta = 1,10) /   &
        0.820143_r8,  0.827223_r8,  0.843166_r8,  0.863002_r8,  0.883905_r8,   &
        0.904128_r8,  0.922585_r8,  0.938687_r8,  0.952222_r8,  0.963255_r8/
      data (bm2ij (  5,  7,ibeta), ibeta = 1,10) /   &
        0.875399_r8,  0.880208_r8,  0.890929_r8,  0.904065_r8,  0.917699_r8,   &
        0.930756_r8,  0.942656_r8,  0.953131_r8,  0.962113_r8,  0.969657_r8/
      data (bm2ij (  5,  8,ibeta), ibeta = 1,10) /   &
        0.934782_r8,  0.937520_r8,  0.943515_r8,  0.950656_r8,  0.957840_r8,   &
        0.964516_r8,  0.970446_r8,  0.975566_r8,  0.979905_r8,  0.983534_r8/
      data (bm2ij (  5,  9,ibeta), ibeta = 1,10) /   &
        0.971369_r8,  0.972679_r8,  0.975505_r8,  0.978797_r8,  0.982029_r8,   &
        0.984964_r8,  0.987518_r8,  0.989685_r8,  0.991496_r8,  0.992994_r8/
      data (bm2ij (  5, 10,ibeta), ibeta = 1,10) /   &
        0.988329_r8,  0.988893_r8,  0.990099_r8,  0.991485_r8,  0.992825_r8,   &
        0.994025_r8,  0.995058_r8,  0.995925_r8,  0.996643_r8,  0.997234_r8/
      data (bm2ij (  6,  1,ibeta), ibeta = 1,10) /   &
        0.933384_r8,  0.937784_r8,  0.947130_r8,  0.957655_r8,  0.967430_r8,   &
        0.975639_r8,  0.982119_r8,  0.987031_r8,  0.990657_r8,  0.993288_r8/
      data (bm2ij (  6,  2,ibeta), ibeta = 1,10) /   &
        0.926445_r8,  0.931227_r8,  0.941426_r8,  0.952975_r8,  0.963754_r8,   &
        0.972845_r8,  0.980044_r8,  0.985514_r8,  0.989558_r8,  0.992498_r8/
      data (bm2ij (  6,  3,ibeta), ibeta = 1,10) /   &
        0.907835_r8,  0.913621_r8,  0.926064_r8,  0.940308_r8,  0.953745_r8,   &
        0.965189_r8,  0.974327_r8,  0.981316_r8,  0.986510_r8,  0.990297_r8/
      data (bm2ij (  6,  4,ibeta), ibeta = 1,10) /   &
        0.879088_r8,  0.886306_r8,  0.901945_r8,  0.920079_r8,  0.937460_r8,   &
        0.952509_r8,  0.964711_r8,  0.974166_r8,  0.981265_r8,  0.986484_r8/
      data (bm2ij (  6,  5,ibeta), ibeta = 1,10) /   &
        0.846500_r8,  0.854862_r8,  0.873189_r8,  0.894891_r8,  0.916264_r8,   &
        0.935315_r8,  0.951197_r8,  0.963812_r8,  0.973484_r8,  0.980715_r8/
      data (bm2ij (  6,  6,ibeta), ibeta = 1,10) /   &
        0.828137_r8,  0.836250_r8,  0.854310_r8,  0.876287_r8,  0.898710_r8,   &
        0.919518_r8,  0.937603_r8,  0.952560_r8,  0.964461_r8,  0.973656_r8/
      data (bm2ij (  6,  7,ibeta), ibeta = 1,10) /   &
        0.848595_r8,  0.854886_r8,  0.868957_r8,  0.886262_r8,  0.904241_r8,   &
        0.921376_r8,  0.936799_r8,  0.950096_r8,  0.961172_r8,  0.970145_r8/
      data (bm2ij (  6,  8,ibeta), ibeta = 1,10) /   &
        0.902919_r8,  0.906922_r8,  0.915760_r8,  0.926427_r8,  0.937312_r8,   &
        0.947561_r8,  0.956758_r8,  0.964747_r8,  0.971525_r8,  0.977175_r8/
      data (bm2ij (  6,  9,ibeta), ibeta = 1,10) /   &
        0.952320_r8,  0.954434_r8,  0.959021_r8,  0.964418_r8,  0.969774_r8,   &
        0.974688_r8,  0.979003_r8,  0.982690_r8,  0.985789_r8,  0.988364_r8/
      data (bm2ij (  6, 10,ibeta), ibeta = 1,10) /   &
        0.979689_r8,  0.980650_r8,  0.982712_r8,  0.985093_r8,  0.987413_r8,   &
        0.989502_r8,  0.991308_r8,  0.992831_r8,  0.994098_r8,  0.995142_r8/
      data (bm2ij (  7,  1,ibeta), ibeta = 1,10) /   &
        0.958611_r8,  0.961598_r8,  0.967817_r8,  0.974620_r8,  0.980752_r8,   &
        0.985771_r8,  0.989650_r8,  0.992543_r8,  0.994653_r8,  0.996171_r8/
      data (bm2ij (  7,  2,ibeta), ibeta = 1,10) /   &
        0.954225_r8,  0.957488_r8,  0.964305_r8,  0.971795_r8,  0.978576_r8,   &
        0.984144_r8,  0.988458_r8,  0.991681_r8,  0.994034_r8,  0.995728_r8/
      data (bm2ij (  7,  3,ibeta), ibeta = 1,10) /   &
        0.942147_r8,  0.946158_r8,  0.954599_r8,  0.963967_r8,  0.972529_r8,   &
        0.979612_r8,  0.985131_r8,  0.989271_r8,  0.992301_r8,  0.994487_r8/
      data (bm2ij (  7,  4,ibeta), ibeta = 1,10) /   &
        0.921821_r8,  0.927048_r8,  0.938140_r8,  0.950598_r8,  0.962118_r8,   &
        0.971752_r8,  0.979326_r8,  0.985046_r8,  0.989254_r8,  0.992299_r8/
      data (bm2ij (  7,  5,ibeta), ibeta = 1,10) /   &
        0.893419_r8,  0.900158_r8,  0.914598_r8,  0.931070_r8,  0.946584_r8,   &
        0.959795_r8,  0.970350_r8,  0.978427_r8,  0.984432_r8,  0.988811_r8/
      data (bm2ij (  7,  6,ibeta), ibeta = 1,10) /   &
        0.863302_r8,  0.871111_r8,  0.888103_r8,  0.907990_r8,  0.927305_r8,   &
        0.944279_r8,  0.958245_r8,  0.969211_r8,  0.977540_r8,  0.983720_r8/
      data (bm2ij (  7,  7,ibeta), ibeta = 1,10) /   &
        0.850182_r8,  0.857560_r8,  0.873890_r8,  0.893568_r8,  0.913408_r8,   &
        0.931591_r8,  0.947216_r8,  0.960014_r8,  0.970121_r8,  0.977886_r8/
      data (bm2ij (  7,  8,ibeta), ibeta = 1,10) /   &
        0.875837_r8,  0.881265_r8,  0.893310_r8,  0.907936_r8,  0.922910_r8,   &
        0.936977_r8,  0.949480_r8,  0.960154_r8,  0.968985_r8,  0.976111_r8/
      data (bm2ij (  7,  9,ibeta), ibeta = 1,10) /   &
        0.926228_r8,  0.929445_r8,  0.936486_r8,  0.944868_r8,  0.953293_r8,   &
        0.961108_r8,  0.968028_r8,  0.973973_r8,  0.978974_r8,  0.983118_r8/
      data (bm2ij (  7, 10,ibeta), ibeta = 1,10) /   &
        0.965533_r8,  0.967125_r8,  0.970558_r8,  0.974557_r8,  0.978484_r8,   &
        0.982050_r8,  0.985153_r8,  0.987785_r8,  0.989982_r8,  0.991798_r8/
      data (bm2ij (  8,  1,ibeta), ibeta = 1,10) /   &
        0.974731_r8,  0.976674_r8,  0.980660_r8,  0.984926_r8,  0.988689_r8,   &
        0.991710_r8,  0.994009_r8,  0.995703_r8,  0.996929_r8,  0.997805_r8/
      data (bm2ij (  8,  2,ibeta), ibeta = 1,10) /   &
        0.972062_r8,  0.974192_r8,  0.978571_r8,  0.983273_r8,  0.987432_r8,   &
        0.990780_r8,  0.993333_r8,  0.995218_r8,  0.996581_r8,  0.997557_r8/
      data (bm2ij (  8,  3,ibeta), ibeta = 1,10) /   &
        0.964662_r8,  0.967300_r8,  0.972755_r8,  0.978659_r8,  0.983921_r8,   &
        0.988181_r8,  0.991444_r8,  0.993859_r8,  0.995610_r8,  0.996863_r8/
      data (bm2ij (  8,  4,ibeta), ibeta = 1,10) /   &
        0.951782_r8,  0.955284_r8,  0.962581_r8,  0.970559_r8,  0.977737_r8,   &
        0.983593_r8,  0.988103_r8,  0.991454_r8,  0.993889_r8,  0.995635_r8/
      data (bm2ij (  8,  5,ibeta), ibeta = 1,10) /   &
        0.931947_r8,  0.936723_r8,  0.946751_r8,  0.957843_r8,  0.967942_r8,   &
        0.976267_r8,  0.982734_r8,  0.987571_r8,  0.991102_r8,  0.993642_r8/
      data (bm2ij (  8,  6,ibeta), ibeta = 1,10) /   &
        0.905410_r8,  0.911665_r8,  0.924950_r8,  0.939908_r8,  0.953798_r8,   &
        0.965469_r8,  0.974684_r8,  0.981669_r8,  0.986821_r8,  0.990556_r8/
      data (bm2ij (  8,  7,ibeta), ibeta = 1,10) /   &
        0.878941_r8,  0.886132_r8,  0.901679_r8,  0.919688_r8,  0.936970_r8,   &
        0.951980_r8,  0.964199_r8,  0.973709_r8,  0.980881_r8,  0.986174_r8/
      data (bm2ij (  8,  8,ibeta), ibeta = 1,10) /   &
        0.871653_r8,  0.878218_r8,  0.892652_r8,  0.909871_r8,  0.927034_r8,   &
        0.942592_r8,  0.955836_r8,  0.966604_r8,  0.975065_r8,  0.981545_r8/
      data (bm2ij (  8,  9,ibeta), ibeta = 1,10) /   &
        0.900693_r8,  0.905239_r8,  0.915242_r8,  0.927232_r8,  0.939335_r8,   &
        0.950555_r8,  0.960420_r8,  0.968774_r8,  0.975651_r8,  0.981188_r8/
      data (bm2ij (  8, 10,ibeta), ibeta = 1,10) /   &
        0.944922_r8,  0.947435_r8,  0.952894_r8,  0.959317_r8,  0.965689_r8,   &
        0.971529_r8,  0.976645_r8,  0.981001_r8,  0.984641_r8,  0.987642_r8/
      data (bm2ij (  9,  1,ibeta), ibeta = 1,10) /   &
        0.984736_r8,  0.985963_r8,  0.988453_r8,  0.991078_r8,  0.993357_r8,   &
        0.995161_r8,  0.996519_r8,  0.997512_r8,  0.998226_r8,  0.998734_r8/
      data (bm2ij (  9,  2,ibeta), ibeta = 1,10) /   &
        0.983141_r8,  0.984488_r8,  0.987227_r8,  0.990119_r8,  0.992636_r8,   &
        0.994632_r8,  0.996137_r8,  0.997238_r8,  0.998030_r8,  0.998595_r8/
      data (bm2ij (  9,  3,ibeta), ibeta = 1,10) /   &
        0.978726_r8,  0.980401_r8,  0.983819_r8,  0.987450_r8,  0.990626_r8,   &
        0.993157_r8,  0.995071_r8,  0.996475_r8,  0.997486_r8,  0.998206_r8/
      data (bm2ij (  9,  4,ibeta), ibeta = 1,10) /   &
        0.970986_r8,  0.973224_r8,  0.977818_r8,  0.982737_r8,  0.987072_r8,   &
        0.990546_r8,  0.993184_r8,  0.995124_r8,  0.996523_r8,  0.997521_r8/
      data (bm2ij (  9,  5,ibeta), ibeta = 1,10) /   &
        0.958579_r8,  0.961700_r8,  0.968149_r8,  0.975116_r8,  0.981307_r8,   &
        0.986301_r8,  0.990112_r8,  0.992923_r8,  0.994954_r8,  0.996404_r8/
      data (bm2ij (  9,  6,ibeta), ibeta = 1,10) /   &
        0.940111_r8,  0.944479_r8,  0.953572_r8,  0.963506_r8,  0.972436_r8,   &
        0.979714_r8,  0.985313_r8,  0.989468_r8,  0.992483_r8,  0.994641_r8/
      data (bm2ij (  9,  7,ibeta), ibeta = 1,10) /   &
        0.916127_r8,  0.921878_r8,  0.934003_r8,  0.947506_r8,  0.959899_r8,   &
        0.970199_r8,  0.978255_r8,  0.984314_r8,  0.988755_r8,  0.991960_r8/
      data (bm2ij (  9,  8,ibeta), ibeta = 1,10) /   &
        0.893848_r8,  0.900364_r8,  0.914368_r8,  0.930438_r8,  0.945700_r8,   &
        0.958824_r8,  0.969416_r8,  0.977603_r8,  0.983746_r8,  0.988262_r8/
      data (bm2ij (  9,  9,ibeta), ibeta = 1,10) /   &
        0.892161_r8,  0.897863_r8,  0.910315_r8,  0.925021_r8,  0.939523_r8,   &
        0.952544_r8,  0.963544_r8,  0.972442_r8,  0.979411_r8,  0.984742_r8/
      data (bm2ij (  9, 10,ibeta), ibeta = 1,10) /   &
        0.922260_r8,  0.925966_r8,  0.934047_r8,  0.943616_r8,  0.953152_r8,   &
        0.961893_r8,  0.969506_r8,  0.975912_r8,  0.981167_r8,  0.985394_r8/
      data (bm2ij ( 10,  1,ibeta), ibeta = 1,10) /   &
        0.990838_r8,  0.991598_r8,  0.993128_r8,  0.994723_r8,  0.996092_r8,   &
        0.997167_r8,  0.997969_r8,  0.998552_r8,  0.998969_r8,  0.999265_r8/
      data (bm2ij ( 10,  2,ibeta), ibeta = 1,10) /   &
        0.989892_r8,  0.990727_r8,  0.992411_r8,  0.994167_r8,  0.995678_r8,   &
        0.996864_r8,  0.997751_r8,  0.998396_r8,  0.998858_r8,  0.999186_r8/
      data (bm2ij ( 10,  3,ibeta), ibeta = 1,10) /   &
        0.987287_r8,  0.988327_r8,  0.990428_r8,  0.992629_r8,  0.994529_r8,   &
        0.996026_r8,  0.997148_r8,  0.997965_r8,  0.998551_r8,  0.998967_r8/
      data (bm2ij ( 10,  4,ibeta), ibeta = 1,10) /   &
        0.982740_r8,  0.984130_r8,  0.986952_r8,  0.989926_r8,  0.992508_r8,   &
        0.994551_r8,  0.996087_r8,  0.997208_r8,  0.998012_r8,  0.998584_r8/
      data (bm2ij ( 10,  5,ibeta), ibeta = 1,10) /   &
        0.975380_r8,  0.977330_r8,  0.981307_r8,  0.985529_r8,  0.989216_r8,   &
        0.992147_r8,  0.994358_r8,  0.995975_r8,  0.997136_r8,  0.997961_r8/
      data (bm2ij ( 10,  6,ibeta), ibeta = 1,10) /   &
        0.963911_r8,  0.966714_r8,  0.972465_r8,  0.978614_r8,  0.984022_r8,   &
        0.988346_r8,  0.991620_r8,  0.994020_r8,  0.995747_r8,  0.996974_r8/
      data (bm2ij ( 10,  7,ibeta), ibeta = 1,10) /   &
        0.947187_r8,  0.951161_r8,  0.959375_r8,  0.968258_r8,  0.976160_r8,   &
        0.982540_r8,  0.987409_r8,  0.991000_r8,  0.993592_r8,  0.995441_r8/
      data (bm2ij ( 10,  8,ibeta), ibeta = 1,10) /   &
        0.926045_r8,  0.931270_r8,  0.942218_r8,  0.954297_r8,  0.965273_r8,   &
        0.974311_r8,  0.981326_r8,  0.986569_r8,  0.990394_r8,  0.993143_r8/
      data (bm2ij ( 10,  9,ibeta), ibeta = 1,10) /   &
        0.908092_r8,  0.913891_r8,  0.926288_r8,  0.940393_r8,  0.953667_r8,   &
        0.964987_r8,  0.974061_r8,  0.981038_r8,  0.986253_r8,  0.990078_r8/
      data (bm2ij ( 10, 10,ibeta), ibeta = 1,10) /   &
        0.911143_r8,  0.915972_r8,  0.926455_r8,  0.938721_r8,  0.950701_r8,   &
        0.961370_r8,  0.970329_r8,  0.977549_r8,  0.983197_r8,  0.987518_r8/


! fsb total correction factor for m2 coagulation j from i

      data  (bm2ji( 1, 1,ibeta), ibeta = 1,10) /   &
        0.753466_r8,  0.756888_r8,  0.761008_r8,  0.759432_r8,  0.748675_r8,   &
        0.726951_r8,  0.693964_r8,  0.650915_r8,  0.600227_r8,  0.545000_r8/
      data  (bm2ji( 1, 2,ibeta), ibeta = 1,10) /   &
        0.824078_r8,  0.828698_r8,  0.835988_r8,  0.838943_r8,  0.833454_r8,   &
        0.817148_r8,  0.789149_r8,  0.750088_r8,  0.701887_r8,  0.647308_r8/
      data  (bm2ji( 1, 3,ibeta), ibeta = 1,10) /   &
        1.007389_r8,  1.014362_r8,  1.028151_r8,  1.041011_r8,  1.047939_r8,   &
        1.045707_r8,  1.032524_r8,  1.007903_r8,  0.972463_r8,  0.927667_r8/
      data  (bm2ji( 1, 4,ibeta), ibeta = 1,10) /   &
        1.246157_r8,  1.255135_r8,  1.274249_r8,  1.295351_r8,  1.313362_r8,   &
        1.325187_r8,  1.329136_r8,  1.324491_r8,  1.311164_r8,  1.289459_r8/
      data  (bm2ji( 1, 5,ibeta), ibeta = 1,10) /   &
        1.450823_r8,  1.459551_r8,  1.478182_r8,  1.499143_r8,  1.518224_r8,   &
        1.533312_r8,  1.543577_r8,  1.548882_r8,  1.549395_r8,  1.545364_r8/
      data  (bm2ji( 1, 6,ibeta), ibeta = 1,10) /   &
        1.575248_r8,  1.581832_r8,  1.595643_r8,  1.610866_r8,  1.624601_r8,   &
        1.635690_r8,  1.643913_r8,  1.649470_r8,  1.652688_r8,  1.653878_r8/
      data  (bm2ji( 1, 7,ibeta), ibeta = 1,10) /   &
        1.638426_r8,  1.642626_r8,  1.651293_r8,  1.660641_r8,  1.668926_r8,   &
        1.675571_r8,  1.680572_r8,  1.684147_r8,  1.686561_r8,  1.688047_r8/
      data  (bm2ji( 1, 8,ibeta), ibeta = 1,10) /   &
        1.669996_r8,  1.672392_r8,  1.677283_r8,  1.682480_r8,  1.687028_r8,   &
        1.690651_r8,  1.693384_r8,  1.695372_r8,  1.696776_r8,  1.697734_r8/
      data  (bm2ji( 1, 9,ibeta), ibeta = 1,10) /   &
        1.686148_r8,  1.687419_r8,  1.689993_r8,  1.692704_r8,  1.695057_r8,   &
        1.696922_r8,  1.698329_r8,  1.699359_r8,  1.700099_r8,  1.700621_r8/
      data  (bm2ji( 1,10,ibeta), ibeta = 1,10) /   &
        1.694364_r8,  1.695010_r8,  1.696313_r8,  1.697676_r8,  1.698853_r8,   &
        1.699782_r8,  1.700482_r8,  1.700996_r8,  1.701366_r8,  1.701631_r8/
      data  (bm2ji( 2, 1,ibeta), ibeta = 1,10) /   &
        0.783166_r8,  0.779369_r8,  0.768044_r8,  0.747572_r8,  0.716709_r8,   &
        0.675422_r8,  0.624981_r8,  0.567811_r8,  0.507057_r8,  0.445975_r8/
      data  (bm2ji( 2, 2,ibeta), ibeta = 1,10) /   &
        0.848390_r8,  0.847100_r8,  0.840874_r8,  0.826065_r8,  0.800296_r8,   &
        0.762625_r8,  0.713655_r8,  0.655545_r8,  0.591603_r8,  0.525571_r8/
      data  (bm2ji( 2, 3,ibeta), ibeta = 1,10) /   &
        1.039894_r8,  1.043786_r8,  1.049445_r8,  1.049664_r8,  1.039407_r8,   &
        1.015322_r8,  0.975983_r8,  0.922180_r8,  0.856713_r8,  0.783634_r8/
      data  (bm2ji( 2, 4,ibeta), ibeta = 1,10) /   &
        1.345995_r8,  1.356064_r8,  1.376947_r8,  1.398304_r8,  1.412685_r8,   &
        1.414611_r8,  1.400652_r8,  1.369595_r8,  1.322261_r8,  1.260993_r8/
      data  (bm2ji( 2, 5,ibeta), ibeta = 1,10) /   &
        1.675575_r8,  1.689859_r8,  1.720957_r8,  1.756659_r8,  1.788976_r8,   &
        1.812679_r8,  1.824773_r8,  1.824024_r8,  1.810412_r8,  1.784630_r8/
      data  (bm2ji( 2, 6,ibeta), ibeta = 1,10) /   &
        1.919835_r8,  1.933483_r8,  1.962973_r8,  1.996810_r8,  2.028377_r8,   &
        2.054172_r8,  2.072763_r8,  2.083963_r8,  2.088190_r8,  2.086052_r8/
      data  (bm2ji( 2, 7,ibeta), ibeta = 1,10) /   &
        2.064139_r8,  2.074105_r8,  2.095233_r8,  2.118909_r8,  2.140688_r8,   &
        2.158661_r8,  2.172373_r8,  2.182087_r8,  2.188330_r8,  2.191650_r8/
      data  (bm2ji( 2, 8,ibeta), ibeta = 1,10) /   &
        2.144871_r8,  2.150990_r8,  2.163748_r8,  2.177731_r8,  2.190364_r8,   &
        2.200712_r8,  2.208687_r8,  2.214563_r8,  2.218716_r8,  2.221502_r8/
      data  (bm2ji( 2, 9,ibeta), ibeta = 1,10) /   &
        2.189223_r8,  2.192595_r8,  2.199540_r8,  2.207033_r8,  2.213706_r8,   &
        2.219125_r8,  2.223297_r8,  2.226403_r8,  2.228660_r8,  2.230265_r8/
      data  (bm2ji( 2,10,ibeta), ibeta = 1,10) /   &
        2.212595_r8,  2.214342_r8,  2.217912_r8,  2.221723_r8,  2.225082_r8,   &
        2.227791_r8,  2.229869_r8,  2.231417_r8,  2.232551_r8,  2.233372_r8/
      data  (bm2ji( 3, 1,ibeta), ibeta = 1,10) /   &
        0.837870_r8,  0.824476_r8,  0.793119_r8,  0.750739_r8,  0.700950_r8,   &
        0.646691_r8,  0.590508_r8,  0.534354_r8,  0.479532_r8,  0.426856_r8/
      data  (bm2ji( 3, 2,ibeta), ibeta = 1,10) /   &
        0.896771_r8,  0.885847_r8,  0.859327_r8,  0.821694_r8,  0.775312_r8,   &
        0.722402_r8,  0.665196_r8,  0.605731_r8,  0.545742_r8,  0.486687_r8/
      data  (bm2ji( 3, 3,ibeta), ibeta = 1,10) /   &
        1.076089_r8,  1.071727_r8,  1.058845_r8,  1.036171_r8,  1.002539_r8,   &
        0.957521_r8,  0.901640_r8,  0.836481_r8,  0.764597_r8,  0.689151_r8/
      data  (bm2ji( 3, 4,ibeta), ibeta = 1,10) /   &
        1.409571_r8,  1.415168_r8,  1.425346_r8,  1.432021_r8,  1.428632_r8,   &
        1.409696_r8,  1.371485_r8,  1.312958_r8,  1.236092_r8,  1.145293_r8/
      data  (bm2ji( 3, 5,ibeta), ibeta = 1,10) /   &
        1.862757_r8,  1.880031_r8,  1.918394_r8,  1.963456_r8,  2.004070_r8,   &
        2.030730_r8,  2.036144_r8,  2.016159_r8,  1.970059_r8,  1.900079_r8/
      data  (bm2ji( 3, 6,ibeta), ibeta = 1,10) /   &
        2.289741_r8,  2.313465_r8,  2.366789_r8,  2.431612_r8,  2.495597_r8,   &
        2.549838_r8,  2.588523_r8,  2.608665_r8,  2.609488_r8,  2.591662_r8/
      data  (bm2ji( 3, 7,ibeta), ibeta = 1,10) /   &
        2.597157_r8,  2.618731_r8,  2.666255_r8,  2.722597_r8,  2.777531_r8,   &
        2.825187_r8,  2.862794_r8,  2.889648_r8,  2.906199_r8,  2.913380_r8/
      data  (bm2ji( 3, 8,ibeta), ibeta = 1,10) /   &
        2.797975_r8,  2.813116_r8,  2.845666_r8,  2.882976_r8,  2.918289_r8,   &
        2.948461_r8,  2.972524_r8,  2.990687_r8,  3.003664_r8,  3.012284_r8/
      data  (bm2ji( 3, 9,ibeta), ibeta = 1,10) /   &
        2.920832_r8,  2.929843_r8,  2.948848_r8,  2.970057_r8,  2.989632_r8,   &
        3.006057_r8,  3.019067_r8,  3.028979_r8,  3.036307_r8,  3.041574_r8/
      data  (bm2ji( 3,10,ibeta), ibeta = 1,10) /   &
        2.989627_r8,  2.994491_r8,  3.004620_r8,  3.015720_r8,  3.025789_r8,   &
        3.034121_r8,  3.040664_r8,  3.045641_r8,  3.049347_r8,  3.052066_r8/
      data  (bm2ji( 4, 1,ibeta), ibeta = 1,10) /   &
        0.893179_r8,  0.870897_r8,  0.820996_r8,  0.759486_r8,  0.695488_r8,   &
        0.634582_r8,  0.579818_r8,  0.532143_r8,  0.490927_r8,  0.454618_r8/
      data  (bm2ji( 4, 2,ibeta), ibeta = 1,10) /   &
        0.948355_r8,  0.927427_r8,  0.880215_r8,  0.821146_r8,  0.758524_r8,   &
        0.697680_r8,  0.641689_r8,  0.591605_r8,  0.546919_r8,  0.506208_r8/
      data  (bm2ji( 4, 3,ibeta), ibeta = 1,10) /   &
        1.109562_r8,  1.093648_r8,  1.056438_r8,  1.007310_r8,  0.951960_r8,   &
        0.894453_r8,  0.837364_r8,  0.781742_r8,  0.727415_r8,  0.673614_r8/
      data  (bm2ji( 4, 4,ibeta), ibeta = 1,10) /   &
        1.423321_r8,  1.417557_r8,  1.402442_r8,  1.379079_r8,  1.347687_r8,   &
        1.308075_r8,  1.259703_r8,  1.201983_r8,  1.134778_r8,  1.058878_r8/
      data  (bm2ji( 4, 5,ibeta), ibeta = 1,10) /   &
        1.933434_r8,  1.944347_r8,  1.968765_r8,  1.997653_r8,  2.023054_r8,   &
        2.036554_r8,  2.029949_r8,  1.996982_r8,  1.934982_r8,  1.845473_r8/
      data  (bm2ji( 4, 6,ibeta), ibeta = 1,10) /   &
        2.547772_r8,  2.577105_r8,  2.645918_r8,  2.735407_r8,  2.830691_r8,   &
        2.917268_r8,  2.981724_r8,  3.013684_r8,  3.007302_r8,  2.961560_r8/
      data  (bm2ji( 4, 7,ibeta), ibeta = 1,10) /   &
        3.101817_r8,  3.139271_r8,  3.225851_r8,  3.336402_r8,  3.453409_r8,   &
        3.563116_r8,  3.655406_r8,  3.724014_r8,  3.766113_r8,  3.781394_r8/
      data  (bm2ji( 4, 8,ibeta), ibeta = 1,10) /   &
        3.540920_r8,  3.573780_r8,  3.647439_r8,  3.737365_r8,  3.828468_r8,   &
        3.911436_r8,  3.981317_r8,  4.036345_r8,  4.076749_r8,  4.103751_r8/
      data  (bm2ji( 4, 9,ibeta), ibeta = 1,10) /   &
        3.856771_r8,  3.879363_r8,  3.928579_r8,  3.986207_r8,  4.042173_r8,   &
        4.091411_r8,  4.132041_r8,  4.164052_r8,  4.188343_r8,  4.206118_r8/
      data  (bm2ji( 4,10,ibeta), ibeta = 1,10) /   &
        4.053923_r8,  4.067191_r8,  4.095509_r8,  4.127698_r8,  4.158037_r8,   &
        4.184055_r8,  4.205135_r8,  4.221592_r8,  4.234115_r8,  4.243463_r8/
      data  (bm2ji( 5, 1,ibeta), ibeta = 1,10) /   &
        0.935846_r8,  0.906814_r8,  0.843358_r8,  0.768710_r8,  0.695885_r8,   &
        0.631742_r8,  0.579166_r8,  0.538471_r8,  0.508410_r8,  0.486863_r8/
      data  (bm2ji( 5, 2,ibeta), ibeta = 1,10) /   &
        0.988308_r8,  0.959524_r8,  0.896482_r8,  0.821986_r8,  0.748887_r8,   &
        0.684168_r8,  0.630908_r8,  0.589516_r8,  0.558676_r8,  0.536056_r8/
      data  (bm2ji( 5, 3,ibeta), ibeta = 1,10) /   &
        1.133795_r8,  1.107139_r8,  1.048168_r8,  0.977258_r8,  0.906341_r8,   &
        0.842477_r8,  0.789093_r8,  0.746731_r8,  0.713822_r8,  0.687495_r8/
      data  (bm2ji( 5, 4,ibeta), ibeta = 1,10) /   &
        1.405692_r8,  1.385781_r8,  1.340706_r8,  1.284776_r8,  1.227085_r8,   &
        1.173532_r8,  1.127008_r8,  1.087509_r8,  1.052712_r8,  1.018960_r8/
      data  (bm2ji( 5, 5,ibeta), ibeta = 1,10) /   &
        1.884992_r8,  1.879859_r8,  1.868463_r8,  1.854995_r8,  1.841946_r8,   &
        1.829867_r8,  1.816972_r8,  1.799319_r8,  1.771754_r8,  1.729406_r8/
      data  (bm2ji( 5, 6,ibeta), ibeta = 1,10) /   &
        2.592275_r8,  2.612268_r8,  2.661698_r8,  2.731803_r8,  2.815139_r8,   &
        2.901659_r8,  2.978389_r8,  3.031259_r8,  3.048045_r8,  3.021122_r8/
      data  (bm2ji( 5, 7,ibeta), ibeta = 1,10) /   &
        3.390321_r8,  3.435519_r8,  3.545615_r8,  3.698419_r8,  3.876958_r8,   &
        4.062790_r8,  4.236125_r8,  4.378488_r8,  4.475619_r8,  4.519170_r8/
      data  (bm2ji( 5, 8,ibeta), ibeta = 1,10) /   &
        4.161376_r8,  4.216558_r8,  4.346896_r8,  4.519451_r8,  4.711107_r8,   &
        4.902416_r8,  5.077701_r8,  5.226048_r8,  5.341423_r8,  5.421764_r8/
      data  (bm2ji( 5, 9,ibeta), ibeta = 1,10) /   &
        4.843961_r8,  4.892035_r8,  5.001492_r8,  5.138515_r8,  5.281684_r8,   &
        5.416805_r8,  5.535493_r8,  5.634050_r8,  5.712063_r8,  5.770996_r8/
      data  (bm2ji( 5,10,ibeta), ibeta = 1,10) /   &
        5.352093_r8,  5.385119_r8,  5.458056_r8,  5.545311_r8,  5.632162_r8,   &
        5.710566_r8,  5.777005_r8,  5.830863_r8,  5.873123_r8,  5.905442_r8/
      data  (bm2ji( 6, 1,ibeta), ibeta = 1,10) /   &
        0.964038_r8,  0.930794_r8,  0.859433_r8,  0.777776_r8,  0.700566_r8,   &
        0.634671_r8,  0.582396_r8,  0.543656_r8,  0.517284_r8,  0.501694_r8/
      data  (bm2ji( 6, 2,ibeta), ibeta = 1,10) /   &
        1.013416_r8,  0.979685_r8,  0.907197_r8,  0.824135_r8,  0.745552_r8,   &
        0.678616_r8,  0.625870_r8,  0.587348_r8,  0.561864_r8,  0.547674_r8/
      data  (bm2ji( 6, 3,ibeta), ibeta = 1,10) /   &
        1.145452_r8,  1.111457_r8,  1.038152_r8,  0.953750_r8,  0.873724_r8,   &
        0.805955_r8,  0.753621_r8,  0.717052_r8,  0.694920_r8,  0.684910_r8/
      data  (bm2ji( 6, 4,ibeta), ibeta = 1,10) /   &
        1.376547_r8,  1.345004_r8,  1.276415_r8,  1.196704_r8,  1.121091_r8,   &
        1.058249_r8,  1.012197_r8,  0.983522_r8,  0.970323_r8,  0.968933_r8/
      data  (bm2ji( 6, 5,ibeta), ibeta = 1,10) /   &
        1.778801_r8,  1.755897_r8,  1.706074_r8,  1.649008_r8,  1.597602_r8,   &
        1.560087_r8,  1.540365_r8,  1.538205_r8,  1.549738_r8,  1.568333_r8/
      data  (bm2ji( 6, 6,ibeta), ibeta = 1,10) /   &
        2.447603_r8,  2.445172_r8,  2.443762_r8,  2.451842_r8,  2.475877_r8,   &
        2.519039_r8,  2.580118_r8,  2.653004_r8,  2.727234_r8,  2.789738_r8/
      data  (bm2ji( 6, 7,ibeta), ibeta = 1,10) /   &
        3.368490_r8,  3.399821_r8,  3.481357_r8,  3.606716_r8,  3.772101_r8,   &
        3.969416_r8,  4.184167_r8,  4.396163_r8,  4.582502_r8,  4.721838_r8/
      data  (bm2ji( 6, 8,ibeta), ibeta = 1,10) /   &
        4.426458_r8,  4.489861_r8,  4.648250_r8,  4.877510_r8,  5.160698_r8,   &
        5.477495_r8,  5.803123_r8,  6.111250_r8,  6.378153_r8,  6.586050_r8/
      data  (bm2ji( 6, 9,ibeta), ibeta = 1,10) /   &
        5.568061_r8,  5.644988_r8,  5.829837_r8,  6.081532_r8,  6.371214_r8,   &
        6.672902_r8,  6.963737_r8,  7.226172_r8,  7.449199_r8,  7.627886_r8/
      data  (bm2ji( 6,10,ibeta), ibeta = 1,10) /   &
        6.639152_r8,  6.707020_r8,  6.863974_r8,  7.065285_r8,  7.281744_r8,   &
        7.492437_r8,  7.683587_r8,  7.847917_r8,  7.983296_r8,  8.090977_r8/
      data  (bm2ji( 7, 1,ibeta), ibeta = 1,10) /   &
        0.980853_r8,  0.945724_r8,  0.871244_r8,  0.787311_r8,  0.708818_r8,   &
        0.641987_r8,  0.588462_r8,  0.547823_r8,  0.518976_r8,  0.500801_r8/
      data  (bm2ji( 7, 2,ibeta), ibeta = 1,10) /   &
        1.026738_r8,  0.990726_r8,  0.914306_r8,  0.828140_r8,  0.747637_r8,   &
        0.679351_r8,  0.625127_r8,  0.584662_r8,  0.556910_r8,  0.540749_r8/
      data  (bm2ji( 7, 3,ibeta), ibeta = 1,10) /   &
        1.146496_r8,  1.108808_r8,  1.028695_r8,  0.938291_r8,  0.854101_r8,   &
        0.783521_r8,  0.728985_r8,  0.690539_r8,  0.667272_r8,  0.657977_r8/
      data  (bm2ji( 7, 4,ibeta), ibeta = 1,10) /   &
        1.344846_r8,  1.306434_r8,  1.224543_r8,  1.132031_r8,  1.046571_r8,   &
        0.976882_r8,  0.926488_r8,  0.896067_r8,  0.884808_r8,  0.891027_r8/
      data  (bm2ji( 7, 5,ibeta), ibeta = 1,10) /   &
        1.670227_r8,  1.634583_r8,  1.558421_r8,  1.472939_r8,  1.396496_r8,   &
        1.339523_r8,  1.307151_r8,  1.300882_r8,  1.319622_r8,  1.360166_r8/
      data  (bm2ji( 7, 6,ibeta), ibeta = 1,10) /   &
        2.224548_r8,  2.199698_r8,  2.148284_r8,  2.095736_r8,  2.059319_r8,   &
        2.050496_r8,  2.075654_r8,  2.136382_r8,  2.229641_r8,  2.347958_r8/
      data  (bm2ji( 7, 7,ibeta), ibeta = 1,10) /   &
        3.104483_r8,  3.105947_r8,  3.118398_r8,  3.155809_r8,  3.230427_r8,   &
        3.350585_r8,  3.519071_r8,  3.731744_r8,  3.976847_r8,  4.235616_r8/
      data  (bm2ji( 7, 8,ibeta), ibeta = 1,10) /   &
        4.288426_r8,  4.331456_r8,  4.447024_r8,  4.633023_r8,  4.891991_r8,   &
        5.221458_r8,  5.610060_r8,  6.036467_r8,  6.471113_r8,  6.880462_r8/
      data  (bm2ji( 7, 9,ibeta), ibeta = 1,10) /   &
        5.753934_r8,  5.837061_r8,  6.048530_r8,  6.363800_r8,  6.768061_r8,   &
        7.241280_r8,  7.755346_r8,  8.276666_r8,  8.771411_r8,  9.210826_r8/
      data  (bm2ji( 7,10,ibeta), ibeta = 1,10) /   &
        7.466219_r8,  7.568810_r8,  7.819032_r8,  8.168340_r8,  8.582973_r8,   &
        9.030174_r8,  9.478159_r8,  9.899834_r8, 10.275940_r8, 10.595910_r8/
      data  (bm2ji( 8, 1,ibeta), ibeta = 1,10) /   &
        0.990036_r8,  0.954782_r8,  0.880531_r8,  0.797334_r8,  0.719410_r8,   &
        0.652220_r8,  0.596923_r8,  0.552910_r8,  0.519101_r8,  0.494529_r8/
      data  (bm2ji( 8, 2,ibeta), ibeta = 1,10) /   &
        1.032428_r8,  0.996125_r8,  0.919613_r8,  0.833853_r8,  0.753611_r8,   &
        0.684644_r8,  0.628260_r8,  0.583924_r8,  0.550611_r8,  0.527407_r8/
      data  (bm2ji( 8, 3,ibeta), ibeta = 1,10) /   &
        1.141145_r8,  1.102521_r8,  1.021017_r8,  0.929667_r8,  0.844515_r8,   &
        0.772075_r8,  0.714086_r8,  0.670280_r8,  0.639824_r8,  0.621970_r8/
      data  (bm2ji( 8, 4,ibeta), ibeta = 1,10) /   &
        1.314164_r8,  1.273087_r8,  1.186318_r8,  1.089208_r8,  0.999476_r8,   &
        0.924856_r8,  0.867948_r8,  0.829085_r8,  0.807854_r8,  0.803759_r8/
      data  (bm2ji( 8, 5,ibeta), ibeta = 1,10) /   &
        1.580611_r8,  1.538518_r8,  1.449529_r8,  1.350459_r8,  1.260910_r8,   &
        1.190526_r8,  1.143502_r8,  1.121328_r8,  1.124274_r8,  1.151974_r8/
      data  (bm2ji( 8, 6,ibeta), ibeta = 1,10) /   &
        2.016773_r8,  1.977721_r8,  1.895727_r8,  1.806974_r8,  1.732891_r8,   &
        1.685937_r8,  1.673026_r8,  1.697656_r8,  1.761039_r8,  1.862391_r8/
      data  (bm2ji( 8, 7,ibeta), ibeta = 1,10) /   &
        2.750093_r8,  2.723940_r8,  2.672854_r8,  2.628264_r8,  2.612250_r8,   &
        2.640406_r8,  2.723211_r8,  2.866599_r8,  3.071893_r8,  3.335217_r8/
      data  (bm2ji( 8, 8,ibeta), ibeta = 1,10) /   &
        3.881905_r8,  3.887143_r8,  3.913667_r8,  3.981912_r8,  4.111099_r8,   &
        4.316575_r8,  4.608146_r8,  4.988157_r8,  5.449592_r8,  5.974848_r8/
      data  (bm2ji( 8, 9,ibeta), ibeta = 1,10) /   &
        5.438870_r8,  5.492742_r8,  5.640910_r8,  5.886999_r8,  6.241641_r8,   &
        6.710609_r8,  7.289480_r8,  7.960725_r8,  8.693495_r8,  9.446644_r8/
      data  (bm2ji( 8,10,ibeta), ibeta = 1,10) /   &
        7.521152_r8,  7.624621_r8,  7.892039_r8,  8.300444_r8,  8.839787_r8,   &
        9.493227_r8, 10.231770_r8, 11.015642_r8, 11.799990_r8, 12.542260_r8/
      data  (bm2ji( 9, 1,ibeta), ibeta = 1,10) /   &
        0.994285_r8,  0.960012_r8,  0.887939_r8,  0.807040_r8,  0.730578_r8,   &
        0.663410_r8,  0.606466_r8,  0.559137_r8,  0.520426_r8,  0.489429_r8/
      data  (bm2ji( 9, 2,ibeta), ibeta = 1,10) /   &
        1.033505_r8,  0.998153_r8,  0.923772_r8,  0.840261_r8,  0.761383_r8,   &
        0.692242_r8,  0.633873_r8,  0.585709_r8,  0.546777_r8,  0.516215_r8/
      data  (bm2ji( 9, 3,ibeta), ibeta = 1,10) /   &
        1.132774_r8,  1.094907_r8,  1.015161_r8,  0.925627_r8,  0.841293_r8,   &
        0.767888_r8,  0.706741_r8,  0.657439_r8,  0.619135_r8,  0.591119_r8/
      data  (bm2ji( 9, 4,ibeta), ibeta = 1,10) /   &
        1.286308_r8,  1.245273_r8,  1.158809_r8,  1.061889_r8,  0.971208_r8,   &
        0.893476_r8,  0.830599_r8,  0.782561_r8,  0.748870_r8,  0.729198_r8/
      data  (bm2ji( 9, 5,ibeta), ibeta = 1,10) /   &
        1.511105_r8,  1.467141_r8,  1.374520_r8,  1.271162_r8,  1.175871_r8,   &
        1.096887_r8,  1.037243_r8,  0.997820_r8,  0.978924_r8,  0.980962_r8/
      data  (bm2ji( 9, 6,ibeta), ibeta = 1,10) /   &
        1.857468_r8,  1.812177_r8,  1.717002_r8,  1.612197_r8,  1.519171_r8,   &
        1.448660_r8,  1.405871_r8,  1.393541_r8,  1.413549_r8,  1.467532_r8/
      data  (bm2ji( 9, 7,ibeta), ibeta = 1,10) /   &
        2.430619_r8,  2.388452_r8,  2.301326_r8,  2.210241_r8,  2.139724_r8,   &
        2.104571_r8,  2.114085_r8,  2.174696_r8,  2.291294_r8,  2.467500_r8/
      data  (bm2ji( 9, 8,ibeta), ibeta = 1,10) /   &
        3.385332_r8,  3.357690_r8,  3.306611_r8,  3.269804_r8,  3.274462_r8,   &
        3.340862_r8,  3.484609_r8,  3.717740_r8,  4.048748_r8,  4.481588_r8/
      data  (bm2ji( 9, 9,ibeta), ibeta = 1,10) /   &
        4.850497_r8,  4.858280_r8,  4.896008_r8,  4.991467_r8,  5.171511_r8,   &
        5.459421_r8,  5.873700_r8,  6.426128_r8,  7.119061_r8,  7.942603_r8/
      data  (bm2ji( 9,10,ibeta), ibeta = 1,10) /   &
        6.957098_r8,  7.020164_r8,  7.197272_r8,  7.499331_r8,  7.946554_r8,   &
        8.555048_r8,  9.330503_r8, 10.263610_r8, 11.327454_r8, 12.478332_r8/
      data  (bm2ji(10, 1,ibeta), ibeta = 1,10) /   &
        0.994567_r8,  0.961842_r8,  0.892854_r8,  0.814874_r8,  0.740198_r8,   &
        0.673303_r8,  0.615105_r8,  0.565139_r8,  0.522558_r8,  0.486556_r8/
      data  (bm2ji(10, 2,ibeta), ibeta = 1,10) /   &
        1.031058_r8,  0.997292_r8,  0.926082_r8,  0.845571_r8,  0.768501_r8,   &
        0.699549_r8,  0.639710_r8,  0.588538_r8,  0.545197_r8,  0.508894_r8/
      data  (bm2ji(10, 3,ibeta), ibeta = 1,10) /   &
        1.122535_r8,  1.086287_r8,  1.009790_r8,  0.923292_r8,  0.840626_r8,   &
        0.766982_r8,  0.703562_r8,  0.650004_r8,  0.605525_r8,  0.569411_r8/
      data  (bm2ji(10, 4,ibeta), ibeta = 1,10) /   &
        1.261142_r8,  1.221555_r8,  1.137979_r8,  1.043576_r8,  0.953745_r8,   &
        0.874456_r8,  0.807292_r8,  0.752109_r8,  0.708326_r8,  0.675477_r8/
      data  (bm2ji(10, 5,ibeta), ibeta = 1,10) /   &
        1.456711_r8,  1.413432_r8,  1.322096_r8,  1.219264_r8,  1.122319_r8,   &
        1.038381_r8,  0.969743_r8,  0.916811_r8,  0.879544_r8,  0.858099_r8/
      data  (bm2ji(10, 6,ibeta), ibeta = 1,10) /   &
        1.741792_r8,  1.695157_r8,  1.596897_r8,  1.487124_r8,  1.385734_r8,   &
        1.301670_r8,  1.238638_r8,  1.198284_r8,  1.181809_r8,  1.190689_r8/
      data  (bm2ji(10, 7,ibeta), ibeta = 1,10) /   &
        2.190197_r8,  2.141721_r8,  2.040226_r8,  1.929245_r8,  1.832051_r8,   &
        1.760702_r8,  1.721723_r8,  1.719436_r8,  1.757705_r8,  1.840677_r8/
      data  (bm2ji(10, 8,ibeta), ibeta = 1,10) /   &
        2.940764_r8,  2.895085_r8,  2.801873_r8,  2.707112_r8,  2.638603_r8,   &
        2.613764_r8,  2.644686_r8,  2.741255_r8,  2.912790_r8,  3.168519_r8/
      data  (bm2ji(10, 9,ibeta), ibeta = 1,10) /   &
        4.186191_r8,  4.155844_r8,  4.101953_r8,  4.069102_r8,  4.089886_r8,   &
        4.189530_r8,  4.389145_r8,  4.707528_r8,  5.161567_r8,  5.765283_r8/
      data  (bm2ji(10,10,ibeta), ibeta = 1,10) /   &
        6.119526_r8,  6.127611_r8,  6.171174_r8,  6.286528_r8,  6.508738_r8,   &
        6.869521_r8,  7.396912_r8,  8.113749_r8,  9.034683_r8, 10.162190_r8/

! *** end of data statements.


! *** start calculations:

      constii = abs( half * ( two ) ** two3rds - one )
      sqrttwo = sqrt(two)
      dlgsqt2 = one / log( sqrttwo )

         esat01   = exp( 0.125_r8 * xxlsgat * xxlsgat )
         esac01   = exp( 0.125_r8 * xxlsgac * xxlsgac )

         esat04  = esat01 ** 4
         esac04  = esac01 ** 4

         esat05  = esat04 * esat01
         esac05  = esac04 * esac01

         esat08  = esat04 * esat04
         esac08  = esac04 * esac04

         esat09  = esat08 * esat01
         esac09  = esac08 * esac01

         esat16  = esat08 * esat08
         esac16  = esac08 * esac08

         esat20  = esat16 * esat04
         esac20  = esac16 * esac04

         esat24  = esat20 * esat04
         esac24  = esac20 * esac04

         esat25  = esat20 * esat05
         esac25  = esac20 * esac05

         esat36  = esat20 * esat16
         esac36  = esac20 * esac16

         esat49  = esat24 * esat25

         esat64  = esat20 * esat20 * esat24
         esac64  = esac20 * esac20 * esac24

         esat100 = esat64 * esat36

         dgat2   = dgatk * dgatk
         dgat3   = dgatk * dgatk * dgatk
         dgac2   = dgacc * dgacc
         dgac3   = dgacc * dgacc * dgacc

         sqdgat  = sqrt( dgatk )
         sqdgac  = sqrt( dgacc )
         sqdgat5 = dgat2 * sqdgat
         sqdgac5 = dgac2 * sqdgac
         sqdgat7 = dgat3 * sqdgat

         xm2at = dgat2 * esat16
         xm3at = dgat3 * esat36

         xm2ac = dgac2 * esac16
         xm3ac = dgac3 * esac36

! *** for the free molecular regime:  page h.3 of whitby et al. (1991)

         r       = sqdgac / sqdgat
         r2      = r * r
         r3      = r2 * r
         rx4     = r2 * r2
         r5      = r3 * r2
         r6      = r3 * r3
         rx8      = rx4 * rx4
         ri1     = one / r
         ri2     = one / r2
         ri3     = one / r3
         ri4     = ri2 * ri2
         kngat   = two * lamda / dgatk
         kngac   = two * lamda / dgacc


! *** calculate ratio of geometric mean diameters
         rat = dgacc / dgatk
! *** trap subscripts for bm0 and bm0i, between 1 and 10
!     see page h.5 of whitby et al. (1991)

      n2n = max( 1, min( 10,   &
            nint( 4.0_r8 * ( sgatk - 0.75_r8 ) ) ) )

      n2a = max( 1, min( 10,   &
            nint( 4.0_r8 * ( sgacc - 0.75_r8 ) ) ) )

      n1  = max( 1, min( 10,   &
             1 + nint( dlgsqt2 * log( rat ) ) ) )

! *** intermodal coagulation


! *** set up for zeroeth moment

! *** near-continuum form:  equation h.10a of whitby et al. (1991)

         coagnc0 = knc * (   &
          two + a * ( kngat * ( esat04 + r2 * esat16 * esac04 )   &
                    + kngac * ( esac04 + ri2 * esac16 * esat04 ) )   &
                    + ( r2 + ri2 ) * esat04 * esac04  )


! *** free-molecular form:  equation h.7a of whitby et al. (1991)

         coagfm0 = kfmatac * sqdgat * bm0ij(n1,n2n,n2a) * (   &
                   esat01 + r * esac01 + two * r2 * esat01 * esac04   &
                 + rx4 * esat09 * esac16 + ri3 * esat16 * esac09   &
                 + two * ri1 * esat04 + esac01  )


! *** loss to accumulation mode

! *** harmonic mean

      coagatac0 = coagnc0 * coagfm0 / ( coagnc0 + coagfm0 )

      qn12 = coagatac0


! *** set up for second moment
!      the second moment equations are new and begin with equations a1
!     through a4 of binkowski and shankar (1995). after some algebraic
!     rearrangement and application of the extended mean value theorem
!     of integral calculus, equations are obtained that can be solved
!     analytically with correction factors as has been done by
!     whitby et al. (1991)

! *** the term ( dp1 + dp2 ) ** (2/3) in equations a3 and a4 of
!     binkowski and shankar (1995) is approximated by
!     (dgat ** 3 + dgac **3 ) ** 2/3

! *** near-continuum form

      i1nc = knc * dgat2 * (   &
             two * esat16   &
           + r2 * esat04 * esac04   &
           + ri2 * esat36 * esac04   &
           + a * kngat * (   &
                 esat04   &
           +     ri2 * esat16 * esac04   &
           +     ri4 * esat36 * esac16   &
           +     r2 * esac04 )  )




! *** free-molecular form

       i1fm =  kfmatac * sqdgat5 * bm2ij(n1,n2n,n2a) * (   &
               esat25   &
            +  two * r2 * esat09 * esac04   &
            +  rx4 * esat01 * esac16   &
            +  ri3 * esat64 * esac09   &
            +  two * ri1 * esat36 * esac01   &
            +  r * esat16 * esac01  )



! *** loss to accumulation mode

! *** harmonic mean

      i1 = ( i1fm * i1nc ) / ( i1fm + i1nc )

      coagatac2 = i1

      qs12 = coagatac2


! *** gain by accumulation mode

      coagacat2 = ( ( one + r6 ) ** two3rds - rx4 ) * i1

      qs21 = coagacat2 * bm2ji(n1,n2n,n2a)

! *** set up for third moment

! *** near-continuum form: equation h.10b of whitby et al. (1991)

      coagnc3 = knc * dgat3 * (   &
                two * esat36   &
              + a * kngat * ( esat16 + r2 * esat04 * esac04 )   &
              + a * kngac * ( esat36 * esac04 + ri2 * esat64 * esac16 )   &
              + r2 * esat16 * esac04 + ri2 * esat64 * esac04 )


! *** free_molecular form: equation h.7b of whitby et al. (1991)

      coagfm3 = kfmatac * sqdgat7 * bm3i( n1, n2n, n2a ) * (   &
               esat49   &
              +  r * esat36  * esac01   &
              + two * r2 * esat25  * esac04   &
              + rx4 * esat09  * esac16   &
              + ri3 * esat100 * esac09   &
              + two * ri1 * esat64  * esac01 )

! *** gain by accumulation mode = loss from aitken mode

! *** harmonic mean

      coagatac3 = coagnc3 * coagfm3 / ( coagnc3 + coagfm3 )

      qv12 = coagatac3

! *** intramodal coagulation

! *** zeroeth moment

! *** aitken mode

! *** near-continuum form: equation h.12a of whitby et al. (1991)

      coagnc_at = knc * (one + esat08 + a * kngat * (esat20 + esat04))

! *** free-molecular form: equation h.11a of whitby et al. (1991)

      coagfm_at = kfmat * sqdgat * bm0(n2n) *   &
                 ( esat01 + esat25 + two * esat05 )


! *** harmonic mean

      coagatat0 = coagfm_at * coagnc_at / ( coagfm_at + coagnc_at )

      qn11 = coagatat0


! *** accumulation mode

! *** near-continuum form: equation h.12a of whitby et al. (1991)

      coagnc_ac = knc * (one + esac08 + a * kngac * (esac20 + esac04))

! *** free-molecular form: equation h.11a of whitby et al. (1991)

      coagfm_ac = kfmac * sqdgac * bm0(n2a) *   &
                   ( esac01 + esac25 + two * esac05 )

! *** harmonic mean

      coagacac0 = coagfm_ac * coagnc_ac / ( coagfm_ac + coagnc_ac )

      qn22 = coagacac0


! *** set up for second moment
!      the second moment equations are new and begin with 3.11a on page
!     45 of whitby et al. (1991). after some algebraic rearrangement and
!     application of the extended mean value theorem of integral calculus
!     equations are obtained that can be solved analytically with
!     correction factors as has been done by whitby et al. (1991)

! *** aitken mode

! *** near-continuum

      i1nc_at = knc * dgat2 * (   &
             two * esat16   &
           + esat04 * esat04   &
           + esat36 * esat04   &
           + a * kngat * (   &
                two * esat04   &
           +     esat16 * esat04   &
           +     esat36 * esat16 )  )

! *** free- molecular form

       i1fm_at =  kfmat * sqdgat5 * bm2ii(n2n) * (   &
               esat25   &
            +  two * esat09 * esat04   &
            +  esat01 * esat16   &
            +  esat64 * esat09   &
            +  two * esat36 * esat01   &
            +  esat16 * esat01  )

      i1_at = ( i1nc_at * i1fm_at ) / ( i1nc_at + i1fm_at  )

      coagatat2 = constii * i1_at

      qs11 = coagatat2 * bm2iitt(n2n)

! *** accumulation mode

! *** near-continuum

      i1nc_ac = knc * dgac2 * (   &
             two * esac16   &
           + esac04 * esac04   &
           + esac36 * esac04   &
           + a * kngac * (   &
                two * esac04   &
           +     esac16 * esac04   &
           +     esac36 * esac16 )  )

! *** free- molecular form

       i1fm_ac =  kfmac * sqdgac5 * bm2ii(n2a) * (   &
               esac25   &
            +  two * esac09 * esac04   &
            +  esac01 * esac16   &
            +  esac64 * esac09   &
            +  two * esac36 * esac01   &
            +  esac16 * esac01  )

      i1_ac = ( i1nc_ac * i1fm_ac ) / ( i1nc_ac + i1fm_ac  )

      coagacac2 = constii * i1_ac

      qs22 = coagacac2 * bm2iitt(n2a)


      return

      end  subroutine getcoags

!----------------------------------------------------------------------
!----------------------------------------------------------------------

   end module modal_aero_coag



