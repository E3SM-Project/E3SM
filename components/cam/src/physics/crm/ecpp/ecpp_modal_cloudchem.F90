module ecpp_modal_cloudchem

!-----------------------------------------------------------------
!  Module interface for cloud chemistry used in the ECPP treatment 
! in the MMF model
!  Adopted the similar one used in the ECPP 
!  for the WRF-chem model written by Dick Easter
!
!  Minghuai Wang, 2009-11
!------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
!==Guangxing Lin
!   use abortutils,   only: endrun
   use cam_abortutils,   only: endrun
!==Guangxing Lin

   implicit none

   public parampollu_tdx_cldchem

contains

!-----------------------------------------------------------------------

subroutine parampollu_tdx_cldchem(               &
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
		laicwpair_of_aerosol,pbuf                          )

!-----------------------------------------------------------------------
! DESCRIPTION
!
! parampollu_tdx_cldchem does cloud chemistry
!    for one main-integ time sub-step
!
! incoming chem_sub_new holds current sub-class mixing ratios
! outgoing chem_sub_new holds updated sub-class mixing ratios
!
! In the beginning of the subroutine, the vertical coordinate (from bottom to top in ECPP)
! is converted into the one used in CAM: from the top to the bottom. And at the end of the
! subroutine, the vertical coordinate is converted back. 
!
!-----------------------------------------------------------------------

!==Guangxing Lin
!	use module_data_ecpp1, only:  p_qv, p_qc
!==Guangxing Lin

	use module_data_radm2, only:  epsilc

	use module_data_mosaic_asect, only:  ai_phase, cw_phase,   &
		massptr_aer, maxd_asize, maxd_atype,   &
		ncomp_aer, nphase_aer, nsize_aer, ntype_aer, numptr_aer

	use module_data_ecpp1

        use mo_setsox,         only : setsox
        use mo_mass_xforms,    only : mmr2vmr, vmr2mmr
        use modal_aero_rename, only : modal_aero_rename_sub
        use physconst,      only: gravit
        use ppgrid,         only: pcols, pver
        use time_manager,  only: get_nstep
        use mo_mean_mass,   only: set_mean_mass
        use chem_mods,      only: gas_pcnst, nfs, indexm, adv_mass !==Guangxing Lin
        use mo_setinv,         only : setinv
        use constituents,    only: pcnst
        use mo_gas_phase_chemdr, only: map2chm
        use chemistry,      only: imozart
        use physics_buffer, only: physics_buffer_desc

	use module_ecpp_util, only:  ecpp_error_fatal, ecpp_message
       
!   arguments
	integer, intent(in) ::                  &
		ktau, ktau_pp, itstep_sub,   &
		it, jt, kts, ktebnd, ktecen, &
                itstep_hybrid
!   ktau - time step number
!   ktau_pp - time step number for "parameterized pollutants" calculations
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!	chem_driver and routines under it do calculations
!	over these spatial indices.

	integer, intent(in) :: idiagaa_ecpp(1:199), ldiagaa_ecpp(1:199)

	real(r8), intent(in) :: dtstep, dtstep_sub
!   dtstep - main model time step (s)
!   dtstep_sub - sub time step (s) currently used in ecpp main-integ routine

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
!	(ug/kg for mass species, #/kg for number species)

	integer, intent(in) :: ncls_ecpp, ncls_use

	integer, intent(in), dimension( 1:2, 1:maxcls_ecpp ) ::   &
		kdraft_bot_use, kdraft_top_use,   &
		mtype_updnenv_use

	real(r8), intent(inout), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:num_chem_ecpp ) ::   &
		chem_sub_new

	real(r8), intent(inout), dimension( 1:num_chem_ecpp ) :: &
		del_chem_clm_cldchem

        real(r8), intent(inout), dimension( 1:num_chem_ecpp ) :: &
                del_chem_clm_rename

        real(r8), intent(inout), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:num_chem_ecpp ) :: &
                del_cldchem3d                 ! 3D change from aqueous chemistry

        real(r8), intent(inout), dimension( kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2, 1:num_chem_ecpp ) :: &
                del_rename3d                 ! 3D change from modal merging 

        real(r8), intent(inout) :: aqso4_h2o2,            &         ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
                                   aqso4_o3                          ! SO4 aqueous phase chemistry due to O3 (kg/m2)

        real(r8), intent(inout), dimension(kts:ktecen, 1:2, 1:maxcls_ecpp, 1:2)  ::  &
                                   xphlwc3d                          ! pH value multiplied by lwc


	real(r8), intent(inout), dimension( kts:ktecen, 0:2, 0:maxcls_ecpp ) :: &
		ardz_cen_old, ardz_cen_new, acen_tavg_use, acen_prec_use

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
	integer :: icc, iccpp, iccpp1, iccpp2, ipp
	integer :: jcls
	integer :: k, kk, l, km
	integer :: numgas_aqfrac
	integer :: p1st
        integer :: m, n
        integer :: im, in, lnumcw
        integer :: ncol

	real(r8) :: tmpa, tmpa1, tmpa2, tmpb1, tmpb2, tmpq, tmpq1, tmpq2, tmpx, tmpx2, tmpy, tmpy2
        real(r8) :: dtmpchem

	real(r8), parameter :: qcldwtr_cutoff = 1.0e-6_r8
	real(r8) :: dt_tmp

	real(r8), allocatable :: p_tmp(:,:,:), t_tmp(:,:,:), rho_tmp(:,:,:), &
		alt_tmp(:,:,:), cldfra_tmp(:,:,:), &
		qlsink_tmp(:,:,:), precr_tmp(:,:,:), &
		precs_tmp(:,:,:), precg_tmp(:,:,:), preci_tmp(:,:,:)
        real(r8), allocatable :: zero_tmp(:, :, :)
	real(r8), allocatable :: moist_tmp(:,:,:,:)
	real(r8), allocatable :: chem_tmpa(:,:,:,:), chem_tmpb(:,:,:,:), chem_tmpc(:,:,:,:)

        real(r8), allocatable :: cwat_tmp(:,:,:)
        real(r8), allocatable :: pdel_tmp(:,:,:)

        real(r8), allocatable :: aqso4_h2o2_tmp(:)
        real(r8), allocatable :: aqso4_o3_tmp(:)
!        real(r8), allocatable :: xphlwc_tmp(:)

        real(r8), allocatable :: mmr(:, :), vmr(:,:), mmrcw(:, :), vmrcw(:, :)
!==Guangxing Lin
        real(r8), allocatable :: vmr_3d(:,:,:), vmrcw_3d(:,:,:)
!==Guangxing Lin
        real(r8), allocatable :: vmr_sv1(:,:), vmrcw_sv1(:,:)
        real(r8), allocatable :: mbar(:)
        real(r8), allocatable :: mmr_3d(:, :, :), mmrcw_3d(:, :, :), mbar_3d(:, :)
!==Guangxing Lin
 !       real(r8), allocatable :: cldnum(:)
        real(r8), allocatable :: cldnum(:,:)
!==Guangxing Lin
        real(r8), allocatable :: invariants(:,:)

        real(r8) :: invariants_full(pcols, pver, nfs)
        real(r8) :: t_full(pcols, pver)
        real(r8) :: pmid_full(pcols, pver)
        real(r8) :: h2ovmr_full(pcols, pver)
        real(r8) :: vmr_full(pcols, pver, gas_pcnst)

        real(r8), allocatable :: qsrflx_full(:, :,:), qqcwsrflx_full(:, :,:)
        integer  :: nsrflx 
        integer  :: nstep
        integer  :: jsrflx_rename
        integer  :: latndx_full(pcols, pver)
        integer  :: lonndx_full(pcols, pver)
        real(r8) :: pdel_full(pcols, pver)
        real(r8) :: dqdt(pver, gas_pcnst)
        real(r8) :: dqdt_other(pver, gas_pcnst)
        real(r8) :: dqqcwdt(pver, gas_pcnst)
        real(r8) :: dqqcwdt_other(pver, gas_pcnst)
        logical  :: dotendrn(gas_pcnst)
        logical  :: dotendqqcwrn(gas_pcnst)
        logical  :: is_dorename_atik
        logical  :: dorename_atik(pver)

	p1st = param_first_ecpp
	numgas_aqfrac = num_chem_ecpp

        nsrflx = 2
        jsrflx_rename = 2
        nstep = get_nstep()
        

!
! load arrays for interfacing with cloud chemistry subroutine 
!
! use the wrfchem "i" index for the ecpp icc & ipp sub-class indices
! use the wrfchem "j" index for the ecpp jcls class index
! all the temporary real*4 arrays must be dimensioned kts:ktebnd
!
	allocate ( p_tmp(         1:4,kts:ktecen,1:ncls_use) )
	allocate ( t_tmp(         1:4,kts:ktecen,1:ncls_use) )
	allocate ( rho_tmp(       1:4,kts:ktecen,1:ncls_use) )
	allocate ( alt_tmp(       1:4,kts:ktecen,1:ncls_use) )
	allocate ( cldfra_tmp(    1:4,kts:ktecen,1:ncls_use) )
	allocate ( qlsink_tmp(    1:4,kts:ktecen,1:ncls_use) )
	allocate ( precr_tmp(     1:4,kts:ktecen,1:ncls_use) )
	allocate ( precs_tmp(     1:4,kts:ktecen,1:ncls_use) )
	allocate ( precg_tmp(     1:4,kts:ktecen,1:ncls_use) )
	allocate ( preci_tmp(     1:4,kts:ktecen,1:ncls_use) )
        allocate ( cwat_tmp(      1:4,kts:ktecen,1:ncls_use) )
        allocate ( pdel_tmp(      1:4,kts:ktecen,1:ncls_use) )
	allocate ( zero_tmp(      1:4,kts:ktecen,1:ncls_use) )
	allocate ( moist_tmp(     1:4,kts:ktecen,1:ncls_use,1:num_moist_ecpp) )
	allocate ( chem_tmpa(     1:4,kts:ktecen,1:ncls_use,1:num_chem_ecpp) )
	allocate ( chem_tmpb(     1:4,kts:ktecen,1:ncls_use,1:num_chem_ecpp) )
        allocate ( chem_tmpc(     1:4,kts:ktecen,1:ncls_use,1:num_chem_ecpp) )

       allocate  ( mmr(kts:ktecen,1:gas_pcnst) )
       allocate  ( vmr(kts:ktecen,1:gas_pcnst) )   
       allocate  ( mmrcw(kts:ktecen,1:gas_pcnst) ) 
       allocate  ( vmrcw(kts:ktecen,1:gas_pcnst) )  
       allocate  ( vmr_sv1(kts:ktecen,1:gas_pcnst) )
       allocate  ( vmrcw_sv1(kts:ktecen,1:gas_pcnst) )
       allocate  ( mbar(kts:ktecen) )
!==Guangxing Lin
       !allocate  ( cldnum(kts:ktecen) )
       allocate  ( cldnum(1,kts:ktecen) )
!==Guangxing Lin
       allocate  ( invariants(kts:ktecen, nfs) )
!==Guangxing Lin
       allocate  ( vmr_3d(1, kts:ktecen,1:gas_pcnst) )
       allocate  ( vmrcw_3d(1, kts:ktecen, 1:gas_pcnst) )
!==Guangxing Lin
       allocate  ( mmr_3d(1, kts:ktecen,1:gas_pcnst) )
       allocate  ( mmrcw_3d(1, kts:ktecen, 1:gas_pcnst) )
       allocate  ( mbar_3d(1, kts:ktecen) )

       allocate  (aqso4_h2o2_tmp(kts:ktecen))
       allocate  (aqso4_o3_tmp(kts:ktecen))
       !allocate  (xphlwc_tmp(kts:ktecen))

       allocate  (qsrflx_full(pcols, gas_pcnst, nsrflx))
       allocate  (qqcwsrflx_full(pcols, gas_pcnst, nsrflx))
     
       zero_tmp(:, :, :) = 0.0_r8 

! chem_tmpa, chem_tmpb and chem_tmpc start from bottom to top, just as chem_sub_new
! But mmr, mmrcw are reordered, starts from top to the bottom for aqueous chemistry at CAM.
        do l = 1, num_chem_ecpp
        do jcls = 1, ncls_use
        do kk = kts, ktecen
            k = min( kk, ktecen )
            do icc = 1, 2
            do ipp = 1, 2
                iccpp = 2*(icc-1) + ipp
                chem_tmpa(iccpp,k,jcls,l) = chem_sub_new(k,icc,jcls,l)
            end do
            end do
        end do
        end do
        end do
        chem_tmpb(:,:,:,:) = chem_tmpa(:,:,:,:)
        chem_tmpc(:,:,:,:) = chem_tmpa(:,:,:,:)

!
! prepare fields for aqueous chemistry at CAM. 
	do kk = kts, ktecen
	  k = min( kk, ktecen )
!
!  vertical coordinate is from bottom to top in the ECPP,
!  so convert it to from top to the bottom for aqueous chemistry at CAM.
          km = ktecen-k+1
!==Guangxing Lin
          do icc=1,4
             do ipp=1, ncls_use 
	  p_tmp(icc,k,ipp) = pcen_bar(km)
	  t_tmp(icc,k,ipp) = tcen_bar(km)
	  rho_tmp(icc,k,ipp) = rhocen_bar(km)
	  alt_tmp(icc,k,ipp) = 1.0_r8/rhocen_bar(km)
          pdel_tmp(icc,k,ipp) = rhocen_bar(km)*dzcen(km)*gravit
	      end do
          end do 
!          p_tmp(1:4,k,1:ncls_use) = pcen_bar(km)
!	  t_tmp(1:4,k,1:ncls_use) = tcen_bar(km)
!	  rho_tmp(1:4,k,1:ncls_use) = rhocen_bar(km)
!	  alt_tmp(1:4,k,1:ncls_use) = 1.0/rhocen_bar(km)
!          pdel_tmp(1:4,k,1:ncls_use) = rhocen_bar(km)*dzcen(km)*gravit
!==Guangxing Lin
	end do

	cldfra_tmp(:,:,:) = 0.0_r8
	moist_tmp(:,:,:,:) = 0.0_r8
	qlsink_tmp(:,:,:) = 0.0_r8
	precr_tmp(:,:,:) = 0.0_r8
	precg_tmp(:,:,:) = 0.0_r8
	precs_tmp(:,:,:) = 0.0_r8
	preci_tmp(:,:,:) = 0.0_r8
        cwat_tmp(:,:,:) = 0.0_r8

	do jcls = 1, ncls_use
	do k = kts, ktecen
!
!  vertical coordinate is from bottom to top in the ECPP,
!  so convert it to from top to the bottom for aqueous chemistry at CAM.
        km = ktecen-k+1
	do icc = 1, 2
	do ipp = 1, 2
	    iccpp = 2*(icc-1) + ipp
	    if (ipp == 1) then
		tmpa = acen_tavg_use(km,icc,jcls) - acen_prec_use(km,icc,jcls)
	    else
		tmpa = acen_prec_use(km,icc,jcls)
	    end if
	    tmpq = qcloud_sub2(km,icc,jcls,ipp)
	    if ((tmpa > afrac_cut_0p5) .and. (tmpq > qcldwtr_cutoff)) then
		moist_tmp(iccpp,k,jcls,p_qc) = tmpq
		qlsink_tmp(iccpp,k,jcls) = qlsink_sub2(km,icc,jcls,ipp)
                cwat_tmp(iccpp,k,jcls) = tmpq
	    end if

	    if (icc == 2) then
               if(tmpa > afrac_cut_0p5) then 
                cldfra_tmp(iccpp,k,jcls) = 1.0_r8
               end if
            end if

	    precr_tmp(iccpp,k,jcls) = precr_sub2(km,icc,jcls,ipp)
	    precs_tmp(iccpp,k,jcls) = precs_sub2(km,icc,jcls,ipp)
	end do
	end do
	end do
	end do


	dt_tmp = dtstep_sub

	if (cldchem_onoff_ecpp > 0) then

       do jcls = 1, ncls_use
       do icc = 2, 2 !  In clear sky, cloud chemistry and renaming are not called. 
       do ipp = 1, 2 
          iccpp = 2*(icc-1) + ipp
          ncol = 1

          !----------------------------------------------------------------------
          !     calculate cldnum from cloud borne aerosol particles
          !     Vertical coordinate is from bottom to top in the ECPP for chem_tempb,
          !     so convert it to from top to the bottom for aqueous chemistry at CAM.
          !----------------------------------------------------------------------
          cldnum(1,:) = 0.0_r8
          do in=1, ntype_aer
            do im=1, nsize_aer(in)
              lnumcw = numptr_aer(im, in, cw_phase)
              do k=kts, ktecen
                km=ktecen-k+1
                cldnum(1,k) = cldnum(1,k)+chem_tmpb(iccpp,k,jcls,lnumcw)           
              end do
            end do
          end do

          !-----------------------------------------------------------------------      
          !        ... map incoming concentrations to working array
          !   Vertical coordinate is from bottom to top in the ECPP for chem_tempb,
          !     so convert it to from top to the bottom for aqueous chemistry at CAM.
          !-----------------------------------------------------------------------      
          mmr(:, :) = 0.0_r8
          mmrcw(:, :) = 0.0_r8
          do m = 1,pcnst
            n = map2chm(m)
            if( n > 0 ) then
              do k = kts, ktecen
                km = ktecen-k+1
                mmr(k,n) = chem_tmpb(iccpp,km,jcls,m)
                mmrcw(k,n) = chem_tmpb(iccpp,km,jcls,m+pcnst)
              end do
            end if
          end do

          !-----------------------------------------------------------------------      
          !        ... Set atmosphere mean mass
          !-----------------------------------------------------------------------      
          call set_mean_mass( ncol, mmr_3d, mbar_3d )
          mbar(:) = mbar_3d(1, :)

          !-----------------------------------------------------------------------      
          !        ... Xform from mmr to vmr
          !-----------------------------------------------------------------------      
	  !vmr_3d(1, :, :) = vmr(:, :) !Guangxing Lin
	  !vmrcw_3d(1, :, :) = vmrcw(:, :)!==Guangxing Lin
	  mmr_3d(1, :, :) = mmr(:, :)
	  mmrcw_3d(1, :, :) = mmrcw(:, :)
!==Guangxing Lin          
          do m = 1,gas_pcnst
           if( adv_mass(m) /= 0._r8 ) then
           do k =kts, ktecen 
             vmr(k,m) = mbar(k) * mmr(k,m) / adv_mass(m)
             vmrcw(k,m) = mbar(k) * mmrcw(k,m) / adv_mass(m)
          end do
          end if
        end do

!          call mmr2vmr( mmr_3d, vmr_3d, mbar_3d, ncol )
!          call mmr2vmr( mmrcw_3d, vmrcw_3d, mbar_3d, ncol )


!          call mmr2vmr( mmr_3d, vmr, mbar, ncol )
!          call mmr2vmr( mmrcw_3d, vmrcw, mbar, ncol )
	  vmr_3d(1, :, :) = vmr(:, :) !Guangxing Lin
	  vmrcw_3d(1, :, :) = vmrcw(:, :)!==Guangxing Lin

          vmr_sv1 = vmr
          vmrcw_sv1 = vmrcw

!            vmr_sv1 = vmr_3d(1,:,:)
!          vmrcw_sv1 = vmrcw_3d(1,:,:)

!          vmr(:,:) = vmr_3d(1,:,:)
!          vmrcw(:,:) = vmrcw_3d(1,:,:)
!==Guangxing Lin          



          !-----------------------------------------------------------------------      
          !        ... Set the "invariants"
          !-----------------------------------------------------------------------  
!==Guangxing Lin          
           !h2ovmr_full(:, :) = 0.0   ! h2ommr is not used in CAM aqueous chemistry, so set it to zero here. 
          h2ovmr_full(:it, :) = 0.0_r8   ! h2ommr is not used in CAM aqueous chemistry, so set it to zero here. 
!==Guangxing Lin          
          do kk = kts, ktecen
            k = min( kk, ktecen)
            t_full(:it, k) = t_tmp(iccpp, k,jcls)  
            pmid_full(:it, k) = p_tmp(iccpp, k, jcls)  
            do n=1, gas_pcnst
              vmr_full(:it, k, n) = vmr(k, n)
            end do
          end do
          call setinv( invariants_full(:it,:,:), t_full, h2ovmr_full(:it,:), vmr_full(:it,:,:), pmid_full, it, jt, pbuf)   ! jt=lchnk
          invariants(:, :) = invariants_full(it, :, :) 

          !--------------------------------------------------------------------------
          !        ... Aqueous chemistry
          !--------------------------------------------------------------------------
!==Guangxing Lin
          call setsox( ncol,   &
               jt,  &
               imozart-1,&   
               dt_tmp,                  &
               p_tmp(iccpp:iccpp, :, jcls),   &
               pdel_tmp(iccpp:iccpp, :, jcls),   &
               t_tmp(iccpp:iccpp, :, jcls),   &
               mbar_3d,   &      
               cwat_tmp(iccpp:iccpp, :, jcls),    &
               cldfra_tmp(iccpp:iccpp, :, jcls),  &
               cldnum, &  
               invariants_full(it:it,:,indexm), &  
               invariants_full(it:it,:,:),  &    
               vmrcw_3d,    &  
               vmr_3d)           
               
!          call setsox( ncol,   &
!               p_tmp(iccpp:iccpp, :, jcls),   &
!               dt_tmp,                  &
!               t_tmp(iccpp:iccpp, :, jcls),   &
!               zero_tmp(iccpp:iccpp, :, jcls),   &
!               cwat_tmp(iccpp, :, jcls),    &
!               invariants(1,indexm), &  
!               vmr, &           
!               invariants,  &    
!               jt,  &
!               pdel_tmp(iccpp:iccpp, :, jcls),   &
!               mbar,   &      
!               zero_tmp(iccpp, :, jcls),  &
!               cldfra_tmp(iccpp, :, jcls),  &
!               cldnum, &  
!               vmrcw,    &  
!               imozart-1,&   
!               aqso4_h2o2_tmp(:), aqso4_o3_tmp(:), xphlwc_tmp(:) )

          !-----------------------------------------------------------------------
          !         ... Xform from vmr to mmr
          !-----------------------------------------------------------------------      
!          call vmr2mmr( vmr, mmr_3d, mbar, ncol )
!          call vmr2mmr( vmrcw, mmrcw_3d, mbar, ncol )
          vmr(:,:)   = vmr_3d(1,:,:)
          vmrcw(:,:) = vmrcw_3d(1,:,:)
          do m = 1,gas_pcnst
           if( adv_mass(m) /= 0._r8 ) then
           do k =kts, ktecen 
             mmr_3d(1,k,m) = adv_mass(m)*vmr(k,m)/mbar(k)
             mmrcw_3d(1,k,m) = adv_mass(m)*vmrcw(k,m)/mbar(k)
          end do
          end if
        end do
!==Guangxing Lin
          mmr(:, :) = mmr_3d(1, :, :)
          mmrcw(:, :) = mmrcw_3d(1, :, :)

          !-----------------------------------------------------------------------      
          !         ... Form the tendencies
          !   Vertical coordinate is from top to bottom in the aqueous chemistry at CAM,  
          !     so convert it to from bottom to the top in the ECPP for chem_tmpb.
          !----------------------------------------------------------------------- 
          do m = 1,pcnst
             n = map2chm(m)
             if( n > 0 ) then
              do k = kts, ktecen
                km = ktecen-k+1
                chem_tmpb(iccpp, k,jcls,m) = mmr(km,n)
                chem_tmpb(iccpp, k,jcls,m+pcnst) = mmrcw(km,n)
              end do
             end if
          end do

          do k = kts, ktecen
              km = ktecen-k+1  ! acen is defined in the ECPP (from bottom to top)
              if (ipp == 1) then
                 tmpa = acen_tavg_use(k,icc,jcls) - acen_prec_use(k,icc,jcls)
              else
                 tmpa = acen_prec_use(k,icc,jcls)
              end if
              if (tmpa > afrac_cut_0p5) then
!==Guangxing Lin
                 !aqso4_h2o2 = aqso4_h2o2+tmpa * aqso4_h2o2_tmp(km)
                 !aqso4_o3 = aqso4_o3 + tmpa * aqso4_o3_tmp(km)
                 aqso4_h2o2 = aqso4_h2o2
                 aqso4_o3 = aqso4_o3 
!==Guangxing Lin
              end if 
!
! xphlwc_tmp is defined in CAM( top to bottom), and xphlwc3d is defined in ECPP (bottom to top)
!==Guangxing Lin
             ! xphlwc3d(k,icc,jcls,ipp) = xphlwc3d(k,icc,jcls,ipp) + xphlwc_tmp(km) * tmpa
              xphlwc3d(k,icc,jcls,ipp) = xphlwc3d(k,icc,jcls,ipp) 
!==Guangxing Lin

          end do

!-----------------------------------------------------------------------------
!        ----- renaming: modal aerosol mode merging ------
!-----------------------------------------------------------------------------
          if(rename_onoff_ecpp > 0) then
            do kk = kts, ktecen
              k = min( kk, ktecen)
              pdel_full(:ncol, k) = p_tmp(iccpp, k, jcls)
            end do
            latndx_full(:ncol,:) = 1 
            lonndx_full(:ncol,:) = 1 
            qsrflx_full(:ncol,:,:) = 0.0_r8
            qqcwsrflx_full(:ncol,:,:) = 0.0_r8
            dotendrn(:) = .false.
            dotendqqcwrn(:) = .false.
            dorename_atik(:) = .true.
            is_dorename_atik = .true.
            dqdt (:,:) = 0.0_r8
            dqqcwdt(:,:) = 0.0_r8
            dqdt_other(:,:)=(vmr-vmr_sv1)/dt_tmp
            dqqcwdt_other(:,:)=(vmrcw-vmrcw_sv1)/dt_tmp

            call modal_aero_rename_sub('ecpp_modal_cloudchem', jt,    &
                                         ncol, nstep,                  &
                                         imozart-1,  dt_tmp,           &
                                         pdel_full,                    &
                                         dotendrn,        vmr,        &
                                         dqdt,        dqdt_other,      &
                                         dotendqqcwrn,    vmrcw,       &
                                         dqqcwdt,    dqqcwdt_other,    &
                                         is_dorename_atik, dorename_atik, &
                                         jsrflx_rename,  nsrflx,       &
                                         qsrflx_full,   qqcwsrflx_full         )
!==Guangxing Lin deleted latndx_full, lonndx_full 
             vmr = vmr + dqdt * dt_tmp
             vmrcw = vmrcw + dqqcwdt * dt_tmp

             !-----------------------------------------------------------------------
             !         ... Xform from vmr to mmr
             !-----------------------------------------------------------------------      
!==Guangxing Lin
!             call vmr2mmr( vmr, mmr_3d, mbar, ncol )
!             call vmr2mmr( vmrcw, mmrcw_3d, mbar, ncol )
          do m = 1,gas_pcnst
           if( adv_mass(m) /= 0._r8 ) then
           do k =kts, ktecen 
             mmr_3d(1,k,m) = adv_mass(m)*vmr(k,m)/mbar(k)
             mmrcw_3d(1,k,m) = adv_mass(m)*vmrcw(k,m)/mbar(k)
          end do
          end if
        end do
!==Guangxing Lin
             mmr(:, :) = mmr_3d(1, :, :)
             mmrcw(:, :) = mmrcw_3d(1, :, :)

             !-----------------------------------------------------------------------      
             !         ... Form the tendencies
             !   Vertical coordinate is from top to bottom in the aqueous chemistry at CAM,  
             !     so convert it to from bottom to the top in the ECPP for chem_tmpb.
             !----------------------------------------------------------------------- 
             do m = 1,pcnst
               n = map2chm(m)
               if( n > 0 ) then
                do k = kts, ktecen
                  km = ktecen-k+1
                  chem_tmpc(iccpp, k,jcls,m) = mmr(km,n)
                  chem_tmpc(iccpp, k,jcls,m+pcnst) = mmrcw(km,n)
                end do
               end if
             end do


          end if  ! (rename_onoff_ecpp > 0)

        end do
        end do
        end do

	do l = p1st, num_chem_ecpp
	    tmpx = 0.0_r8
            tmpx2 = 0.0_r8
	    do k = kts, ktecen
		tmpy = 0.0_r8
                tmpy2 = 0.0_r8
		do jcls = 1, ncls_use
		do icc = 1, 2
		do ipp = 1, 2
		    iccpp = 2*(icc-1) + ipp
		    if (ipp == 1) then
			tmpa = acen_tavg_use(k,icc,jcls) - acen_prec_use(k,icc,jcls)
		    else
			tmpa = acen_prec_use(k,icc,jcls)
		    end if

		    if (tmpa > afrac_cut_0p5) then
			tmpq = (chem_tmpb(iccpp,k,jcls,l) - chem_tmpa(iccpp,k,jcls,l))
			tmpy = tmpy + tmpa*tmpq
                        del_cldchem3d(k,icc,jcls,ipp,l)=del_cldchem3d(k,icc,jcls,ipp,l)+tmpa*tmpq
                    else 
                        del_cldchem3d(k,icc,jcls,ipp,l)=del_cldchem3d(k,icc,jcls,ipp,l)+0.0
		    end if

                    if(rename_onoff_ecpp > 0 ) then
                      if (tmpa > afrac_cut_0p5) then
                        tmpq = (chem_tmpc(iccpp,k,jcls,l) - chem_tmpb(iccpp,k,jcls,l))
                        tmpy2 = tmpy2 + tmpa*tmpq
                        del_rename3d(k,icc,jcls,ipp,l)=del_rename3d(k,icc,jcls,ipp,l)+tmpa*tmpq
                      else
                        del_rename3d(k,icc,jcls,ipp,l)=del_rename3d(k,icc,jcls,ipp,l)+0.0
                      end if
                    end if  ! (rename_onoff_ecpp > 0.)

		end do ! ipp
		end do ! icc
		end do ! jcls
		tmpx = tmpx + tmpy*rhodz_cen(k)
                if(rename_onoff_ecpp > 0 ) tmpx2 = tmpx2+tmpy2 * rhodz_cen(k)
	    end do ! k

	    del_chem_clm_cldchem(l) = del_chem_clm_cldchem(l) + tmpx
            if(rename_onoff_ecpp > 0 )  & 
              del_chem_clm_rename(l) = del_chem_clm_rename(l) + tmpx2
	end do ! l

	end if ! (cldchem_onoff_ecpp > 0) 

	if ((cldchem_onoff_ecpp > 0)) then 

	do l = p1st, num_chem_ecpp
	do k = kts, ktecen
	do jcls = 1, ncls_use
	do icc = 1, 2
	    tmpa1 = acen_tavg_use(k,icc,jcls) - acen_prec_use(k,icc,jcls)
	    tmpa2 = acen_prec_use(k,icc,jcls)
	    if ((tmpa1 <= afrac_cut_0p5) .and. (tmpa2 <= afrac_cut_0p5)) cycle

	    iccpp1 = 2*(icc-1) + 1
	    iccpp2 = 2*(icc-1) + 2
            
            if(rename_onoff_ecpp > 0 ) then
              if ((tmpa1 > afrac_cut_0p5) .and. (tmpa2 > afrac_cut_0p5)) then
		tmpb1 = max( 0.0_r8, min( 1.0_r8, (tmpa1/(tmpa1+tmpa2)) ) )
		tmpb2 = 1.0_r8 - tmpb1
		tmpq1 = chem_tmpa(iccpp1,k,jcls,l)*tmpb1 &
		      + chem_tmpa(iccpp2,k,jcls,l)*tmpb2
		tmpq2 = chem_tmpc(iccpp1,k,jcls,l)*tmpb1 &
		      + chem_tmpc(iccpp2,k,jcls,l)*tmpb2
	      else if (tmpa1 > afrac_cut_0p5) then
		tmpq1 = chem_tmpa(iccpp1,k,jcls,l)
		tmpq2 = chem_tmpc(iccpp1,k,jcls,l)
	      else
		tmpq1 = chem_tmpa(iccpp2,k,jcls,l)
		tmpq2 = chem_tmpc(iccpp2,k,jcls,l)
	      end if
            else   ! no renaming
              if ((tmpa1 > afrac_cut_0p5) .and. (tmpa2 > afrac_cut_0p5)) then
                tmpb1 = max( 0.0_r8, min( 1.0_r8, (tmpa1/(tmpa1+tmpa2)) ) )
                tmpb2 = 1.0_r8 - tmpb1
                tmpq1 = chem_tmpa(iccpp1,k,jcls,l)*tmpb1 &
                      + chem_tmpa(iccpp2,k,jcls,l)*tmpb2
                tmpq2 = chem_tmpb(iccpp1,k,jcls,l)*tmpb1 &
                      + chem_tmpb(iccpp2,k,jcls,l)*tmpb2
              else if (tmpa1 > afrac_cut_0p5) then
                tmpq1 = chem_tmpa(iccpp1,k,jcls,l)
                tmpq2 = chem_tmpb(iccpp1,k,jcls,l)
              else
                tmpq1 = chem_tmpa(iccpp2,k,jcls,l)
                tmpq2 = chem_tmpb(iccpp2,k,jcls,l)
              end if
            end if ! (rename_onoff_ecpp > 0)
	    if (tmpq1 /= tmpq2) chem_sub_new(k,icc,jcls,l) = tmpq2

	end do ! icc
	end do ! jcls
	end do ! k
	end do ! l

	end if ! ((cldchem_onoff_ecpp > 0)) 


	deallocate ( p_tmp, t_tmp, rho_tmp, alt_tmp, &
	             cldfra_tmp, &
	             qlsink_tmp, &
	             precr_tmp, precs_tmp, precg_tmp, preci_tmp )
	deallocate ( moist_tmp, zero_tmp,&
	             chem_tmpa, chem_tmpb, chem_tmpc)
        deallocate(cwat_tmp, pdel_tmp,aqso4_h2o2_tmp,aqso4_o3_tmp)
        deallocate ( mmr, mmrcw, vmr, vmrcw, vmr_3d, vmrcw_3d,vmr_sv1, vmrcw_sv1, &
                    mbar, cldnum, invariants, mmr_3d, mmrcw_3d, mbar_3d,  &
                    qsrflx_full, qqcwsrflx_full)

	return
	end subroutine parampollu_tdx_cldchem

end module ecpp_modal_cloudchem
