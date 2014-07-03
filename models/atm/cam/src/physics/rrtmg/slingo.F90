module slingo

!------------------------------------------------------------------------------------------------
!  Implements Slingo Optics for MG/RRTMG for liquid clouds and
!  a copy of the old cloud routine for reference
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
use radconstants,     only: nswbands, nlwbands, idx_sw_diag, ot_length, idx_lw_diag, get_sw_spectral_boundaries
use abortutils,       only: endrun
use cam_history,      only: outfld

implicit none
private
save

public :: &
   slingo_rad_props_init,        &
   cloud_rad_props_get_sw, & ! return SW optical props of total bulk aerosols
   cloud_rad_props_get_lw,  & ! return LW optical props of total bulk aerosols
   slingo_liq_get_rad_props_lw, &
   slingo_liq_optics_sw

! Minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
! 
   real(r8) cldmin
   parameter (cldmin = 1.0e-80_r8)
!
! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
! 
   real(r8) cldeps
   parameter (cldeps = 0.0_r8)

! 
! indexes into pbuf for optical parameters of MG clouds
! 
   integer :: iclwp_idx  = 0 
   integer :: iciwp_idx  = 0
   integer :: cld_idx    = 0 
   integer :: rel_idx  = 0
   integer :: rei_idx  = 0

! indexes into constituents for old optics
   integer :: &
        ixcldliq,   &         ! cloud liquid water index
        ixcldice              ! cloud liquid water index


!==============================================================================
contains
!==============================================================================

subroutine slingo_rad_props_init()

   use cam_history, only: addfld, phys_decomp
   use netcdf
   use spmd_utils,     only: masterproc
   use ioFileMod,      only: getfil
   use cam_logfile,    only: iulog
   use error_messages, only: handle_ncerr
#if ( defined SPMD )
   use mpishorthand
#endif
   use constituents,   only: cnst_get_ind

   integer :: err

   iciwp_idx  = pbuf_get_index('ICIWP',errcode=err)
   iclwp_idx  = pbuf_get_index('ICLWP',errcode=err)
   cld_idx    = pbuf_get_index('CLD')
   rel_idx    = pbuf_get_index('REL')
   rei_idx    = pbuf_get_index('REI')

   ! old optics
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   !call addfld ('CLWPTH_OLD','Kg/m2   ',pver, 'I','old In Cloud Liquid Water Path',phys_decomp, sampling_seq='rad_lwsw')
   !call addfld ('KEXT_OLD','m^2/kg',pver,'I','old extinction',phys_decomp)
   !call addfld ('CLDOD_OLD','1',pver,'I','old liquid OD',phys_decomp)
   !call addfld ('REL_OLD','1',pver,'I','old liquid effective radius (liquid)',phys_decomp)

   !call addfld ('CLWPTH_NEW','Kg/m2   ',pver, 'I','In Cloud Liquid Water Path',phys_decomp, sampling_seq='rad_lwsw')
   !call addfld ('KEXT_NEW','m^2/kg',pver,'I','extinction',phys_decomp)
   !call addfld ('CLDOD_NEW','1',pver,'I','liquid OD',phys_decomp)

   !call addfld('CIWPTH_NEW','Kg/m2   ',pver, 'I','In Cloud Ice Water Path',phys_decomp, sampling_seq='rad_lwsw')
   !call addfld('CIWPTH_OLD','Kg/m2   ',pver, 'I','In Cloud Ice Water Path (old)',phys_decomp, sampling_seq='rad_lwsw')

   return

end subroutine slingo_rad_props_init

!==============================================================================

subroutine cloud_rad_props_get_sw(state, pbuf, &
                                  tau, tau_w, tau_w_g, tau_w_f,&
                                  diagnosticindex)

! return totaled (across all species) layer tau, omega, g, f 
! for all spectral interval for aerosols affecting the climate

   ! Arguments
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, optional,   intent(in) :: diagnosticindex      ! index (if present) to radiation diagnostic information

   real(r8), intent(out) :: tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
   real(r8), intent(out) :: tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
   real(r8), intent(out) :: tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
   real(r8), intent(out) :: tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w

   ! Local variables

   integer :: ncol
   integer :: lchnk
   integer :: k, i    ! lev and daycolumn indices
   integer :: iswband ! sw band indices

   real(r8) :: liq_tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
   real(r8) :: liq_tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
   real(r8) :: liq_tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
   real(r8) :: liq_tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w


   !-----------------------------------------------------------------------------

   ncol  = state%ncol
   lchnk = state%lchnk

   call slingo_liq_optics_sw(state, pbuf, tau, tau_w, tau_w_g, tau_w_f, oldliqwp=.true. )

end subroutine cloud_rad_props_get_sw
!==============================================================================

subroutine cloud_rad_props_get_lw(state, pbuf, cld_abs_od, diagnosticindex, oldliq, oldice, oldcloud)

! Purpose: Compute cloud longwave absorption optical depth
!    cloud_rad_props_get_lw() is called by radlw() 

   ! Arguments
   type(physics_state), intent(in)  :: state
   type(physics_buffer_desc), pointer  :: pbuf(:)
   real(r8),            intent(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
   integer, optional,   intent(in)  :: diagnosticindex
   logical, optional,   intent(in)  :: oldliq  ! use old liquid optics
   logical, optional,   intent(in)  :: oldice  ! use old ice optics
   logical, optional,   intent(in)  :: oldcloud  ! use old optics for both (b4b)

   ! Local variables

   integer :: bnd_idx     ! LW band index
   integer :: i           ! column index
   integer :: k           ! lev index
   integer :: ncol        ! number of columns
   integer :: lchnk

   ! rad properties for liquid clouds
   real(r8) :: liq_tau_abs_od(nlwbands,pcols,pver) ! liquid cloud absorption optical depth

   !-----------------------------------------------------------------------------

   ncol = state%ncol
   lchnk = state%lchnk

   ! compute optical depths cld_absod 
   cld_abs_od = 0._r8

   call slingo_liq_get_rad_props_lw(state, pbuf, liq_tau_abs_od, oldliqwp=.true.)
      
   cld_abs_od(:,1:ncol,:) = liq_tau_abs_od(:,1:ncol,:) 

end subroutine cloud_rad_props_get_lw

!==============================================================================
! Private methods
!==============================================================================


subroutine slingo_liq_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp)

   use physconst, only: gravit

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8),intent(out) :: liq_tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: liq_tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: liq_tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
   real(r8),intent(out) :: liq_tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w
   logical, intent(in) :: oldliqwp

   real(r8), pointer, dimension(:,:) :: rel
   real(r8), pointer, dimension(:,:) :: cldn
   real(r8), pointer, dimension(:,:) :: tmpptr
   real(r8), dimension(pcols,pver) :: cliqwp
   real(r8), dimension(nswbands) :: wavmin
   real(r8), dimension(nswbands) :: wavmax

   ! Minimum cloud amount (as a fraction of the grid-box area) to 
   ! distinguish from clear sky
   real(r8), parameter :: cldmin = 1.0e-80_r8

   ! Decimal precision of cloud amount (0 -> preserve full resolution;
   ! 10^-n -> preserve n digits of cloud amount)
   real(r8), parameter :: cldeps = 0.0_r8

   ! A. Slingo's data for cloud particle radiative properties (from 'A GCM
   ! Parameterization for the Shortwave Properties of Water Clouds' JAS
   ! vol. 46 may 1989 pp 1419-1427)
   real(r8) :: abarl(4) = &  ! A coefficient for extinction optical depth
      (/ 2.817e-02_r8, 2.682e-02_r8,2.264e-02_r8,1.281e-02_r8/)
   real(r8) :: bbarl(4) = &  ! B coefficient for extinction optical depth
      (/ 1.305_r8    , 1.346_r8    ,1.454_r8    ,1.641_r8    /)
   real(r8) :: cbarl(4) = &  ! C coefficient for single scat albedo
      (/-5.62e-08_r8 ,-6.94e-06_r8 ,4.64e-04_r8 ,0.201_r8    /)
   real(r8) :: dbarl(4) = &  ! D coefficient for single  scat albedo
      (/ 1.63e-07_r8 , 2.35e-05_r8 ,1.24e-03_r8 ,7.56e-03_r8 /)
   real(r8) :: ebarl(4) = &  ! E coefficient for asymmetry parameter
      (/ 0.829_r8    , 0.794_r8    ,0.754_r8    ,0.826_r8    /)
   real(r8) :: fbarl(4) = &  ! F coefficient for asymmetry parameter
      (/ 2.482e-03_r8, 4.226e-03_r8,6.560e-03_r8,4.353e-03_r8/)

   real(r8) :: abarli        ! A coefficient for current spectral band
   real(r8) :: bbarli        ! B coefficient for current spectral band
   real(r8) :: cbarli        ! C coefficient for current spectral band
   real(r8) :: dbarli        ! D coefficient for current spectral band
   real(r8) :: ebarli        ! E coefficient for current spectral band
   real(r8) :: fbarli        ! F coefficient for current spectral band

   ! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
   ! greater than 20 micro-meters

   integer :: ns, i, k, indxsl, Nday
   integer :: i_rel, lchnk, icld, itim
   real(r8) :: tmp1l, tmp2l, tmp3l, g
   real(r8) :: kext(pcols,pver)
   real(r8), pointer, dimension(:,:) :: iclwpth

   Nday = state%ncol
   lchnk = state%lchnk

   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, rel_idx, rel)

   if (oldliqwp) then
     do k=1,pver
        do i = 1,Nday
           cliqwp(i,k) = state%q(i,k,ixcldliq)*state%pdel(i,k)/(gravit*max(0.01_r8,cldn(i,k)))
        end do
     end do
   else
     if (iclwp_idx<=0) then 
        call endrun('slingo_liq_optics_sw: oldliqwp must be set to true since ICLWP was not found in pbuf')
     endif
     ! The following is the eventual target specification for in cloud liquid water path.
     call pbuf_get_field(pbuf, iclwp_idx, tmpptr)
     cliqwp = tmpptr
   endif
   
   call get_sw_spectral_boundaries(wavmin,wavmax,'microns')
  
   do ns = 1, nswbands
      ! Set index for cloud particle properties based on the wavelength,
      ! according to A. Slingo (1989) equations 1-3:
      ! Use index 1 (0.25 to 0.69 micrometers) for visible
      ! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
      ! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
      ! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
      if(wavmax(ns) <= 0.7_r8) then
         indxsl = 1
      else if(wavmax(ns) <= 1.25_r8) then
         indxsl = 2
      else if(wavmax(ns) <= 2.38_r8) then
         indxsl = 3
      else if(wavmax(ns) > 2.38_r8) then
         indxsl = 4
      end if

      ! Set cloud extinction optical depth, single scatter albedo,
      ! asymmetry parameter, and forward scattered fraction:
      abarli = abarl(indxsl)
      bbarli = bbarl(indxsl)
      cbarli = cbarl(indxsl)
      dbarli = dbarl(indxsl)
      ebarli = ebarl(indxsl)
      fbarli = fbarl(indxsl)

      do k=1,pver
         do i=1,Nday

            ! note that optical properties for liquid valid only
            ! in range of 4.2 > rel > 16 micron (Slingo 89)
            if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
               tmp1l = abarli + bbarli/min(max(4.2_r8,rel(i,k)),16._r8)
               liq_tau(ns,i,k) = 1000._r8*cliqwp(i,k)*tmp1l
            else
               liq_tau(ns,i,k) = 0.0_r8
            endif

            tmp2l = 1._r8 - cbarli - dbarli*min(max(4.2_r8,rel(i,k)),16._r8)
            tmp3l = fbarli*min(max(4.2_r8,rel(i,k)),16._r8)
            ! Do not let single scatter albedo be 1.  Delta-eddington solution
            ! for non-conservative case has different analytic form from solution
            ! for conservative case, and raddedmx is written for non-conservative case.
            liq_tau_w(ns,i,k) = liq_tau(ns,i,k) * min(tmp2l,.999999_r8)
            g = ebarli + tmp3l
            liq_tau_w_g(ns,i,k) = liq_tau_w(ns,i,k) * g
            liq_tau_w_f(ns,i,k) = liq_tau_w(ns,i,k) * g * g

         end do ! End do i=1,Nday
      end do    ! End do k=1,pver
   end do ! nswbands

   !call outfld('CL_OD_SW_OLD',liq_tau(idx_sw_diag,:,:), pcols, lchnk)
   !call outfld('REL_OLD',rel(:,:), pcols, lchnk)
   !call outfld('CLWPTH_OLD',cliqwp(:,:), pcols, lchnk)
   !call outfld('KEXT_OLD',kext(:,:), pcols, lchnk)


end subroutine slingo_liq_optics_sw

subroutine slingo_liq_get_rad_props_lw(state, pbuf, abs_od, oldliqwp)
   use physconst, only: gravit

   type(physics_state), intent(in)  :: state
   type(physics_buffer_desc),pointer  :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)
   logical, intent(in) :: oldliqwp

   real(r8) :: gicewp(pcols,pver)
   real(r8) :: gliqwp(pcols,pver)
   real(r8) :: cicewp(pcols,pver)
   real(r8) :: cliqwp(pcols,pver)
   real(r8) :: ficemr(pcols,pver)
   real(r8) :: cwp(pcols,pver)
   real(r8) :: cldtau(pcols,pver)

   real(r8), pointer, dimension(:,:) :: cldn
   real(r8), pointer, dimension(:,:) :: rei
   integer :: ncol, icld, itim, i_rei, lwband, i, k, lchnk 

    real(r8) :: kabs, kabsi
    real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
    parameter (kabsl = 0.090361_r8)

   real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth

    ncol=state%ncol
   lchnk = state%lchnk

   itim  =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, rei_idx,   rei)
   call pbuf_get_field(pbuf, cld_idx,   cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/))

   if (oldliqwp) then
     do k=1,pver
         do i = 1,ncol
            gicewp(i,k) = state%q(i,k,ixcldice)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box ice water path.
            gliqwp(i,k) = state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box liquid water path.
            cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
            cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
            ficemr(i,k) = state%q(i,k,ixcldice) /                 &
                 max(1.e-10_r8,(state%q(i,k,ixcldice)+state%q(i,k,ixcldliq)))
         end do
     end do
     cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
   else
     if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
        call endrun('slingo_liq_get_rad_props_lw: oldliqwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
     endif
     call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
     call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
     do k=1,pver
         do i = 1,ncol
              cwp   (i,k) = 1000.0_r8 * iclwpth(i,k) + 1000.0_r8 * iciwpth(i, k)
              ficemr(i,k) = 1000.0_r8 * iciwpth(i,k)/(max(1.e-18_r8, cwp(i,k)))
         end do
     end do
   endif


   do k=1,pver
       do i=1,ncol

          ! Note from Andrew Conley:
          !  Optics for RK no longer supported, This is constructed to get
          !  close to bit for bit.  Otherwise we could simply use liquid water path
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabs = kabsl*(1._r8-ficemr(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       end do
   end do
!
   do lwband = 1,nlwbands
      abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
   enddo


end subroutine slingo_liq_get_rad_props_lw

end module slingo
