
module ebert_curry

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
use radconstants,     only: nswbands, nlwbands, idx_sw_diag, ot_length, idx_lw_diag, get_sw_spectral_boundaries
use cam_abortutils,       only: endrun
use cam_history,      only: outfld

implicit none
private
save

public :: &
   ec_rad_props_init,        &
   cloud_rad_props_get_sw, & ! return SW optical props of total bulk aerosols
   cloud_rad_props_get_lw, & ! return LW optical props of total bulk aerosols
   ec_ice_optics_sw,       &
   ec_ice_get_rad_props_lw


real, public, parameter:: scalefactor = 1._r8 !500._r8/917._r8

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
   integer :: dei_idx     = 0 
   integer :: mu_idx      = 0
   integer :: lambda_idx  = 0 
   integer :: iciwp_idx   = 0 
   integer :: iclwp_idx   = 0 
   integer :: cld_idx     = 0  
   integer :: rei_idx     = 0  

! indexes into constituents for old optics
   integer :: &
        ixcldice,           & ! cloud ice water index
        ixcldliq              ! cloud liquid water index


!==============================================================================
contains
!==============================================================================

subroutine ec_rad_props_init()

   use cam_history, only: addfld
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
   rei_idx    = pbuf_get_index('REI')

   ! old optics
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)

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

end subroutine ec_rad_props_init

!==============================================================================

subroutine cloud_rad_props_get_sw(state, pbuf, &
                                  tau, tau_w, tau_w_g, tau_w_f,&
                                  diagnosticindex, oldliq, oldice)

! return totaled (across all species) layer tau, omega, g, f 
! for all spectral interval for aerosols affecting the climate

   ! Arguments
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc),  pointer :: pbuf(:)
   integer, optional,   intent(in) :: diagnosticindex      ! index (if present) to radiation diagnostic information

   real(r8), intent(out) :: tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
   real(r8), intent(out) :: tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
   real(r8), intent(out) :: tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
   real(r8), intent(out) :: tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w

   logical, optional, intent(in) :: oldliq,oldice

   ! Local variables

   integer :: ncol
   integer :: lchnk
   integer :: k, i    ! lev and daycolumn indices
   integer :: iswband ! sw band indices

   !-----------------------------------------------------------------------------

   ncol  = state%ncol
   lchnk = state%lchnk

   ! initialize to conditions that would cause failure
   tau     (:,:,:) = -100._r8
   tau_w   (:,:,:) = -100._r8
   tau_w_g (:,:,:) = -100._r8
   tau_w_f (:,:,:) = -100._r8

   ! initialize layers to accumulate od's
   tau    (:,1:ncol,:) = 0._r8
   tau_w  (:,1:ncol,:) = 0._r8
   tau_w_g(:,1:ncol,:) = 0._r8
   tau_w_f(:,1:ncol,:) = 0._r8


   call ec_ice_optics_sw   (state, pbuf, tau, tau_w, tau_w_g, tau_w_f, oldicewp=.true.)
!  call outfld ('CI_OD_SW_OLD', ice_tau(idx_sw_diag,:,:), pcols, lchnk)


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

   !-----------------------------------------------------------------------------

   ncol = state%ncol
   lchnk = state%lchnk

   ! compute optical depths cld_absod 
   cld_abs_od = 0._r8

   call ec_ice_get_rad_props_lw(state, pbuf, cld_abs_od, oldicewp=.true.)
   !call outfld('CI_OD_LW_OLD', ice_tau_abs_od(idx_lw_diag ,:,:), pcols, lchnk)
      
end subroutine cloud_rad_props_get_lw

!==============================================================================
! Private methods
!==============================================================================

subroutine ec_ice_optics_sw   (state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp)

   use physconst,      only: gravit

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8),intent(out) :: ice_tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: ice_tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: ice_tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
   real(r8),intent(out) :: ice_tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w
   logical, intent(in) :: oldicewp

   real(r8), pointer, dimension(:,:) :: rei
   real(r8), pointer, dimension(:,:) :: cldn
   real(r8), pointer, dimension(:,:) :: tmpptr
   real(r8), dimension(pcols,pver) :: cicewp
   real(r8), dimension(nswbands) :: wavmin
   real(r8), dimension(nswbands) :: wavmax
   !
   ! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
   real(r8) :: abari(4) = &     ! a coefficient for extinction optical depth
      (/ 3.448e-03_r8, 3.448e-03_r8,3.448e-03_r8,3.448e-03_r8/)
   real(r8) :: bbari(4) = &     ! b coefficient for extinction optical depth
      (/ 2.431_r8    , 2.431_r8    ,2.431_r8    ,2.431_r8    /)
   real(r8) :: cbari(4) = &     ! c coefficient for single scat albedo
      (/ 1.00e-05_r8 , 1.10e-04_r8 ,1.861e-02_r8,.46658_r8   /)
   real(r8) :: dbari(4) = &     ! d coefficient for single scat albedo
      (/ 0.0_r8      , 1.405e-05_r8,8.328e-04_r8,2.05e-05_r8 /)
   real(r8) :: ebari(4) = &     ! e coefficient for asymmetry parameter
      (/ 0.7661_r8   , 0.7730_r8   ,0.794_r8    ,0.9595_r8   /)
   real(r8) :: fbari(4) = &     ! f coefficient for asymmetry parameter
      (/ 5.851e-04_r8, 5.665e-04_r8,7.267e-04_r8,1.076e-04_r8/)

   real(r8) :: abarii           ! A coefficient for current spectral band
   real(r8) :: bbarii           ! B coefficient for current spectral band
   real(r8) :: cbarii           ! C coefficient for current spectral band
   real(r8) :: dbarii           ! D coefficient for current spectral band
   real(r8) :: ebarii           ! E coefficient for current spectral band
   real(r8) :: fbarii           ! F coefficient for current spectral band

   ! Minimum cloud amount (as a fraction of the grid-box area) to 
   ! distinguish from clear sky
   real(r8), parameter :: cldmin = 1.0e-80_r8

   ! Decimal precision of cloud amount (0 -> preserve full resolution;
   ! 10^-n -> preserve n digits of cloud amount)
   real(r8), parameter :: cldeps = 0.0_r8

   integer :: ns, i, k, indxsl, lchnk, Nday
   integer :: itim_old
   real(r8) :: tmp1i, tmp2i, tmp3i, g

   Nday = state%ncol
   lchnk = state%lchnk

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, rei_idx,rei)

   if(oldicewp) then
     do k=1,pver
        do i = 1,Nday
           cicewp(i,k) = 1000.0_r8*state%q(i,k,ixcldice)*state%pdel(i,k) /(gravit* max(0.01_r8,cldn(i,k)))
        end do
     end do
   else
     if (iciwp_idx<=0) then 
        call endrun('ec_ice_optics_sw: oldicewp must be set to true since ICIWP was not found in pbuf')
     endif
     call pbuf_get_field(pbuf, iciwp_idx, tmpptr)
     cicewp(1:pcols,1:pver) =  1000.0_r8*tmpptr(1:pcols,1:pver)
   endif
   
   call get_sw_spectral_boundaries(wavmin,wavmax,'microns')

   do ns = 1, nswbands

      if(wavmax(ns) <= 0.7_r8) then
         indxsl = 1
      else if(wavmax(ns) <= 1.25_r8) then
         indxsl = 2
      else if(wavmax(ns) <= 2.38_r8) then
         indxsl = 3
      else if(wavmax(ns) > 2.38_r8) then
         indxsl = 4
      end if

      abarii = abari(indxsl)
      bbarii = bbari(indxsl)
      cbarii = cbari(indxsl)
      dbarii = dbari(indxsl)
      ebarii = ebari(indxsl)
      fbarii = fbari(indxsl)

      do k=1,pver
         do i=1,Nday

            ! note that optical properties for ice valid only
            ! in range of 13 > rei > 130 micron (Ebert and Curry 92)
            if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
               tmp1i = abarii + bbarii/max(13._r8,min(scalefactor*rei(i,k),130._r8))
               ice_tau(ns,i,k) = cicewp(i,k)*tmp1i
            else
               ice_tau(ns,i,k) = 0.0_r8
            endif

            tmp2i = 1._r8 - cbarii - dbarii*min(max(13._r8,scalefactor*rei(i,k)),130._r8)
            tmp3i = fbarii*min(max(13._r8,scalefactor*rei(i,k)),130._r8)
            ! Do not let single scatter albedo be 1.  Delta-eddington solution
            ! for non-conservative case has different analytic form from solution
            ! for conservative case, and raddedmx is written for non-conservative case.
            ice_tau_w(ns,i,k) = ice_tau(ns,i,k) * min(tmp2i,.999999_r8)
            g = ebarii + tmp3i
            ice_tau_w_g(ns,i,k) = ice_tau_w(ns,i,k) * g
            ice_tau_w_f(ns,i,k) = ice_tau_w(ns,i,k) * g * g

         end do ! End do i=1,Nday
      end do    ! End do k=1,pver
   end do ! nswbands

end subroutine ec_ice_optics_sw
!==============================================================================

subroutine ec_ice_get_rad_props_lw(state, pbuf, abs_od, oldicewp)
   use physconst,      only: gravit

   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)
    real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)
   logical, intent(in) :: oldicewp

   real(r8) :: gicewp(pcols,pver)
   real(r8) :: gliqwp(pcols,pver)
   real(r8) :: cicewp(pcols,pver)
   real(r8) :: cliqwp(pcols,pver)
   real(r8) :: ficemr(pcols,pver)
   real(r8) :: cwp(pcols,pver)
   real(r8) :: cldtau(pcols,pver)

   real(r8), pointer, dimension(:,:) :: cldn
   real(r8), pointer, dimension(:,:) :: rei
   integer :: ncol, itim_old, lwband, i, k, lchnk

    real(r8) :: kabs, kabsi

    real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
    parameter (kabsl = 0.090361_r8)

   real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth


   ncol = state%ncol
   lchnk = state%lchnk

   itim_old  =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, rei_idx,   rei)
   call pbuf_get_field(pbuf, cld_idx,   cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))


   if(oldicewp) then
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
         call endrun('ec_ice_get_rad_props_lw: oldicewp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
      endif
      call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
      call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
      do k=1,pver
         do i = 1,ncol
            cwp(i,k) = 1000.0_r8 *iciwpth(i,k) + 1000.0_r8 *iclwpth(i,k)
            ficemr(i,k) = 1000.0_r8*iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
         end do
      end do
   endif

   do k=1,pver
       do i=1,ncol

          ! Note from Andrew Conley:
          !  Optics for RK no longer supported, This is constructed to get
          !  close to bit for bit.  Otherwise we could simply use ice water path
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs =  kabsi*ficemr(i,k) ! kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
          !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       end do
   end do
!
   do lwband = 1,nlwbands
      abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
   enddo

   !if(oldicewp) then
   !  call outfld('CIWPTH_OLD',cicewp(:,:)/1000,pcols,lchnk)
   !else
   !  call outfld('CIWPTH_OLD',iciwpth(:,:),pcols,lchnk)
   !endif
   !call outfld('CI_OD_LW_OLD',cldtau(:,:),pcols,lchnk)

end subroutine ec_ice_get_rad_props_lw
!==============================================================================

end module ebert_curry
