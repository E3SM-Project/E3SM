module oldcloud

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use radconstants,     only: nswbands, nlwbands, idx_sw_diag, ot_length, idx_lw_diag, get_sw_spectral_boundaries
use abortutils,       only: endrun
use cam_history,      only: outfld
use rad_constituents, only: iceopticsfile, liqopticsfile
use ebert_curry,      only: scalefactor

implicit none
private
save

public :: &
 oldcloud_init, oldcloud_lw, old_liq_get_rad_props_lw, old_ice_get_rad_props_lw

integer :: nmu, nlambda
real(r8), allocatable :: g_mu(:)           ! mu samples on grid
real(r8), allocatable :: g_lambda(:,:)     ! lambda scale samples on grid
real(r8), allocatable :: ext_sw_liq(:,:,:)
real(r8), allocatable :: ssa_sw_liq(:,:,:)
real(r8), allocatable :: asm_sw_liq(:,:,:)
real(r8), allocatable :: abs_lw_liq(:,:,:)

integer :: n_g_d
real(r8), allocatable :: g_d_eff(:)        ! radiative effective diameter samples on grid
real(r8), allocatable :: ext_sw_ice(:,:)
real(r8), allocatable :: ssa_sw_ice(:,:)
real(r8), allocatable :: asm_sw_ice(:,:)
real(r8), allocatable :: abs_lw_ice(:,:)

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
   integer ::   iciwp_idx   = 0 
   integer ::   iclwp_idx   = 0 
   integer ::   cld_idx     = 0 
   integer ::   rel_idx     = 0 
   integer ::   rei_idx     = 0 

! indexes into constituents for old optics
   integer :: &
        ixcldice,           & ! cloud ice water index
        ixcldliq              ! cloud liquid water index


!==============================================================================
contains
!==============================================================================

subroutine oldcloud_init()

   use constituents,   only: cnst_get_ind

   integer :: err

   iciwp_idx  = pbuf_get_index('ICIWP',errcode=err)
   iclwp_idx  = pbuf_get_index('ICLWP',errcode=err)
   cld_idx    = pbuf_get_index('CLD')
   rel_idx    = pbuf_get_index('REL')
   rei_idx    = pbuf_get_index('REI')

   ! old optics
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)

   return

end subroutine oldcloud_init

!==============================================================================
! Private methods
!==============================================================================

subroutine old_liquid_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp)

   use physconst, only: gravit

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

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
   integer :: lchnk, itim_old
   real(r8) :: tmp1l, tmp2l, tmp3l, g
   real(r8) :: kext(pcols,pver)
   real(r8), pointer, dimension(:,:) :: iclwpth

   Nday = state%ncol
   lchnk = state%lchnk

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, rel_idx,rel)

   if (oldliqwp) then
     do k=1,pver
        do i = 1,Nday
           cliqwp(i,k) = state%q(i,k,ixcldliq)*state%pdel(i,k)/(gravit*max(0.01_r8,cldn(i,k)))
        end do
     end do
   else
     if (iclwp_idx<0) then 
        call endrun('old_liquid_optics_sw: oldliqwp must be set to true since ICLWP was not found in pbuf')
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
      else if(wavmin(ns) > 2.38_r8) then
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


end subroutine old_liquid_optics_sw
!==============================================================================

subroutine old_ice_optics_sw   (state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp)

   use physconst, only: gravit

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

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
        call endrun('old_ice_optics_sw: oldicewp must be set to true since ICIWP was not found in pbuf')
     endif
     call pbuf_get_field(pbuf, iciwp_idx, tmpptr)
     cicewp(1:pcols,1:pver) =  1000.0_r8*tmpptr
   endif
   
   call get_sw_spectral_boundaries(wavmin,wavmax,'microns')

   do ns = 1, nswbands

      if(wavmax(ns) <= 0.7_r8) then
         indxsl = 1
      else if(wavmax(ns) <= 1.25_r8) then
         indxsl = 2
      else if(wavmax(ns) <= 2.38_r8) then
         indxsl = 3
      else if(wavmin(ns) > 2.38_r8) then
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

end subroutine old_ice_optics_sw
!==============================================================================

subroutine oldcloud_lw(state,pbuf,cld_abs_od,oldwp)
   use physconst, only: gravit
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc),pointer :: pbuf(:)
   real(r8),            intent(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
   logical,intent(in)              :: oldwp ! use old definition of waterpath


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
   real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth

    real(r8) :: kabs, kabsi
    real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
    parameter (kabsl = 0.090361_r8)


 
   ncol = state%ncol
   lchnk = state%lchnk

   itim_old  =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, rei_idx,   rei)
   call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   if (oldwp) then
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
         call endrun('oldcloud_lw: oldwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
      endif
      call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
      call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
      do k=1,pver
         do i = 1,ncol
            cwp(i,k) = 1000.0_r8 *iclwpth(i,k) + 1000.0_r8 *iciwpth(i, k)
            ficemr(i,k) = 1000.0_r8 * iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
         end do
      end do
   endif

   do k=1,pver
       do i=1,ncol

          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs = kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
          !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       end do
   end do
!
   do lwband = 1,nlwbands
      cld_abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
   enddo

end subroutine oldcloud_lw

!==============================================================================
subroutine old_liq_get_rad_props_lw(state, pbuf, abs_od, oldliqwp)
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
   integer :: ncol, itim_old, lwband, i, k, lchnk 

   real(r8) :: kabs, kabsi
   real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
   parameter (kabsl = 0.090361_r8)

   real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth

   ncol=state%ncol
   lchnk = state%lchnk

   itim_old  =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, rei_idx,   rei)
   call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

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
         call endrun('old_liq_get_rad_props_lw: oldliqwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
      endif
      call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
      call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
      do k=1,pver
         do i = 1,ncol
            cwp(i,k) = 1000.0_r8 *iclwpth(i,k) + 1000.0_r8 *iciwpth(i, k)
            ficemr(i,k) = 1000.0 * iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
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
          kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs = kabsl*(1._r8-ficemr(i,k)) ! + kabsi*ficemr(i,k)
          !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       end do
   end do
!
   do lwband = 1,nlwbands
      abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
   enddo


end subroutine old_liq_get_rad_props_lw
!==============================================================================

subroutine old_ice_get_rad_props_lw(state, pbuf, abs_od, oldicewp)
   use physconst, only: gravit
   type(physics_state), intent(in)  :: state
   type(physics_buffer_desc),pointer  :: pbuf(:)
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
   call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

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
         call endrun('old_ice_get_rad_props_lw: oldicewp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
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

end subroutine old_ice_get_rad_props_lw
!==============================================================================

subroutine cloud_total_vis_diag_out(lchnk, nnite, idxnite, tau, radsuffix)

   ! output total aerosol optical depth for the visible band

   use cam_history, only: outfld
   use cam_history_support, only : fillvalue

   integer,          intent(in) :: lchnk
   integer,          intent(in) :: nnite          ! number of night columns
   integer,          intent(in) :: idxnite(nnite) ! local column indices of night columns
   real(r8),         intent(in) :: tau(:,:)
   character(len=*), intent(in) :: radsuffix ! identifies whether the radiation call
                                             ! is for the climate calc or a diagnostic calc
 
   ! Local variables
   integer  :: i
   real(r8) :: tmp(pcols)
   !-----------------------------------------------------------------------------

   ! compute total aerosol optical depth output where only daylight columns
   tmp(:) = sum(tau(:,:), 2)
   do i = 1, nnite
      tmp(idxnite(i)) = fillvalue
   end do
   !call outfld('cloudOD_v'//trim(radsuffix), tmp, pcols, lchnk)

end subroutine cloud_total_vis_diag_out

!==============================================================================

end module oldcloud
