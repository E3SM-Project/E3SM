module ebert_curry

   use shr_kind_mod,     only: r8 => shr_kind_r8
   use ppgrid,           only: pcols, pver
   use physics_types,    only: physics_state
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
   use radconstants,     only: nswbands, nlwbands, get_sw_spectral_boundaries

   implicit none
   private
   save

   public :: &
      ec_rad_props_init, &
      ec_ice_optics_sw,  &
      ec_ice_optics_lw


   real, public, parameter:: scalefactor = 1._r8 !500._r8/917._r8

   ! indexes into pbuf for optical parameters of MG clouds
   integer :: iciwp_idx   = 0 
   integer :: iclwp_idx   = 0 
   integer :: cld_idx     = 0  
   integer :: rei_idx     = 0  

contains

   !===========================================================================

   subroutine ec_rad_props_init()

      iciwp_idx  = pbuf_get_index('ICIWP')
      iclwp_idx  = pbuf_get_index('ICLWP')
      cld_idx    = pbuf_get_index('CLD')
      rei_idx    = pbuf_get_index('REI')

      return

   end subroutine ec_rad_props_init

   !===========================================================================

   subroutine ec_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)

      real(r8),intent(out) :: ice_tau    (nswbands,pcols,pver) ! extinction optical depth
      real(r8),intent(out) :: ice_tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
      real(r8),intent(out) :: ice_tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
      real(r8),intent(out) :: ice_tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

      real(r8), pointer, dimension(:,:) :: rei
      real(r8), pointer, dimension(:,:) :: cldn
      real(r8), pointer, dimension(:,:) :: tmpptr
      real(r8), dimension(pcols,pver) :: cicewp
      real(r8), dimension(nswbands) :: wavmin
      real(r8), dimension(nswbands) :: wavmax

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

      ! Optical properties for ice are valid only in the range of
      ! 13 < rei < 130 micron (Ebert and Curry 92)
      real(r8), parameter :: rei_min = 13._r8
      real(r8), parameter :: rei_max = 130._r8

      integer :: ns, i, k, indxsl, lchnk, Nday
      integer :: itim_old
      real(r8) :: tmp1i, tmp2i, tmp3i, g

      Nday = state%ncol
      lchnk = state%lchnk

      itim_old = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, cld_idx,cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, rei_idx,rei)
      call pbuf_get_field(pbuf, iciwp_idx, tmpptr)
      cicewp(1:pcols,1:pver) =  1000.0_r8*tmpptr(1:pcols,1:pver)
      
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
                  tmp1i = abarii + bbarii/max(rei_min,min(scalefactor*rei(i,k),rei_max))
                  ice_tau(ns,i,k) = cicewp(i,k)*tmp1i
               else
                  ice_tau(ns,i,k) = 0.0_r8
               endif

               tmp2i = 1._r8 - cbarii - dbarii*min(max(rei_min,scalefactor*rei(i,k)),rei_max)
               tmp3i = fbarii*min(max(rei_min,scalefactor*rei(i,k)),rei_max)
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

   !===========================================================================

   subroutine ec_ice_optics_lw(state, pbuf, abs_od)

      type(physics_state), intent(in)   :: state
      type(physics_buffer_desc),pointer :: pbuf(:)
      real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

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

      ! Optical properties for ice are valid only in the range of
      ! 13 < rei < 130 micron (Ebert and Curry 92)
      real(r8), parameter :: rei_min = 13._r8
      real(r8), parameter :: rei_max = 130._r8

      real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth


      ncol = state%ncol
      lchnk = state%lchnk

      itim_old  =  pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, rei_idx,   rei)
      call pbuf_get_field(pbuf, cld_idx,   cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
      call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
      do k=1,pver
         do i = 1,ncol
            cwp(i,k) = 1000.0_r8 *iciwpth(i,k) + 1000.0_r8 *iclwpth(i,k)
            ficemr(i,k) = 1000.0_r8*iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
         end do
      end do

      do k=1,pver
         do i=1,ncol
            ! Note from Andrew Conley:
            !  Optics for RK no longer supported, This is constructed to get
            !  close to bit for bit.  Otherwise we could simply use ice water path
            !note that optical properties for ice valid only
            !in range of 13 > rei > 130 micron (Ebert and Curry 92)
            kabsi = 0.005_r8 + 1._r8/min(max(rei_min,scalefactor*rei(i,k)),rei_max)
            kabs =  kabsi*ficemr(i,k) ! kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
            cldtau(i,k) = kabs*cwp(i,k)
         end do
      end do

      do lwband = 1,nlwbands
         abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
      enddo

   end subroutine ec_ice_optics_lw

   !===========================================================================

end module ebert_curry
