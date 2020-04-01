module cam_optics

   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range
   use radiation_state, only: nlev_rad, ktop, kbot
   use radiation_utils, only: handle_error
   use radconstants, only: nswbands, nlwbands
   use rad_constituents, only: liqcldoptics, icecldoptics

   implicit none
   private

   public ::                  &
      do_snow_optics,         &
      free_optics_sw,         &
      free_optics_lw,         &
      set_cloud_optics_sw,    &
      set_cloud_optics_lw,    &
      get_cloud_optics_sw,    &
      get_cloud_optics_lw,    &
      sample_cloud_optics_sw, &
      sample_cloud_optics_lw, &
      set_aerosol_optics_sw,  &
      set_aerosol_optics_lw

   ! Mapping from old RRTMG sw bands to new band ordering in RRTMGP
   integer, dimension(14) :: map_rrtmg_to_rrtmgp_swbands = (/ &
      14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 &
   /)
   integer, dimension(14) :: map_rrtmgp_to_rrtmg_swbands = (/ &
      2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1 &
   /)

contains

   subroutine get_cloud_optics_sw( &
         ncol, nlev, nbnd, do_snow, cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rel, rei, &
         tau_out, ssa_out, asm_out)

      use ppgrid, only: pcols
      use cam_abortutils, only: endrun
      use physics_types, only: physics_state
      use cloud_rad_props, only: mitchell_ice_optics_sw, gammadist_liq_optics_sw
      use ebert_curry, only: ec_ice_optics_sw
      use slingo, only: slingo_liq_optics_sw

      integer, intent(in) :: ncol, nlev, nbnd
      logical, intent(in) :: do_snow
      real(r8), intent(in), dimension(:,:) :: &
         cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rel, rei

      ! Outputs are shortwave cloud optical properties *by band*. Dimensions should
      ! be nswbands,ncol,nlev. Generally, this should be able to handle cases were
      ! ncol might be something like nday, and nlev could be arbitrary so long as
      ! corresponding fields were defined for all indices of nlev. This
      ! isn't the case right now I don't think, as cloud_rad_props makes explicit
      ! assumptions about array sizes.
      real(r8), intent(out), dimension(:,:,:) :: tau_out, ssa_out, asm_out

      ! Temporary variables to hold cloud optical properties before combining into
      ! output arrays. Same shape as output arrays, so get shapes from output.
      real(r8), dimension(nbnd,pcols,nlev) :: &
            liq_tau, liq_tau_ssa, liq_tau_ssa_g, liq_tau_ssa_f, &
            ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f, &
            cld_tau, cld_tau_ssa, cld_tau_ssa_g, cld_tau_ssa_f, &
            snow_tau, snow_tau_ssa, snow_tau_ssa_g, snow_tau_ssa_f, &
            combined_tau, combined_tau_ssa, combined_tau_ssa_g, combined_tau_ssa_f

      integer :: iband, ilev, icol

      ! Initialize
      ice_tau = 0
      ice_tau_ssa = 0
      ice_tau_ssa_g = 0
      ice_tau_ssa_f = 0
      liq_tau = 0
      liq_tau_ssa = 0
      liq_tau_ssa_g = 0
      liq_tau_ssa_f = 0
      snow_tau = 0
      snow_tau_ssa = 0
      snow_tau_ssa_g = 0
      snow_tau_ssa_f = 0
      combined_tau = 0
      combined_tau_ssa = 0
      combined_tau_ssa_g = 0
      combined_tau_ssa_f = 0

      ! Get ice cloud optics
      if (trim(icecldoptics) == 'mitchell') then
         call mitchell_ice_optics_sw( &
            ncol, nlev, iciwp, dei, &
            ice_tau, ice_tau_ssa, &
            ice_tau_ssa_g, ice_tau_ssa_f &
         )
         ! We need to fix band ordering because the old input files assume RRTMG band
         ! ordering, but this has changed in RRTMGP.
         ! TODO: fix the input files themselves!
         do ilev = 1,nlev
            do icol = 1,ncol
               ice_tau      (:,icol,ilev) = reordered(ice_tau      (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa  (:,icol,ilev) = reordered(ice_tau_ssa  (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa_g(:,icol,ilev) = reordered(ice_tau_ssa_g(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa_f(:,icol,ilev) = reordered(ice_tau_ssa_f(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
            end do
         end do
      else if (trim(icecldoptics) == 'ebertcurry') then
         call ec_ice_optics_sw(ncol, nlev, cld, iciwp, rei, ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f)
      else
         call endrun('icecldoptics ' // trim(icecldoptics) // ' not supported.')
      end if
      call assert_range(ice_tau(1:nbnd,1:ncol,1:nlev), 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: ice_tau')
      
      ! Get liquid cloud optics
      if (trim(liqcldoptics) == 'gammadist') then
         call gammadist_liq_optics_sw( &
            ncol, nlev, iclwp, lambdac, mu, &
            liq_tau, liq_tau_ssa, &
            liq_tau_ssa_g, liq_tau_ssa_f &
         )
         do ilev = 1,nlev
            do icol = 1,ncol
               liq_tau      (:,icol,ilev) = reordered(liq_tau      (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa  (:,icol,ilev) = reordered(liq_tau_ssa  (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa_g(:,icol,ilev) = reordered(liq_tau_ssa_g(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa_f(:,icol,ilev) = reordered(liq_tau_ssa_f(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
            end do
         end do
      else if (trim(liqcldoptics) == 'slingo') then
         call slingo_liq_optics_sw( &
            ncol, nlev, cld, iclwp, rel, &
            liq_tau, liq_tau_ssa, &
            liq_tau_ssa_g, liq_tau_ssa_f &
         )
      else
         call endrun('liqcldoptics ' // trim(liqcldoptics) // ' not supported.')
      end if

      call assert_range(liq_tau(1:nbnd,1:ncol,1:nlev), 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: liq_tau')

      ! Get snow cloud optics
      if (do_snow) then
         call mitchell_ice_optics_sw( &
            ncol, nlev, icswp, des, &
            snow_tau, snow_tau_ssa, &
            snow_tau_ssa_g, snow_tau_ssa_f &
         )
         do ilev = 1,nlev
            do icol = 1,ncol
               snow_tau      (:,icol,ilev) = reordered(snow_tau      (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               snow_tau_ssa  (:,icol,ilev) = reordered(snow_tau_ssa  (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               snow_tau_ssa_g(:,icol,ilev) = reordered(snow_tau_ssa_g(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               snow_tau_ssa_f(:,icol,ilev) = reordered(snow_tau_ssa_f(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
            end do
         end do
      else
         ! We are not doing snow optics, so set these to zero so we can still use 
         ! the arrays without additional logic
         snow_tau(:,:,:) = 0.0
         snow_tau_ssa(:,:,:) = 0.0
         snow_tau_ssa_g(:,:,:) = 0.0
         snow_tau_ssa_f(:,:,:) = 0.0
      end if

      ! Combine all cloud optics from CAM routines
      cld_tau = ice_tau + liq_tau
      cld_tau_ssa = ice_tau_ssa + liq_tau_ssa
      cld_tau_ssa_g = ice_tau_ssa_g + liq_tau_ssa_g
      if (do_snow) then
         call combine_properties( &
            nbnd, ncol, nlev, &
            cld(1:ncol,1:nlev), cld_tau(1:nbnd,1:ncol,1:nlev), &
            cldfsnow(1:ncol,1:nlev), snow_tau(1:nbnd,1:ncol,1:nlev), &
            combined_tau(1:nbnd,1:ncol,1:nlev) &
         )
         call combine_properties( &
            nbnd, ncol, nlev, &
            cld(1:ncol,1:nlev), cld_tau_ssa(1:nbnd,1:ncol,1:nlev), &
            cldfsnow(1:ncol,1:nlev), snow_tau_ssa(1:nbnd,1:ncol,1:nlev), &
            combined_tau_ssa(1:nbnd,1:ncol,1:nlev) &
         )
         call combine_properties( &
            nbnd, ncol, nlev, &
            cld(1:ncol,1:nlev), cld_tau_ssa_g(1:nbnd,1:ncol,1:nlev), &
            cldfsnow(1:ncol,1:nlev), snow_tau_ssa_g(1:nbnd,1:ncol,1:nlev), &
            combined_tau_ssa_g(1:nbnd,1:ncol,1:nlev) &
         )
      else
         combined_tau = cld_tau
         combined_tau_ssa = cld_tau_ssa
         combined_tau_ssa_g = cld_tau_ssa_g
      end if
      
      ! Copy to output arrays, converting to optical depth, single scattering
      ! albedo, and assymmetry parameter from the products that the CAM routines
      ! return. Make sure we do not try to divide by zero...
      do iband = 1,nbnd
         tau_out(:ncol,:nlev,iband) = combined_tau(iband,:ncol,:nlev)
         where (combined_tau(iband,:ncol,:nlev) > 0)
            ssa_out(:ncol,:nlev,iband) &
               = combined_tau_ssa(iband,:ncol,:nlev) / combined_tau(iband,:ncol,:nlev)
         elsewhere
            ssa_out(:ncol,:nlev,iband) = 1.0
         endwhere
         where (combined_tau_ssa(iband,:ncol,:nlev) > 0)
            asm_out(:ncol,:nlev,iband) &
               = combined_tau_ssa_g(iband,:ncol,:nlev) / combined_tau_ssa(iband,:ncol,:nlev)
         elsewhere
            asm_out(:ncol,:nlev,iband) = 0.0
         end where
      end do
!     ! Copy to output arrays, converting to optical depth, single scattering
!     ! albedo, and assymmetry parameter from the products that the CAM routines
!     ! return. Make sure we do not try to divide by zero.
!     ncol = state%ncol
!     optics_out%tau = 0
!     optics_out%ssa = 0
!     optics_out%g   = 0
!     do iband = 1,nswbands
!        do ilev = 1,pver
!           do icol = 1,ncol
!              optics_out%tau(icol,ilev,iband) = combined_tau(iband,icol,ilev)
!              if (combined_tau(iband,icol,ilev) > 0) then
!                 optics_out%ssa(icol,ilev,iband) = combined_tau_ssa(iband,icol,ilev) / combined_tau(iband,icol,ilev)
!              end if
!              if (combined_tau_ssa(iband,icol,ilev) > 0) then
!                 optics_out%g(icol,ilev,iband) = combined_tau_ssa_g(iband,icol,ilev) / combined_tau_ssa(iband,icol,ilev)
!              end if
!           end do
!        end do
!     end do

   end subroutine get_cloud_optics_sw
    
   !----------------------------------------------------------------------------

   subroutine get_cloud_optics_lw( &
         ncol, nlev, nbnd, do_snow, cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rei, &
         tau_out)

      use ppgrid, only: pcols
      use cam_abortutils, only: endrun
      use cloud_rad_props, only: gammadist_liq_optics_lw, mitchell_ice_optics_lw
      use ebert_curry, only: ec_ice_optics_lw
      use slingo, only: slingo_liq_optics_lw
      use radconstants, only: nlwbands

      integer, intent(in) :: ncol, nlev, nbnd
      logical, intent(in) :: do_snow
      real(r8), intent(in), dimension(:,:) :: &
         cld, cldfsnow, &
         iclwp, iciwp, icswp, &
         mu, lambdac, dei, des, rei
      real(r8), intent(out), dimension(:,:,:) :: tau_out

      ! Temporary variables to hold absorption optical depth
      real(r8), dimension(nbnd,pcols,nlev) :: &
            ice_tau, liq_tau, snow_tau, cld_tau, combined_tau

      integer :: iband

      ! initialize
      ice_tau(:,:,:) = 0.0
      liq_tau(:,:,:) = 0.0
      snow_tau(:,:,:) = 0.0
      cld_tau(:,:,:) = 0.0
      combined_tau(:,:,:) = 0.0

      ! Get ice optics
      if (trim(icecldoptics) == 'mitchell') then
         call mitchell_ice_optics_lw(ncol, nlev, iciwp, dei, ice_tau)
      else if (trim(icecldoptics) == 'ebertcurry') then
         call ec_ice_optics_lw(ncol, nlev, cld, iclwp, iciwp, rei, ice_tau)
      else
         call endrun('icecldoptics ' // trim(icecldoptics) // ' not supported.')
      end if

      ! Get liquid optics
      if (trim(liqcldoptics) == 'gammadist') then
         call gammadist_liq_optics_lw(ncol, nlev, iclwp, lambdac, mu, liq_tau)
      else if (trim(liqcldoptics) == 'slingo') then
         call slingo_liq_optics_lw(ncol, nlev, cld, iclwp, iciwp, liq_tau)
      else
         call endrun('liqcldoptics ' // trim(liqcldoptics) // ' not supported.')
      end if

      ! Combined cloud optics
      cld_tau = liq_tau + ice_tau

      ! Get snow optics?
      if (do_snow) then
         call mitchell_ice_optics_lw(ncol, nlev, icswp, des, snow_tau)
         call combine_properties(nbnd, ncol, nlev, &
            cld(1:ncol,1:nlev), cld_tau(1:nbnd,1:ncol,1:nlev), &
            cldfsnow(1:ncol,1:nlev), snow_tau(1:nbnd,1:ncol,1:nlev), &
            combined_tau(1:nbnd,1:ncol,1:nlev) &
         )
      else
         combined_tau(1:nbnd,1:ncol,1:nlev) = cld_tau(1:nbnd,1:ncol,1:nlev)
      end if

      ! Set output optics
      do iband = 1,nbnd
         tau_out(1:ncol,1:nlev,iband) = combined_tau(iband,1:ncol,1:nlev)
      end do

   end subroutine get_cloud_optics_lw

   !----------------------------------------------------------------------------

   ! Provide a procedure to combine cloud optical properties by weighting
   ! contributions by fraction present. I.e., for combining cloud and snow
   ! optical properties, we weight the cloud and snow properties by the fraction
   ! of cloud and snow present.
   subroutine combine_properties(nbands, ncols, nlevs, &
                                 fraction1, property1, &
                                 fraction2, property2, &
                                 combined_property)
      
      ! Input dimensions for automatic error checking
      integer, intent(in) :: nbands, ncols, nlevs

      ! Input fractions of each type of constituent
      real(r8), intent(in) :: fraction1(ncols,nlevs), fraction2(ncols,nlevs)

      ! Individual optical properties for each constituent
      real(r8), intent(in) :: property1(nbands,ncols,nlevs), property2(nbands,ncols,nlevs)

      ! Combined optical property
      real(r8), intent(out) :: combined_property(nbands,ncols,nlevs)

      ! Combined fraction (max of property 1 and 2)
      real(r8) :: combined_fraction(ncols,nlevs)

      ! Loop variables
      integer :: iband, icol, ilev

      ! Combined fraction
      combined_fraction = max(fraction1, fraction2)

      ! Combine optical properties by weighting by amount of cloud and snow
      combined_property = 0
      do ilev = 1,nlevs
         do icol = 1,ncols
            do iband = 1,nbands
               if (combined_fraction(icol,ilev) > 0) then
                  combined_property(iband,icol,ilev) = ( &
                     fraction1(icol,ilev) * property1(iband,icol,ilev) &
                   + fraction2(icol,ilev) * property2(iband,icol,ilev) &
                  ) / combined_fraction(icol,ilev)
               else
                  combined_property(iband,icol,ilev) = 0
               end if
            end do
         end do
      end do

   end subroutine combine_properties

   !----------------------------------------------------------------------------

   subroutine set_cloud_optics_sw( &
         ncol, nlev, kdist, pmid, cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rel, rei, &
         optics_out)
      
      use ppgrid, only: pcols, pver, pverp
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, &
                                pbuf_get_index, pbuf_old_tim_idx
      use mo_optical_props, only: ty_optical_props_2str
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mcica_subcol_gen, only: mcica_subcol_mask

      integer, intent(in) :: ncol, nlev
      type(ty_gas_optics_rrtmgp), intent(in) :: kdist
      real(r8), intent(in), dimension(:,:) :: &
         pmid, cld, cldfsnow, &
         iclwp, iciwp, icswp, dei, des, lambdac, mu, &
         rei, rel
      type(ty_optical_props_2str), intent(inout) :: optics_out

      real(r8), dimension(ncol,nlev,nswbands) :: cld_tau_bnd, cld_ssa_bnd, cld_asm_bnd

      integer :: ngpt, icol, ilev, igpt, iband, ilev_cam, ilev_rad

      ! Set a name for this subroutine to write to error messages
      character(len=32) :: subname = 'set_cloud_optics_sw'

      ngpt = kdist%get_ngpt()

      ! Get optics by band
      call get_cloud_optics_sw( &
         ncol, nlev, nswbands, do_snow_optics(), cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rel, rei, &
         cld_tau_bnd, cld_ssa_bnd, cld_asm_bnd &
      )
      ! Initialize (or reset) output cloud optics object
      optics_out%tau = 0.0
      optics_out%ssa = 1.0
      optics_out%g = 0.0
      ! Do MCICA sampling
      call sample_cloud_optics_sw( &
         ncol, nlev, ngpt, kdist%get_gpoint_bands(), &
         pmid, cld, cldfsnow, &
         cld_tau_bnd, cld_ssa_bnd, cld_asm_bnd, &
         optics_out%tau(1:ncol,ktop:kbot,1:ngpt), &
         optics_out%ssa(1:ncol,ktop:kbot,1:ngpt), &
         optics_out%g  (1:ncol,ktop:kbot,1:ngpt)  &
      )
      ! Apply delta scaling to account for forward-scattering
      ! TODO: delta_scale takes the forward scattering fraction as an optional
      ! parameter. In the current cloud optics_sw scheme, forward scattering is taken
      ! just as g^2, which delta_scale assumes if forward scattering fraction is
      ! omitted in the function call. In the future, we should explicitly pass
      ! this. This just requires modifying the get_cloud_optics_sw procedures to also
      ! pass the foward scattering fraction that the CAM cloud optics_sw assumes.
      call handle_error(optics_out%delta_scale())

   end subroutine set_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine set_cloud_optics_lw( &
         ncol, nlev, kdist, pmid, cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rei, &
         optics_out)
      
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, pbuf_get_index, pbuf_old_tim_idx
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mo_optical_props, only: ty_optical_props_1scl
      use mcica_subcol_gen, only: mcica_subcol_mask

      integer, intent(in) :: ncol, nlev
      type(ty_gas_optics_rrtmgp), intent(in) :: kdist
      real(r8), intent(in), dimension(:,:) :: &
         pmid, cld, cldfsnow, &
         iclwp, iciwp, icswp, dei, des, lambdac, mu, &
         rei
      type(ty_optical_props_1scl), intent(inout) :: optics_out

      real(r8) :: combined_cld(pcols,pver)
      real(r8), dimension(pcols,pver,nlwbands) :: tau_bnd

      ! For MCICA sampling routine
      integer, parameter :: changeseed = 1

      ! Dimension sizes
      integer :: ngpt

      ! Temporary arrays to hold mcica-sampled cloud optics
      logical, allocatable :: iscloudy(:,:,:)

      ! Loop variables
      integer :: icol, ilev_rad, igpt, iband

      ! Set dimension size working variables
      ngpt = kdist%get_ngpt()
      ! Get optics by band
      call get_cloud_optics_lw( &
         ncol, nlev, nlwbands, do_snow_optics(), cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rei, &
         tau_bnd &
      )
      ! Do MCICA sampling
      optics_out%tau = 0
      call sample_cloud_optics_lw( &
         ncol, nlev, ngpt, kdist%get_gpoint_bands(), &
         pmid, cld, cldfsnow, &
         tau_bnd, optics_out%tau(1:ncol,ktop:kbot,1:ngpt) &
      )
   end subroutine set_cloud_optics_lw

   !----------------------------------------------------------------------------
   ! Do MCICA sampling of optics here. This will map bands to gpoints,
   ! while doing stochastic sampling of cloud state
   subroutine sample_cloud_optics_sw( &
         ncol, nlev, ngpt, gpt2bnd, &
         pmid, cld, cldfsnow, &
         tau_bnd, ssa_bnd, asm_bnd, &
         tau_gpt, ssa_gpt, asm_gpt)
      use mcica_subcol_gen, only: mcica_subcol_mask
      integer, intent(in) :: ncol, nlev, ngpt
      integer, intent(in), dimension(:) :: gpt2bnd
      real(r8), intent(in), dimension(:,:) :: pmid, cld, cldfsnow
      real(r8), intent(in), dimension(:,:,:) :: tau_bnd, ssa_bnd, asm_bnd
      real(r8), intent(out), dimension(:,:,:) :: tau_gpt, ssa_gpt, asm_gpt
      real(r8), dimension(ncol,nlev) :: combined_cld
      logical, dimension(ngpt,ncol,nlev) :: iscloudy

      ! For MCICA sampling routine, how many times to permute random seed
      integer, parameter :: changeseed = 1

      ! Loop variables
      integer :: icol, ilev, igpt

      ! Combined snow and cloud fraction
      combined_cld(1:ncol,1:nlev) = max(cld(1:ncol,1:nlev), &
                                        cldfsnow(1:ncol,1:nlev))
      ! Get stochastic subcolumn cloud mask
      call mcica_subcol_mask(ngpt, ncol, nlev, changeseed, &
                             pmid(1:ncol,1:nlev), &
                             combined_cld(1:ncol,1:nlev), &
                             iscloudy(1:ngpt,1:ncol,1:nlev))
      ! Generate subcolumns for homogeneous clouds
      do igpt = 1,ngpt
         do ilev = 1,nlev
            do icol = 1,ncol
               if (iscloudy(igpt,icol,ilev) .and. combined_cld(icol,ilev) > 0._r8) then
                  tau_gpt(icol,ilev,igpt) = tau_bnd(icol,ilev,gpt2bnd(igpt))
                  ssa_gpt(icol,ilev,igpt) = ssa_bnd(icol,ilev,gpt2bnd(igpt))
                  asm_gpt(icol,ilev,igpt) = asm_bnd(icol,ilev,gpt2bnd(igpt))
               else
                  tau_gpt(icol,ilev,igpt) = 0
                  ssa_gpt(icol,ilev,igpt) = 1
                  asm_gpt(icol,ilev,igpt) = 0
               end if
            end do
         end do
      end do
   end subroutine sample_cloud_optics_sw

   !----------------------------------------------------------------------------
   ! Do MCICA sampling of optics here. This will map bands to gpoints,
   ! while doing stochastic sampling of cloud state
   subroutine sample_cloud_optics_lw( &
         ncol, nlev, ngpt, gpt2bnd, &
         pmid, cld, cldfsnow, &
         tau_bnd, tau_gpt)
      use mcica_subcol_gen, only: mcica_subcol_mask

      integer, intent(in) :: ncol, nlev, ngpt
      integer, intent(in), dimension(:) :: gpt2bnd
      real(r8), intent(in), dimension(:,:) :: pmid, cld, cldfsnow
      real(r8), intent(in), dimension(:,:,:) :: tau_bnd
      real(r8), intent(out), dimension(:,:,:) :: tau_gpt
      real(r8), dimension(ncol,nlev) :: combined_cld
      logical, dimension(ngpt,ncol,nlev) :: iscloudy

      ! For MCICA sampling routine, how many times to permute random seed
      integer, parameter :: changeseed = 1

      ! Loop variables
      integer :: icol, ilev, igpt

      ! Combine cloud and snow fractions for MCICA sampling
      combined_cld(1:ncol,1:nlev) = max(cld(1:ncol,1:nlev), cldfsnow(1:ncol,1:nlev))

      ! Get the stochastic subcolumn cloudy mask
      call mcica_subcol_mask(ngpt, ncol, nlev, changeseed, &
                             pmid(1:ncol,1:nlev), &
                             combined_cld(1:ncol,1:nlev), &
                             iscloudy(1:ngpt,1:ncol,1:nlev))

      ! Map optics to g-points, selecting a single subcolumn for each
      ! g-point. This implementation generates homogeneous clouds, but it would be
      ! straightforward to extend this to handle horizontally heterogeneous clouds
      ! as well.
      do igpt = 1,ngpt
         do ilev = 1,nlev
            do icol = 1,ncol
               if (iscloudy(igpt,icol,ilev) .and. combined_cld(icol,ilev) > 0._r8) then
                  tau_gpt(icol,ilev,igpt) = tau_bnd(icol,ilev,gpt2bnd(igpt))
               else
                  tau_gpt(icol,ilev,igpt) = 0
               end if
            end do
         end do
      end do
   end subroutine sample_cloud_optics_lw

   !----------------------------------------------------------------------------

   subroutine set_aerosol_optics_lw(icall, state, pbuf, is_cmip6_volc, tau)
     
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, &
                                pbuf_get_field, pbuf_old_tim_idx
      use aer_rad_props, only: aer_rad_props_lw
      use radconstants, only: nlwbands

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      logical, intent(in) :: is_cmip6_volc
      real(r8), intent(out), dimension(:,:,:) :: tau(pcols,pver,nlwbands)

      ! Subroutine name for error messages
      character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_lw'

      ! Get aerosol absorption optical depth from CAM routine
      tau = 0.0
      call aer_rad_props_lw(is_cmip6_volc, icall, state, pbuf, tau)

   end subroutine set_aerosol_optics_lw

   !----------------------------------------------------------------------------

   subroutine set_aerosol_optics_sw( &
         icall, state, pbuf, night_indices, is_cmip6_volc, &
         tau_out, ssa_out, asm_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use aer_rad_props, only: aer_rad_props_sw

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      integer, intent(in) :: night_indices(:)
      logical, intent(in) :: is_cmip6_volc
      real(r8), intent(out) :: tau_out(pcols,pver,nswbands)
      real(r8), intent(out) :: ssa_out(pcols,pver,nswbands)
      real(r8), intent(out) :: asm_out(pcols,pver,nswbands)

      ! NOTE: aer_rad_props expects 0:pver indexing on these! It appears this is to
      ! account for the extra layer added above model top, but it is not entirely
      ! clear. This is not done for the longwave, and it is not really documented
      ! anywhere that I can find. Regardless, optical properties for the zero index
      ! are set to zero in aer_rad_props_sw as far as I can tell.
      !
      ! NOTE: dimension ordering is different than for cloud optics!
      real(r8), dimension(pcols,0:pver,nswbands) :: tau, tau_w, tau_w_g, tau_w_f

      ! Loop indices
      integer :: ncol, icol, ilev, ibnd

      ncol = state%ncol

      ! Get aerosol absorption optical depth from CAM routine
      tau = 0.0
      tau_w = 0.0
      tau_w_g = 0.0
      tau_w_f = 0.0
      call aer_rad_props_sw( &
         icall, state, pbuf, &
         count(night_indices > 0), night_indices, is_cmip6_volc, &
         tau, tau_w, tau_w_g, tau_w_f &
      )
      ! Convert from products to optical properties
      tau_out = 0
      ssa_out = 0
      asm_out = 0
      do ibnd = 1,nswbands
         do ilev = 1,pver
            do icol = 1,ncol 
               tau_out(icol,ilev,ibnd) = tau(icol,ilev,ibnd)
               if (tau(icol,ilev,ibnd) > 0) then
                  ssa_out(icol,ilev,ibnd) = tau_w(icol,ilev,ibnd) / tau(icol,ilev,ibnd)
               end if
               if (tau_w(icol,ilev,ibnd) > 0) then
                  asm_out(icol,ilev,ibnd) = tau_w_g(icol,ilev,ibnd) / tau_w(icol,ilev,ibnd)
               end if
            end do
         end do
      end do
      ! Set values
      ! We need to fix band ordering because the old input files assume RRTMG band
      ! ordering, but this has changed in RRTMGP.
      ! TODO: fix the input files themselves!
      do icol = 1,ncol
         do ilev = 1,pver
            tau_out(icol,ilev,:) = reordered(tau_out(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
            ssa_out(icol,ilev,:) = reordered(ssa_out(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
            asm_out(icol,ilev,:) = reordered(asm_out(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
         end do
      end do
   end subroutine set_aerosol_optics_sw

   !----------------------------------------------------------------------------

   subroutine free_optics_sw(optics)
      use mo_optical_props, only: ty_optical_props_2str
      type(ty_optical_props_2str), intent(inout) :: optics
      if (allocated(optics%tau)) deallocate(optics%tau)
      if (allocated(optics%ssa)) deallocate(optics%ssa)
      if (allocated(optics%g)) deallocate(optics%g)
      call optics%finalize()
   end subroutine free_optics_sw

   !----------------------------------------------------------------------------

   subroutine free_optics_lw(optics)
      use mo_optical_props, only: ty_optical_props_1scl
      type(ty_optical_props_1scl), intent(inout) :: optics
      if (allocated(optics%tau)) deallocate(optics%tau)
      call optics%finalize()
   end subroutine free_optics_lw

   !----------------------------------------------------------------------------

   ! Utility function to reorder an array given a new indexing
   function reordered(array_in, new_indexing) result(array_out)

      ! Inputs
      real(r8), intent(in) :: array_in(:)
      integer, intent(in) :: new_indexing(:)

      ! Output, reordered array
      real(r8), dimension(size(array_in)) :: array_out

      ! Loop index
      integer :: ii

      ! Check inputs
      call assert(size(array_in) == size(new_indexing), 'reorder_array: sizes inconsistent')

      ! Reorder array based on input index mapping, which maps old indices to new
      do ii = 1,size(new_indexing)
         array_out(ii) = array_in(new_indexing(ii))
      end do

   end function reordered

   !----------------------------------------------------------------------------

   subroutine output_cloud_optics_sw(state, tau, ssa, asm)

      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld
      use radconstants, only: idx_sw_diag

      type(physics_state), intent(in) :: state
      real(r8), intent(in), dimension(:,:,:) :: tau, ssa, asm
      character(len=*), parameter :: subname = 'output_cloud_optics_sw'

      ! Send outputs to history buffer
      call outfld('CLOUD_TAU_SW', &
                  tau(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('CLOUD_SSA_SW', &
                  ssa(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('CLOUD_G_SW', &
                  asm(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('TOT_ICLD_VISTAU', &
                  tau(1:state%ncol,1:pver,idx_sw_diag), &
                  state%ncol, state%lchnk)

   end subroutine output_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine output_cloud_optics_lw(state, tau)

      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld

      type(physics_state), intent(in) :: state
      real(r8), intent(in), dimension(:,:,:) :: tau

      ! Output
      call outfld('CLOUD_TAU_LW', &
                  tau(1:state%ncol,1:pver,1:nlwbands), &
                  state%ncol, state%lchnk)

   end subroutine output_cloud_optics_lw

   !----------------------------------------------------------------------------

   ! Should we do snow optics? Check for existence of "cldfsnow" variable
   logical function do_snow_optics()
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index
      use phys_control, only: phys_getopts
      use cam_abortutils, only: endrun
      real(r8), pointer :: pbuf(:)
      integer :: err, idx
      logical :: use_MMF
      character(len=16) :: MMF_microphysics_scheme

      idx = pbuf_get_index('CLDFSNOW', errcode=err)
      if (idx > 0) then
         do_snow_optics = .true.
      else
         do_snow_optics = .false.
      end if

      ! Reset to false if using MMF with 1-mom scheme
      call phys_getopts(use_MMF_out           = use_MMF          )
      call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)
      if (use_MMF .and. (trim(MMF_microphysics_scheme) == 'sam1mom')) then
         do_snow_optics = .false.
      end if

      return
   end function do_snow_optics

end module cam_optics
