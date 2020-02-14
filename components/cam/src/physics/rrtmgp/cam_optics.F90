module cam_optics

   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range
   use radiation_state, only: nlev_rad, ktop, kbot
   use radiation_utils, only: handle_error
   use radconstants, only: nswbands, nlwbands
   use rad_constituents, only: icecldoptics, liqcldoptics
   use ebert_curry, only: ec_ice_optics_sw, ec_ice_optics_lw
   use slingo, only: slingo_liq_optics_sw, slingo_liq_optics_lw
   use cam_abortutils, only: endrun

   implicit none
   private

   public set_cloud_optics_sw, &
          set_cloud_optics_lw, &
          set_aerosol_optics_sw, &
          set_aerosol_optics_lw

   ! Mapping from old RRTMG sw bands to new band ordering in RRTMGP
   integer, dimension(14) :: map_rrtmg_to_rrtmgp_swbands = (/ &
      14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 &
   /)

contains

   !-------------------------------------------------------------------------------

   subroutine get_cloud_optics_sw( &
         ncol, cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rel, rei, &
         tau_out, ssa_out, asm_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use cloud_rad_props, only: mitchell_ice_optics_sw, &
                                 gammadist_liq_optics_sw

      integer, intent(in) :: ncol
      real(r8), intent(in), dimension(:,:) :: &
         cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rel, rei

      ! Outputs are shortwave cloud optical properties *by band*. Dimensions should
      ! be nswbands,ncol,pver. Generally, this should be able to handle cases were
      ! ncol might be something like nday, and pver could be arbitrary so long as
      ! corresponding fields were defined for all indices of pver. This
      ! isn't the case right now I don't think, as cloud_rad_props makes explicit
      ! assumptions about array sizes.
      real(r8), intent(out), dimension(:,:,:) :: tau_out, ssa_out, asm_out

      ! Temporary variables to hold cloud optical properties before combining into
      ! output arrays. Same shape as output arrays, so get shapes from output.
      real(r8), dimension(nswbands,pcols,pver) :: &
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
            ncol, iciwp, dei, &
            ice_tau, ice_tau_ssa, &
            ice_tau_ssa_g, ice_tau_ssa_f &
         )
         ! We need to fix band ordering because the old input files assume RRTMG band
         ! ordering, but this has changed in RRTMGP.
         ! TODO: fix the input files themselves!
         do ilev = 1,pver
            do icol = 1,ncol
               ice_tau      (:,icol,ilev) = reordered(ice_tau      (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa  (:,icol,ilev) = reordered(ice_tau_ssa  (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa_g(:,icol,ilev) = reordered(ice_tau_ssa_g(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa_f(:,icol,ilev) = reordered(ice_tau_ssa_f(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
            end do
         end do
      else if (trim(icecldoptics) == 'ebertcurry') then
         call ec_ice_optics_sw(ncol, cld, iciwp, rei, ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f)
      else
         call endrun('icecldoptics ' // trim(icecldoptics) // ' not supported.')
      end if
      call assert_range(ice_tau(1:nswbands,1:ncol,1:pver), 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: ice_tau')
      
      ! Get liquid cloud optics
      if (trim(liqcldoptics) == 'gammadist') then
         call gammadist_liq_optics_sw( &
            ncol, iclwp, lambdac, mu, &
            liq_tau, liq_tau_ssa, &
            liq_tau_ssa_g, liq_tau_ssa_f &
         )
         do ilev = 1,pver
            do icol = 1,ncol
               liq_tau      (:,icol,ilev) = reordered(liq_tau      (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa  (:,icol,ilev) = reordered(liq_tau_ssa  (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa_g(:,icol,ilev) = reordered(liq_tau_ssa_g(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa_f(:,icol,ilev) = reordered(liq_tau_ssa_f(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
            end do
         end do
      else if (trim(liqcldoptics) == 'slingo') then
         call slingo_liq_optics_sw( &
            ncol, cld, iclwp, rel, &
            liq_tau, liq_tau_ssa, &
            liq_tau_ssa_g, liq_tau_ssa_f &
         )
      else
         call endrun('liqcldoptics ' // trim(liqcldoptics) // ' not supported.')
      end if

      call assert_range(liq_tau(1:nswbands,1:ncol,1:pver), 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: liq_tau')

      ! Get snow cloud optics
      if (do_snow_optics()) then
         call mitchell_ice_optics_sw( &
            ncol, icswp, des, &
            snow_tau, snow_tau_ssa, &
            snow_tau_ssa_g, snow_tau_ssa_f &
         )
         do ilev = 1,pver
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
      call combine_properties( &
         nswbands, ncol, pver, &
         cld(1:ncol,1:pver), cld_tau(1:nswbands,1:ncol,1:pver), &
         cldfsnow(1:ncol,1:pver), snow_tau(1:nswbands,1:ncol,1:pver), &
         combined_tau(1:nswbands,1:ncol,1:pver) &
      )
      call combine_properties( &
         nswbands, ncol, pver, &
         cld(1:ncol,1:pver), cld_tau_ssa(1:nswbands,1:ncol,1:pver), &
         cldfsnow(1:ncol,1:pver), snow_tau_ssa(1:nswbands,1:ncol,1:pver), &
         combined_tau_ssa(1:nswbands,1:ncol,1:pver) &
      )
      call combine_properties( &
         nswbands, ncol, pver, &
         cld(1:ncol,1:pver), cld_tau_ssa_g(1:nswbands,1:ncol,1:pver), &
         cldfsnow(1:ncol,1:pver), snow_tau_ssa_g(1:nswbands,1:ncol,1:pver), &
         combined_tau_ssa_g(1:nswbands,1:ncol,1:pver) &
      )
      
      ! Copy to output arrays, converting to optical depth, single scattering
      ! albedo, and assymmetry parameter from the products that the CAM routines
      ! return. Make sure we do not try to divide by zero...
      do iband = 1,nswbands
         tau_out(:ncol,:pver,iband) = combined_tau(iband,:ncol,:pver)
         where (combined_tau(iband,:ncol,:pver) > 0)
            ssa_out(:ncol,:pver,iband) &
               = combined_tau_ssa(iband,:ncol,:pver) / combined_tau(iband,:ncol,:pver)
         elsewhere
            ssa_out(:ncol,:pver,iband) = 1.0
         endwhere
         where (combined_tau_ssa(iband,:ncol,:pver) > 0)
            asm_out(:ncol,:pver,iband) &
               = combined_tau_ssa_g(iband,:ncol,:pver) / combined_tau_ssa(iband,:ncol,:pver)
         elsewhere
            asm_out(:ncol,:pver,iband) = 0.0
         end where
      end do

      ! Check values
      call assert_range(tau_out, 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: tau_out')
      call assert_range(ssa_out, 0._r8, 1._r8, &
                        'get_cloud_optics_sw: ssa_out')
      call assert_range(asm_out, -1._r8, 1._r8, &
                        'get_cloud_optics_sw: asm_out')
   end subroutine get_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine get_cloud_optics_lw( &
         ncol, cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rei, &
         tau_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index
      use cloud_rad_props, only: gammadist_liq_optics_lw, &
                                 mitchell_ice_optics_lw
      use radconstants, only: nlwbands

      integer, intent(in) :: ncol
      real(r8), intent(in), dimension(:,:) :: &
         cld, cldfsnow, &
         iclwp, iciwp, icswp, &
         mu, lambdac, dei, des, rei
      real(r8), intent(out), dimension(:,:,:) :: tau_out

      ! Temporary variables to hold absorption optical depth
      real(r8), dimension(nlwbands,pcols,pver) :: &
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
         call mitchell_ice_optics_lw(ncol, iciwp, dei, ice_tau)
      else if (trim(icecldoptics) == 'ebertcurry') then
         call ec_ice_optics_lw(ncol, cld, iclwp, iciwp, rei, ice_tau)
      else
         call endrun('icecldoptics ' // trim(icecldoptics) // ' not supported.')
      end if

      ! Get liquid optics
      if (trim(liqcldoptics) == 'gammadist') then
         call gammadist_liq_optics_lw(ncol, iclwp, lambdac, mu, liq_tau)
      else if (trim(liqcldoptics) == 'slingo') then
         call slingo_liq_optics_lw(ncol, cld, iclwp, iciwp, liq_tau)
      else
         call endrun('liqcldoptics ' // trim(liqcldoptics) // ' not supported.')
      end if

      ! Get snow optics?
      if (do_snow_optics()) then
         call mitchell_ice_optics_lw(ncol, icswp, des, snow_tau)

         ! Combined cloud optics
         cld_tau = liq_tau + ice_tau
         call combine_properties(nlwbands, ncol, pver, &
            cld(1:ncol,1:pver), cld_tau(1:nlwbands,1:ncol,1:pver), &
            cldfsnow(1:ncol,1:pver), snow_tau(1:nlwbands,1:ncol,1:pver), &
            combined_tau(1:nlwbands,1:ncol,1:pver) &
         )
      else
         combined_tau(1:nlwbands,1:ncol,1:pver) = cld_tau(1:nlwbands,1:ncol,1:pver)
      end if

      ! Set output optics
      do iband = 1,nlwbands
         tau_out(1:ncol,1:pver,iband) = combined_tau(iband,1:ncol,1:pver)
      end do

      ! Check values
      call assert_range(tau_out, 0._r8, 1e20_r8, &
                        'get_cloud_optics_lw: tau_out')

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

   subroutine set_cloud_optics_sw(state, pbuf, day_indices, kdist, optics_out)
      
      use ppgrid, only: pcols, pver, pverp
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, &
                                pbuf_get_index
      use mo_optical_props, only: ty_optical_props_2str
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mcica_subcol_gen, only: mcica_subcol_mask

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      integer, intent(in) :: day_indices(:)
      type(ty_gas_optics_rrtmgp), intent(in) :: kdist
      type(ty_optical_props_2str), intent(inout) :: optics_out

      ! Cloud optics by band
      real(r8), dimension(pcols, pver, nswbands) :: tau_bnd, ssa_bnd, asm_bnd

      ! Pointers to fields on the physics buffer
      real(r8), pointer, dimension(:,:) :: &
         cld, cldfsnow, &
         iclwp, iciwp, icswp, dei, des, lambdac, mu, &
         rei, rel

      ! Combined cloud and snow fraction
      real(r8) :: combined_cld(pcols,pver)

      ! For MCICA sampling routine
      integer, parameter :: changeseed = 1

      ! Dimension sizes
      integer :: ngpt, ncol, nday, nlay

      ! McICA subcolumn cloud flag
      logical, allocatable :: iscloudy(:,:,:)

      ! Loop variables
      integer :: icol, ilev, igpt, iband, iday, ilev_cam, ilev_rad

      ! Set a name for this subroutine to write to error messages
      character(len=32) :: subname = 'set_cloud_optics_sw'

      ncol = state%ncol
      ngpt = kdist%get_ngpt()

      ! Allocate array to hold subcolumn cloud flag
      allocate(iscloudy(ngpt,ncol,pver))

      ! Initialize output cloud optics object
      nday = count(day_indices > 0)

      ! Get fields from pbuf
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cld)
      call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), cldfsnow)
      call pbuf_get_field(pbuf, pbuf_get_index('ICLWP'), iclwp)
      call pbuf_get_field(pbuf, pbuf_get_index('ICIWP'), iciwp)
      call pbuf_get_field(pbuf, pbuf_get_index('ICSWP'), icswp)
      call pbuf_get_field(pbuf, pbuf_get_index('DEI'), dei)
      call pbuf_get_field(pbuf, pbuf_get_index('DES'), des)
      call pbuf_get_field(pbuf, pbuf_get_index('REL'), rel)
      call pbuf_get_field(pbuf, pbuf_get_index('REI'), rei)
      call pbuf_get_field(pbuf, pbuf_get_index('LAMBDAC'), lambdac)
      call pbuf_get_field(pbuf, pbuf_get_index('MU'), mu)

      ! Retrieve the mean in-cloud optical properties via CAM cloud radiative
      ! properties interface (cloud_rad_props). This retrieves cloud optical
      ! properties by *band* -- these will be mapped to g-points when doing
      ! the subcolumn sampling to account for cloud overlap.
      call get_cloud_optics_sw( &
         ncol, cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rel, rei, &
         tau_bnd, ssa_bnd, asm_bnd &
      )

      ! Send in-cloud optical depth for visible band to history buffer
      call output_cloud_optics_sw(state, tau_bnd, ssa_bnd, asm_bnd)

      ! Initialize (or reset) output cloud optics object
      optics_out%tau = 0.0
      optics_out%ssa = 1.0
      optics_out%g = 0.0

      ! Get cloud and snow fractions, and combine
      combined_cld(1:ncol,1:pver) = max(cld(1:ncol,1:pver), &
                                                   cldfsnow(1:ncol,1:pver))

      ! Do MCICA sampling of optics here. This will map bands to gpoints,
      ! while doing stochastic sampling of cloud state
      call mcica_subcol_mask(ngpt, ncol, pver, changeseed, &
                             state%pmid(1:ncol,1:pver), &
                             combined_cld(1:ncol,1:pver), &
                             iscloudy(1:ngpt,1:ncol,1:pver))
      
      ! -- generate subcolumns for homogeneous clouds -----
      ! where there is a cloud, set the subcolumn cloud properties;
      optics_out%tau(:,:,:) = 0
      optics_out%ssa(:,:,:) = 1
      optics_out%g(:,:,:) = 0
      do ilev_cam = 1,pver  ! Loop over indices on the CAM grid

         ! Index to radiation grid
         ilev_rad = ilev_cam + (nlev_rad - pver)

         ! Loop over columns and map CAM columns to those on radiation grid
         ! (daytime-only columns)
         do iday = 1,size(day_indices)

            ! Map daytime to column indices
            icol = day_indices(iday)

            ! Loop over g-points and map bands to g-points; each subcolumn
            ! corresponds to a single g-point. This is how this code implements the
            ! McICA assumptions: simultaneously sampling over cloud state and
            ! g-point.
            do igpt = 1,ngpt
               if (iscloudy(igpt,icol,ilev_cam) .and. &
                   combined_cld(icol,ilev_cam) > 0._r8) then
               
                  iband = kdist%convert_gpt2band(igpt)
                  optics_out%tau(iday,ilev_rad,igpt) = tau_bnd(icol,ilev_cam,iband)
                  optics_out%ssa(iday,ilev_rad,igpt) = ssa_bnd(icol,ilev_cam,iband)
                  optics_out%g(iday,ilev_rad,igpt) = asm_bnd(icol,ilev_cam,iband)
               else
                  optics_out%tau(iday,ilev_rad,igpt) = 0._r8
                  optics_out%ssa(iday,ilev_rad,igpt) = 1._r8
                  optics_out%g(iday,ilev_rad,igpt) = 0._r8
               end if
            end do  ! igpt
         end do  ! iday
      end do  ! ilev_cam

      ! Apply delta scaling to account for forward-scattering
      ! TODO: delta_scale takes the forward scattering fraction as an optional
      ! parameter. In the current cloud optics_sw scheme, forward scattering is taken
      ! just as g^2, which delta_scale assumes if forward scattering fraction is
      ! omitted in the function call. In the future, we should explicitly pass
      ! this. This just requires modifying the get_cloud_optics_sw procedures to also
      ! pass the foward scattering fraction that the CAM cloud optics_sw assumes.
      call handle_error(optics_out%delta_scale())

      ! Check cloud optics_sw
      call handle_error(optics_out%validate())

      deallocate(iscloudy)

   end subroutine set_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine set_cloud_optics_lw(state, pbuf, kdist, optics_out)
      
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, pbuf_get_index
      use mo_optical_props, only: ty_optical_props_1scl
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mcica_subcol_gen, only: mcica_subcol_mask

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(ty_gas_optics_rrtmgp), intent(in) :: kdist
      type(ty_optical_props_1scl), intent(inout) :: optics_out

      real(r8), dimension(pcols,pver,nlwbands) :: tau_bnd
      real(r8) :: combined_cld(pcols,pver)

      ! Pointers for pbuf fields
      real(r8), pointer, dimension(:,:) :: &
         cld, cldfsnow, &
         iclwp, iciwp, icswp, mu, lambdac, dei, des, &
         rel, rei

      ! For MCICA sampling routine
      integer, parameter :: changeseed = 1

      ! Dimension sizes
      integer :: ngpt, ncol

      ! Temporary arrays to hold mcica-sampled cloud optics (ngpt,ncol,pver)
      logical, allocatable :: iscloudy(:,:,:)

      ! Loop variables
      integer :: icol, ilev_rad, igpt, iband, ilev_cam

      ! Initialize (or reset) output cloud optics object
      optics_out%tau = 0.0

      ! Set dimension size working variables
      ngpt = kdist%get_ngpt()
      ncol = state%ncol

      ! Allocate array to hold subcolumn cloudy flag
      allocate(iscloudy(ngpt,ncol,pver))

      ! Get fields from pbuf
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cld)
      call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), cldfsnow)
      call pbuf_get_field(pbuf, pbuf_get_index('ICLWP'), iclwp)
      call pbuf_get_field(pbuf, pbuf_get_index('ICIWP'), iciwp)
      call pbuf_get_field(pbuf, pbuf_get_index('ICSWP'), icswp)
      call pbuf_get_field(pbuf, pbuf_get_index('DEI'), dei)
      call pbuf_get_field(pbuf, pbuf_get_index('DES'), des)
      call pbuf_get_field(pbuf, pbuf_get_index('REI'), rei)
      call pbuf_get_field(pbuf, pbuf_get_index('LAMBDAC'), lambdac)
      call pbuf_get_field(pbuf, pbuf_get_index('MU'), mu)

      ! Initialize cloud optics object; cloud_optics_lw will be indexed by
      ! g-point, rather than by band, and subcolumn routines will associate each
      ! g-point with a stochastically-sampled cloud state
      call handle_error(optics_out%alloc_1scl(ncol, nlev_rad, kdist))
      call optics_out%set_name('longwave cloud optics')

      ! Get cloud optics using CAM routines. This should combine cloud with snow
      ! optics, if "snow clouds" are being considered
      call get_cloud_optics_lw( &
         ncol, cld, cldfsnow, iclwp, iciwp, icswp, &
         lambdac, mu, dei, des, rei, &
         tau_bnd &
      )

      ! Check values
      call assert_range(tau_bnd, 0._r8, 1e20_r8, &
                        'set_cloud_optics_lw: tau_bnd')

      ! Send cloud optics to history buffer
      call output_cloud_optics_lw(state, tau_bnd)

      ! Get cloud and snow fractions, and combine
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cld)
      call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), cldfsnow)

      ! Combine cloud and snow fractions for MCICA sampling
      combined_cld(1:ncol,1:pver) = max(cld(1:ncol,1:pver), &
                                        cldfsnow(1:ncol,1:pver))

      ! Do MCICA sampling of optics here. This will map bands to gpoints,
      ! while doing stochastic sampling of cloud state
      !
      ! First, just get the stochastic subcolumn cloudy mask...
      call mcica_subcol_mask(ngpt, ncol, pver, changeseed, &
                             state%pmid(1:ncol,1:pver), &
                             combined_cld(1:ncol,1:pver), &
                             iscloudy(1:ngpt,1:ncol,1:pver))

      ! ... and now map optics to g-points, selecting a single subcolumn for each
      ! g-point. This implementation generates homogeneous clouds, but it would be
      ! straightforward to extend this to handle horizontally heterogeneous clouds
      ! as well.
      ! NOTE: incoming optics should be in-cloud quantites and not grid-averaged 
      ! quantities!
      optics_out%tau = 0
      do ilev_cam = 1,pver

         ! Get level index on CAM grid (i.e., the index that this rad level
         ! corresponds to in CAM fields). If this index is above the model top
         ! (index less than 0) then skip setting the optical properties for this
         ! level (leave set to zero)
         ilev_rad = ilev_cam + (nlev_rad - pver)

         do icol = 1,ncol
            do igpt = 1,ngpt
               if (iscloudy(igpt,icol,ilev_cam) .and. (combined_cld(icol,ilev_cam) > 0._r8) ) then
                  iband = kdist%convert_gpt2band(igpt)
                  optics_out%tau(icol,ilev_rad,igpt) = tau_bnd(icol,ilev_cam,iband)
               else
                  optics_out%tau(icol,ilev_rad,igpt) = 0._r8
               end if
            end do
         end do
      end do

      ! Apply delta scaling to account for forward-scattering
      ! TODO: delta_scale takes the forward scattering fraction as an optional
      ! parameter. In the current cloud optics_lw scheme, forward scattering is taken
      ! just as g^2, which delta_scale assumes if forward scattering fraction is
      ! omitted in the function call. In the future, we should explicitly pass
      ! this. This just requires modifying the get_cloud_optics_lw procedures to also
      ! pass the foward scattering fraction that the CAM cloud optics_lw assumes.
      call handle_error(optics_out%delta_scale())

      ! Check values
      call assert_range(optics_out%tau, 0._r8, 1e20_r8, &
                        'set_cloud_optics_lw: optics_out%tau')

      ! Check cloud optics
      call handle_error(optics_out%validate())

      deallocate(iscloudy)

   end subroutine set_cloud_optics_lw

   !----------------------------------------------------------------------------

   subroutine set_aerosol_optics_lw(icall, state, pbuf, is_cmip6_volc, optics_out)
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, &
                                pbuf_get_field
      use mo_optical_props, only: ty_optical_props_1scl
      use aer_rad_props, only: aer_rad_props_lw
      use radconstants, only: nlwbands

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      logical, intent(in) :: is_cmip6_volc
      type(ty_optical_props_1scl), intent(inout) :: optics_out

      ! Subroutine name for error messages
      character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_lw'

      ! Output from CAM routines; expected to have dimension pcols
      real(r8) :: absorption_tau(pcols,pver,nlwbands)

      ! Loop variables
      integer :: ilev
      integer :: ncol

      ! Get aerosol absorption optical depth from CAM routine
      absorption_tau = 0.0
      call aer_rad_props_lw(is_cmip6_volc, icall, state, pbuf, absorption_tau)

      ! Populate the RRTMGP optical properties object with CAM optical depth
      optics_out%tau(:,:,:) = 0.0
      ncol = state%ncol
      optics_out%tau(1:ncol,ktop:kbot,1:nlwbands) = absorption_tau(1:ncol,1:pver,1:nlwbands)

   end subroutine set_aerosol_optics_lw

   !----------------------------------------------------------------------------

   subroutine set_aerosol_optics_sw(icall, state, pbuf, &
                                    day_indices, night_indices, &
                                    is_cmip6_volc, &
                                    optics_out)
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use aer_rad_props, only: aer_rad_props_sw
      use radconstants, only: nswbands
      use mo_optical_props, only: ty_optical_props_2str
      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      integer, intent(in) :: day_indices(:), night_indices(:)
      logical, intent(in) :: is_cmip6_volc
      type(ty_optical_props_2str), intent(inout) :: optics_out

      ! NOTE: aer_rad_props expects 0:pver indexing on these! It appears this is to
      ! account for the extra layer added above model top, but it is not entirely
      ! clear. This is not done for the longwave, and it is not really documented
      ! anywhere that I can find. Regardless, optical properties for the zero index
      ! are set to zero in aer_rad_props_sw as far as I can tell.
      !
      ! NOTE: dimension ordering is different than for cloud optics!
      real(r8), dimension(pcols,0:pver,nswbands) :: tau, tau_w, tau_w_g, tau_w_f

      integer :: ncol
      integer :: iday, icol, ilay

      ! Everyone needs a name
      character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_sw'

      ncol = state%ncol

      ! Get aerosol absorption optical depth from CAM routine
      tau = 0.0
      tau_w = 0.0
      tau_w_g = 0.0
      tau_w_f = 0.0
      call aer_rad_props_sw(icall, state, pbuf, &
                            count(night_indices > 0), night_indices, is_cmip6_volc, &
                            tau, tau_w, tau_w_g, tau_w_f)

      ! Reset outputs (also handles case where radiation grid contains an extra
      ! layer above CAM grid)
      optics_out%tau = 0
      optics_out%ssa = 1
      optics_out%g = 0

      ! Assign daytime columns
      do iday = 1,count(day_indices > 0)

         ! Get index into full chunk-wide array for this daytime index
         icol = day_indices(iday)

         ! Copy cloud optical depth over directly
         optics_out%tau(iday,ktop:kbot,1:nswbands) = tau(icol,1:pver,1:nswbands)

         ! Extract single scattering albedo from the product-defined fields
         where (tau(icol,1:pver,1:nswbands) > 0)
            optics_out%ssa(iday,ktop:kbot,1:nswbands) &
               = tau_w(icol,1:pver,1:nswbands) / tau(icol,1:pver,1:nswbands)
         elsewhere
            optics_out%ssa(iday,ktop:kbot,1:nswbands) = 1
         endwhere

         ! Extract assymmetry parameter from the product-defined fields
         where (tau_w(icol,1:pver,1:nswbands) > 0)
            optics_out%g(iday,ktop:kbot,1:nswbands) &
               = tau_w_g(icol,1:pver,1:nswbands) / tau_w(icol,1:pver,1:nswbands)
         elsewhere
            optics_out%g(iday,ktop:kbot,1:nswbands) = 0
         endwhere

      end do

      ! We need to fix band ordering because the old input files assume RRTMG band
      ! ordering, but this has changed in RRTMGP.
      ! TODO: fix the input files themselves!
      do icol = 1,size(optics_out%tau,1)
         do ilay = 1,size(optics_out%tau,2)
            optics_out%tau(icol,ilay,:) = reordered(optics_out%tau(icol,ilay,:), map_rrtmg_to_rrtmgp_swbands)
            optics_out%ssa(icol,ilay,:) = reordered(optics_out%ssa(icol,ilay,:), map_rrtmg_to_rrtmgp_swbands)
            optics_out%g(icol,ilay,:) = reordered(optics_out%g(icol,ilay,:), map_rrtmg_to_rrtmgp_swbands)
         end do
      end do

      ! Check values
      call handle_error(optics_out%validate())
      
   end subroutine set_aerosol_optics_sw

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

      ! Check values
      call assert_valid(tau(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%optical_depth')
      call assert_valid(ssa(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%single_scattering_albedo')
      call assert_valid(asm(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%assymmetry_parameter')

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

      ! Check values
      call assert_valid(tau(1:state%ncol,1:pver,1:nlwbands), 'cld_tau_lw')

      ! Output
      call outfld('CLOUD_TAU_LW', &
                  tau(1:state%ncol,1:pver,1:nlwbands), &
                  state%ncol, state%lchnk)

   end subroutine output_cloud_optics_lw

   !----------------------------------------------------------------------------

   ! Should we do snow optics? Check for existence of "cldfsnow" variable
   logical function do_snow_optics()
      use physics_buffer, only: pbuf_get_index
      use cam_abortutils, only: endrun
      real(r8), pointer :: pbuf(:)
      integer :: err, idx

      idx = pbuf_get_index('CLDFSNOW', errcode=err)
      if (idx > 0) then
         do_snow_optics = .true.
      else
         do_snow_optics = .false.
      end if

      return
   end function do_snow_optics 

end module cam_optics
