module cam_optics

   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range
   use radiation_state, only: nlev_rad, ktop, kbot
   use radiation_utils, only: handle_error
   use radconstants, only: nswbands, nlwbands
   use rad_constituents, only: liqcldoptics, icecldoptics

   implicit none
   private

   public ::                 &
      free_optics_sw,        &
      free_optics_lw,        &
      set_cloud_optics_sw,   &
      set_cloud_optics_lw,   &
      set_aerosol_optics_sw, &
      set_aerosol_optics_lw

   ! Mapping from old RRTMG sw bands to new band ordering in RRTMGP
   integer, dimension(14) :: map_rrtmg_to_rrtmgp_swbands = (/ &
      14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 &
   /)
   integer, dimension(14) :: map_rrtmgp_to_rrtmg_swbands = (/ &
      2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1 &
   /)

contains

   subroutine get_cloud_optics_sw(state, pbuf, optics_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index, pbuf_old_tim_idx
      use cloud_rad_props, only: get_mitchell_ice_optics_sw, &
                                 get_conley_liq_optics_sw, &
                                 get_snow_optics_sw
      use ebert_curry, only: ec_ice_optics_sw
      use slingo, only: slingo_liq_optics_sw
      use phys_control, only: phys_getopts
      use cam_abortutils, only: endrun
      use mo_optical_props, only: ty_optical_props_2str

      ! Inputs. Right now, this uses state and pbuf, and passes these along to the
      ! individual get_*_optics routines from cloud_rad_props. This is not very
      ! flexible, in that it would be difficult to use this to get optical
      ! properties from an SP configuration if wanted. So rather, it might be
      ! better to adapt this to pass in things like liquid/ice water paths,
      ! effective drop sizes, etc.
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)

      ! Outputs are shortwave cloud optical properties *by band*. Dimensions should
      ! be nswbands,ncol,pver. Generally, this should be able to handle cases were
      ! ncol might be something like nday, and pver could be arbitrary so long as
      ! corresponding pbuf/state fields were defined for all indices of pver. This
      ! isn't the case right now I don't think, as cloud_rad_props makes explicit
      ! assumptions about array sizes.
      type(ty_optical_props_2str), intent(inout) :: optics_out

      ! Temporary variables to hold cloud optical properties before combining into
      ! output arrays. Same shape as output arrays, so get shapes from output.
      real(r8), dimension(nswbands,pcols,pver) :: &
            liq_tau, liq_tau_ssa, liq_tau_ssa_g, liq_tau_ssa_f, &
            ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f, &
            cloud_tau, cloud_tau_ssa, cloud_tau_ssa_g, cloud_tau_ssa_f, &
            snow_tau, snow_tau_ssa, snow_tau_ssa_g, snow_tau_ssa_f, &
            combined_tau, combined_tau_ssa, combined_tau_ssa_g, combined_tau_ssa_f

      ! Pointers to fields on the physics buffer
      real(r8), pointer :: cloud_fraction(:,:), snow_fraction(:,:)

      integer :: ncol, iband, icol, ilev


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

      ncol = state%ncol

      ! Get ice cloud optics
      if (trim(icecldoptics) == 'mitchell') then
         call get_mitchell_ice_optics_sw(state, pbuf, &
                                ice_tau, ice_tau_ssa, &
                                ice_tau_ssa_g, ice_tau_ssa_f)

         ! Conley optics hard-coded for RRTMG band-ordering
         do ilev = 1,pver
            do icol = 1,ncol
               ice_tau      (:,icol,ilev) = reordered(ice_tau      (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa  (:,icol,ilev) = reordered(ice_tau_ssa  (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa_g(:,icol,ilev) = reordered(ice_tau_ssa_g(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               ice_tau_ssa_f(:,icol,ilev) = reordered(ice_tau_ssa_f(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
            end do
         end do
      else if (trim(icecldoptics) == 'ebertcurry') then
         call ec_ice_optics_sw (state, pbuf, &
                                ice_tau, ice_tau_ssa, &
                                ice_tau_ssa_g, ice_tau_ssa_f)
      else
         call endrun('Ice optics scheme ' // trim(icecldoptics) // ' not recognized.')
      end if
      
      ! Get liquid cloud optics
      if (trim(liqcldoptics) == 'gammadist') then
         call get_conley_liq_optics_sw(state, pbuf, &
                                       liq_tau, liq_tau_ssa, &
                                       liq_tau_ssa_g, liq_tau_ssa_f)

         ! Mitchell optics hard-coded for RRTMG band-ordering
         do ilev = 1,pver
            do icol = 1,ncol
               liq_tau      (:,icol,ilev) = reordered(liq_tau      (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa  (:,icol,ilev) = reordered(liq_tau_ssa  (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa_g(:,icol,ilev) = reordered(liq_tau_ssa_g(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               liq_tau_ssa_f(:,icol,ilev) = reordered(liq_tau_ssa_f(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
            end do
         end do
      else if (trim(liqcldoptics) == 'slingo') then
         call slingo_liq_optics_sw(state, pbuf, &
                                   liq_tau, liq_tau_ssa, &
                                   liq_tau_ssa_g, liq_tau_ssa_f)
      else
         call endrun('Liquid optics scheme ' // trim(liqcldoptics) // ' not recognized.')
      end if 

      ! Combine all cloud optics from CAM routines
      cloud_tau = ice_tau + liq_tau
      cloud_tau_ssa = ice_tau_ssa + liq_tau_ssa
      cloud_tau_ssa_g = ice_tau_ssa_g + liq_tau_ssa_g

      ! Get snow cloud optics?
      if (do_snow_optics()) then
         ! Doing snow optics; call procedure to get these from CAM state and pbuf
         call get_snow_optics_sw(state, pbuf, &
                                 snow_tau, snow_tau_ssa, &
                                 snow_tau_ssa_g, snow_tau_ssa_f)
         ! Conley optics hard-coded for RRTMG band-ordering
         do ilev = 1,pver
            do icol = 1,ncol
               snow_tau      (:,icol,ilev) = reordered(snow_tau      (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               snow_tau_ssa  (:,icol,ilev) = reordered(snow_tau_ssa  (:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               snow_tau_ssa_g(:,icol,ilev) = reordered(snow_tau_ssa_g(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
               snow_tau_ssa_f(:,icol,ilev) = reordered(snow_tau_ssa_f(:,icol,ilev), map_rrtmg_to_rrtmgp_swbands)
            end do
         end do

         ! Get cloud and snow fractions. This is used to weight the contribution to
         ! the total lw absorption by the fraction of the column that contains
         ! cloud vs snow. TODO: is this the right thing to do here?
         call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction, &
                             start=(/1,1,pbuf_old_tim_idx()/), &
                             kount=(/pcols,pver,1/))
         call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction, &
                             start=(/1,1,pbuf_old_tim_idx()/), &
                             kount=(/pcols,pver,1/))
         call combine_properties( &
            nswbands, ncol, pver, &
            cloud_fraction(1:ncol,1:pver), cloud_tau(1:nswbands,1:ncol,1:pver), &
            snow_fraction(1:ncol,1:pver), snow_tau(1:nswbands,1:ncol,1:pver), &
            combined_tau(1:nswbands,1:ncol,1:pver) &
         )
         call combine_properties( &
            nswbands, ncol, pver, &
            cloud_fraction(1:ncol,1:pver), cloud_tau_ssa(1:nswbands,1:ncol,1:pver), &
            snow_fraction(1:ncol,1:pver), snow_tau_ssa(1:nswbands,1:ncol,1:pver), &
            combined_tau_ssa(1:nswbands,1:ncol,1:pver) &
         )
         call combine_properties( &
            nswbands, ncol, pver, &
            cloud_fraction(1:ncol,1:pver), cloud_tau_ssa_g(1:nswbands,1:ncol,1:pver), &
            snow_fraction(1:ncol,1:pver), snow_tau_ssa_g(1:nswbands,1:ncol,1:pver), &
            combined_tau_ssa_g(1:nswbands,1:ncol,1:pver) &
         )
      else
         combined_tau = cloud_tau
         combined_tau_ssa = cloud_tau_ssa
         combined_tau_ssa_g = cloud_tau_ssa_g
      end if
     
      ! Copy to output arrays, converting to optical depth, single scattering
      ! albedo, and assymmetry parameter from the products that the CAM routines
      ! return. Make sure we do not try to divide by zero.
      ncol = state%ncol
      optics_out%tau = 0
      optics_out%ssa = 0
      optics_out%g   = 0
      do iband = 1,nswbands
         do ilev = 1,pver
            do icol = 1,ncol
               optics_out%tau(icol,ilev,iband) = combined_tau(iband,icol,ilev)
               if (combined_tau(iband,icol,ilev) > 0) then
                  optics_out%ssa(icol,ilev,iband) = combined_tau_ssa(iband,icol,ilev) / combined_tau(iband,icol,ilev)
               end if
               if (combined_tau_ssa(iband,icol,ilev) > 0) then
                  optics_out%g(icol,ilev,iband) = combined_tau_ssa_g(iband,icol,ilev) / combined_tau_ssa(iband,icol,ilev)
               end if
            end do
         end do
      end do
                 
   end subroutine get_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine get_cloud_optics_lw(state, pbuf, optics_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index, pbuf_old_tim_idx
      use phys_control, only: phys_getopts
      use cloud_rad_props, only: get_conley_liq_optics_lw, &
                                 get_mitchell_ice_optics_lw, &
                                 get_snow_optics_lw
      use ebert_curry, only: ec_ice_optics_lw
      use slingo, only: slingo_liq_optics_lw
      use radconstants, only: nlwbands
      use cam_abortutils, only: endrun
      use mo_optical_props, only: ty_optical_props_1scl

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(ty_optical_props_1scl), intent(inout) :: optics_out

      ! Cloud and snow fractions, used to weight optical properties by
      ! contributions due to cloud vs snow
      real(r8), pointer :: cloud_fraction(:,:), snow_fraction(:,:)
      real(r8), pointer, dimension(:,:) :: rei, iclwp, iciwp

      ! Temporary variables to hold absorption optical depth
      real(r8), dimension(nlwbands,pcols,pver) :: &
            ice_tau, liq_tau, snow_tau, cloud_tau, combined_tau

      integer :: iband, ncol


      ! Number of columns in this chunk
      ncol = state%ncol

      ! initialize
      ice_tau(:,:,:) = 0.0
      liq_tau(:,:,:) = 0.0
      snow_tau(:,:,:) = 0.0
      cloud_tau(:,:,:) = 0.0
      combined_tau(:,:,:) = 0.0

      ! Get ice optics
      if (trim(icecldoptics) == 'mitchell') then
         call get_mitchell_ice_optics_lw(state, pbuf, ice_tau)
      else if (trim(icecldoptics) == 'ebertcurry') then
         call ec_ice_optics_lw(state, pbuf, ice_tau)
      else
         call endrun('Ice optics scheme ' // trim(icecldoptics) // ' not recognized.')
      end if

      ! Get liquid optics
      if (trim(liqcldoptics) == 'gammadist') then
         call get_conley_liq_optics_lw(state, pbuf, liq_tau)
      else if (trim(liqcldoptics) == 'slingo') then
         call slingo_liq_optics_lw(state, pbuf, liq_tau)
      else
         call endrun('Ice optics scheme ' // trim(liqcldoptics) // ' not recognized.')
      end if

      ! Combine liquid and ice
      cloud_tau = liq_tau + ice_tau

      ! Get snow optics?
      if (do_snow_optics()) then
         call get_snow_optics_lw(state, pbuf, snow_tau)

         ! Get cloud and snow fractions. This is used to weight the contribution to
         ! the total lw absorption by the fraction of the column that contains
         ! cloud vs snow. TODO: is this the right thing to do here?
         call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction, &
                             start=(/1,1,pbuf_old_tim_idx()/), &
                             kount=(/pcols,pver,1/))
         call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction, &
                             start=(/1,1,pbuf_old_tim_idx()/), &
                             kount=(/pcols,pver,1/))

         ! Combined cloud optics
         call combine_properties(nlwbands, ncol, pver, &
            cloud_fraction(1:ncol,1:pver), cloud_tau(1:nlwbands,1:ncol,1:pver), &
            snow_fraction(1:ncol,1:pver), snow_tau(1:nlwbands,1:ncol,1:pver), &
            combined_tau(1:nlwbands,1:ncol,1:pver) &
         )
      else
         combined_tau = cloud_tau
      end if

      ! Set optics_out
      optics_out%tau = 0
      do iband = 1,nlwbands
         optics_out%tau(1:ncol,1:pver,iband) = combined_tau(iband,1:ncol,1:pver)
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

   subroutine set_cloud_optics_sw(state, pbuf, kdist, optics_out)
      
      use ppgrid, only: pcols, pver, pverp
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, &
                                pbuf_get_index, pbuf_old_tim_idx
      use mo_optical_props, only: ty_optical_props_2str
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mcica_subcol_gen, only: mcica_subcol_mask

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(ty_gas_optics_rrtmgp), intent(in) :: kdist
      type(ty_optical_props_2str), intent(inout) :: optics_out

      ! Temporary optics object to hold optical properties by band
      type(ty_optical_props_2str) :: optics_bnd

      ! Pointer to cloud fraction on physics buffer
      real(r8), pointer :: cloud_fraction(:,:), snow_fraction(:,:)

      ! Combined cloud and snow fraction
      real(r8) :: combined_cloud_fraction(pcols,pver)

      ! For MCICA sampling routine
      integer, parameter :: changeseed = 1

      ! Dimension sizes
      integer :: ncol, nlay, ngpt

      ! McICA subcolumn cloud flag
      logical, allocatable :: iscloudy(:,:,:)

      ! Loop variables
      integer :: icol, ilev, igpt, iband, ilev_cam, ilev_rad

      ! Set a name for this subroutine to write to error messages
      character(len=32) :: subname = 'set_cloud_optics_sw'

      ncol = state%ncol
      ngpt = kdist%get_ngpt()

      ! Allocate array to hold subcolumn cloud flag
      allocate(iscloudy(ngpt,ncol,pver))

      ! Get optics by band
      call handle_error(optics_bnd%alloc_2str(ncol, pver, kdist%get_band_lims_wavenumber()))
      call get_cloud_optics_sw(state, pbuf, optics_bnd)

      ! Initialize (or reset) output cloud optics object
      optics_out%tau = 0.0
      optics_out%ssa = 1.0
      optics_out%g = 0.0

      ! Get cloud and snow fractions, and combine
      ! Set pointer to cloud fraction; this is used by McICA routines
      ! TODO: why the extra arguments to pbuf_get_field here? Are these necessary?
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction, &
                          start=(/1,1,pbuf_old_tim_idx()/), &
                          kount=(/pcols,pver,1/))
      !call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
      if (do_snow_optics()) then
         call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction, &
                             start=(/1,1,pbuf_old_tim_idx()/), &
                             kount=(/pcols,pver,1/))
         combined_cloud_fraction(1:ncol,1:pver) = max(cloud_fraction(1:ncol,1:pver), &
                                                      snow_fraction(1:ncol,1:pver))
      else
         combined_cloud_fraction = cloud_fraction
      end if

      ! Do MCICA sampling of optics here. This will map bands to gpoints,
      ! while doing stochastic sampling of cloud state
      call mcica_subcol_mask(ngpt, ncol, pver, changeseed, &
                             state%pmid(1:ncol,1:pver), &
                             combined_cloud_fraction(1:ncol,1:pver), &
                             iscloudy(1:ngpt,1:ncol,1:pver))

      ! -- generate subcolumns for homogeneous clouds -----
      ! where there is a cloud, set the subcolumn cloud properties;
      optics_out%tau(:,:,:) = 0
      optics_out%ssa(:,:,:) = 1
      optics_out%g(:,:,:) = 0
      do ilev_cam = 1,pver  ! Loop over indices on the CAM grid

         ! Index to radiation grid
         ilev_rad = ilev_cam + (nlev_rad - pver)

         ! Loop over columns
         do icol = 1,ncol

            ! Loop over g-points and map bands to g-points; each subcolumn
            ! corresponds to a single g-point. This is how this code implements the
            ! McICA assumptions: simultaneously sampling over cloud state and
            ! g-point.
            do igpt = 1,ngpt
               if (iscloudy(igpt,icol,ilev_cam) .and. &
                   combined_cloud_fraction(icol,ilev_cam) > 0._r8) then
                  iband = kdist%convert_gpt2band(igpt)
                  optics_out%tau(icol,ilev_rad,igpt) = optics_bnd%tau(icol,ilev_cam,iband)
                  optics_out%ssa(icol,ilev_rad,igpt) = optics_bnd%ssa(icol,ilev_cam,iband)
                  optics_out%g  (icol,ilev_rad,igpt) = optics_bnd%g  (icol,ilev_cam,iband)
               else
                  optics_out%tau(icol,ilev_rad,igpt) = 0._r8
                  optics_out%ssa(icol,ilev_rad,igpt) = 1._r8
                  optics_out%g(icol,ilev_rad,igpt) = 0._r8
               end if
            end do  ! igpt
         end do  ! icol
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

      ! Free memory for optics by band
      call free_optics_sw(optics_bnd)

      deallocate(iscloudy)

   end subroutine set_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine set_cloud_optics_lw(state, pbuf, kdist, optics_out)
      
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, pbuf_get_index, pbuf_old_tim_idx
      use mo_optical_props, only: ty_optical_props_1scl
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mcica_subcol_gen, only: mcica_subcol_mask

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(ty_gas_optics_rrtmgp), intent(in) :: kdist
      type(ty_optical_props_1scl), intent(inout) :: optics_out

      type(ty_optical_props_1scl) :: optics_bnd

      real(r8), pointer :: cloud_fraction(:,:)
      real(r8), pointer :: snow_fraction(:,:)
      real(r8) :: combined_cloud_fraction(pcols,pver)

      ! For MCICA sampling routine
      integer, parameter :: changeseed = 1

      ! Dimension sizes
      integer :: ncol, ngpt

      ! Temporary arrays to hold mcica-sampled cloud optics
      logical, allocatable :: iscloudy(:,:,:)

      ! Loop variables
      integer :: icol, ilev_rad, igpt, iband, ilev_cam

      ! Set dimension size working variables
      ncol = state%ncol
      ngpt = kdist%get_ngpt()

      ! Allocate array to hold subcolumn cloudy flag
      allocate(iscloudy(ngpt,ncol,pver))

      ! Get optics by band
      call handle_error(optics_bnd%alloc_1scl(ncol, pver, kdist%get_band_lims_wavenumber()))
      call get_cloud_optics_lw(state, pbuf, optics_bnd)

      ! Combine cloud and snow fractions for MCICA sampling
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction, &
                             start=(/1,1,pbuf_old_tim_idx()/), &
                             kount=(/pcols,pver,1/))
      if (do_snow_optics()) then
         call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction, &
                             start=(/1,1,pbuf_old_tim_idx()/), &
                             kount=(/pcols,pver,1/))
         combined_cloud_fraction(1:ncol,1:pver) = max(cloud_fraction(1:ncol,1:pver), &
                                                      snow_fraction(1:ncol,1:pver))
      else
         combined_cloud_fraction = cloud_fraction
      end if

      ! Do MCICA sampling of optics here. This will map bands to gpoints,
      ! while doing stochastic sampling of cloud state
      !
      ! First, just get the stochastic subcolumn cloudy mask...
      call mcica_subcol_mask(ngpt, ncol, pver, changeseed, &
                             state%pmid(1:ncol,1:pver), &
                             combined_cloud_fraction(1:ncol,1:pver), &
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
               if (iscloudy(igpt,icol,ilev_cam) .and. (combined_cloud_fraction(icol,ilev_cam) > 0._r8) ) then
                  iband = kdist%convert_gpt2band(igpt)
                  optics_out%tau(icol,ilev_rad,igpt) = optics_bnd%tau(icol,ilev_cam,iband)
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

      ! Check cloud optics
      call handle_error(optics_out%validate())

      ! Free memory for optics by band
      call free_optics_lw(optics_bnd)

      deallocate(iscloudy)

   end subroutine set_cloud_optics_lw

   !----------------------------------------------------------------------------

   subroutine set_aerosol_optics_lw(icall, state, pbuf, is_cmip6_volc, optics_out)
     
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, &
                                pbuf_get_field, pbuf_old_tim_idx
      use aer_rad_props, only: aer_rad_props_lw
      use radconstants, only: nlwbands
      use mo_optical_props, only: ty_optical_props_1scl

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      logical, intent(in) :: is_cmip6_volc
      type(ty_optical_props_1scl), intent(inout) :: optics_out
      real(r8) :: tau(pcols,pver,nlwbands)
      integer :: ncol

      ! Subroutine name for error messages
      character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_lw'

      ! Get aerosol absorption optical depth from CAM routine
      tau = 0.0
      call aer_rad_props_lw(is_cmip6_volc, icall, state, pbuf, tau)

      ! Populate the RRTMGP optical properties object with CAM optical depth
      ncol = size(optics_out%tau, 1)
      optics_out%tau(:,:,:) = 0.0
      optics_out%tau(1:ncol,ktop:kbot,1:nlwbands) = tau(1:ncol,1:pver,1:nlwbands)

   end subroutine set_aerosol_optics_lw

   !----------------------------------------------------------------------------

   subroutine set_aerosol_optics_sw(icall, state, pbuf, &
                                    night_indices, &
                                    is_cmip6_volc, &
                                    optics_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use aer_rad_props, only: aer_rad_props_sw
      use mo_optical_props, only: ty_optical_props_2str

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      integer, intent(in) :: night_indices(:)
      logical, intent(in) :: is_cmip6_volc
      type(ty_optical_props_2str), intent(inout) :: optics_out

      real(r8) :: tau_out(pcols,pver,nswbands)
      real(r8) :: ssa_out(pcols,pver,nswbands)
      real(r8) :: asm_out(pcols,pver,nswbands)

      ! NOTE: aer_rad_props expects 0:pver indexing on these! It appears this is to
      ! account for the extra layer added above model top, but it is not entirely
      ! clear. This is not done for the longwave, and it is not really documented
      ! anywhere that I can find. Regardless, optical properties for the zero index
      ! are set to zero in aer_rad_props_sw as far as I can tell.
      !
      ! NOTE: dimension ordering is different than for cloud optics!
      real(r8), dimension(pcols,0:pver,nswbands) :: tau, tau_w, tau_w_g, tau_w_f

      ! Loop indices
      integer :: ncol, icol, ilev, ibnd, ilev_rad

      ncol = state%ncol

      ! Get aerosol absorption optical depth from CAM routine
      tau = 0.0
      tau_w = 0.0
      tau_w_g = 0.0
      tau_w_f = 0.0
      call aer_rad_props_sw(icall, state, pbuf, &
                            count(night_indices > 0), night_indices, is_cmip6_volc, &
                            tau, tau_w, tau_w_g, tau_w_f)

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

      ! Reset outputs (also handles case where radiation grid contains an extra
      ! layer above CAM grid)
      optics_out%tau = 0
      optics_out%ssa = 1
      optics_out%g = 0

      ! Set values
      ! We need to fix band ordering because the old input files assume RRTMG band
      ! ordering, but this has changed in RRTMGP.
      ! TODO: fix the input files themselves!
      do icol = 1,ncol
         do ilev = 1,pver
            ilev_rad = ilev + (nlev_rad - pver)
            optics_out%tau(icol,ilev_rad,:) = reordered(tau_out(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
            optics_out%ssa(icol,ilev_rad,:) = reordered(ssa_out(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
            optics_out%g  (icol,ilev_rad,:) = reordered(asm_out(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
         end do
      end do

      ! Check values
      call handle_error(optics_out%validate())

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

      ! Reset to false if using SPCAM with 1-mom scheme
      call phys_getopts(use_MMF_out           = use_MMF          )
      call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)
      if (use_MMF .and. (trim(MMF_microphysics_scheme) == 'sam1mom')) then
         do_snow_optics = .false.
      end if

      return
   end function do_snow_optics

end module cam_optics
