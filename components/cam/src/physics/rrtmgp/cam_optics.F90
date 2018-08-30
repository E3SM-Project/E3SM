module cam_optics

   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range
   use radiation_state, only: nlev_rad, ktop, kbot
   use radiation_utils, only: handle_error
   use radconstants, only: nswbands, nlwbands

   implicit none
   private

   public cam_optics_type, &
          set_cloud_optics_sw, &
          set_cloud_optics_lw, &
          set_aerosol_optics_sw, &
          set_aerosol_optics_lw

   type cam_optics_type
      integer :: nbands, ncolumns, nlevels
      real(r8), allocatable :: optical_depth(:,:,:)
      real(r8), allocatable :: single_scattering_albedo(:,:,:)
      real(r8), allocatable :: assymmetry_parameter(:,:,:)
      real(r8), allocatable :: forward_scattering_fraction(:,:,:)
   contains
      procedure :: initialize => cam_optics_initialize
      procedure :: finalize => cam_optics_finalize
   end type cam_optics_type

   ! Mapping from old RRTMG sw bands to new band ordering in RRTMGP
   integer, dimension(14) :: map_rrtmg_to_rrtmgp_swbands = (/ &
      14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 &
   /)

contains
   !-------------------------------------------------------------------------------
   ! Type-bound procedures for cam_optics_type
   subroutine cam_optics_initialize(this, nbands, ncolumns, nlevels)
      class(cam_optics_type), intent(inout) :: this
      integer, intent(in) :: nbands, ncolumns, nlevels

      this%nbands = nbands
      this%ncolumns = ncolumns
      this%nlevels = nlevels

      allocate(this%optical_depth(ncolumns,nlevels,nbands), &
               this%single_scattering_albedo(ncolumns,nlevels,nbands), &
               this%assymmetry_parameter(ncolumns,nlevels,nbands), &
               this%forward_scattering_fraction(ncolumns,nlevels,nbands))

      this%optical_depth = 0
      this%single_scattering_albedo = 1
      this%assymmetry_parameter = 0
      this%forward_scattering_fraction = 0
   end subroutine cam_optics_initialize
   !-------------------------------------------------------------------------------
   subroutine cam_optics_finalize(this)
      class(cam_optics_type), intent(inout) :: this
      deallocate(this%optical_depth, &
                 this%single_scattering_albedo, &
                 this%assymmetry_parameter, &
                 this%forward_scattering_fraction)
   end subroutine cam_optics_finalize
   !-------------------------------------------------------------------------------

   subroutine get_cloud_optics_sw(state, pbuf, optics_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
      use cloud_rad_props, only: get_ice_optics_sw, &
                                 get_liquid_optics_sw, &
                                 get_snow_optics_sw

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
      type(cam_optics_type), intent(inout) :: optics_out

      ! Temporary variables to hold cloud optical properties before combining into
      ! output arrays. Same shape as output arrays, so get shapes from output.
      !real(r8), dimension(size(optics_out%optical_depth,1), &
      !                    size(optics_out%optical_depth,2), &
      !                    size(optics_out%optical_depth,3)) :: &
      real(r8), dimension(nswbands,pcols,pver) :: &
            liquid_tau, liquid_tau_ssa, liquid_tau_ssa_g, liquid_tau_ssa_f, &
            ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f, &
            cloud_tau, cloud_tau_ssa, cloud_tau_ssa_g, cloud_tau_ssa_f, &
            snow_tau, snow_tau_ssa, snow_tau_ssa_g, snow_tau_ssa_f, &
            combined_tau, combined_tau_ssa, combined_tau_ssa_g, combined_tau_ssa_f

      ! Pointers to fields on the physics buffer
      real(r8), pointer :: iciwp(:,:), dei(:,:)
      real(r8), pointer :: cloud_fraction(:,:), snow_fraction(:,:)

      ! Flag to see if we should be doing snow optics. To set this, we can look for
      ! the "snow cloud fraction" on the physics buffer, and if found then set this
      ! to .true.
      logical :: do_snow_optics = .true.
      integer :: err

      integer :: ncol, iband

      ! Initialize
      ice_tau = 0
      ice_tau_ssa = 0
      ice_tau_ssa_g = 0
      ice_tau_ssa_f = 0
      liquid_tau = 0
      liquid_tau_ssa = 0
      liquid_tau_ssa_g = 0
      liquid_tau_ssa_f = 0
      snow_tau = 0
      snow_tau_ssa = 0
      snow_tau_ssa_g = 0
      snow_tau_ssa_f = 0
      combined_tau = 0
      combined_tau_ssa = 0
      combined_tau_ssa_g = 0
      combined_tau_ssa_f = 0

      ! Get ice cloud optics
      !call pbuf_get_field(pbuf, pbuf_get_index('ICIWP'), iciwp)
      !call pbuf_get_field(pbuf, pbuf_get_index('DEI'), dei)
      ncol = state%ncol
      call get_ice_optics_sw(state, pbuf, &
                             ice_tau, ice_tau_ssa, &
                             ice_tau_ssa_g, ice_tau_ssa_f)
      call assert_range(ice_tau(1:nswbands,1:ncol,1:pver), 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: ice_tau')
      
      ! Get liquid cloud optics
      call get_liquid_optics_sw(state, pbuf, &
                                liquid_tau, liquid_tau_ssa, &
                                liquid_tau_ssa_g, liquid_tau_ssa_f)
      call assert_range(liquid_tau(1:nswbands,1:ncol,1:pver), 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: liquid_tau')

      ! Should we do snow optics? Check for existence of "cldfsnow" variable
      ! NOTE: turned off for now...we need to figure out how to adjust the cloud
      ! fraction seen by the mcica sampling as well when we are doing snow optics.
      ! The thing to do then is probably to set this at the module level.
      !call pbuf_get_index(pbuf, 'CLDFSNOW', err=err)
      !if (err > 0) then
      !   do_snow_optics = .true.
      !else
      !   do_snow_optics = .false.
      !end if

      ! Get snow cloud optics
      if (do_snow_optics) then
         ! Doing snow optics; call procedure to get these from CAM state and pbuf
         call get_snow_optics_sw(state, pbuf, &
                                 snow_tau, snow_tau_ssa, &
                                 snow_tau_ssa_g, snow_tau_ssa_f)
      else
         ! We are not doing snow optics, so set these to zero so we can still use 
         ! the arrays without additional logic
         snow_tau(:,:,:) = 0.0
         snow_tau_ssa(:,:,:) = 0.0
         snow_tau_ssa_g(:,:,:) = 0.0
         snow_tau_ssa_f(:,:,:) = 0.0
      end if

      ! Get cloud and snow fractions. This is used to weight the contribution to
      ! the total lw absorption by the fraction of the column that contains
      ! cloud vs snow. TODO: is this the right thing to do here?
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
      call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction)

      ! Combine all cloud optics from CAM routines
      cloud_tau = ice_tau + liquid_tau
      cloud_tau_ssa = ice_tau_ssa + liquid_tau_ssa
      cloud_tau_ssa_g = ice_tau_ssa_g + liquid_tau_ssa_g
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
      
      ! Copy to output arrays, converting to optical depth, single scattering
      ! albedo, and assymmetry parameter from the products that the CAM routines
      ! return. Make sure we do not try to divide by zero...
      ncol = state%ncol
      do iband = 1,nswbands
         optics_out%optical_depth(:ncol,:pver,iband) = combined_tau(iband,:ncol,:pver)
         where (combined_tau(iband,:ncol,:pver) > 0)
            optics_out%single_scattering_albedo(:ncol,:pver,iband) &
               = combined_tau_ssa(iband,:ncol,:pver) / combined_tau(iband,:ncol,:pver)
         elsewhere
            optics_out%single_scattering_albedo(:ncol,:pver,iband) = 1.0
         endwhere
         where (combined_tau_ssa(iband,:ncol,:pver) > 0)
            optics_out%assymmetry_parameter(:ncol,:pver,iband) &
               = combined_tau_ssa_g(iband,:ncol,:pver) / combined_tau_ssa(iband,:ncol,:pver)
         elsewhere
            optics_out%assymmetry_parameter(:ncol,:pver,iband) = 0.0
         end where
      end do

      ! Check values
      call assert_range(optics_out%optical_depth, 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: optics_out%optical_depth')
      call assert_range(optics_out%single_scattering_albedo, 0._r8, 1._r8, &
                        'get_cloud_optics_sw: optics_out%single_scattering_albedo')
      call assert_range(optics_out%assymmetry_parameter, -1._r8, 1._r8, &
                        'get_cloud_optics_sw: optics_out%assymmetry_parameter')
   end subroutine get_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine get_cloud_optics_lw(state, pbuf, optics_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index
      use cloud_rad_props, only: get_liquid_optics_lw, &
                                 get_ice_optics_lw, &
                                 get_snow_optics_lw
      use radconstants, only: nlwbands

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(cam_optics_type), intent(inout) :: optics_out

      ! Cloud and snow fractions, used to weight optical properties by
      ! contributions due to cloud vs snow
      real(r8), pointer :: cloud_fraction(:,:), snow_fraction(:,:)

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
      call get_ice_optics_lw(state, pbuf, ice_tau)

      ! Get liquid optics
      call get_liquid_optics_lw(state, pbuf, liq_tau)

      ! Get snow optics?
      call get_snow_optics_lw(state, pbuf, snow_tau)

      ! Get cloud and snow fractions. This is used to weight the contribution to
      ! the total lw absorption by the fraction of the column that contains
      ! cloud vs snow. TODO: is this the right thing to do here?
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
      call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction)

      ! Combined cloud optics
      cloud_tau = liq_tau + ice_tau
      call combine_properties(nlwbands, ncol, pver, &
         cloud_fraction(1:ncol,1:pver), cloud_tau(1:nlwbands,1:ncol,1:pver), &
         snow_fraction(1:ncol,1:pver), snow_tau(1:nlwbands,1:ncol,1:pver), &
         combined_tau(1:nlwbands,1:ncol,1:pver) &
      )

      ! Set optics_out
      do iband = 1,nlwbands
         optics_out%optical_depth(1:ncol,1:pver,iband) &
            = combined_tau(iband,1:ncol,1:pver)
      end do

      ! Check values
      call assert_range(optics_out%optical_depth, 0._r8, 1e20_r8, &
                        'get_cloud_optics_lw: optics_out%optical_depth')

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
      use mo_gas_optics, only: ty_gas_optics
      use mcica_subcol_gen, only: mcica_subcol_mask

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      integer, intent(in) :: day_indices(:)
      type(ty_gas_optics), intent(in) :: kdist
      type(ty_optical_props_2str), intent(inout) :: optics_out

      ! Type to hold optics on CAM grid
      type(cam_optics_type) :: optics_cam

      ! Pointer to cloud fraction on physics buffer
      real(r8), pointer :: cloud_fraction(:,:), snow_fraction(:,:)

      ! Combined cloud and snow fraction
      real(r8) :: combined_cloud_fraction(pcols,pver)

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
      call handle_error(optics_out%alloc_2str(nday, nlev_rad, kdist))
      call optics_out%set_name('shortwave cloud optics')

      ! Retrieve the mean in-cloud optical properties via CAM cloud radiative
      ! properties interface (cloud_rad_props). This retrieves cloud optical
      ! properties by *band* -- these will be mapped to g-points when doing
      ! the subcolumn sampling to account for cloud overlap.
      call optics_cam%initialize(nswbands, ncol, pver)
      call get_cloud_optics_sw(state, pbuf, optics_cam)

      ! We need to fix band ordering because the old input files assume RRTMG band
      ! ordering, but this has changed in RRTMGP.
      ! TODO: fix the input files themselves!
      do icol = 1,size(optics_cam%optical_depth,1)
         do ilev = 1,size(optics_cam%optical_depth,2)
            optics_cam%optical_depth(icol,ilev,:) = reordered( &
               optics_cam%optical_depth(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands &
            )
            optics_cam%single_scattering_albedo(icol,ilev,:) = reordered( &
               optics_cam%single_scattering_albedo(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands &
            )
            optics_cam%assymmetry_parameter(icol,ilev,:) = reordered( &
               optics_cam%assymmetry_parameter(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands &
            )
         end do
      end do

      ! Send in-cloud optical depth for visible band to history buffer
      call output_cloud_optics_sw(state, optics_cam)

      ! Initialize (or reset) output cloud optics object
      optics_out%tau = 0.0
      optics_out%ssa = 1.0
      optics_out%g = 0.0

      ! Set pointer to cloud fraction; this is used by McICA routines
      ! TODO: why the extra arguments to pbuf_get_field here? Are these necessary?
      !call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction, &
      !                    start=(/1,1,pbuf_old_tim_idx()/), &
      !                    kount=(/pcols,pver,1/))

      ! Get cloud and snow fractions, and combine
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
      call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction)
      combined_cloud_fraction(1:ncol,1:pver) = max(cloud_fraction(1:ncol,1:pver), &
                                                   snow_fraction(1:ncol,1:pver))

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
                   combined_cloud_fraction(icol,ilev_cam) > 0._r8) then
               
                  iband = kdist%convert_gpt2band(igpt)
                  optics_out%tau(iday,ilev_rad,igpt) = optics_cam%optical_depth(icol,ilev_cam,iband)
                  optics_out%ssa(iday,ilev_rad,igpt) = optics_cam%single_scattering_albedo(icol,ilev_cam,iband)
                  optics_out%g(iday,ilev_rad,igpt) = optics_cam%assymmetry_parameter(icol,ilev_cam,iband)
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

      call optics_cam%finalize()

      deallocate(iscloudy)

   end subroutine set_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine set_cloud_optics_lw(state, pbuf, kdist, optics_out)
      
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, pbuf_get_index
      use mo_optical_props, only: ty_optical_props_1scl
      use mo_gas_optics, only: ty_gas_optics
      use mcica_subcol_gen, only: mcica_subcol_mask

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(ty_gas_optics), intent(in) :: kdist
      type(ty_optical_props_1scl), intent(inout) :: optics_out

      type(cam_optics_type) :: optics_cam
      real(r8), pointer :: cloud_fraction(:,:)
      real(r8), pointer :: snow_fraction(:,:)
      real(r8) :: combined_cloud_fraction(pcols,pver)

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

      ! Initialize cloud optics object; cloud_optics_lw will be indexed by
      ! g-point, rather than by band, and subcolumn routines will associate each
      ! g-point with a stochastically-sampled cloud state
      call handle_error(optics_out%alloc_1scl(ncol, nlev_rad, kdist))
      call optics_out%set_name('longwave cloud optics')

      ! Get cloud optics using CAM routines. This should combine cloud with snow
      ! optics, if "snow clouds" are being considered
      call optics_cam%initialize(nlwbands, ncol, pver)
      call get_cloud_optics_lw(state, pbuf, optics_cam)

      ! Check values
      call assert_range(optics_cam%optical_depth, 0._r8, 1e20_r8, &
                        'set_cloud_optics_lw: optics_cam%optical_depth')

      ! Send cloud optics to history buffer
      call output_cloud_optics_lw(state, optics_cam)

      ! Get cloud and snow fractions, and combine
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
      call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction)

      ! Combine cloud and snow fractions for MCICA sampling
      combined_cloud_fraction(1:ncol,1:pver) = max(cloud_fraction(1:ncol,1:pver), &
                                                   snow_fraction(1:ncol,1:pver))

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
                  optics_out%tau(icol,ilev_rad,igpt) = optics_cam%optical_depth(icol,ilev_cam,iband)
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

      call optics_cam%finalize()

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

   subroutine output_cloud_optics_sw(state, optics)
      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld
      use radconstants, only: idx_sw_diag

      type(physics_state), intent(in) :: state
      type(cam_optics_type), intent(in) :: optics
      character(len=*), parameter :: subname = 'output_cloud_optics_sw'

      ! Check values
      call assert_valid(optics%optical_depth(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%optical_depth')
      call assert_valid(optics%single_scattering_albedo(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%single_scattering_albedo')
      call assert_valid(optics%assymmetry_parameter(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%assymmetry_parameter')

      ! Send outputs to history buffer
      call outfld('CLOUD_TAU_SW', &
                  optics%optical_depth(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('CLOUD_SSA_SW', &
                  optics%single_scattering_albedo(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('CLOUD_G_SW', &
                  optics%assymmetry_parameter(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('TOT_ICLD_VISTAU', &
                  optics%optical_depth(1:state%ncol,1:pver,idx_sw_diag), &
                  state%ncol, state%lchnk)
   end subroutine output_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine output_cloud_optics_lw(state, optics)

      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld

      type(physics_state), intent(in) :: state
      type(cam_optics_type), intent(in) :: optics

      ! Check values
      call assert_valid(optics%optical_depth(1:state%ncol,1:pver,1:nlwbands), 'cloud_tau_lw')

      ! Output
      call outfld('CLOUD_TAU_LW', &
                  optics%optical_depth(1:state%ncol,1:pver,1:nlwbands), &
                  state%ncol, state%lchnk)

   end subroutine output_cloud_optics_lw

   !----------------------------------------------------------------------------

end module cam_optics
