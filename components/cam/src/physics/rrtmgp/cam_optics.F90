module cam_optics

   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range

   implicit none
   private

   public cam_optics_type, &
          get_cloud_optics_sw, &
          get_cloud_optics_lw

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
      use radconstants, only: nswbands

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

end module cam_optics
