module cam_optics

   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range

   implicit none
   private

   public cam_optics_type, &
          get_aerosol_optics_sw, &
          get_aerosol_optics_lw, &
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
   !subroutine cam_optics_get_aerosol_optics(this, icall, state, pbuf, night_indices)

   subroutine get_aerosol_optics_sw(icall, state, pbuf, night_indices, optics_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use aer_rad_props, only: aer_rad_props_sw
      use radconstants, only: nswbands

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      integer, intent(in) :: night_indices(:)

      ! Outputs (nbands,ncol,pver)
      type(cam_optics_type), intent(inout) :: optics_out

      ! NOTE: aer_rad_props expects 0:pver indexing on these!
      ! NOTE: dimension ordering is different than for cloud optics!
      real(r8), dimension(pcols,0:pver,nswbands) :: tau, tau_w, tau_w_g, tau_w_f

      integer :: ncol, nbnd
      integer :: iband

      ! Everyone needs a name
      character(len=*), parameter :: subroutine_name = 'get_aerosol_optics_sw'

      ncol = state%ncol
      nbnd = optics_out%nbands

      ! Check dimensions
      call assert(nbnd == nswbands, 'nbnd /= nswbands')

      ! Get aerosol absorption optical depth from CAM routine
      tau = 0.0
      tau_w = 0.0
      tau_w_g = 0.0
      tau_w_f = 0.0
      call aer_rad_props_sw(icall, state, pbuf, &
                            count(night_indices > 0), night_indices, &
                            tau, tau_w, tau_w_g, tau_w_f)

      ! Convert outputs from the product definitions to actual optical depth,
      ! single scattering albedo, and assymmetry parameter.

      ! Copy cloud optical depth over directly
      optics_out%optical_depth(1:ncol,1:pver,1:nbnd) = tau(1:ncol,1:pver,1:nbnd)

      ! Extract single scattering albedo from the product-defined fields
      where (tau(1:ncol,1:pver,1:nbnd) > 0)
         optics_out%single_scattering_albedo(1:ncol,1:pver,1:nbnd) &
            = tau_w(1:ncol,1:pver,1:nbnd) / tau(1:ncol,1:pver,1:nbnd)
      elsewhere
         optics_out%single_scattering_albedo(1:ncol,1:pver,1:nbnd) = 1
      endwhere

      ! Extract assymmetry parameter from the product-defined fields
      where (tau_w(1:ncol,1:pver,1:nbnd) > 0)
         optics_out%assymmetry_parameter(1:ncol,1:pver,1:nbnd) &
            = tau_w_g(1:ncol,1:pver,1:nbnd) / tau_w(1:ncol,1:pver,1:nbnd)
      elsewhere
         optics_out%assymmetry_parameter(1:ncol,1:pver,1:nbnd) = 0
      endwhere

      ! Extract forward scattering fraction from the product-defined fields
      where (tau_w(1:ncol,1:pver,1:nbnd) > 0)
         optics_out%forward_scattering_fraction(1:ncol,1:pver,1:nbnd) &
            = tau_w_f(1:ncol,1:pver,1:nbnd) / tau_w(1:ncol,1:pver,1:nbnd)
      elsewhere
         optics_out%forward_scattering_fraction(1:ncol,1:pver,1:nbnd) = 0
      endwhere

      ! Check values
      call assert_valid(optics_out%optical_depth, 'aerosol optical_depth')
      call assert_valid(optics_out%single_scattering_albedo, &
                        'aerosol single_scattering_albedo')
      call assert_valid(optics_out%assymmetry_parameter, &
                        'aerosol assymmetry_parameter')

      call assert_range(optics_out%assymmetry_parameter, -1._r8, 1._r8, &
                        'assymmetry parameter')
      
   end subroutine get_aerosol_optics_sw


   subroutine get_aerosol_optics_lw(icall, state, pbuf, optics_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use aer_rad_props, only: aer_rad_props_lw
      use radconstants, only: nlwbands

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(cam_optics_type), intent(inout) :: optics_out

      ! Output from CAM routines; expected to have dimension pcols
      real(r8) :: absorption_tau(pcols,pver,nlwbands)

      ! Subroutine name for error messages
      character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_lw'

      ! Make sure aerosol optical properties have been initialized
      call assert(allocated(optics_out%optical_depth), &
                  subroutine_name // ': optical_depth not initialized.')

      ! Get aerosol absorption optical depth from CAM routine
      absorption_tau = 0.0
      call aer_rad_props_lw(icall, state, pbuf, absorption_tau)

      ! Make sure all optical depths are within range
      call assert(all(absorption_tau >= 0), &
                  subroutine_name // ': tau has negative values.')

      ! Copy to output
      optics_out%optical_depth(:,:,:) = 0.0
      optics_out%optical_depth(:state%ncol,:pver,:) = absorption_tau(:state%ncol,:pver,:)

      ! Validate
      call assert(all(optics_out%optical_depth >= 0), &
                  subroutine_name // ': tau has negative values')

   end subroutine get_aerosol_optics_lw


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
            snow_tau, snow_tau_ssa, snow_tau_ssa_g, snow_tau_ssa_f, &
            combined_tau, combined_tau_ssa, combined_tau_ssa_g, combined_tau_ssa_f

      ! Pointers to fields on the physics buffer
      real(r8), pointer :: iciwp(:,:), dei(:,:)

      ! Flag to see if we should be doing snow optics. To set this, we can look for
      ! the "snow cloud fraction" on the physics buffer, and if found then set this
      ! to .true.
      logical :: do_snow_optics = .false.
      integer :: err

      integer :: nbnd, ncol, i

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
      call get_ice_optics_sw(state, pbuf, &
                             ice_tau, ice_tau_ssa, &
                             ice_tau_ssa_g, ice_tau_ssa_f)
      call assert_range(ice_tau(:nswbands,:ncol,:pver), 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: ice_tau')
      
      ! Get liquid cloud optics
      call get_liquid_optics_sw(state, pbuf, &
                                liquid_tau, liquid_tau_ssa, &
                                liquid_tau_ssa_g, liquid_tau_ssa_f)
      call assert_range(liquid_tau(:nswbands,:ncol,:pver), 0._r8, 1e20_r8, &
                        'get_cloud_optics_sw: liquid_tau')

      ! Should we do snow optics? Check for existence of "cldfsnow" variable
      ! NOTE: turned off for now...we need to figure out how to adjust the cloud
      ! fraction seen by the mcica sampling as well when we are doing snow optics.
      ! The thing to do then is probably to set this at the module level.
      !call pbuf_get_index(pbuf, 'cldfsnow', err=err)
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

      ! Combine all cloud optics from CAM routines
      combined_tau = ice_tau + liquid_tau + snow_tau
      combined_tau_ssa = ice_tau_ssa + liquid_tau_ssa + snow_tau_ssa
      combined_tau_ssa_g = ice_tau_ssa_g + liquid_tau_ssa_g + snow_tau_ssa_g

      ! Copy to output arrays, converting to optical depth, single scattering
      ! albedo, and assymmetry parameter from the products that the CAM routines
      ! return. Make sure we do not try to divide by zero...
      ncol = state%ncol
      nbnd = nswbands
      do i = 1,nbnd
         optics_out%optical_depth(:ncol,:pver,i) = combined_tau(i,:ncol,:pver)
         where (combined_tau(i,:ncol,:pver) > 0)
            optics_out%single_scattering_albedo(:ncol,:pver,i) &
               = combined_tau_ssa(i,:ncol,:pver) / combined_tau(i,:ncol,:pver)
         elsewhere
            optics_out%single_scattering_albedo(:ncol,:pver,i) = 1.0
         endwhere
         where (combined_tau_ssa(i,:ncol,:pver) > 0)
            optics_out%assymmetry_parameter(:ncol,:pver,i) &
               = combined_tau_ssa_g(i,:ncol,:pver) / combined_tau_ssa(i,:ncol,:pver)
         elsewhere
            optics_out%assymmetry_parameter(:ncol,:pver,i) = 0.0
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
      use physics_buffer, only: physics_buffer_desc
      use cloud_rad_props, only: get_liquid_optics_lw, &
                                 get_ice_optics_lw, &
                                 get_snow_optics_lw
      use radconstants, only: nlwbands

      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(cam_optics_type), intent(inout) :: optics_out

      ! Temporary variables to hold absorption optical depth
      real(r8), dimension(nlwbands,pcols,pver) :: &
            ice_tau, liq_tau, snow_tau, combined_tau

      integer :: i

      ! initialize
      ice_tau(:,:,:) = 0.0
      liq_tau(:,:,:) = 0.0
      snow_tau(:,:,:) = 0.0
      combined_tau(:,:,:) = 0.0

      ! Get ice optics
      call get_ice_optics_lw(state, pbuf, ice_tau)

      ! Get liquid optics
      call get_liquid_optics_lw(state, pbuf, liq_tau)

      ! Get snow optics?

      combined_tau = ice_tau + liq_tau + snow_tau

      ! Set optics_out
      do i = 1,nlwbands
         optics_out%optical_depth(:state%ncol,:pver,i) &
            = combined_tau(i,:state%ncol,:pver)
      end do

   end subroutine get_cloud_optics_lw


end module cam_optics
