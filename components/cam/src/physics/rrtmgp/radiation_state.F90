module radiation_state

   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid, only: pcols, pver

   implicit none
   private
   save

   public :: set_rad_state

   ! Indices to the top and bottom of the model vertical grid. This is to handle
   ! cases where we want to add an extra layer ABOVE the model grid. In the case of
   ! 1 level above the model top, ktop would be 2 and kbot would be pver, since the
   ! top model level would be the second level in the radiation grid. The bottom
   ! level should always be kbot = pver
   integer, public :: ktop, kbot

   ! Number of vertical levels in radiation calculations. This will generally
   ! include an extra layer above the model top (nlev_rad = pver + 1).
   integer, public :: nlev_rad

contains

   !----------------------------------------------------------------------------
   ! Create copies of state variables for radiation calculations that (probably)
   ! include an extra level above model top
   subroutine set_rad_state(state, cam_in, tmid, tint, pmid, pint)

      use physics_types, only: physics_state
      use camsrfexch, only: cam_in_t
      use physconst, only: stebol

      type(physics_state), intent(in) :: state
      type(cam_in_t), intent(in) :: cam_in
      real(r8), intent(out) :: tmid(:,:)
      real(r8), intent(out) :: tint(:,:)
      real(r8), intent(out) :: pmid(:,:)
      real(r8), intent(out) :: pint(:,:)
      integer :: ncol

      ncol = state%ncol

      ! Map CAM to rad fields
      tmid(1:ncol,ktop:kbot) = state%t(1:ncol,1:pver)
      pmid(1:ncol,ktop:kbot) = state%pmid(1:ncol,1:pver)
      pint(1:ncol,ktop:kbot+1) = state%pint(1:ncol,1:pver+1)

      ! Calculate interface temperature explicitly
      call set_interface_temperature(state, cam_in, tint(1:ncol,ktop:kbot+1))

      ! Fill layer above model top; this is done 
      ! consistent with the RRTMG implementation
      if (nlev_rad > pver) then
         tmid(:,1) = tmid(:,ktop) !state%t(:ncol,1)
         tint(:,1) = tint(:,ktop)
         pmid(:,1) = 0.5_r8 * pint(:,ktop) !state%pint(:ncol,1)
         pint(:,1) = 1.01_r8
      end if

   end subroutine set_rad_state

   !----------------------------------------------------------------------------
   ! Calculate temperatures at model level interfaces
   subroutine set_interface_temperature(state, cam_in, tint)

      use physics_types, only: physics_state
      use camsrfexch, only: cam_in_t
      use physconst, only: stebol

      type(physics_state), intent(in) :: state
      type(cam_in_t), intent(in) :: cam_in
      real(r8), intent(inout) :: tint(:,:)
      real(r8) :: dy
      integer :: i, k

      ! Calculate interface temperatures (following method used in radtpl for the
      ! longwave), using surface upward flux and stebol constant in mks units.
      !
      ! NOTE: this code copied from RRTMG implementation, and DIFFERS from what is
      ! done in RRTMGP if interface temperatures are omitted! Letting RRTMGP handle
      ! this leads to large differences in near-surface fluxes with RRTMG. TODO:
      ! why is this? Are we doing something wrong here? Maybe it's just the surface
      ! temperature (tint(:,pverp)) that is to blame?
      !
      ! TODO: physics_state provides interface temperatures, we should use those
      ! instead, but this is consistent with RRTMG implementation
      do i = 1,state%ncol

         ! Set top level interface temperature
         tint(i,1) = state%t(i,1)

         ! Calculate interface temperatures between top and bottom
         do k = 2,pver
            dy = (state%lnpint(i,k) - state%lnpmid(i,k)) / (state%lnpmid(i,k-1) - state%lnpmid(i,k))
            tint(i,k) = state%t(i,k) - dy * (state%t(i,k) - state%t(i,k-1))
         end do

         ! Surface temperature
         ! TODO: we should *not* be doing this like this here! Land model ought to
         ! say *something* about emissivity or surface temperature, but this
         ! assumes that the emissivity is 1, all the time, regardless! This is bad!
         tint(i,pver+1) = sqrt(sqrt(cam_in%lwup(i)/stebol))
      end do

   end subroutine set_interface_temperature

   !----------------------------------------------------------------------------

end module radiation_state
