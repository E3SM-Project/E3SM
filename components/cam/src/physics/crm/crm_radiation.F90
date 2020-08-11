module crm_radiation

   use params, only: rknd => crm_rknd
   use crmdims, only: crm_nx_rad, crm_ny_rad, crm_nz
   use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
   use mo_optical_props, only: ty_optical_props_2str, ty_optical_props_1scl

   ! Use k-distribution objects from the main RRTMGP interface in E3SM. This way
   ! we do not have to initialize these ourselves, and we can be sure we are
   ! using the same absorption coefficient data as the clear-sky calc, which
   ! will probably be outside this interface. Also, this saves us some trouble
   ! for now.
   use radiation, only: k_dist_sw, k_dist_lw, nswgpts

   implicit none


contains

   ! Initialize this interface, including loading absorption coefficient data
   subroutine radiation_init()
      ! Do nothing here for now
   end subroutine radiation_init

   subroutine radiation_tend(ncol, mu0, sfc_alb_dir, sfc_alb_dif)
      use grid, only: pres, presi
      use vars, only: tabs
      integer, intent(in) :: ncol
      real(rknd), intent(in), dimension(:) :: mu0
      real(rknd), intent(in), dimension(:,:) :: sfc_alb_dir
      real(rknd), intent(in), dimension(:,:) :: sfc_alb_dif
      integer :: ix, iy

      ! Loop over CRM columns and call driver for each; alternatively, we could
      ! pack columnns together but that might  mean more data movement
      do ix = 1,crm_nx_rad
         do iy = 1,crm_ny_rad
            call radiation_sw(ncol, nswgpts, pres, presi, tabs(:,ix,iy,:), mu0, sfc_alb_dir, sfc_alb_dif)
         end do
      end do
      call radiation_lw()
   end subroutine radiation_tend

   subroutine radiation_finalize()
   end subroutine radiation_finalize

   subroutine radiation_sw(ncol, ngpt, p_lay, p_lev, t_lay, mu0, sfc_alb_dir, sfc_alb_dif)
      integer, intent(in) :: ncol, ngpt
      real(rknd), intent(in), dimension(:,:) :: p_lay, p_lev, t_lay
      real(rknd), intent(in), dimension(:) :: mu0
      real(rknd), intent(in), dimension(:,:) :: sfc_alb_dir
      real(rknd), intent(in), dimension(:,:) :: sfc_alb_dif
      real(rknd), dimension(ncol,ngpt) :: toa_flux
      type(ty_optical_props_2str) :: cloud_optics, combined_optics
      logical :: top_at_1

      ! Determine if profiles are top down or bottom up
      top_at_1 = p_lay(1,1) < p_lay(1,2)

!     ! Compute optical properties; this is going to require gas concentrations
!     call handle_error(k_dist%gas_optics( &
!        p_lay, p_lev, t_lay, gas_concs, atmos, toa_flux &
!     ))


!     call handle_error(rte_sw( &
!        combined_optics, top_at_1, &
!        mu0,   toa_flux, &
!        sfc_alb_dir, sfc_alb_dif, &
!        fluxes))
   end subroutine radiation_sw

   subroutine radiation_lw()
   end subroutine radiation_lw

   subroutine handle_error(errmsg)
      character(len=*), intent(in) :: errmsg
      if (trim(errmsg) /= '') then
         print *, trim(errmsg)
         stop
      end if
   end subroutine handle_error

end module crm_radiation
