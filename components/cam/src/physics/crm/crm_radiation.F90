module crm_radiation

   use params, only: rknd => crm_rknd

   use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
   use mo_optical_props, only: ty_optical_props_2str, ty_optical_props_1scl

   ! Use k-distribution objects from the main RRTMGP interface in E3SM. This way
   ! we do not have to initialize these ourselves, and we can be sure we are
   ! using the same absorption coefficient data as the clear-sky calc, which
   ! will probably be outside this interface. Also, this saves us some trouble
   ! for now.
   use radiation, only: k_dist_sw, k_dist_lw

   implicit none


contains

   ! Initialize this interface, including loading absorption coefficient data
   subroutine radiation_init()
      ! Do nothing here for now
   end subroutine radiation_init

   subroutine radiation_tend(mu0, sfc_alb_dir, sfc_alb_dif)
      use grid, only: pres
      real(rknd), intent(in), dimension(:) :: mu0
      real(rknd), intent(in), dimension(:,:) :: sfc_alb_dir
      real(rknd), intent(in), dimension(:,:) :: sfc_alb_dif
      call radiation_sw(pres, mu0, sfc_alb_dir, sfc_alb_dif)
      call radiation_lw()
   end subroutine radiation_tend

   subroutine radiation_finalize()
   end subroutine radiation_finalize

   subroutine radiation_sw(p_lay, mu0, sfc_alb_dir, sfc_alb_dif)
      real(rknd), intent(in), dimension(:,:) :: p_lay
      real(rknd), intent(in), dimension(:) :: mu0
      real(rknd), intent(in), dimension(:,:) :: sfc_alb_dir
      real(rknd), intent(in), dimension(:,:) :: sfc_alb_dif
      type(ty_optical_props_2str) :: cloud_optics, combined_optics
      logical :: top_at_1

      ! Determine if profiles are top down or bottom up
      top_at_1 = p_lay(1,1) < p_lay(1,2)

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
