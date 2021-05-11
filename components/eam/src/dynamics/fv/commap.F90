module commap

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plon, splon, plat

   real(r8), target :: w(plat)                  ! integration weights for physics grid
   real(r8), target :: w_staggered(plat-1)      ! integration weights for the staggered wind arrays
!
   real(r8), target :: clat(plat)               ! model latitudes (radians)
   real(r8), target :: clat_staggered (plat-1)  ! model latitudes on staggered grid (radians)
   real(r8), target :: clon(plon,plat)          ! model longitudes (radians)
   real(r8), target :: latdeg(plat)             ! model latitudes (degrees)
   real(r8), target :: londeg(plon,plat)        ! model longitudes (degrees)
   real(r8), target :: latdeg_st(plat-1)        ! model staggered latitudes (degrees)
   real(r8), target :: londeg_st(splon,plat)    ! model staggered longitudes (degrees)
end module commap

