! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine setups up parameters related to the atmospheric state. It assumes that the
!! pressure, temperature, and dimensional fields (xc, dx, yc, dy, zc, zl) have already been
!! specified and all state arrays allocated via CARMASTATE_Create().
!! 
!! @author Chuck Bardeen
!! @ version Feb-1995
!! @see CARMASTATE_Create
subroutine setupatm(carma, cstate, rescale, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  logical, intent(in)                  :: rescale     !! rescale the fall velocity for zmet change, this is instead of realculating
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  ! Local declarations
  !--
  ! For air viscosity calculations
  ! Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
  !   Meteorological Tables, in which there is a misprint -- T is deg_K, not
  !   deg_C.
  real(kind=f), parameter :: rmu_0 = 1.8325e-4_f  
  real(kind=f), parameter :: rmu_t0 = 296.16_f     
  real(kind=f), parameter :: rmu_c = 120._f        
  real(kind=f), parameter :: rmu_const = rmu_0 * (rmu_t0 + rmu_c)
    
  integer :: ielem, ibin, i, j, ix, iy, iz, ie, ig, ip, igrp, jgrp, igroup


  ! Calculate the dry air density at each level, using the ideal gas
  ! law. This will be used to calculate zmet.
  rhoa(:) = p(:) / (R_AIR * t(:))

  ! Calculate the dimensions and the dimensional metrics.
  dz(:) = abs(zl(2:NZP1) - zl(1:NZ))
  
  ! Horizontal Metrics
  select case(igridh)
    ! Cartesian
    case (I_CART)
      xmet(:) = 1._f
      ymet(:) = 1._f
    
    ! Latitude/Longitude
    case (I_LL)
      xmet(:) = REARTH * DEG2RAD * cos(DEG2RAD * yc(:))
      ymet(:) = REARTH * DEG2RAD
      
    case default
      if (do_print) write(LUNOPRT,*) "setupatm:: ERROR - The specified horizontal grid type (", igridh, &
        ") is not supported."
      rc = -1
  end select

  
  ! Put the fall velocity back into cgs units, so that we can determine
  ! new metrics and then scale it back. This is optional and is done instead
  ! of recalculating everything from scratch to improve performance.
  if (rescale .and. (igridv /= I_CART)) then
    do ibin = 1, NBIN
      do igroup = 1, NGROUP
        vf(:, ibin, igroup)   = vf(:, ibin, igroup)  * zmetl(:)
        dkz(:, ibin, igroup)  = dkz(:, ibin, igroup) * (zmetl(:)**2)
      end do
    end do
  end if

 
  ! Vertical Metrics
  select case(igridv)
    ! Cartesian
    case (I_CART)
      zmet = 1._f
      
    ! Sigma
    case (I_SIG)
      zmet(:) = abs(((pl(1:NZ) - pl(2:NZP1)) / (zl(1:NZ) - zl(2:NZP1))) / &
        (GRAV * rhoa(:)))
        
    ! Hybrid
    case (I_HYBRID)
      zmet(:) = abs(((pl(1:NZ) - pl(2:NZP1)) / (zl(1:NZ) - zl(2:NZP1))) / &
        (GRAV * rhoa(:)))
        
    case default
      if (do_print) write(LUNOPRT,*) "setupatm:: ERROR - The specified vertical grid type (", igridv, &
        ") is not supported."
      rc = -1
  end select
 
  ! Interpolate the z metric to the grid box edges.
  if (NZ == 1) then
    zmetl(:) = zmet(1)
  else
    
    ! Extrpolate the top and bottom.
    zmetl(1)    = zmet(1)  + (zmet(2) - zmet(1))     / (zc(2) - zc(1))     * (zl(1) - zc(1))
    zmetl(NZP1) = zmet(NZ) + (zmet(NZ) - zmet(NZ-1)) / (zc(NZ) - zc(NZ-1)) * (zl(NZP1) - zc(NZ))
    
    ! Interpolate the middles.
    if (NZ > 2) then
      do iz = 2, NZ
        zmetl(iz) = zmet(iz-1) + (zmet(iz) - zmet(iz-1)) / (zc(iz) - zc(iz-1)) * (zl(iz) - zc(iz-1))
      end do
    end if
  end if

 
  ! Determine the z metrics at the grid box edges and then use this to put the
  ! fall velocity back into /x/y/z units.
  if (rescale .and. (igridv /= I_CART)) then
    do ibin = 1, NBIN
      do igroup = 1, NGROUP
        vf(:, ibin, igroup)   = vf(:, ibin, igroup)  / zmetl(:)
        dkz(:, ibin, igroup)  = dkz(:, ibin, igroup) / (zmetl(:)**2)
      end do
    end do
  end if
  
  
  ! Scale the density into the units carma wants (i.e. /x/y/z)     
  rhoa(:) = rhoa(:) * xmet(:) * ymet(:) * zmet(:)

  ! Use the pressure difference across the cell and the fact that the 
  ! atmosphere is hydrostatic to caclulate an average density in the
  ! grid box.
  rhoa_wet(:) = abs((pl(2:NZP1) - pl(1:NZ))) / (GRAV)
  rhoa_wet(:) = (rhoa_wet(:) * xmet(:) * ymet(:)) / dz(:)

  ! Calculate the thermal properties of the atmosphere.
  rmu(:)     = rmu_const / ( t(:) + rmu_c ) * (t(:) / rmu_t0 )**1.5_f
  thcond(:)  = (5.69_f + .017_f*(t(:) - T0)) * 4.186e2_f
end subroutine
