!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2016-2017,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:  Heating rate calculation

module mo_heating_rates
  use mo_rte_kind,      only: wp, wl
  use mo_rrtmgp_constants, only: cp_dry, grav ! Only needed for heating rate calculation
  implicit none
  private
  public ::  compute_heating_rate
contains
  ! Compute heating rate from fluxes
  ! heating rate H [K/sec] = 1/(rho cp) d f_net/d z
  ! Here use hydrostatic equation for density and heat capacity of dry air
  function compute_heating_rate(flux_up, flux_dn, plev, heating_rate) result(error_msg)
    real(wp), dimension(:,:), intent(in ) :: flux_up, flux_dn, & !< fluxes at interfaces [W/m2]
                                             plev                !< pressure at interfaces [Pa]
    real(wp), dimension(:,:), intent(out) :: heating_rate        !< heating rate within layer [K/sec]
    character(len=128)                    :: error_msg
    ! ---------
    integer :: ncol, nlay, ilay
    ! ---------
    error_msg = ""
    ncol = size(flux_up, 1)
    nlay = size(flux_up, 2) - 1

    if(any([size(flux_dn, 1), size(flux_dn, 2)] /= [ncol, nlay+1])) &
      error_msg = "heating_rate: flux_dn array inconsistently sized."
    if(any([size(plev,    1), size(plev,     2)] /= [ncol, nlay+1])) &
      error_msg = "heating_rate: plev array inconsistently sized."
    if(any([size(heating_rate, 1), size(heating_rate, 2)] /= [ncol, nlay])) &
      error_msg = "heating_rate: heating_rate array inconsistently sized."
    if(error_msg /= "") return

    do ilay = 1, nlay
      heating_rate(1:ncol,ilay) = (flux_up(1:ncol,ilay+1) - flux_up(1:ncol,ilay) - &
                                   flux_dn(1:ncol,ilay+1) + flux_dn(1:ncol,ilay)) * &
                                  grav / (cp_dry * (plev(1:ncol,ilay+1) - plev(1:ncol,ilay)))
    end do

  end function compute_heating_rate
end module mo_heating_rates
