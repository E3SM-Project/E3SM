module clmgridMod
!
! DESCRIPTION
!
! initialize computation grid
  use shr_kind_mod , only: r8 => shr_kind_r8
implicit none

  private
  real(r8), public, pointer :: zsoi(:)  => null() !soil depth, node center 1 : nlevsoi
  real(r8), public, pointer :: zisoi(:) => null() !soil depth, interface,  0 : nlevsoi
  real(r8), public, pointer :: dzsoi(:) => null() !soil layer thickness

  public :: init_clm_vertgrid
contains


  subroutine init_clm_vertgrid(nlevgrnd)
  !
  ! DESCRIPTION
  ! initialize the vertical grid for computation
  !
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none
  integer, intent(in) :: nlevgrnd

  real(r8)              :: scalez = 0.025_r8        ! Soil layer thickness discretization (m)
  integer :: j

  allocate(zsoi(1:nlevgrnd))
  allocate(zisoi(0:nlevgrnd))
  allocate(dzsoi(1:nlevgrnd))

  do j = 1, nlevgrnd
    zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
  enddo

  dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces

  do j = 2,nlevgrnd-1
    dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
  enddo
  dzsoi(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)

  zisoi(0) = 0._r8
  do j = 1, nlevgrnd-1
    zisoi(j) = 0.5_r8*(zsoi(j)+zsoi(j+1))         !interface depths
  enddo
  zisoi(nlevgrnd) = zsoi(nlevgrnd) + 0.5_r8*dzsoi(nlevgrnd)

  end subroutine init_clm_vertgrid

end module clmgridMod
