module mjo_mod
  ! ====================
  use kinds
  ! ====================
  use physical_constants
  ! ====================
  use dimensions_mod, only : nlev
  ! ====================
  use quadrature_mod
  ! ====================
  use vertical_mod
  ! ====================
  use global_norms_mod
  ! ====================
  use hybrid_mod
  ! ====================
  use control_mod, only : integration
  ! ====================
  use hybvcoord_mod, only : hvcoord_t
  ! ====================
implicit none
private

  public :: mjo_cooling

contains

  subroutine mjo_cooling(cooling,hvcoord)

    type (hvcoord_t), intent(in) :: hvcoord
    real (kind=real_kind), dimension(nlev) :: cooling

    ! Local variables

    integer, parameter :: nhl  = nlev+1, nfl = nlev
    real (kind=real_kind), parameter :: href = 7.34
    real (kind=real_kind), parameter :: hmx = 12.0
    real (kind=real_kind), parameter :: hmn = 15.0
    real (kind=real_kind) :: zcoef,func32,a,b,c

    real (kind=real_kind), dimension(nlev+1) :: shalf
    real (kind=real_kind), dimension(nlev)   :: sfull, zfull
    real (kind=real_kind) :: func23
    integer :: k

    func23(a,b,c)=max(0.0_real_kind,min(1.0_real_kind,(b-a)/(b-c)))

    do k=1,nhl  ! nlev+1
       shalf(k) = (k-1)*1.0_real_kind/REAL(nlev,kind=real_kind)
    end do

    sfull = ( shalf(2:nhl) + shalf(1:nhl-1) ) / 2.0_real_kind
    zfull = -href * log(sfull)

    do k=1,nfl
     ZCOEF=func23(zfull(k),HMX,HMN)
     cooling(k) = -1.5_real_kind/(24.0_real_kind*3600.0_real_kind)*ZCOEF  ! [K/sec]
    end do

  end subroutine mjo_cooling

end module mjo_mod
