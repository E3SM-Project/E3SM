
module zero_mod
  use params, only: asyncid
  implicit none

contains

  subroutine zero(ncrms)
    use vars
    implicit none
    integer, intent(in) :: ncrms
    integer k,icrm, j, i
    !dudt(ncrms,nxp1, ny  , nzm, 3)
    !dvdt(ncrms,nx  , nyp1, nzm, 3)
    !dwdt(ncrms,nx  , ny  , nz , 3)
    !misc(ncrms,nx  , ny  , nz )
    
    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1 , nz
      do j = 1 , nyp1
        do i = 1 , nxp1
          do icrm = 1 , ncrms
            if (i <= nxp1 .and. j <= ny   .and. k <= nzm) dudt(icrm,i,j,k,na) = 0.
            if (i <= nx   .and. j <= nyp1 .and. k <= nzm) dvdt(icrm,i,j,k,na) = 0.
            if (i <= nx   .and. j <= ny   .and. k <= nz ) dwdt(icrm,i,j,k,na) = 0.
            if (i <= nx   .and. j <= ny   .and. k <= nz ) misc(icrm,i,j,k) = 0.
          enddo
        enddo
      enddo
    enddo
  end

end module zero_mod
