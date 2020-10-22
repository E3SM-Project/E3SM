module uvw_mod
  implicit none

contains

  subroutine uvw(ncrms)
    ! update the velocity field
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    integer :: i, j, k, icrm

    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1 , nzm
      do j = 1 , ny
        do i = 1 , nx
          do icrm = 1 , ncrms
            u(icrm,i,j,k) = dudt(icrm,i,j,k,nc)
            v(icrm,i,j,k) = dvdt(icrm,i,j,k,nc)
            w(icrm,i,j,k) = dwdt(icrm,i,j,k,nc)
          enddo
        enddo
      enddo
    enddo

  end subroutine uvw

end module uvw_mod
