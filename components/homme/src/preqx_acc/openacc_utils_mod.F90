
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module openacc_utils_mod
  use kinds, only: real_kind
  use dimensions_mod, only: nelemd
  implicit none
  private

  public :: copy_qdp_h2d
  public :: copy_qdp_d2h
  public :: update_host_async
  public :: update_device_async
  public :: copy_ondev
  public :: copy_ondev_async
  public :: acc_async_test_wrap

contains

  function acc_async_test_wrap( asyncid )  result(rslt)
#   ifdef _OPENACC
      use openacc, only: acc_async_test
#   endif
    implicit none
    integer, intent(in) :: asyncid
    logical             :: rslt
#   ifdef _OPENACC
      rslt = .false.
      rslt = acc_async_test(asyncid)
#   else
      rslt = .true.
#   endif
  end function acc_async_test_wrap

  subroutine copy_qdp_h2d( elem , tl )
    use element_mod, only: element_t, state_qdp
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer        , intent(in) :: tl
    integer :: ie
    !$omp barrier
    !$omp master
    do ie = 1 , nelemd
      !$acc update device(state_qdp(:,:,:,:,tl,ie))
    enddo
    !$omp end master
    !$omp barrier
  end subroutine copy_qdp_h2d

  subroutine copy_qdp_d2h( elem , tl )
    use element_mod, only: element_t, state_qdp
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer        , intent(in) :: tl
    integer :: ie
    !$omp barrier
    !$omp master
    do ie = 1 , nelemd
      !$acc update host(state_qdp(:,:,:,:,tl,ie))
    enddo
    !$omp end master
    !$omp barrier
  end subroutine copy_qdp_d2h

  subroutine copy_ondev(dest,src,len)
    implicit none
    real(kind=real_kind), intent(  out) :: dest(len)
    real(kind=real_kind), intent(in   ) :: src (len)
    integer             , intent(in   ) :: len
    integer :: i
    !$acc parallel loop gang vector present(src,dest)
    do i = 1 , len
      dest(i) = src(i)
    enddo
  end subroutine copy_ondev

  subroutine update_host_async(dat,len,id)
    implicit none
    real(kind=real_kind) :: dat(len)
    integer              :: len, id
    !$acc update host(dat(1:len)) async(id)
  end subroutine update_host_async

  subroutine update_device_async(dat,len,id)
    implicit none
    real(kind=real_kind) :: dat(len)
    integer              :: len, id
    !$acc update device(dat(1:len)) async(id)
  end subroutine update_device_async

  subroutine copy_ondev_async(dest,src,len,id)
    implicit none
    real(kind=real_kind) :: dest(len)
    real(kind=real_kind) :: src (len)
    integer              :: len, id
    integer :: i
    !$acc parallel loop gang vector present(dest,src) async(id)
    do i = 1 , len
      dest(i) = src(i)
    enddo
  end subroutine copy_ondev_async

end module openacc_utils_mod

