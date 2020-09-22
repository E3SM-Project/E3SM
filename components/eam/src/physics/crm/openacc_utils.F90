

module openacc_utils
#ifdef _OPENACC
  use openacc
#endif
#ifdef _CUDA
  use cudafor
#endif
  implicit none
  integer :: ierr
  integer :: asyncid_loc = 1

  interface memset_async
    module procedure memset_r4_1d
    module procedure memset_r4_2d
    module procedure memset_r4_3d
    module procedure memset_r4_4d
    module procedure memset_r4_5d
    module procedure memset_r4_6d
    module procedure memset_r4_7d

    module procedure memset_r8_1d
    module procedure memset_r8_2d
    module procedure memset_r8_3d
    module procedure memset_r8_4d
    module procedure memset_r8_5d
    module procedure memset_r8_6d
    module procedure memset_r8_7d

    module procedure memset_i4_1d
    module procedure memset_i4_2d
    module procedure memset_i4_3d
    module procedure memset_i4_4d
    module procedure memset_i4_5d
    module procedure memset_i4_6d
    module procedure memset_i4_7d

    module procedure memset_i8_1d
    module procedure memset_i8_2d
    module procedure memset_i8_3d
    module procedure memset_i8_4d
    module procedure memset_i8_5d
    module procedure memset_i8_6d
    module procedure memset_i8_7d

    module procedure memset_log_1d
    module procedure memset_log_2d
    module procedure memset_log_3d
    module procedure memset_log_4d
    module procedure memset_log_5d
    module procedure memset_log_6d
    module procedure memset_log_7d
  end interface memset_async

  interface prefetch
    module procedure prefetch_r4_1d
    module procedure prefetch_r4_2d
    module procedure prefetch_r4_3d
    module procedure prefetch_r4_4d
    module procedure prefetch_r4_5d
    module procedure prefetch_r4_6d
    module procedure prefetch_r4_7d

    module procedure prefetch_r8_1d
    module procedure prefetch_r8_2d
    module procedure prefetch_r8_3d
    module procedure prefetch_r8_4d
    module procedure prefetch_r8_5d
    module procedure prefetch_r8_6d
    module procedure prefetch_r8_7d

    module procedure prefetch_i4_1d
    module procedure prefetch_i4_2d
    module procedure prefetch_i4_3d
    module procedure prefetch_i4_4d
    module procedure prefetch_i4_5d
    module procedure prefetch_i4_6d
    module procedure prefetch_i4_7d

    module procedure prefetch_i8_1d
    module procedure prefetch_i8_2d
    module procedure prefetch_i8_3d
    module procedure prefetch_i8_4d
    module procedure prefetch_i8_5d
    module procedure prefetch_i8_6d
    module procedure prefetch_i8_7d

    module procedure prefetch_log_1d
    module procedure prefetch_log_2d
    module procedure prefetch_log_3d
    module procedure prefetch_log_4d
    module procedure prefetch_log_5d
    module procedure prefetch_log_6d
    module procedure prefetch_log_7d
  end interface prefetch

contains


  subroutine memset_r4_flat(a,n,v,asyncid)
    implicit none
    real(4) :: a(n)
    real(4) :: v
    integer :: n, asyncid, i
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemsetAsync( a , v , n , acc_get_cuda_stream(asyncid) )
    !$acc end host_data
#else
    !$acc parallel loop async(asyncid)
    do i=1,n
      a(i) = v
    enddo
#endif
  end subroutine memset_r4_flat


  subroutine memset_r8_flat(a,n,v,asyncid)
    implicit none
    real(8) :: a(n)
    real(8) :: v
    integer :: n, asyncid, i
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemsetAsync( a , v , n , acc_get_cuda_stream(asyncid) )
    !$acc end host_data
#else
    !$acc parallel loop async(asyncid)
    do i=1,n
      a(i) = v
    enddo
#endif
  end subroutine memset_r8_flat


  subroutine memset_i4_flat(a,n,v,asyncid)
    implicit none
    integer(4) :: a(n)
    integer(4) :: v
    integer :: n, asyncid, i
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemsetAsync( a , v , n , acc_get_cuda_stream(asyncid) )
    !$acc end host_data
#else
    !$acc parallel loop async(asyncid)
    do i=1,n
      a(i) = v
    enddo
#endif
  end subroutine memset_i4_flat


  subroutine memset_i8_flat(a,n,v,asyncid)
    implicit none
    integer(8) :: a(n)
    integer(8) :: v
    integer :: n, asyncid, i
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemsetAsync( a , v , n , acc_get_cuda_stream(asyncid) )
    !$acc end host_data
#else
    !$acc parallel loop async(asyncid)
    do i=1,n
      a(i) = v
    enddo
#endif
  end subroutine memset_i8_flat


  subroutine memset_log_flat(a,n,v,asyncid)
    implicit none
    logical :: a(n)
    logical :: v
    integer :: n, asyncid, i
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemsetAsync( a , v , n , acc_get_cuda_stream(asyncid) )
    !$acc end host_data
#else
    !$acc parallel loop async(asyncid)
    do i=1,n
      a(i) = v
    enddo
#endif
  end subroutine memset_log_flat


  subroutine prefetch_r4_flat(a,n)
    implicit none
    real(4) :: a(n)
    integer :: n
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemPrefetchAsync( a , n , acc_get_device_num(acc_device_nvidia) , acc_get_cuda_stream(asyncid_loc+1) )
    !$acc end host_data
#endif
  end subroutine prefetch_r4_flat


  subroutine prefetch_r8_flat(a,n)
    implicit none
    real(8) :: a(n)
    integer :: n
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemPrefetchAsync( a , n , acc_get_device_num(acc_device_nvidia) , acc_get_cuda_stream(asyncid_loc+1) )
    !$acc end host_data
#endif
  end subroutine prefetch_r8_flat


  subroutine prefetch_i4_flat(a,n)
    implicit none
    integer(4) :: a(n)
    integer :: n
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemPrefetchAsync( a , n , acc_get_device_num(acc_device_nvidia) , acc_get_cuda_stream(asyncid_loc+1) )
    !$acc end host_data
#endif
  end subroutine prefetch_i4_flat


  subroutine prefetch_i8_flat(a,n)
    implicit none
    integer(8) :: a(n)
    integer :: n
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemPrefetchAsync( a , n , acc_get_device_num(acc_device_nvidia) , acc_get_cuda_stream(asyncid_loc+1) )
    !$acc end host_data
#endif
  end subroutine prefetch_i8_flat


  subroutine prefetch_log_flat(a,n)
    implicit none
    logical :: a(n)
    integer :: n
#if defined(_OPENACC) && defined (_CUDA)
    !$acc host_data use_device(a)
    ierr = cudaMemPrefetchAsync( a , n , acc_get_device_num(acc_device_nvidia) , acc_get_cuda_stream(asyncid_loc+1) )
    !$acc end host_data
#endif
  end subroutine prefetch_log_flat


  subroutine memset_r4_1d(a,v,asyncid)
    implicit none
    real(4) :: a(:)
    real(4) :: v
    integer :: asyncid
    call memset_r4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r4_1d
  subroutine memset_r4_2d(a,v,asyncid)
    implicit none
    real(4) :: a(:,:)
    real(4) :: v
    integer :: asyncid
    call memset_r4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r4_2d
  subroutine memset_r4_3d(a,v,asyncid)
    implicit none
    real(4) :: a(:,:,:)
    real(4) :: v
    integer :: asyncid
    call memset_r4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r4_3d
  subroutine memset_r4_4d(a,v,asyncid)
    implicit none
    real(4) :: a(:,:,:,:)
    real(4) :: v
    integer :: asyncid
    call memset_r4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r4_4d
  subroutine memset_r4_5d(a,v,asyncid)
    implicit none
    real(4) :: a(:,:,:,:,:)
    real(4) :: v
    integer :: asyncid
    call memset_r4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r4_5d
  subroutine memset_r4_6d(a,v,asyncid)
    implicit none
    real(4) :: a(:,:,:,:,:,:)
    real(4) :: v
    integer :: asyncid
    call memset_r4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r4_6d
  subroutine memset_r4_7d(a,v,asyncid)
    implicit none
    real(4) :: a(:,:,:,:,:,:,:)
    real(4) :: v
    integer :: asyncid
    call memset_r4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r4_7d


  subroutine memset_r8_1d(a,v,asyncid)
    implicit none
    real(8) :: a(:)
    real(8) :: v
    integer :: asyncid
    call memset_r8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r8_1d
  subroutine memset_r8_2d(a,v,asyncid)
    implicit none
    real(8) :: a(:,:)
    real(8) :: v
    integer :: asyncid
    call memset_r8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r8_2d
  subroutine memset_r8_3d(a,v,asyncid)
    implicit none
    real(8) :: a(:,:,:)
    real(8) :: v
    integer :: asyncid
    call memset_r8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r8_3d
  subroutine memset_r8_4d(a,v,asyncid)
    implicit none
    real(8) :: a(:,:,:,:)
    real(8) :: v
    integer :: asyncid
    call memset_r8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r8_4d
  subroutine memset_r8_5d(a,v,asyncid)
    implicit none
    real(8) :: a(:,:,:,:,:)
    real(8) :: v
    integer :: asyncid
    call memset_r8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r8_5d
  subroutine memset_r8_6d(a,v,asyncid)
    implicit none
    real(8) :: a(:,:,:,:,:,:)
    real(8) :: v
    integer :: asyncid
    call memset_r8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r8_6d
  subroutine memset_r8_7d(a,v,asyncid)
    implicit none
    real(8) :: a(:,:,:,:,:,:,:)
    real(8) :: v
    integer :: asyncid
    call memset_r8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_r8_7d


  subroutine memset_i4_1d(a,v,asyncid)
    implicit none
    integer(4) :: a(:)
    integer(4) :: v
    integer :: asyncid
    call memset_i4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i4_1d
  subroutine memset_i4_2d(a,v,asyncid)
    implicit none
    integer(4) :: a(:,:)
    integer(4) :: v
    integer :: asyncid
    call memset_i4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i4_2d
  subroutine memset_i4_3d(a,v,asyncid)
    implicit none
    integer(4) :: a(:,:,:)
    integer(4) :: v
    integer :: asyncid
    call memset_i4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i4_3d
  subroutine memset_i4_4d(a,v,asyncid)
    implicit none
    integer(4) :: a(:,:,:,:)
    integer(4) :: v
    integer :: asyncid
    call memset_i4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i4_4d
  subroutine memset_i4_5d(a,v,asyncid)
    implicit none
    integer(4) :: a(:,:,:,:,:)
    integer(4) :: v
    integer :: asyncid
    call memset_i4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i4_5d
  subroutine memset_i4_6d(a,v,asyncid)
    implicit none
    integer(4) :: a(:,:,:,:,:,:)
    integer(4) :: v
    integer :: asyncid
    call memset_i4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i4_6d
  subroutine memset_i4_7d(a,v,asyncid)
    implicit none
    integer(4) :: a(:,:,:,:,:,:,:)
    integer(4) :: v
    integer :: asyncid
    call memset_i4_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i4_7d


  subroutine memset_i8_1d(a,v,asyncid)
    implicit none
    integer(8) :: a(:)
    integer(8) :: v
    integer :: asyncid
    call memset_i8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i8_1d
  subroutine memset_i8_2d(a,v,asyncid)
    implicit none
    integer(8) :: a(:,:)
    integer(8) :: v
    integer :: asyncid
    call memset_i8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i8_2d
  subroutine memset_i8_3d(a,v,asyncid)
    implicit none
    integer(8) :: a(:,:,:)
    integer(8) :: v
    integer :: asyncid
    call memset_i8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i8_3d
  subroutine memset_i8_4d(a,v,asyncid)
    implicit none
    integer(8) :: a(:,:,:,:)
    integer(8) :: v
    integer :: asyncid
    call memset_i8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i8_4d
  subroutine memset_i8_5d(a,v,asyncid)
    implicit none
    integer(8) :: a(:,:,:,:,:)
    integer(8) :: v
    integer :: asyncid
    call memset_i8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i8_5d
  subroutine memset_i8_6d(a,v,asyncid)
    implicit none
    integer(8) :: a(:,:,:,:,:,:)
    integer(8) :: v
    integer :: asyncid
    call memset_i8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i8_6d
  subroutine memset_i8_7d(a,v,asyncid)
    implicit none
    integer(8) :: a(:,:,:,:,:,:,:)
    integer(8) :: v
    integer :: asyncid
    call memset_i8_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_i8_7d


  subroutine memset_log_1d(a,v,asyncid)
    implicit none
    logical :: a(:)
    logical :: v
    integer :: asyncid
    call memset_log_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_log_1d
  subroutine memset_log_2d(a,v,asyncid)
    implicit none
    logical :: a(:,:)
    logical :: v
    integer :: asyncid
    call memset_log_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_log_2d
  subroutine memset_log_3d(a,v,asyncid)
    implicit none
    logical :: a(:,:,:)
    logical :: v
    integer :: asyncid
    call memset_log_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_log_3d
  subroutine memset_log_4d(a,v,asyncid)
    implicit none
    logical :: a(:,:,:,:)
    logical :: v
    integer :: asyncid
    call memset_log_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_log_4d
  subroutine memset_log_5d(a,v,asyncid)
    implicit none
    logical :: a(:,:,:,:,:)
    logical :: v
    integer :: asyncid
    call memset_log_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_log_5d
  subroutine memset_log_6d(a,v,asyncid)
    implicit none
    logical :: a(:,:,:,:,:,:)
    logical :: v
    integer :: asyncid
    call memset_log_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_log_6d
  subroutine memset_log_7d(a,v,asyncid)
    implicit none
    logical :: a(:,:,:,:,:,:,:)
    logical :: v
    integer :: asyncid
    call memset_log_flat(a,product(shape(a)),v,asyncid)
  end subroutine memset_log_7d


  subroutine prefetch_r4_1d(a)
    implicit none
    real(4) :: a(:)
    call prefetch_r4_flat(a,product(shape(a)))
  end subroutine prefetch_r4_1d
  subroutine prefetch_r4_2d(a)
    implicit none
    real(4) :: a(:,:)
    call prefetch_r4_flat(a,product(shape(a)))
  end subroutine prefetch_r4_2d
  subroutine prefetch_r4_3d(a)
    implicit none
    real(4) :: a(:,:,:)
    call prefetch_r4_flat(a,product(shape(a)))
  end subroutine prefetch_r4_3d
  subroutine prefetch_r4_4d(a)
    implicit none
    real(4) :: a(:,:,:,:)
    call prefetch_r4_flat(a,product(shape(a)))
  end subroutine prefetch_r4_4d
  subroutine prefetch_r4_5d(a)
    implicit none
    real(4) :: a(:,:,:,:,:)
    call prefetch_r4_flat(a,product(shape(a)))
  end subroutine prefetch_r4_5d
  subroutine prefetch_r4_6d(a)
    implicit none
    real(4) :: a(:,:,:,:,:,:)
    call prefetch_r4_flat(a,product(shape(a)))
  end subroutine prefetch_r4_6d
  subroutine prefetch_r4_7d(a)
    implicit none
    real(4) :: a(:,:,:,:,:,:,:)
    call prefetch_r4_flat(a,product(shape(a)))
  end subroutine prefetch_r4_7d


  subroutine prefetch_r8_1d(a)
    implicit none
    real(8) :: a(:)
    call prefetch_r8_flat(a,product(shape(a)))
  end subroutine prefetch_r8_1d
  subroutine prefetch_r8_2d(a)
    implicit none
    real(8) :: a(:,:)
    call prefetch_r8_flat(a,product(shape(a)))
  end subroutine prefetch_r8_2d
  subroutine prefetch_r8_3d(a)
    implicit none
    real(8) :: a(:,:,:)
    call prefetch_r8_flat(a,product(shape(a)))
  end subroutine prefetch_r8_3d
  subroutine prefetch_r8_4d(a)
    implicit none
    real(8) :: a(:,:,:,:)
    call prefetch_r8_flat(a,product(shape(a)))
  end subroutine prefetch_r8_4d
  subroutine prefetch_r8_5d(a)
    implicit none
    real(8) :: a(:,:,:,:,:)
    call prefetch_r8_flat(a,product(shape(a)))
  end subroutine prefetch_r8_5d
  subroutine prefetch_r8_6d(a)
    implicit none
    real(8) :: a(:,:,:,:,:,:)
    call prefetch_r8_flat(a,product(shape(a)))
  end subroutine prefetch_r8_6d
  subroutine prefetch_r8_7d(a)
    implicit none
    real(8) :: a(:,:,:,:,:,:,:)
    call prefetch_r8_flat(a,product(shape(a)))
  end subroutine prefetch_r8_7d


  subroutine prefetch_i4_1d(a)
    implicit none
    integer(4) :: a(:)
    call prefetch_i4_flat(a,product(shape(a)))
  end subroutine prefetch_i4_1d
  subroutine prefetch_i4_2d(a)
    implicit none
    integer(4) :: a(:,:)
    call prefetch_i4_flat(a,product(shape(a)))
  end subroutine prefetch_i4_2d
  subroutine prefetch_i4_3d(a)
    implicit none
    integer(4) :: a(:,:,:)
    call prefetch_i4_flat(a,product(shape(a)))
  end subroutine prefetch_i4_3d
  subroutine prefetch_i4_4d(a)
    implicit none
    integer(4) :: a(:,:,:,:)
    call prefetch_i4_flat(a,product(shape(a)))
  end subroutine prefetch_i4_4d
  subroutine prefetch_i4_5d(a)
    implicit none
    integer(4) :: a(:,:,:,:,:)
    call prefetch_i4_flat(a,product(shape(a)))
  end subroutine prefetch_i4_5d
  subroutine prefetch_i4_6d(a)
    implicit none
    integer(4) :: a(:,:,:,:,:,:)
    call prefetch_i4_flat(a,product(shape(a)))
  end subroutine prefetch_i4_6d
  subroutine prefetch_i4_7d(a)
    implicit none
    integer(4) :: a(:,:,:,:,:,:,:)
    call prefetch_i4_flat(a,product(shape(a)))
  end subroutine prefetch_i4_7d


  subroutine prefetch_i8_1d(a)
    implicit none
    integer(8) :: a(:)
    call prefetch_i8_flat(a,product(shape(a)))
  end subroutine prefetch_i8_1d
  subroutine prefetch_i8_2d(a)
    implicit none
    integer(8) :: a(:,:)
    call prefetch_i8_flat(a,product(shape(a)))
  end subroutine prefetch_i8_2d
  subroutine prefetch_i8_3d(a)
    implicit none
    integer(8) :: a(:,:,:)
    call prefetch_i8_flat(a,product(shape(a)))
  end subroutine prefetch_i8_3d
  subroutine prefetch_i8_4d(a)
    implicit none
    integer(8) :: a(:,:,:,:)
    call prefetch_i8_flat(a,product(shape(a)))
  end subroutine prefetch_i8_4d
  subroutine prefetch_i8_5d(a)
    implicit none
    integer(8) :: a(:,:,:,:,:)
    call prefetch_i8_flat(a,product(shape(a)))
  end subroutine prefetch_i8_5d
  subroutine prefetch_i8_6d(a)
    implicit none
    integer(8) :: a(:,:,:,:,:,:)
    call prefetch_i8_flat(a,product(shape(a)))
  end subroutine prefetch_i8_6d
  subroutine prefetch_i8_7d(a)
    implicit none
    integer(8) :: a(:,:,:,:,:,:,:)
    call prefetch_i8_flat(a,product(shape(a)))
  end subroutine prefetch_i8_7d


  subroutine prefetch_log_1d(a)
    implicit none
    logical :: a(:)
    call prefetch_log_flat(a,product(shape(a)))
  end subroutine prefetch_log_1d
  subroutine prefetch_log_2d(a)
    implicit none
    logical :: a(:,:)
    call prefetch_log_flat(a,product(shape(a)))
  end subroutine prefetch_log_2d
  subroutine prefetch_log_3d(a)
    implicit none
    logical :: a(:,:,:)
    call prefetch_log_flat(a,product(shape(a)))
  end subroutine prefetch_log_3d
  subroutine prefetch_log_4d(a)
    implicit none
    logical :: a(:,:,:,:)
    call prefetch_log_flat(a,product(shape(a)))
  end subroutine prefetch_log_4d
  subroutine prefetch_log_5d(a)
    implicit none
    logical :: a(:,:,:,:,:)
    call prefetch_log_flat(a,product(shape(a)))
  end subroutine prefetch_log_5d
  subroutine prefetch_log_6d(a)
    implicit none
    logical :: a(:,:,:,:,:,:)
    call prefetch_log_flat(a,product(shape(a)))
  end subroutine prefetch_log_6d
  subroutine prefetch_log_7d(a)
    implicit none
    logical :: a(:,:,:,:,:,:,:)
    call prefetch_log_flat(a,product(shape(a)))
  end subroutine prefetch_log_7d


end module openacc_utils
