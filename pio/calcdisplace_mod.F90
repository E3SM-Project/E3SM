MODULE calcdisplace_mod

  use pio_kinds, only: i4
  use pio_support, only : piodie

  private
  public :: GCDblocksize,gcd
  public :: calcdisplace

  interface gcd
     module procedure gcd_array
     module procedure gcd_pair
  end interface

CONTAINS

  subroutine calcdisplace(bsize,dof,displace)
    implicit none
    integer(i4), intent(in) :: bsize    ! length of contigious blocks of numbers
    integer(i4), intent(in) :: dof(:)   ! degree of freedom on which to setup the displacement array
    integer(i4), intent(inout) :: displace(:)  ! array of mpi displacments

    integer(i4) :: numblocks,lenblocks,i,ii,dis

    numblocks = size(displace)
    lenblocks = bsize
    do i=1,numblocks
       ii = (i-1)*lenblocks+1
       dis = dof(ii)-1
       dis = dis/lenblocks   
       displace(i) = dis
    enddo
!    if (numblocks > 1) then
!       do i=1,numblocks-1	
!          if(displace(i+1) .lt. displace(i)) then
!             print *,'calcdisplace: monotonic error: displace(i:i+1)',displace(i:i+1)
!             print *,'calcdisplace: monotonic error: i,lenblocks,numblocks,size(dof): ',i,lenblocks,numblocks,size(dof)
!             !         call piodie( _FILE_,__LINE__)
!          endif
!       enddo
!    endif

  end subroutine calcdisplace

  SUBROUTINE GCDblocksize(arr_in,bsize)
    implicit none

    integer(i4),intent(in) ,dimension(:)  :: arr_in  !arr_in = rindex array from box_rearrange
    integer(i4),intent(out)               :: bsize   ! the gcd of the block length array

    ! Locals
    integer(i4),dimension(:),allocatable   :: del_arr,loc_arr,blk_len
    integer(i4) :: i,j,k,n,numblks,numtimes,tloc


    n = size(arr_in)

    allocate(del_arr(n-1))

    del_arr = 0    ! forward diff of the input array to deterine where contiguous blocks end.


    do  i = 1,n-1  ! compute foward diff of the elements; if =1, still in a contiguous block,
       ! if /= 1 , the end of a block has been reached. 

       del_arr(i) = (arr_in(i+1) - arr_in(i))


    end do

    numblks  = count( del_arr /= 1) + 1   ! the number of contiguous blocks.
    numtimes = count( del_arr /= 1) 

    if ( numtimes == 0 ) then    ! new logic to account for the case that there is only
       allocate(loc_arr(numblks))  ! one contigious block in which case numtimes=0 and the 
    else                         ! error from the assignment in line 87 goes away
       allocate(loc_arr(numtimes))
    end if
    loc_arr = 1

    j=0

    do i = 1, n-1

       if ( del_arr(i) == 1 ) cycle

       tloc = i

       j = j+1

       loc_arr(j) = tloc

    end do

    allocate(blk_len(numblks))
    blk_len(1) = loc_arr(1)

    do k = 2,numblks-1           ! computes the the length of each block by differencing the 
       ! eleemts of the res array. 

       blk_len(k)  = loc_arr(k) - loc_arr(k-1)

    end do

    blk_len(numblks) = n - sum(blk_len(1:numblks-1)) ! computes the length of the last block


    bsize = gcd_array(blk_len) ! call to compute the gcd of the blk_len array.    
    deallocate(del_arr,loc_arr,blk_len)

  end SUBROUTINE GCDblocksize


  integer function gcd_array(ain) result(bsize)
    implicit none


    integer(i4), intent(in),dimension(:) :: ain

    ! locals
    integer(i4),allocatable :: gcr(:) 
    integer(i4) :: i,j,a_min,x,index,n

    bsize=1
    n = size(ain) 
    ! First check, if an element is 1, then 1 is the gcd (i.e bsize)
    if(n==0 .or. any(ain == 1)) return

    allocate(gcr(n))
    gcr = 0
    !

    ! Find and set the min value.

    a_min = minval(ain)


    !print*,"amin = ", a_min

    ! Now compute the gcd between a_min and the rest of the array elements as long as a_min /= 1
    ! otherwise gcd = a_min = 1.
    ! Done by calling the external function that is below.

    do i = 1,n

       gcr(i) = gcd_pair(a_min,ain(i))

    end do

    !
    ! Now look for the smallest value in the gcr array and and assign it to bsize.
    !

    bsize = minval(gcr)

    deallocate(gcr)

  end function gcd_array

  integer FUNCTION gcd_pair(u,v) result(gcd)
    implicit none

    integer(i4),intent(in) :: u,v

    ! locals
    integer(i4) :: x,a,b

    a = u
    b = v

    do 

       x = mod(a,b)
       if ( x == 0 ) EXIT
       gcd = b 

       a = b 
       b = x 

       gcd = x 

    end do

    gcd = b 

  end FUNCTION gcd_pair

END MODULE calcdisplace_mod
