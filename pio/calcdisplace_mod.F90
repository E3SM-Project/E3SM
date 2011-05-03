#define __PIO_FILE__ "calcdisplace_mod.F90"
MODULE calcdisplace_mod

  use pio_kinds, only: i4, PIO_OFFSET
  use pio_support, only : piodie

  private
  public :: GCDblocksize,gcd
  public :: calcdisplace, calcdisplace_box

  interface gcd
     module procedure gcd_array
     module procedure gcd_pair
  end interface

CONTAINS
  !*****************************
  ! calcdisplace
  !

  subroutine calcdisplace(bsize,dof,displace)

    integer(i4), intent(in) :: bsize    ! length of contigious blocks of numbers
    integer(i4), intent(in) :: dof(:)   ! degree of freedom on which to setup the displacement array
    integer(i4), intent(inout) :: displace(:)  ! array of mpi displacments

    integer :: numblocks,lenblocks,i,ii,dis

    numblocks = size(displace)
    lenblocks = bsize
    do i=1,numblocks
       ii = (i-1)*lenblocks+1
       dis = dof(ii)-1
       dis = dis/lenblocks   
       displace(i) = dis
    enddo

  end subroutine calcdisplace



  subroutine calcdisplace_box(gsize,start,count,ndim,displace)

    integer(i4),intent(in) :: gsize(:)   ! global size of output domain
    integer(kind=PIO_offset),intent(in) :: start(:), count(:)
    integer(i4), intent(in) :: ndim
    integer(i4),intent(inout) :: displace(:)  ! mpi displacments

    !!

    integer ndisp
    integer(i4) :: gstride(ndim)
    integer i,j
    integer iosize
    integer(i4) :: myloc(ndim)
    integer(i4) :: ub(ndim)
    integer idim
    logical done
    integer gindex

    gstride(1)=gsize(1)
    do i=2,ndim
       gstride(i)=gsize(i)*gstride(i-1)
    end do

    iosize=min(int(count(1)),1)
    do i=2,ndim
       iosize=iosize*count(i)
    end do

    ndisp=size(displace)


    if (iosize<1 .or. ndisp<1) return

    if (ndisp/=iosize) then
       call piodie(__PIO_FILE__,__LINE__,'ndisp=',ndisp,' /= iosize=',iosize)
    endif

    do i=1,ndim
       ub(i)=start(i)+count(i)-1
    end do

    ! skip x dimension (start of each block)
    ! generate displacement for every 1,y,z
    !  i.e. loop over y,z,...
    !       compute corresponding global index
    !       divide by lenblocks

    displace(1)=1
    myloc=start

    do i=1,iosize
       ! go from myloc() to 1-based global index
       gindex=myloc(1)
       do j=2,ndim
          gindex=gindex+(myloc(j)-1)*gstride(j-1)
       end do

       ! rml
       ! following original but is that right???
       ! seems like the 'if' is erroneous

       gindex=gindex-1

       gindex=gindex/count(1)    ! gindex/lenblock

       displace(i)=gindex

       ! increment myloc to next position


       idim=2                    ! dimension to increment
       done=.false.

       if (i<iosize) then
          do while (.not. done)
             if (myloc(idim)<ub(idim)) then
                myloc(idim)=myloc(idim)+1
                done=.true.
             else
                myloc(idim)=start(idim)
                idim=idim+1
                if (idim>ndim) call piodie(__PIO_FILE__,__LINE__,'dim overflow')
             endif
          end do
       endif

    end do

    do i=2,ndim
       if (myloc(i) /= ub(i)) then
          print *,'myloc=',myloc
          print *,'ub=',ub
          call piodie( __PIO_FILE__,__LINE__,'myloc/=ub')
       endif
    end do


    ! check for strictly increasing

    do i=1,ndisp-1	
       if(displace(i) .gt. displace(i+1)) then
          call piodie(__PIO_FILE__,__LINE__,'displace is not increasing')
       endif
    enddo

  end subroutine calcdisplace_box


  SUBROUTINE GCDblocksize(arr_in,bsize)
    implicit none

    integer(i4),intent(in) ,dimension(:)  :: arr_in  !arr_in = rindex array from box_rearrange
    integer(i4),intent(out)               :: bsize   ! the gcd of the block length array

    ! Locals
    integer(i4),dimension(:),allocatable   :: del_arr,loc_arr,blk_len
    integer(i4),dimension(:),allocatable :: gaps
    integer(i4) :: i,j,k,n,numblks,numtimes,tloc,bsizeg,ii


    n = size(arr_in)

    allocate(del_arr(n-1))

    del_arr = 0    ! forward diff of the input array to deterine where contiguous blocks end.


    do  i = 1,n-1  ! compute foward diff of the elements; if =1, still in a contiguous block,
       ! if /= 1 , the end of a block has been reached. 

       del_arr(i) = (arr_in(i+1) - arr_in(i))


    end do
!    print *,'del_arr: ',del_arr

    numtimes = count( del_arr /= 1) 
    numblks  = numtimes + 1   ! the number of contiguous blocks.
     

    if ( numtimes == 0 ) then    ! new logic to account for the case that there is only
       allocate(loc_arr(numblks))  ! one contigious block in which case numtimes=0 and the 
    else                         ! error from the assignment in line 87 goes away
       allocate(loc_arr(numtimes))
       allocate(gaps(numtimes))
    end if
    loc_arr = 1

    j=0

    do i = 1, n-1

       if ( del_arr(i) == 1 ) cycle

       tloc = i

       j = j+1

       loc_arr(j) = tloc

    end do
    
    if(numtimes>0) then 
       ii=1
       do i=1,n-1
         if(del_arr(i) .gt. 1) then
            gaps(ii) = del_arr(i) -1
	    ii=ii+1
         endif
       enddo
    endif

    allocate(blk_len(numblks))
    blk_len(1) = loc_arr(1)

    do k = 2,numblks-1           ! computes the the length of each block by differencing the 
       ! eleemts of the res array. 

       blk_len(k)  = loc_arr(k) - loc_arr(k-1)

    end do

    blk_len(numblks) = n - sum(blk_len(1:numblks-1)) ! computes the length of the last block


    bsize = gcd_array(blk_len) ! call to compute the gcd of the blk_len array.    
    if(numtimes>0) then 
       bsizeg = gcd_array(gaps(1:numtimes)) 
       bsize = gcd_pair(bsize,bsizeg)
       if(arr_in(1)>0) then    ! account for an initial gap
          bsize = gcd_pair(bsize,arr_in(1))
       end if
       deallocate(gaps)
    endif
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
    if(n==0 .or. any(ain <= 1)) return

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

    if(a < b) then
       x = a
       a = b
       b = x
    end if

    do 
       x = mod(a,b)
       if ( x == 0 ) EXIT
       a = b 
       b = x 
    end do

    gcd = b 

  end FUNCTION gcd_pair

END MODULE calcdisplace_mod
