#define __PIO_FILE__ "calcdecomp.F90"
module calcdecomp
#ifdef TESTCALCDECOMP
  implicit none 
  integer, parameter :: i4=selected_int_kind(6), &
       i8=selected_int_kind(13), pio_offset=i8, r8=selected_real_kind(13)
  logical, parameter :: debug=.false.
  integer, parameter :: pio_real=1,pio_int=2, pio_double=3
#else
  use pio_kinds, only: i4, r4,r8,i4,i8, PIO_offset
  use pio_types, only: PIO_int, PIO_real, PIO_double
  use pio_support, only : debug, piodie
  implicit none 
#endif

  public :: CalcStartandCount

contains 

  subroutine CalcStartandCount(basetype, ndims, gdims, num_io_procs, iorank, start, kount, use_io_procs)
    integer(i4), intent(in) :: ndims, num_io_procs, iorank, basetype
    integer(i4), intent(in) :: gdims(ndims)
    integer(kind=PIO_OFFSET), intent(out) :: start(ndims), kount(ndims)
    integer, intent(out) :: use_io_procs
    integer :: i, p, dims(ndims), lb, ub, inc
    integer :: extras, subrank

    integer, parameter :: stripeSize = 864*1024
    integer :: minbytes = stripeSize-256   ! minimum number of contigous blocks in bytes to put on a IO task
    integer :: maxbytes = stripeSize+256   ! maximum length of contigous block in bytes to put on a IO task

    integer :: minblocksize, basesize, maxiosize


    select case(basetype)
    case(PIO_int)
       basesize = 4
    case(PIO_real)
       basesize = 4 
    case(PIO_double)
       basesize = 8
    end select

    minblocksize = minbytes/basesize

    p=1
    do i=1,ndims
       p=p*gdims(i)
    end do
    use_io_procs = max(1,min(ceiling(real(p)/real(minblocksize)),num_io_procs))

    start(:)=1
    kount(:)=0
    if(iorank>=use_io_procs) return 

    dims(:) = 0

    call outerdimfirst(use_io_procs, ndims, gdims,  minblocksize, maxbytes/basesize, dims)

    !    print *,__LINE__,dims

    kount(:)=gdims(:)
    p=1
    maxiosize = 1
    do i=1,ndims
       p=p*dims(i)
       if(dims(i)>1) then
          kount(i)=gdims(i)/dims(i)
          subrank = (p*iorank+1)/use_io_procs

          start(i)= mod(kount(i)*subrank ,kount(i)*dims(i)) +1

          extras = gdims(i)-kount(i)*dims(i)

          ! We return the maxio size so that buffers for serial netcdf are always big enough
          !
          ! We couldnt divide into equally sized domains, we have extras
          ! So we add those to some of the io tasks
          !
	  if(i==4) print *, iorank,extras,subrank,p,use_io_procs           
          if(extras>0 .and. subrank >= p-extras) then
             start(i)=start(i)+ (subrank+extras-p) 
             kount(i)=kount(i)+1
          end if
          
          if(mod(kount(i)*(p*(iorank+1)/use_io_procs),kount(i)*dims(i))+1<start(i)) then
             kount(i) = gdims(i)-start(i)+1
          end if
               

       end if
    end do

    

  end subroutine Calcstartandcount

  subroutine outerdimfirst(iocnt,ndims,extent,minblocksize,maxblocksize,dims)
    implicit none
    integer, intent(in) :: iocnt,ndims
    integer, intent(in) :: extent(ndims),minblocksize,maxblocksize
    integer, intent(out) :: dims(ndims)

    integer :: sub(0:ndims)
    integer :: lextent(ndims)
    integer :: i, p, j(1), k, tsize, myiocnt
    logical :: done

    tsize = product(extent)
    dims = 1

    myiocnt=iocnt

    i=ndims

    do while(tsize>maxblocksize .and. i>0)
       myiocnt=iocnt
       if(extent(i)>1) then
          if(extent(i)<myiocnt) then
             dims(i)=max(2,gcd(extent(i),myiocnt))
             myiocnt=max(2,myiocnt/dims(i))
          else
             dims(i)=iocnt/product(dims(i:ndims))
          endif
          tsize=product(extent)/product(dims)
       end if
       !        print *,i,tsize,myiocnt,dims(i)
       i=i-1
    end do

    !     print *, __LINE__,i,iocnt,myiocnt,tsize,maxblocksize,minblocksize, dims
    !     stop
    return

  end subroutine outerdimfirst
  recursive function gcd( a, b ) result( divisor )

    ! Note that the function itself is not declared when the result
    ! variable is present.  The type of the function is the type of 
    ! the result variable.  Thus, only the result variable may be 
    ! declared. 

    integer divisor
    integer, intent(in) :: a, b

    integer m, n

    ! Multiple statements may be written on a single source line
    ! provided they are delimited with semicolons.

    m = abs(a); n = abs(b)
    if ( m > n ) call swap( m, n )  ! Insure that m <= n.

    ! When the function invokes itself recursively, the result variable
    ! should be used to store the result of the function.  The function
    ! name is used to invoke the function.  Thus, the function name should
    ! not appear on the left-hand side of an assignment statement.

    if ( m == 0 ) then 
       divisor = n
    else
       divisor = gcd( mod( n, m ), m )
    end if

  end function gcd
  subroutine swap( x, y )
    integer, intent(inout) :: x, y

    integer tmp

    tmp = x; x = y; y = tmp 
  end subroutine swap

end module calcdecomp

#ifdef TESTCALCDECOMP
program sandctest
  use calcdecomp
  implicit none
  
!  integer, parameter :: ndims=4
!  integer, parameter :: gdims(ndims) = (/144,96,30,3/)
  integer, parameter :: ndims=3
  integer, parameter :: gdims(ndims) = (/576,384,1/)
  integer, parameter :: num_io_procs=16

  integer :: psize, n

  integer(kind=pio_offset) :: start(ndims), count(ndims)
  integer :: iorank, numaiotasks

  do iorank=0,num_io_procs-1
     call Calcstartandcount(PIO_double, ndims, gdims, num_io_procs, iorank, start, count, numaiotasks)

     psize=1
     do n=1,ndims
        psize=psize*count(n)
     end do
     if(ndims==3) then
        write(*,'(i2,a,3i5,a,3i5,2i12)') iorank,' start =',start,' count=', count, product(gdims), psize
     else if(ndims==4) then
        write(*,'(i2,a,4i5,a,4i5,2i12)') iorank,' start =',start,' count=', count, product(gdims), psize
     end if
  end do
end program sandctest
#endif
