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


#define REALLYOLD
  public :: CalcStartandCount

contains 
#ifdef REALLYOLD
subroutine CalcStartandCount(basetype, ndims, gdims, num_io_procs, iorank, start, count, use_io_procs)
    integer(i4), intent(in) :: ndims, num_io_procs, iorank, basetype
    integer(i4), intent(in) :: gdims(ndims)
    integer(kind=PIO_OFFSET), intent(out) :: start(ndims), count(ndims)
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
    use_io_procs = num_io_procs
    do while(p/minblocksize < use_io_procs .and. use_io_procs>1)
       use_io_procs=use_io_procs-1
    end do
    start(:)=1
    count(:)=0
    if(iorank>=use_io_procs) return 

    dims(:) = 0

    call outerdimfirst(use_io_procs, ndims, gdims,  minblocksize, maxbytes/basesize, dims)
    lb=1
    ub=ndims
    inc=1

!    print *,__LINE__,dims

    count(:)=gdims(:)
    p=1
    maxiosize = 1
    do i=lb,ub,inc
       p=p*dims(i)
       if(dims(i)>1) then
          count(i)=gdims(i)/dims(i)
          subrank = (p*iorank)/use_io_procs

          start(i)= mod(count(i)*subrank ,count(i)*dims(i)) +1
          
          extras = gdims(i)-count(i)*dims(i)
 
! We return the maxio size so that buffers for serial netcdf are always big enough
          if(extras>0) then
             maxiosize=maxiosize*(count(i)+1)
          else
             maxiosize=maxiosize*count(i)
          end if
!
! We couldnt divide into equally sized domains, we have extras
! So we add those to some of the io tasks
!           
          if(extras>0 .and. subrank >= p-extras) then
             start(i)=start(i)+ (subrank+extras-p) 
             count(i)=count(i)+1
          end if
       else
          maxiosize=maxiosize*count(i)
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
        if(extent(i)<myiocnt) then
           dims(i)=max(2,gcd(extent(i),myiocnt))
           myiocnt=max(2,myiocnt/dims(i))
        else
           dims(i)=iocnt/product(dims(i:ndims))
       endif
        tsize=product(extent)/product(dims)
!        print *,i,tsize,myiocnt,dims(i)
        i=i-1
     end do

!     print *, __LINE__,i,iocnt,myiocnt,tsize,maxblocksize,minblocksize, dims
!     stop
return

     
     done = .false.
     lextent(:)=extent(:)
     sub(ndims)=iocnt
     dims=1
     k=1
     do while(.not. done)
        p=1
        k=k+1
        do i=ndims,1,-1
           dims(i) = gcd(lextent(i),sub(i))
           p=p*dims(i)
           sub(i-1)=sub(i)/dims(i)
           if(p>=iocnt) exit
        end do
        if(i==0) then
           j = maxloc(lextent)
           lextent(j(1))=lextent(j(1))-1        
        else
          done=.true.
        end if
        if(k>100 .or. minval(lextent)< 1) then
#ifndef TESTCALCDECOMP           
           call piodie( _FILE_,__LINE__, &
                   'OuterDimFirst: Failed to converge')
#else
           stop    'OuterDimFirst: Failed to converge'
#endif
        end if

     end do

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


#else
subroutine CalcStartandCount(basetype,ndims,gdims,numiotasks,iorank,start,kount,numaiotasks)


    integer, intent(in) :: basetype 
    integer, intent(in) :: ndims
    integer, intent(in) :: gdims(:)
    integer, intent(in) :: numiotasks
    integer, intent(in) :: iorank
    integer(kind=PIO_Offset), intent(out) :: start(:)
    integer(kind=PIO_Offset), intent(out) :: kount(:)
    integer(i4),intent(out) :: numaiotasks

    logical, allocatable :: decompose_dim(:)
    integer, allocatable :: npes_per_dim(:)
    integer, allocatable :: indx(:)

    integer :: nbtmp
    integer :: it
    integer :: basesize,n
    integer :: idim
    integer :: totActive

    integer, parameter :: stripeSize = 864*1024
    integer :: minbytes = stripeSize-256   ! minimum number of contigous blocks in bytes to put on a IO task
    integer :: maxbytes = stripeSize+256   ! maximum length of contigous block in bytes to put on a IO task
    integer :: maxiter = 20
    integer(kind=PIO_offset) :: nbytes, totalSize
    integer :: numOPS
    logical :: finished_block
    integer(i8) :: blocksize
    integer :: itmp,itmp_pes
    logical :: firstdim
    integer :: lastdim
    integer :: per_dim, rem
    logical :: ltest,ltest0,ltest1

    !    print *,'point #1'
    !     print *,'gdims: ',gdims
    !    print *,'basetype: ',basetype
    !    print *,'PIO_double: ',PIO_double
    select case(basetype)
    case(PIO_int)
       basesize = 4
    case(PIO_real)
       basesize = 4 
    case(PIO_double)
       basesize = 8
    end select

    allocate(decompose_dim(ndims))
    allocate(npes_per_dim(ndims))
    allocate(indx(ndims))
    start=1
    kount=0

    !    print *,'point #2'
    !-----------------------------------------------------
    ! determine which dimensions to potentiall decomposing
    !-----------------------------------------------------
    decompose_dim(:) = .false.
    nbytes = basesize
    do n=1,ndims
       nbytes = nbytes*gdims(n)
       if(debug) then
          print *,'n: ',n
          print *,'basesize: ',basesize
          print *,'gdims(n): ',gdims(n)
          print *,'nbytes: ',nbytes
          print *,'minbytes: ',minbytes
          print *,''
       end if
       if(nbytes > minbytes) then 
          decompose_dim(n) = .true.
          !	   if(debug) print *,'Need to decompse dimension #:',n
       endif
    enddo
    totalSize=nbytes
    npes_per_dim=1
!    print *,totalsize,totalsize/minbytes,totalsize/maxbytes
    numaiotasks=min(numiotasks,totalsize/minbytes)

    do n=ndims,1,-1
       nbytes = nbytes/min(gdims(n),numaiotasks)
       if(nbytes>maxbytes) then
          npes_per_dim(n)=min(gdims(n),numaiotasks)
       else if(nbytes<minbytes) then
          itmp=1
          do while(nbytes<minbytes .and. itmp<min(numaiotasks,gdims(n)))
             nbytes = nbytes*2
             itmp=itmp+1
          end do
          npes_per_dim(n)=itmp
          exit
       end if
       if(npes_per_dim(n)>=numaiotasks) exit
    end do
    per_dim=1
    do n=ndims,1,-1
       per_dim=per_dim*gdims(n)
       if(npes_per_dim(n)>1) then
          kount(n) = gdims(n)/npes_per_dim(n)

          start(n) = kount(n)*MOD(iorank,npes_per_dim(n)) + 1
          rem = gdims(n)-kount(n)*npes_per_dim(n)
          print *,' rem = ',rem
          if(rem>0) then
             kount(n)=kount(n)+1
             start(n)=start(n)+ (iorank+rem-per_dim)
          end if
       else
          start(n) = 1
          kount(n) = gdims(n)   
       end if
    end do
!    print *,start,kount

#ifdef DOTHIS


    nbytes = basesize
    kount = gdims
    npes_per_dim(:) = 1 
    !    print *,'point #3'
    if(debug) print *,'npes_per_dim: ',npes_per_dim
    n = 1
    totActive = 1
    finished_block=.false.  ! indicates if we have finished contigous block
    blocksize=-1
    do n=1,ndims
       !       print *,'n: ',n
       if(decompose_dim(n)) then 
	  nbtmp = maxbytes+1  ! size of contigous block
          it=1 
	  itmp=totActive      ! total number of active iotasks
          if(.not. finished_block) then 
             !-------------------------------------
             ! first form a contigous block where: 
             !     minbytes < size_block < max_bytes
             !-------------------------------------
             ltest0 = (nbtmp < minbytes)
             ltest = ((nbtmp > maxbytes) .or. ltest0)
          else
             ltest0 = .false.
             ltest = .true.
          endif
          do while ((it<=maxiter) .and. ltest .and. (itmp<numiotasks) )
	     if(ltest0) then 
                !---------------------------------------
                ! This routine determines a core count 
                ! that is 10% smaller than the previous
                !---------------------------------------
		itmp_pes = nextsmaller(npes_per_dim(n),gdims(n))
             else
                !------------------------------------------
		! This subroutine tries to calculate the 
                ! next largest evenly divisible core count
                !------------------------------------------
                itmp_pes = nextlarger(npes_per_dim(n),gdims(n))
             endif
             !---------------------------------------------
             ! total number of iotasks based on this decomp
             !---------------------------------------------
             itmp = totActive*itmp_pes
             !             print *,'itmp_pes: ',itmp_pes
             !	     print *,'itmp: ',itmp
             if(itmp<=numiotasks) then 
		!------------------------------------------------
		! if we have not exceeded our iotask count 
                ! then set the values for count and npes_per_dim
		!------------------------------------------------
	        kount(n) = gdims(n)/itmp_pes
	        nbtmp = nbytes*kount(n)
                if(debug) then
                   print *,'it: ',it
                   print *,'nbtmp: ',nbtmp
                   print *,'nbytes: ',nbytes
                   print *,'maxbytes: ',maxbytes
                   print *,'' 
                end if
                npes_per_dim(n)=itmp_pes
                !                print *,'npes_per_dim: ',npes_per_dim(n)
             endif
	     !----------------------------
             ! increment iteration counter
	     !----------------------------
	     it=it+1
             if(.not. finished_block) then 
                !-------------------------------------
	        ! first form a contigous block where: 
	        !     minbytes < size_block < max_bytes
                !-------------------------------------
		ltest0 = (nbtmp < minbytes)
                ltest = ((nbtmp > maxbytes) .or. ltest0)
             else
                if(kount(n) == 1) then 
                   !---------------------------------
	           ! We can't decompose this dimension 
	           ! any further so exit.
                   !---------------------------------
		   ltest = .false.
		else 
		   !--------------------------
		   ! keep trying to decompose
                   ! dimension further 
		   !--------------------------
                   ltest = .true.
		endif
             endif
          enddo
          if(debug) print *,'npes_per_dim(n): ',npes_per_dim(n)
          !-------------------------------------
          ! we have finished the contigous block
          !-------------------------------------
          finished_block = .true.
	  idim = n
          if(blocksize<0) then 
             blocksize=nbtmp
          endif
          !---------------------------------------------
          ! calculate the number of total active iotasks
          !---------------------------------------------
          totActive = totActive*npes_per_dim(n)
       else
 	  !------------------------------------
	  ! this is a non-decomposed dimension
 	  !------------------------------------
	  nbytes = nbytes*gdims(n)
          totActive = totActive*npes_per_dim(n)
       endif
    enddo
    !JMD    if(iorank == 0) then 
    !JMD       numOPS = NINT(real(totalSize,kind=8)/real(totActive*blocksize,kind=8))
    !JMD       write(*,101) 'PIO: calcdecomp: IO tasks:= ',totActive,' # of ops:= ', numOPS,' size:= ',blocksize
    !JMD       write(*,100) 'PIO: calcdecomp: global dimensions: ',gdims
    !JMD       write(*,*) 'PIO: calcdecomp: decompse dimensions: ',decompose_dim
    !JMD       write(*,100) 'PIO: calcdecomp: count: ',kount
    !JMD       write(*,100) 'PIO: calcdecomp: start: ',start
    !JMD    endif


    !----------------------------------------
    ! correct decompose_dim variable based on
    ! limitations imposed by numiotask value
    !----------------------------------------
    do n=1,ndims
       if(npes_per_dim(n) == 1) then 
	  decompose_dim(n) = .false.
       endif
    enddo
    !    print *,'point #5'
    !    stop 'after point #5'
    !    if(debug) print *,'idim: ',idim
    !    if(debug) print *,'gdims: ',gdims(:)
    !    print *,'Dimensions to decompose: ',decompose_dim(:)
    !    print *,'numiotasks: ',numiotasks
    !    print *,'npes_per_dim: ',npes_per_dim
    !    print *,'point #6'
    !    print *,'Total number of active iotasks: ',totActive

    !    stop 'point #6'

    !-------------------------------------------
    ! figure out which is the last decompsed 
    ! dimension this gets treated specially when 
    ! calculating multidimensional iotask index.
    !-------------------------------------------
    lastdim = -1
    do n=ndims,1,-1
       if(decompose_dim(n) .and. lastdim<0) then 
          lastdim=n
       endif
    enddo

    !---------------------------------------------------
    ! calculate the multi-dimensional iotask index (indx)
    !---------------------------------------------------
    per_dim=1
    do n=ndims,1,-1
!       if(n == lastdim) then 
          !-------------------------------
          ! the last decomposed dimension 
          !-------------------------------
!          indx(n) = floor(real(iorank,kind=i8)/real(per_dim,kind=i8))+1
!       else
          !---------------------------------
          ! Not the last decompsed dimension
          !---------------------------------
          if(decompose_dim(n)) then 
             !--------------------------------------------
             ! if a decompsed dimension calculate the indx
             !--------------------------------------------
             per_dim = per_dim*npes_per_dim(n)
!             indx(n) = MOD(iorank,per_dim) + 1
             indx(n) = (iorank*per_dim)/totactive + 1
             print *,__LINE__,n,indx(n), per_dim, iorank/(totactive)/per_dim
          else
             !---------------------------
             ! a non-decomposed dimension
             !---------------------------
             indx(n)=1
!          endif
       endif
    enddo
    !   print *,'indx: ',indx
    !   print *,'point #7'
    !   stop 'point #7'
    ! Rewrite the number of active IO tasks

    ! print *,'totActive: ',totActive
    !-----------------------------------------
    ! If I am an active io task then calculate
    ! the start index and correct the count 
    ! if necessary
    !-----------------------------------------
    if(iorank < totActive) then 
       do n=1,ndims
          start(n) = (indx(n)-1)*kount(n) + 1
          !---------------------------------------
          ! if I am on the edge of all the iotasks
          ! grid check to make sure that all 
          ! gridpoints are accounted for.  
          !---------------------------------------
          if(indx(n) == npes_per_dim(n)) then 
             if((start(n)+kount(n)-1) .ne. gdims(n)) then 
                !-------------------------------------
                ! looks like the edges need a bit of 
                ! fixing up so that all values of the 
                ! array are included
                !-------------------------------------
                kount(n) = gdims(n)-start(n)+1
             endif
          endif
       enddo
    else 
       !----------------
       ! This iotask is
       ! not active
       !----------------
       kount=0
    endif
    ! stop 'point #8'
    !    print *,'point #8'
    !    print *,'iorank: ',iorank,' indx: ',indx
    !   print *,'count: ',kount
    !  stop 'point #9'

    numaiotasks=totActive
    if(debug .and. iorank == 0) then 
       numOPS = NINT(real(totalSize,kind=8)/real(totActive*blocksize,kind=8))
       write(*,101) 'PIO: calcdecomp: IO tasks:= ',totActive,' # of ops:= ', numOPS,' size:= ',blocksize
       write(*,100) 'PIO: calcdecomp: global dimensions: ',gdims
       write(*,*) 'PIO: calcdecomp: decompse dimensions: ',decompose_dim
       write(*,100) 'PIO: calcdecomp: count: ',kount
       write(*,100) 'PIO: calcdecomp: start: ',start
    endif
100 format (a,6(i8))
101 format (a,i4,a,i5,a,i10)

#endif
  end subroutine CalcStartandCount

  integer function nextsmaller(current,value) result(res)
    integer :: current, value, rem

    res = floor(real(current,kind=i8)/1.1)

  end function nextsmaller

  integer function nextlarger(current,value) result(res)

    integer :: current,value,rem
    !-----------------------------------------
    ! This function finds a value res that is
    !  larger than current that divides value evenly
    !-----------------------------------------
    res=current
    rem = 1
    do while ((rem .ne. 0) .and. (res .lt. value))
       res = res + 1
       rem = MOD(value,res)
       !     if(debug) print *,'res,rem: ',res,rem
    enddo
    if(res .eq. value) then 
       res = 2*current;
    endif
    !  print *,'rem: ',rem
    !  print *,'res: ',res
    !  print *,'current: ',current
    !  print *,'value: ',value
    !  stop 'end of nextlarger'

  end function nextlarger
#endif
end module calcdecomp

#ifdef TESTCALCDECOMP
program sandctest
  use calcdecomp
  implicit none
  
!  integer, parameter :: ndims=4
!  integer, parameter :: gdims(ndims) = (/576,384,2,7/)
  integer, parameter :: ndims=3
  integer, parameter :: gdims(ndims) = (/576,384,26/)
  integer, parameter :: num_io_procs=12

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
