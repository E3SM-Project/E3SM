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
    nbytes = basesize
    kount = gdims
    npes_per_dim(:) = 1 
    !    print *,'point #3'
    if(debug) print *,'npes_per_dim: ',npes_per_dim
    n = 1
    numaiotasks = 1
    finished_block=.false.  ! indicates if we have finished contigous block
    blocksize=-1
    do n=1,ndims
       !       print *,'n: ',n
       if(decompose_dim(n)) then 
	  nbtmp = maxbytes+1  ! size of contigous block
          it=1 
	  itmp=numaiotasks      ! total number of active iotasks
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
             itmp = numaiotasks*itmp_pes
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
!                               print *,'npes_per_dim: ',npes_per_dim(n)
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
       else
 	  !------------------------------------
	  ! this is a non-decomposed dimension
 	  !------------------------------------
	  nbytes = nbytes*gdims(n)
       endif
       numaiotasks = numaiotasks*npes_per_dim(n)
    enddo

    !----------------------------------------
    ! correct decompose_dim variable based on
    ! limitations imposed by numiotask value
    !----------------------------------------
    do n=1,ndims
       if(npes_per_dim(n) == 1) then 
	  decompose_dim(n) = .false.
       endif
    enddo

    !-------------------------------------------
    ! figure out which is the last decompsed 
    ! dimension this gets treated specially when 
    ! calculating multidimensional iotask index.
    !-------------------------------------------
    per_dim=1
    do n=1,ndims
       !---------------------------------------------------
       ! calculate the multi-dimensional iotask index (indx)
       !---------------------------------------------------       
       if(iorank<numaiotasks) then
          if(decompose_dim(n)) then 
             !--------------------------------------------
             ! if a decompsed dimension calculate the indx
             !--------------------------------------------
             per_dim = per_dim*npes_per_dim(n)
             
             start(n) = MOD(iorank,npes_per_dim(n))*kount(n)+1
             rem = gdims(n) - kount(n)*(mod(iorank, npes_per_dim(n))+1)
             if(rem>0 .and. rem < kount(n)) then
                kount(n)=kount(n)+1
             end if

          else 
             !---------------------------
             ! a non-decomposed dimension
             !---------------------------
             start(n)=1
             kount(n)=gdims(n)
          endif
       else
          start(n)=1
          kount(n)=0
       end if
    enddo


    do n=1,ndims
       if((start(n)+ kount(n) -1) > gdims(n)) then
          print *,__FILE__,__LINE__,' start=',start, ' kount=',kount, ' gdims=',gdims, ' numaiotasks ',&
               numaiotasks, ' npes_per_dim ',npes_per_dim
#ifdef TESTCALCDECOMP
          stop 'bad start or count'
#else
          call piodie(__PIO_FILE__,__LINE__,'bad start or count')
#endif
       end if
    end do


    if(debug .and. iorank == 0) then 
       numOPS = NINT(real(totalSize,kind=8)/real(numaiotasks*blocksize,kind=8))
       write(*,101) 'PIO: calcdecomp: IO tasks:= ',numaiotasks,' # of ops:= ', numOPS,' size:= ',blocksize
       write(*,100) 'PIO: calcdecomp: global dimensions: ',gdims
       write(*,*) 'PIO: calcdecomp: decompse dimensions: ',decompose_dim
       write(*,100) 'PIO: calcdecomp: count: ',kount
       write(*,100) 'PIO: calcdecomp: start: ',start
    endif
100 format (a,6(i8))
101 format (a,i4,a,i5,a,i10)


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

end module calcdecomp

#ifdef TESTCALCDECOMP
program sandctest
  use calcdecomp
  implicit none
  
  integer, parameter :: ntasks=10, ndims=4
  integer, parameter :: gdims(ndims) = (/576,384,2,2/)
  integer, parameter :: num_io_procs=9

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
