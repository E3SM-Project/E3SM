module calcdecomp
  use pio_kinds, only: i4, r4,r8,i4,i8, PIO_offset
  use pio_types, only: PIO_int, PIO_real, PIO_double

 public :: CalcStartandCount

contains 

  subroutine CalcStartandCount(numtask,basetype,ndims,gdims,numiotasks,iorank,start,count,numaiotasks)

    implicit none 


    integer :: numtask
    integer :: basetype 
    integer :: ndims
    integer :: gdims(:)
    integer :: numiotasks
    integer :: iorank
    integer(kind=PIO_Offset) :: start(:)
    integer(kind=PIO_Offset) :: count(:)
    integer(i4) :: numaiotasks
    
    logical, allocatable :: decompose_dim(:)
    integer, allocatable :: npes_per_dim(:)
    integer, allocatable :: indx(:)
   
    integer :: nbtmp
    integer :: it
    integer :: basesize,n
    integer :: idim
    integer :: totActive

    integer :: minbytes = (1024-256)*1024      ! minimum number of contigous blocks in bytes to put on a IO task
    integer :: maxbytes = (1024+256)*1024   ! maximum length of contigous block in bytes to put on a IO task
    integer :: maxiter = 10
    integer(kind=PIO_offset) :: nbytes
    logical :: finished_block
    logical, parameter :: verbose = .false.
    integer :: itmp,itmp_pes
    logical :: firstdim
    integer :: lastdim
    integer :: per_dim
    logical :: ltest,ltest0,ltest1

    iorank = iorank+1 
    print *,'numtask: ',numtask
    print *,'gdims: ',gdims
    print *,'numiotasks: ',numiotasks
    print *,'iorank: ',iorank
!    print *,'point #1'
!     print *,'gdims: ',gdims
    print *,'basetype: ',basetype
    print *,'PIO_double: ',PIO_double
     select case(basetype)
	case(PIO_int)
	    basesize = 4
        case(PIO_real)
           basesize = 4 
        case(PIO_double)
           basesize = 8
    end select
    ! KLUDGE this for now 
!    basesize=8


    allocate(decompose_dim(ndims))
    allocate(npes_per_dim(ndims))
    allocate(indx(ndims))
    start=1
    count=0
     
!    print *,'point #2'
    !-----------------------------------------------------
    ! determine which dimensions to potentiall decomposing
    !-----------------------------------------------------
    decompose_dim(:) = .false.
    nbytes = basesize
    do n=1,ndims
       nbytes = nbytes*gdims(n)
       if(verbose) print *,'n: ',n
       if(verbose) print *,'basesize: ',basesize
       if(verbose) print *,'gdims(n): ',gdims(n)
       if(verbose) print *,'nbytes: ',nbytes
       if(verbose) print *,'minbytes: ',minbytes
       if(verbose) print *,''
       if(nbytes > minbytes) then 
	   decompose_dim(n) = .true.
!	   if(verbose) print *,'Need to decompse dimension #:',n
       endif
    enddo
    nbytes = basesize
    count = gdims
    npes_per_dim(:) = 1 
!    print *,'point #3'
    if(verbose) print *,'npes_per_dim: ',npes_per_dim
    n = 1
    totActive = 1
    finished_block=.false.  ! indicates if we have finished contigous block
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
	        count(n) = gdims(n)/itmp_pes
	        nbtmp = nbytes*count(n)
	        if(verbose) print *,'it: ',it
	        if(verbose) print *,'nbtmp: ',nbtmp
                if(verbose) print *,'nbytes: ',nbytes
                if(verbose) print *,'maxbytes: ',maxbytes
	        if(verbose) print *,'' 
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
                if(count(n) == 1) then 
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
          if(verbose) print *,'npes_per_dim(n): ',npes_per_dim(n)
          !-------------------------------------
          ! we have finished the contigous block
          !-------------------------------------
          finished_block = .true.
	  idim = n
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
!     print *,''
    enddo
    print *,'point #4: totActive',totActive
    print *,'decompse_dim: ',decompose_dim
    print *,'basesize: ',basesize
    print *,'count: ',count
!    stop 'point #4'

!    print *,'count: ',count
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
!    if(verbose) print *,'idim: ',idim
!    if(verbose) print *,'gdims: ',gdims(:)
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
   do n=1,ndims
      if(n == lastdim) then 
	 !-------------------------------
	 ! the last decomposed dimension 
	 !-------------------------------
         indx(n) = floor(real(iorank-1,kind=i8)/real(per_dim,kind=i8))+1
      else
	 !---------------------------------
         ! Not the last decompsed dimension
	 !---------------------------------
         if(decompose_dim(n)) then 
           !--------------------------------------------
	   ! if a decompsed dimension calculate the indx
           !--------------------------------------------
           per_dim = per_dim*npes_per_dim(n)
	   indx(n) = MOD(iorank-1,per_dim) + 1
         else
	   !---------------------------
           ! a non-decomposed dimension
	   !---------------------------
           indx(n)=1
         endif
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
 if(iorank <= totActive) then 
    do n=1,ndims
       start(n) = (indx(n)-1)*count(n) + 1
       !---------------------------------------
       ! if I am on the edge of all the iotasks
       ! grid check to make sure that all 
       ! gridpoints are accounted for.  
       !---------------------------------------
       if(indx(n) == npes_per_dim(n)) then 
          if((start(n)+count(n)-1) .ne. gdims(n)) then 
	     !-------------------------------------
	     ! looks like the edges need a bit of 
	     ! fixing up so that all values of the 
             ! array are included
	     !-------------------------------------
	     count(n) = gdims(n)-start(n)+1
          endif
      endif
    enddo
 else 
    !----------------
    ! This iotask is
    ! not active
    !----------------
    count=0
 endif
! stop 'point #8'
!    print *,'point #8'
!    print *,'iorank: ',iorank,' indx: ',indx
!   print *,'count: ',count
!  stop 'point #9'
  
  numaiotasks=totActive
  iorank=iorank-1
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
!     if(verbose) print *,'res,rem: ',res,rem
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
