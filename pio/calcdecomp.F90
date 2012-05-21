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

  public :: CalcStartandCount, pio_set_blocksize
  integer, parameter :: default_blocksize=884736
  integer :: blocksize=default_blocksize
contains 
  subroutine pio_set_blocksize(newsize)
    integer, intent(in) :: newsize
#ifndef TESTCALCDECOMP
    if(newsize<0) then
       call piodie(__PIO_FILE__,__LINE__,'bad value to blocksize: ',newsize)
    end if
#endif
    blocksize=newsize


  end subroutine pio_set_blocksize


!
!  Determine start and kount values for an array of global size gdims over at most num_io_procs tasks.
!  The algorythm creates contigous blocks of approximate size stripesize.  Blocksize should be adjusted
!  to be optimal for the filesystem being used.   The actual number of io tasks used is output in variable
!  use_io_procs
!
  subroutine CalcStartandCount(basetype, ndims, gdims, num_io_procs, iorank, start, kount, use_io_procs, innermostdecomposed)
    integer(i4), intent(in) :: ndims, num_io_procs, iorank, basetype
    integer(i4), intent(in) :: gdims(ndims)
    integer(kind=PIO_OFFSET), intent(out) :: start(ndims), kount(ndims)
    integer, intent(out) :: use_io_procs
    integer, intent(out), optional :: innermostdecomposed
    integer :: i,  dims(ndims), lb, ub, inc
    integer(kind=pio_offset) :: p
    integer :: extras, subrank, tioprocs, rem
    integer :: minbytes, maxbytes

    integer :: minblocksize, basesize, maxiosize, ioprocs, tiorank
    integer :: ldims


    minbytes = blocksize-256   ! minimum number of contigous blocks in bytes to put on a IO task
    maxbytes = blocksize+256   ! maximum length of contigous block in bytes to put on a IO task

    select case(basetype)
    case(PIO_int)
       basesize = 4
    case(PIO_real)
       basesize = 4 
    case(PIO_double)
       basesize = 8
    end select

    minblocksize = minbytes/basesize

    p=product(gdims)
    use_io_procs = max(1, min(int(real(p)/real(minblocksize)+0.5),num_io_procs))




    start(:)=1
    kount(:)=0

    ldims=ndims
    p=basesize
    do i=1,ndims
       p=p*gdims(i)
       if(p/use_io_procs > maxbytes) then
          ldims=i
          exit
       end if
    end do

! Things work best if use_io_procs is a multiple of gdims(ndims)
! this adjustment makes it so, potentially increasing the blocksize a bit
    if (gdims(ldims)<use_io_procs) then
       do while(mod(use_io_procs,gdims(ldims))>0)
          use_io_procs=use_io_procs-1
       end do
    end if
    if(iorank>=use_io_procs) return 


    kount=gdims

    ioprocs=use_io_procs
    tiorank=iorank

    do i=ldims,1,-1
       if(gdims(i)>1) then
          if(gdims(i)>=ioprocs) then
             
             call computestartandcount(gdims(i),ioprocs,tiorank,start(i),kount(i))
             if(start(i)+kount(i)>gdims(i)+1) then
                print *,__PIO_FILE__,__LINE__,i,ioprocs,gdims(i),start(i),kount(i)
#if TESTCALCDECOMP
                stop
#else
                call piodie(__PIO_FILE__,__LINE__,'Start plus count exceeds dimension bound')
#endif
             endif
             exit  ! Decomposition is complete
          else
             ! The current dimension cannot complete the decomposition.   Decompose this 
             ! dimension in groups then go on to decompose the next dimesion in each of those
             ! groups.
             tioprocs=gdims(i)
             tiorank = (iorank*tioprocs)/ioprocs
             
             call computestartandcount(gdims(i),tioprocs, tiorank  , start(i),kount(i))
             ioprocs=ioprocs/tioprocs
             tiorank = mod(iorank,ioprocs)
          end if
       end if
    end do
    if(present(innermostdecomposed)) then
       innermostdecomposed=i
    end if
    if(Debug) print *,__PIO_FILE__,__LINE__,'iorank: ',iorank,'start: ',start,' kount: ',kount


  end subroutine Calcstartandcount

!
! Compute start and kount values to distribute gdim over ioprocs
! as evenly as possible.  gdim must be >= ioprocs
!
!
  subroutine computestartandcount(gdim,ioprocs,rank,start,kount)
    implicit none
    integer,intent(in) :: gdim,ioprocs,rank
    integer(kind=pio_offset),intent(out) :: start,kount
    integer :: remainder, irank

    if(gdim<ioprocs) then
       stop 'bad arguments'
    end if
!    print *,__LINE__,gdim,ioprocs,rank

    irank = mod(rank,ioprocs)

    kount = gdim/ioprocs
    start = kount*irank+1
    remainder = gdim-kount*ioprocs
    if(remainder>=ioprocs-irank) then
       kount=kount+1
       start=start+max(0,(irank+remainder-ioprocs))
    end if
    write(99,*) __LINE__,gdim,ioprocs,rank,start,kount,remainder
  end subroutine computestartandcount


end module calcdecomp

#ifdef TESTCALCDECOMP
program sandctest
  use calcdecomp  !_EXTERNAL
  implicit none
  

!  integer, parameter :: ndims=4
!  integer, parameter :: gdims(ndims) = (/66,199,10,8/)
  integer, parameter :: ndims=3
  integer, parameter :: gdims(ndims) = (/777602,2,2/)
  integer, parameter :: num_io_procs=32
!  integer :: gdims(ndims)
  integer :: psize, n, i,j,k,m
  integer, parameter :: imax=200,jmax=200,kmax=30,mmax=7
  integer(kind=pio_offset) :: start(ndims), count(ndims)
  integer :: iorank, numaiotasks, tpsize
#ifdef DOTHIS
  do i=1,imax
     gdims(1)=i
     do j=1,jmax
        gdims(2)=j
        do k=20,kmax
           gdims(3)=k
           do m=1,mmax
              gdims(4)=m
#endif
              tpsize = 0
              numaiotasks=0
              do while(tpsize==0)
                 do iorank=0,num_io_procs-1
                    call Calcstartandcount(PIO_double, ndims, gdims, num_io_procs, iorank, start, count, numaiotasks)
                    
                    psize=1
                    do n=1,ndims
                       psize=psize*count(n)
                    end do
                    tpsize = tpsize+psize
                    !                 if(ndims==3) then
!                                        write(*,'(i2,a,3i5,a,3i5,2i12)') iorank,' start =',start,' count=', count, product(gdims), psize
                    !                 else if(ndims==4) then
                    if(sum(count)>0) then
                       write(*,'(i2,a,3i8,a,3i8,2i12)') iorank,' start =',start,' count=', count, product(gdims), psize
                       if(any(start<0)) then
                          print *, gdims
                          stop 
                       endif
                    end if
                    
                 end do
                 if(tpsize/=product(gdims)) then
                    print *,'Failed to converge: ',tpsize,product(gdims),gdims
                    stop
                 end if
              end do
#ifdef DOTHIS
           end do
        end do
     end do
  end do
#endif
end program sandctest
#endif
