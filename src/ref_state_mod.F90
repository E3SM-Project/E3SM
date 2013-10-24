#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module ref_state_mod
#ifdef _REFSOLN
  ! ------------------
  use kinds, only : real_kind
  ! ------------------
  use parallel_mod, only : parallel_t, iam, syncmp, abortmp
  ! ------------------
  use schedtype_mod, only : Schedule
  ! ------------------
!  use dimensions_mod
  ! ------------------
implicit none
private

  interface ref_state_write
      module procedure ref_state_write_2d
      module procedure ref_state_write_jr
  end interface

  interface ref_state_read
      module procedure ref_state_read_2d
      module procedure ref_state_read_jr
  end interface

  public :: ref_state_write
  public :: ref_state_read

contains

  subroutine ref_state_write_2d(phi,v,fstub,timetag,nets,nete)
     use dimensions_mod, only : np
     implicit none
     integer              , intent(in) :: nets,nete
     real (kind=real_kind), intent(in) :: phi(np,np,nets:nete)
     real (kind=real_kind), intent(in) :: v(np,np,2,nets:nete)
     integer              , intent(in) :: timetag         ! time tag (day)
     character(len=*)     , intent(in) :: fstub           ! file stub

     ! =========================
     ! Local variables...
     ! =========================

     integer iunit
     integer ie
     integer reclen

     character(len=6) :: chartag
     character(len=80):: fname

     reclen = np*np + 2*np*np 
     iunit = 66                 ! hardwire hack the unit number (should use navu)
     
     write(chartag,'(i6)') timetag
     fname=TRIM(ADJUSTL(fstub))//"."//TRIM(ADJUSTL(chartag))
     
      open(unit=iunit,file=fname,access="DIRECT",recl=reclen*8,form="UNFORMATTED",status="UNKNOWN")

      do ie=nets,nete
        write(iunit,rec=ie)phi(:,:,ie),v(:,:,:,ie)
        write(6,*) 'writing elem=',ie
      end do
     
      close(iunit)

  end subroutine ref_state_write_2d

  subroutine ref_state_write_jr(phi,v,fstub,timetag,nets,nete,par)
     use dimensions_mod, only : np
     implicit none
     integer              , intent(in) :: nets,nete
     real (kind=real_kind), intent(in) :: phi(np,np,nets:nete)
     real (kind=real_kind), intent(in) :: v(np,np,2,nets:nete)
     integer              , intent(in) :: timetag         ! time tag (day)
     character(len=*)     , intent(in) :: fstub           ! file stub
     type(parallel_t)     , intent(in) :: par             ! communicator

     ! =========================
     ! Local variables...
     ! =========================

     integer handle, ret
     integer ie, ig
     integer reclen
     integer iunit

     character(len=6) :: chartag
     character(len=20):: fname

     reclen = np*np
!     reclen = np*np + 2*np*np 
     iunit = 66                 ! hardwire hack the unit number (should use navu)
     
      write(chartag,'(i6)') timetag
      fname=TRIM(ADJUSTL(fstub))//"."//TRIM(ADJUSTL(chartag))

     call jropen_direct(TRIM(fname),'replace',reclen*8,handle,ret)
     if (ret/= 0) then
      write(6,*) 'file,handle=',fname,handle
      call abortmp('error opening binary reference file')
     end if
     if (par%masterproc) write(6,*) 'opening binary reference file=',TRIM(fname),handle

     call syncmp(par)
 
     do ie=nets,nete
      ig = Schedule(1)%Local2Global(ie)
      call jrwrite_direct(handle, ig, phi(1,1,ie), 0, ret)
      if (ret /= 0) then
        write(6,*) 'file,handle:',fname,handle
        call abortmp('failure to write reference file')
      end if
     end do
     
     call syncmp(par)

     call jrclose_direct(handle)

  end subroutine ref_state_write_jr

  subroutine ref_state_read_2d(phi,v,fstub,timetag,nets,nete)
     use dimensions_mod, only : np
     implicit none
     integer              , intent(in)    :: nets,nete
     real (kind=real_kind), intent(out)   :: phi(np,np,nets:nete)
     real (kind=real_kind), intent(out)   :: v(np,np,2,nets:nete)
     integer              , intent(in)    :: timetag         ! time tag (day)
     character(len=*)     , intent(in)    :: fstub           ! file stub

     ! =========================
     ! Local variables...
     ! =========================

     integer iunit
     integer ie
     integer reclen

     character(len=6) :: chartag
     character(len=20):: fname
    
     reclen = np*np + 2*np*np 
     iunit = 66                 ! hardwire hack the unit number (should use navu)
     write(chartag,'(i6)') timetag
     fname=TRIM(ADJUSTL(fstub))//"."//TRIM(ADJUSTL(chartag))
          
!$OMP CRITICAL (IOCRIT)
     print *,"opening file ",fname

     open(unit=iunit,         &
          file=fname,         & 
          access="DIRECT",    &
          recl=reclen*8,      &
          form="UNFORMATTED", &
          STATUS='OLD')

     print *,"reading elements:",LBOUND(phi,3),UBOUND(phi,3)
     do ie=nets,nete
        read(iunit,rec=ie)phi(:,:,ie),v(:,:,:,ie)
#if 0
        print *,ie,"read v min,max =",MAXVAL(phi(:,:,ie))
#endif 
     end do
     close(iunit)
!$OMP END CRITICAL (IOCRIT)

  end subroutine ref_state_read_2d

  subroutine ref_state_read_jr(phi,v,fstub,timetag,nets,nete,par)
     use dimensions_mod, only : np
     implicit none
     integer              , intent(in)    :: nets,nete
     real (kind=real_kind), intent(out)   :: phi(np,np,nets:nete)
     real (kind=real_kind), intent(out)   :: v(np,np,2,nets:nete)
     integer              , intent(in)    :: timetag         ! time tag (day)
     character(len=*)     , intent(in)    :: fstub           ! file stub
     type(parallel_t)     , intent(in)    :: par           ! file stub

     ! =========================
     ! Local variables...
     ! =========================

     integer handle, ret
     integer ie,ig
     integer reclen
     integer iunit

     character(len=6) :: chartag
     character(len=80):: fname
    
     reclen = np*np
!     reclen = np*np + 2*np*np 
     iunit = 66                 ! hardwire hack the unit number (should use navu)
     write(chartag,'(i6)') timetag
     fname=TRIM(ADJUSTL(fstub))//"."//TRIM(ADJUSTL(chartag))
          
     if (par%masterproc) print *,"opening file ",fname

     call jropen_direct(TRIM(fname),'old',reclen*8,handle,ret)
     if (ret/= 0) then
      write(6,*) 'filename,handle=',fname,handle
      call abortmp('error opening existing binary reference file=')
     end if
     if (par%masterproc) write(6,*) 'opening existing reference file=',TRIM(fname),handle
 
     do ie=nets,nete
      ig = Schedule(1)%Local2Global(ie)
      call jrread_direct(handle, ig, phi(1,1,ie), 0, ret)
      if (ret /= 0) then
       write(6,*) 'fname,handle:',fname,handle
       call abortmp('failure to read existing binary reference file=')
      end if
     end do
     
     call syncmp(par)
 
     call jrclose_direct(handle)

  end subroutine ref_state_read_jr
#endif
end module ref_state_mod

