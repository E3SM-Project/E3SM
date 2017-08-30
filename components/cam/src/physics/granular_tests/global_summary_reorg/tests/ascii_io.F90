module ascii_io

  use shr_kind_mod,  only: r8=>shr_kind_r8, longchar=>SHR_KIND_CL
  use constituents,  only: pcnst, cnst_name 
  use ppgrid,        only: pcols, pver
  use physics_types, only: physics_state
  use cam_abortutils,only: endrun
  use cam_logfile,   only : iulog

  implicit none

CONTAINS

!---------------------------------------------------------------------------------
! Write out components of model state vector to ASCII file
!---------------------------------------------------------------------------------
subroutine read_state_ascii( ic_filepath, nstep, state )

   character(len=*),   intent(in)    :: ic_filepath
   integer,            intent(in)    :: nstep
   type(physics_state),intent(inout) :: state       ! Physics state variables

   integer :: ncol, lchnk, ii, kk , m, idummy1, idummy2, idummy3
   integer :: funit, ierr
   character(len=longchar) :: ic_filename
   character(len=longchar) :: cdummy1, cdummy2, cdummy3

   !--------
   ! chunk index will be part of the output file name

   lchnk = state%lchnk

   ! output will include all active columns in this chunk

   ncol  = state%ncol

   ! open a file

   funit = 1234
   write(ic_filename,'(2(a,i5.5),a)') trim(ic_filepath)//'/nstep_',nstep,'_chunk_',lchnk,'_state.asc'

   write(iulog,*) 'Attempting to open file '//trim(ic_filename)//' for IC conditions'
   open(unit=funit,file=trim(ic_filename),access='sequential',action='read',form='formatted',iostat=ierr)
   if (ierr/=0) then
      write(iulog,*) 'Failed to open file '//trim(ic_filename)//' for ASCII output of model state.' 
      call endrun
   end if

   ! write out components of state

   do ii = 1,ncol

      read(funit,*) cdummy1, state%lat(ii),cdummy2,state%lon(ii)

      read(funit,*) cdummy1
      do kk = 1,pver
         read(funit,*) idummy1, state%u(ii,kk)
      end do

      read(funit,*) cdummy1
      do kk = 1,pver
         read(funit,*) idummy1, state%v(ii,kk)
      end do

      read(funit,*) cdummy1
      do kk = 1,pver
         read(funit,*) idummy1, state%s(ii,kk)
      end do

      read(funit,*) cdummy1
      do kk = 1,pver
         read(funit,*) idummy1, state%pdel(ii,kk)
      end do

      ! tracers
     !do m=1,pcnst
      do m=1,3

         read(funit,*) cdummy1
         do kk = 1,pver
            read(funit,*) idummy1, state%q(ii,kk,m)
         end do

      end do

   end do  ! column loop

   close(unit=funit)

end subroutine read_state_ascii

end module ascii_io
