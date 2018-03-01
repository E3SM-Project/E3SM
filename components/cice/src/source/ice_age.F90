!=======================================================================
!
!BOP
!
! !MODULE: ice_age - Age tracer for sea ice
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_age
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size
      use ice_constants
      use ice_fileunits
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file, &
                             runtype
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
!
!EOP
!
      implicit none

      logical (kind=log_kind) :: & 
         restart_age      ! if .true., read age tracer restart file

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_age
!
! !DESCRIPTION:
!
!  Initialize ice age tracer (call prior to reading restart data)
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_age 
!
! !USES:
!
      use ice_state, only: filename_iage
!
!EOP
!

      if (trim(filename_iage) /= 'none') restart_age = .true.

      if (restart_age) then
         if (trim(runtype) == 'continue') then
            call read_restart_age
         else
            call read_restart_age(filename_iage)
         endif
      endif

      end subroutine init_age

!=======================================================================

!BOP
!
! !ROUTINE: increment_age 
!
! !DESCRIPTION:
!
!  Increase ice age tracer by timestep length.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine increment_age (nx_block, ny_block, &
                                dt,       icells,   &
                                indxi,    indxj,    &
                                iage)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         iage
!
!  local variables
!
      integer (kind=int_kind) :: i, j, ij
!
!EOP
!
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         iage(i,j) = iage(i,j) + dt 
      enddo

      end subroutine increment_age

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_age - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_age(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for restarting
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
         
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.age.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_age,filename,0)

      if (my_task == master_task) then
        write(nu_dump_age) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         call ice_write(nu_dump_age,0,trcrn(:,:,nt_iage,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_age)

      end subroutine write_restart_age

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_age - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_age(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for an ice age restart
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag

      if (my_task == master_task) then
         ! reconstruct path/file
         if (present(filename_spec)) then
            filename = filename_spec
         else
            open(nu_rst_pointer,file=pointer_file)
            read(nu_rst_pointer,'(a)') filename0
            filename = trim(filename0)
            close(nu_rst_pointer)

            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice('iage restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.age', &
               string2(1:lenstr(string2))
         endif
      endif ! master_task

      call ice_open(nu_restart_age,filename,0)

      if (my_task == master_task) then
        read(nu_restart_age) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         call ice_read(nu_restart_age,0,trcrn(:,:,nt_iage,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_restart_age)

      end subroutine read_restart_age

!=======================================================================

      end module ice_age

!=======================================================================
