!=======================================================================
!
!BOP
!
! !MODULE: ice_lvl - Ridged ice tracers for sea ice
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
      module ice_lvl
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size
      use ice_constants
      use ice_fileunits
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, &
                             pointer_file, runtype
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
!
!EOP
!
      implicit none

      logical (kind=log_kind) :: & 
         restart_lvl      ! if .true., read lvl tracer restart file

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_lvl
!
! !DESCRIPTION:
!
!  Initialize ice lvl tracers (call prior to reading restart data)
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_lvl 
!
! !USES:
!
      use ice_state, only: filename_lvl
!
!EOP
!

      if (trim(filename_lvl) /= 'none') restart_lvl = .true.

      if (restart_lvl) then
         if (trim(runtype) == 'continue') then
            call read_restart_lvl
         else
            call read_restart_lvl(filename_lvl)
         endif
      endif

      end subroutine init_lvl

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_lvl - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_lvl(filename_spec)
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
      use ice_state, only: nt_alvl, nt_vlvl, trcrn
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
              restart_file(1:lenstr(restart_file)),'.lvl.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_lvl,filename,0)

      if (my_task == master_task) then
        write(nu_dump_lvl) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         call ice_write(nu_dump_lvl,0,trcrn(:,:,nt_alvl,n,:),'ruf8',diag)
         call ice_write(nu_dump_lvl,0,trcrn(:,:,nt_vlvl,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_lvl)

      end subroutine write_restart_lvl

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_lvl - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_lvl(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for an ice lvl restart
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
      use ice_state, only: nt_alvl, nt_vlvl, trcrn
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

            ! reconstruct path/file
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice('ilvl restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.lvl', &
               string2(1:lenstr(string2))
         endif
      endif ! master_task

      call ice_open(nu_restart_lvl,filename,0)

      if (my_task == master_task) then
        read(nu_restart_lvl) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         call ice_read(nu_restart_lvl,0,trcrn(:,:,nt_alvl,n,:),'ruf8',diag)
         call ice_read(nu_restart_lvl,0,trcrn(:,:,nt_vlvl,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_restart_lvl)

      end subroutine read_restart_lvl

!=======================================================================

      end module ice_lvl

!=======================================================================
