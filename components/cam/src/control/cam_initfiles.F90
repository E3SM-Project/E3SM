module cam_initfiles
!----------------------------------------------------------------------- 
! 
! Open, close, and provide access to the initial conditions and topography files.
! 
!-----------------------------------------------------------------------

!### MBRANSON: q3d additions
use pio,              only: file_desc_t, pio_offset_kind, pio_global, &
                            pio_inq_att, pio_get_att, pio_nowrite,    &
                            pio_closefile
use spmd_utils,       only: masterproc

!use cam_control_mod,  only: initial_run, restart_run, branch_run, caseid, brnch_retain_casename
use cam_control_mod,  only: nsrest
use filenames,        only: caseid, brnch_retain_casename
use ioFileMod,        only: getfil, opnfil
use cam_pio_utils,    only: cam_pio_openfile
use shr_kind_mod,     only: r8=>shr_kind_r8, cl=>shr_kind_cl
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

use time_manager,     only: dtime
!!use rayleigh_friction, only: raytau0

implicit none
private
save

! Public methods

public :: &
!### MBRANSON: q3d additions
   cam_initfiles_readnl,      &! read namelist
   cam_initfiles_open,   &! open initial and topo files
   initial_file_get_id,  &! returns filehandle for initial file
   topo_file_get_id,     &! returns filehandle for topo file
   cam_initfiles_close     ! close initial and topo files

type(file_desc_t), pointer :: fh_ini, fh_topo

!### MBRANSON: q3d additions
! Namelist inputs
!### mdb hack since E3SM is currently not accepting use_topo_file in user_nl_cam
!!!logical :: use_topo_file = .true.
logical :: use_topo_file = .false.
character(len=cl), public, protected :: ncdata = 'ncdata'     ! full pathname for initial dataset
character(len=cl), public, protected :: bnd_topo = 'bnd_topo' ! full pathname for topography dataset

real(r8), public, protected :: pertlim = 0.0_r8 ! maximum abs value of scale factor used to perturb
character(len=cl) :: cam_branch_file = ' '      ! Filepath of primary restart file for a branch run

! The restart pointer file contains name of most recently written primary restart file.
! The contents of this file are updated by cam_write_restart as new restart files are written.
character(len=cl), public, protected :: rest_pfile
                                                   
! Filename for initial restart file.
character(len=cl) :: restart_file = ' '

! case name read from initial restart file.  This case name matches the caseid
! which is embedded in the filename.
character(len=cl) :: caseid_prev = ' '

type(file_desc_t), target  :: fh_restart

!======================================================================= 
contains
!======================================================================= 
subroutine cam_initfiles_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpir8=>mpi_real8, &
                              mpichar=>mpi_character, mpi_logical
   use cam_instance,    only: inst_suffix

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr

   character(len=cl)        :: locfn
   logical                  :: filefound
   integer                  :: xtype
   integer(pio_offset_kind) :: slen

   character(len=*), parameter :: sub = 'cam_initfiles_readnl'
   real(r8) :: raytau0    ! this won't be used but put it here for now
   !!namelist /cam_initfiles_nl/ ncdata, use_topo_file, bnd_topo, pertlim, &
   namelist /cam_inparm/ ncdata, use_topo_file, bnd_topo, pertlim, &
                         cam_branch_file, dtime, raytau0
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      !write(*,*) '### cam_initfiles_readnl: unitn = ',unitn,' nlfile = ',trim(nlfile)
      open( unitn, file=trim(nlfile), status='old' )
      !!!call find_group_name(unitn, 'cam_initfiles_nl', status=ierr)
      call find_group_name(unitn, 'cam_inparm', status=ierr)
      if (ierr == 0) then
         read(unitn, cam_inparm, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub // ': ERROR: reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(ncdata, len(ncdata), mpichar, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": ERROR: mpi_bcast: ncdata")
   call mpi_bcast(use_topo_file, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": ERROR: mpi_bcast: use_topo_file")
   call mpi_bcast(bnd_topo, len(bnd_topo), mpichar, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": ERROR: mpi_bcast: bnd_topo")
   call mpi_bcast(pertlim, 1, mpir8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": ERROR: mpi_bcast: pertlim")
   call mpi_bcast(cam_branch_file, len(cam_branch_file), mpichar, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": ERROR: mpi_bcast: cam_branch_file")

   ! Set pointer file name based on instance suffix
   rest_pfile = './rpointer.atm' // trim(inst_suffix)

   ! Set name of primary restart file
   !if (restart_run) then
   if (nsrest==1) then
      ! Read name of restart file from pointer file
      if (masterproc) then
         unitn = getunit()
         call opnfil(rest_pfile, unitn, 'f', status="old")
         read (unitn, '(a)', iostat=ierr) restart_file
         if (ierr /= 0) then
            call endrun(sub // ': ERROR: reading rpointer file')
         end if
         close(unitn)
         call freeunit(unitn)
      end if

      call mpi_bcast(restart_file, len(restart_file), mpichar, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": ERROR: mpi_bcast: restart_file")
      
   !else if (branch_run) then
   else if (nsrest==3) then
   
      ! use namelist input
      restart_file = trim(cam_branch_file)
   end if

   ! Get caseid from restart or branch file.
   !if (restart_run .or. branch_run) then
   if (nsrest .ge. 1) then

      call getfil(restart_file, locfn)
      inquire(file=trim(locfn), exist=filefound)
      if (.not.filefound) then
         call endrun(sub//': ERROR: could not find restart file '//trim(locfn))
      end if

      call cam_pio_openfile(fh_restart, trim(locfn), pio_nowrite)

      ierr = pio_inq_att(fh_restart, pio_global, 'caseid', xtype, slen)
      ierr = pio_get_att(fh_restart, pio_global, 'caseid', caseid_prev)
      caseid_prev(slen+1:len(caseid_prev)) = ' '

      !if (branch_run .and. caseid_prev==caseid .and. .not.brnch_retain_casename) then
      if (nsrest==3 .and. caseid_prev==caseid .and. .not.brnch_retain_casename) then
         write(iulog,*) sub//': Must change case name on branch run'
         write(iulog,*) 'Prev case = ',caseid_prev,' current case = ',caseid
         call endrun(sub//': ERROR: Must change case name on branch run')
      end if
   end if

   if (masterproc) then
      write(iulog,*) sub//' options:'

      !if (initial_run) then
      if (nsrest==0) then
      
         write(iulog,*)'  Initial run will start from: ', trim(ncdata)

         if (use_topo_file) then
            write(iulog,*) '  Topography dataset is: ', trim(bnd_topo)
         else
            write(iulog,*) '  Topography dataset not used: PHIS, SGH, SGH30, LANDM_COSLAT set to zero'
         end if

      !else if (restart_run) then
      else if (nsrest==1) then
         write(iulog,*)'  Continuation of case:             ', trim(caseid_prev)
         write(iulog,*)'  Restart run will start from file: ', trim(restart_file)
      !else if (branch_run) then
      else if (nsrest==3) then
         write(iulog,*)'  Continuation of case:             ', trim(caseid_prev)
         write(iulog,*)'  Branch run will start from file:  ', trim(restart_file)
      end if

      write(iulog,*) &
         '  Maximum abs value of scale factor used to perturb initial conditions, pertlim= ', pertlim

#ifdef PERGRO
      write(iulog,*)'  The PERGRO CPP token is defined.'
#endif

   end if

end subroutine cam_initfiles_readnl

!======================================================================= 

function initial_file_get_id()
  type(file_desc_t), pointer :: initial_file_get_id
  initial_file_get_id => fh_ini
end function initial_file_get_id

function topo_file_get_id()
  type(file_desc_t), pointer :: topo_file_get_id
  topo_file_get_id => fh_topo
end function topo_file_get_id

!======================================================================= 

subroutine cam_initfiles_open()

   ! Open the initial conditions and topography files.

   !!!use filenames,        only: ncdata, bnd_topo
   use filenames,        only: bnd_topo
   use ioFileMod,        only: getfil

   use cam_pio_utils,    only: cam_pio_openfile
   use pio,              only: pio_nowrite

   use readinitial,      only: read_initial

   character(len=256) :: ncdata_loc     ! filepath of initial file on local disk
   character(len=256) :: bnd_topo_loc   ! filepath of topo file on local disk

   integer :: allocate_status
   !----------------------------------------------------------------------- 
   
   ! Open initial, topography, and landfrac datasets
   call getfil (ncdata, ncdata_loc)

   !allocate(fh_ini)
   allocate(fh_ini,stat=allocate_status)
   if (allocate_status /= 0) write(*,*) '### allocate fh_ini FAILED ###'

   call cam_pio_openfile(fh_ini, ncdata_loc, PIO_NOWRITE)
   ! Backward compatibility: look for topography data on initial file if topo file name not provided.
   if (trim(bnd_topo) /= 'bnd_topo' .and. len_trim(bnd_topo) > 0) then
      allocate(fh_topo)
      call getfil(bnd_topo, bnd_topo_loc)
      call cam_pio_openfile(fh_topo, bnd_topo_loc, PIO_NOWRITE)
   else
      fh_topo => fh_ini
   end if

   ! Check for consistent settings on initial dataset -- this is dycore
   ! dependent -- should move to dycore interface
   call read_initial (fh_ini)

   if (use_topo_file) then
      if (trim(bnd_topo) /= 'bnd_topo' .and. len_trim(bnd_topo) > 0) then
         allocate(fh_topo)
         call getfil(bnd_topo, bnd_topo_loc)
         call cam_pio_openfile(fh_topo, bnd_topo_loc, pio_nowrite)
      else
         ! Allow topography data to be read from the initial file if topo file name
         ! is not provided.
         fh_topo => fh_ini
      end if
   else
      nullify(fh_topo)
   end if

end subroutine cam_initfiles_open

!======================================================================= 

subroutine cam_initfiles_close()

  use pio,          only: pio_closefile

  if(associated(fh_ini)) then
     if(.not. associated(fh_ini, target=fh_topo)) then
        call pio_closefile(fh_topo)
        deallocate(fh_topo)
     end if
     
     call pio_closefile(fh_ini)
     deallocate(fh_ini)
     nullify(fh_ini)
     nullify(fh_topo)
  end if
end subroutine cam_initfiles_close

!======================================================================= 

end module cam_initfiles
