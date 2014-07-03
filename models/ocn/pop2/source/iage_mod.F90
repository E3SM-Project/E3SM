!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module iage_mod

!BOP
! !MODULE: iage_mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: iage_mod.F90 26603 2011-01-28 23:09:02Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use blocks, only: nx_block, ny_block
   use domain_size, only: max_blocks_clinic, km
   use domain, only: nblocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use prognostic, only: tracer_field
   use kinds_mod
   use constants, only: c0, c1, char_blank, delim_fmt
   use io, only: data_set
   use io_types, only: stdout, nml_in, nml_filename
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
   use passive_tracer_tools, only: ind_name_pair, tracer_read, &
       rest_read_tracer_block, file_read_tracer_block
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: iage_tracer_cnt,        &
             iage_init,              &
             iage_set_interior,      &
             iage_reset

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracer
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      iage_tracer_cnt = 1

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      iage_ind = 1     ! iage index

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(iage_tracer_cnt) :: &
      ind_name_table = (/ ind_name_pair(iage_ind, 'IAGE') /)

!-----------------------------------------------------------------------
!  tavg ids for non-standard tavg variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_IAGE_RESET_TEND       ! tavg id for surface reset tendency of IAGE

!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: iage_init
! !INTERFACE:

 subroutine iage_init(init_ts_file_fmt, read_restart_filename, &
                      tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize iage tracer module. This involves setting metadata, reading
!  the module namelist and setting initial conditions.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use broadcast, only: broadcast_scalar
   use prognostic, only: curtime, oldtime
   use grid, only: KMT, n_topo_smooth, fill_points

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(iage_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,iage_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'iage_mod:iage_init'

   character(char_len) :: &
      init_iage_option,           & ! option for initialization of iage
      init_iage_init_file,        & ! filename for option 'file'
      init_iage_init_file_fmt       ! file format for option 'file'

   logical(log_kind) :: &
      lnml_found             ! Was iage_nml found ?

   integer(int_kind) :: &
      n,                   & ! index for looping over tracers
      k,                   & ! index for looping over depth levels
      iblock,              & ! index for looping over blocks
      nml_error              ! namelist i/o error flag

!     l,                   & ! index for looping over time levels

   type(tracer_read), dimension(iage_tracer_cnt) :: &
      tracer_init_ext        ! namelist variable for initializing tracers

   namelist /iage_nml/ &
      init_iage_option, init_iage_init_file, tracer_init_ext, &
      init_iage_init_file_fmt

   character (char_len) ::  &
      iage_restart_filename  ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   tracer_d_module(iage_ind)%short_name = 'IAGE'
   tracer_d_module(iage_ind)%long_name  = 'Ideal Age'
   tracer_d_module(iage_ind)%units      = 'years'
   tracer_d_module(iage_ind)%tend_units = 'years/s'
   tracer_d_module(iage_ind)%flux_units = 'cm years/s'

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_iage_option = 'unknown'
   init_iage_init_file = 'unknown'
   init_iage_init_file_fmt = 'bin'

   do n = 1,iage_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then  
         nml_error = -1
      else
         nml_error =  1      
      endif
      do while (nml_error > 0)
         read(nml_in, nml=iage_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'iage_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_iage_option , master_task)
   call broadcast_scalar(init_iage_init_file, master_task)
   call broadcast_scalar(init_iage_init_file_fmt, master_task)

   do n = 1,iage_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   select case (init_iage_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d Ideal Age set to all zeros' 
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
       
   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      iage_restart_filename = char_blank

      if (init_iage_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read iage from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         iage_restart_filename = read_restart_filename
         init_iage_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         iage_restart_filename = trim(init_iage_init_file)

      endif

      call rest_read_tracer_block(init_iage_init_file_fmt, &
                                  iage_restart_filename,   &
                                  tracer_d_module,         &
                                  TRACER_MODULE)

   case ('file')
      call document(subname, 'iage being read from separate file')

      call file_read_tracer_block(init_iage_init_file_fmt, &
                                  init_iage_init_file,     &
                                  tracer_d_module,         &
                                  ind_name_table,          &
                                  tracer_init_ext,         &
                                  TRACER_MODULE)
 
      if (n_topo_smooth > 0) then
         do k=1,km
            call fill_points(k,TRACER_MODULE(:,:,k,1,curtime,:), errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'iage_init: error in fill_points for tracer 1')
               return
            endif
         enddo
      endif

   case default
      call document(subname, 'unknown init_iage_option = ', init_iage_option)
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
      do n = 1,iage_tracer_cnt
         do k = 1,km
            where (k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do
   enddo

!-----------------------------------------------------------------------

   call define_tavg_field(tavg_IAGE_RESET_TEND, 'IAGE_RESET_TEND',2,  &
                          long_name='surface reset tendency of IAGE', &
                          units='years/s', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!EOC

 end subroutine iage_init

!***********************************************************************
!BOP
! !IROUTINE: iage_set_interior
! !INTERFACE:

 subroutine iage_set_interior(k, DTRACER_MODULE)

! !DESCRIPTION:
!  set interior source/sink term for ideal age tracer
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use time_management, only: seconds_in_year

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: &
      k                   ! vertical level index

! !OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,iage_tracer_cnt), intent(out) :: &
      DTRACER_MODULE      ! computed source/sink term

!EOP
!BOC

!-----------------------------------------------------------------------

    if (k > 1) then
       DTRACER_MODULE = c1 / seconds_in_year
    else
       DTRACER_MODULE = c0
    endif

!-----------------------------------------------------------------------
!EOC

 end subroutine iage_set_interior

!***********************************************************************
!BOP
! !IROUTINE: iage_reset
! !INTERFACE:

 subroutine iage_reset(TRACER_MODULE, bid)

! !DESCRIPTION:
!  reset surface value for ideal age tracer
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use time_management, only: mix_pass, c2dtt
   use prognostic, only: PSURF, newtime
   use constants, only: grav
   use grid, only: dz

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: bid

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,km,iage_tracer_cnt), intent(inout) :: &
      TRACER_MODULE      ! ideal age tracer

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block) :: WORK

!-----------------------------------------------------------------------

   if (mix_pass /= 1) then
      if (accumulate_tavg_now(tavg_IAGE_RESET_TEND)) then
            WORK = -TRACER_MODULE(:,:,1,1) / c2dtt(1)
            WORK = WORK * (c1 + PSURF(:,:,newtime,bid)/grav/dz(1))
            call accumulate_tavg_field(WORK,tavg_IAGE_RESET_TEND,bid,1)
      endif
   endif

   TRACER_MODULE(:,:,1,:) = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine iage_reset

!***********************************************************************

end module iage_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
