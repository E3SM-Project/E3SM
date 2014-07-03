module cam_history
!-------------------------------------------------------------------------------------------
! 
! The cam_history module provides the user interface for CAM's history output capabilities.
! It maintains the lists of fields that are written to each history file, and the associated
! metadata for those fields such as descriptive names, physical units, time axis properties,
! etc.  It also contains the programmer interface which provides routines that are called from
! the physics and dynamics initialization routines to define the fields that are produced by
! the model and are available for output, and the routine that is called from the corresponding
! run method to add the field values into a history buffer so that they may be output to disk.
!
! There are two special history files.  The initial file and the satellite track file.
!
! Public functions/subroutines:
!   addfld, add_default
!   intht
!   write_restart_history
!   read_restart_history
!   outfld
!   wshist
!   initialize_iop_history
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use spmd_utils,   only: masterproc, npes
   use ppgrid,       only: pcols
   use cam_instance, only: inst_suffix
   use filenames,    only: interpret_filename_spec, ncdata, bnd_topo
   use abortutils,   only: endrun
   use pmgrid,       only : dyndecomp_set

   use units,        only: getunit

   use hycoef,       only: hyai, hybi, hyam, hybm, ps0

   use dyn_grid,     only: get_horiz_grid_dim_d, get_horiz_grid_d, get_dyn_grid_parm

   use pio,          only: file_desc_t, var_desc_t, io_desc_t, &
                           iotype_pnetcdf, iotype_netcdf, &
                           pio_noerr, pio_bcast_error, pio_internal_error, &
                           pio_seterrorhandling, pio_setdebuglevel, pio_setframe, &
                           pio_rearr_box, pio_rearr_none,  &
                           pio_nofill, pio_clobber,pio_offset, &
                           pio_int, pio_real, pio_double, pio_char, pio_offset, pio_unlimited, pio_global, &
                           pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, pio_inq_varname, pio_inq_varndims, &
                           pio_def_dim, pio_def_var, pio_enddef, &
                           pio_put_att, pio_put_var, pio_get_att, pio_get_var, &
                           pio_write, pio_write_darray, &
                           pio_read_darray, pio_closefile, pio_freedecomp, pio_syncfile

   use scamMod,       only: scmlon,single_column
   use perf_mod,      only: t_startf, t_stopf
   use cam_logfile,   only: iulog	
   use cam_pio_utils, only: phys_decomp, dyn_decomp, dyn_stagger_decomp
   use cam_history_buffers, only : dim_index_3d,hbuf_accum_inst, hbuf_accum_add, hbuf_accum_min, hbuf_accum_max, &
        hbuf_accum_add00z, hbuf_accum_addlcltime
   use cam_history_support, only : max_fieldname_len, fieldname_suffix_len, max_chars, ptapes, fieldname_len, max_string_len, &
        date2yyyymmdd, registeredmdims, fieldname_lenp2, sec2hms, pflds, column_info, field_info, active_entry, hentry, &
        hist_mdims, register_hist_mdim, lookup_hist_mdim_indicies
   use sat_hist,    only: is_satfile
   use phys_control,only: cam_chempkg_is

   implicit none
   private
   save

   public :: phys_decomp, dyn_decomp, dyn_stagger_decomp, fieldname_len ! forward from cam_pio_utils


!
! master_entry: elements of an entry in the master field list
!
   type master_entry
      type (field_info)     :: field            ! field information
      character(len=1)           :: avgflag(ptapes)  ! averaging flag
      character(len=max_chars) :: time_op(ptapes)  ! time operator (e.g. max, min, avg)
      logical               :: act_sometape     ! Field is active on some tape
      logical               :: actflag(ptapes)  ! Per tape active/inactive flag
      integer               :: htapeindx(ptapes)! This field's index on particular history tape
      type(master_entry), pointer :: next_entry => null() ! The next master entry
   end type master_entry

   type (master_entry), pointer :: masterlinkedlist => null()   ! master field linkedlist top

   type master_list
      type(master_entry), pointer :: thisentry => null()
   end type master_list

   type (master_list), pointer :: masterlist(:) => null() ! master field array for hash lookup of field

!  column_info structure has moved to cam_pio_utils



!


   type (active_entry), pointer :: tape(:) => null()         ! history tapes
   type (active_entry), target,allocatable :: history_tape(:)          ! history tapes
   type (active_entry), target, allocatable :: restarthistory_tape(:)          ! restart history tapes


   type rvar_id
      type(var_desc_t), pointer :: vdesc => null()
      integer :: type
      integer :: ndims
      integer :: dims(4)
      character(len=fieldname_lenp2) :: name
   end type rvar_id
   type rdim_id
      integer :: len
      integer :: dimid
      character(len=fieldname_lenp2) :: name      
   end type rdim_id
!
!   The size of these parameters should match the assignments in restart_vars_setnames and restart_dims_setnames below
!   
   integer, parameter :: restartvarcnt=29
   integer, parameter :: restartdimcnt=8
   type(rvar_id) :: restartvars(restartvarcnt)
   type(rdim_id) :: restartdims(restartdimcnt)



   integer :: nfmaster = 0             ! number of fields in master field list
   integer :: nflds(ptapes)            ! number of fields per tape

! per tape sampling frequency (0=monthly avg)

   integer :: i                        ! index for nhtfrq initialization
   integer :: nhtfrq(ptapes) = (/0, (-24, i=2,ptapes)/)  ! history write frequency (0 = monthly)
   integer :: mfilt(ptapes) = 30       ! number of time samples per tape
   integer :: nfils(ptapes)            ! Array of no. of files on current h-file
   integer :: ngroup(ptapes)           ! Array of no. of contiguous columns on current h-file
   integer :: nexcl(ptapes)            ! Actual number of excluded fields
   integer :: nincl(ptapes)            ! Actual number of included primary file fields
   integer :: ndens(ptapes) = 2        ! packing density (double (1) or real (2))
   integer :: ncprec(ptapes) = -999    ! netcdf packing parameter based on ndens
   real(r8) :: beg_time(ptapes)        ! time at beginning of an averaging interval
!
! Netcdf ids
!

   type(var_desc_t) :: gwid             ! var id for gaussian weights


   integer :: nscurf(ptapes)           ! First "current" second of day for each h-file
   integer :: ncsecf(ptapes)           ! First "current" second of date for each h-file

   logical :: rgnht(ptapes) = .false.  ! flag array indicating regeneration volumes
   logical :: hstwr(ptapes) = .false.  ! Flag for history writes
   logical :: empty_htapes  = .false.  ! Namelist flag indicates no default history fields
   logical :: htapes_defined = .false. ! flag indicates history contents have been defined

   integer :: luhrest = -1             ! history restart file unit number
   character(len=max_string_len) :: hrestpath(ptapes) = (/(' ',i=1,ptapes)/) ! Full history restart pathnames
   character(len=max_string_len) :: nfpath(ptapes) = (/(' ',i=1,ptapes)/) ! Array of first pathnames, for header
   character(len=max_string_len) :: cpath(ptapes)                   ! Array of current pathnames
   character(len=max_string_len) :: nhfil(ptapes)                   ! Array of current file names
   character(len=1)  :: avgflag_pertape(ptapes) = (/(' ',i=1,ptapes)/) ! per tape averaging flag
   character(len=16)  :: logname             ! user name
   character(len=16) :: host                ! host name
   character(len=max_string_len) :: ctitle = ' '      ! Case title
   character(len=8)  :: inithist = 'YEARLY' ! If set to '6-HOURLY, 'DAILY', 'MONTHLY' or
                                            ! 'YEARLY' then write IC file 
   character(len=fieldname_lenp2) :: fincl(pflds,ptapes) ! List of fields to add to primary h-file
   character(len=max_chars)       :: fincllonlat(pflds,ptapes) ! List of fields to add to primary h-file
   character(len=fieldname_lenp2) :: fexcl(pflds,ptapes) ! List of fields to rm from primary h-file
   character(len=fieldname_lenp2) :: fwrtpr(pflds,ptapes) ! List of fields to change default history output prec
   character(len=fieldname_suffix_len ) :: fieldname_suffix = '&IC' ! Suffix appended to field names for IC file

   character(len=22) ,parameter :: LT_DESC = 'mean (over local time)' ! local time description
   logical :: collect_column_output(ptapes)

   integer :: maxvarmdims=1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Hashing.
!
!  Accelerate outfld processing by using a hash function of the field name
!  to index masterlist and determine whehter the particular field is to
!  be written to any history tape.
!
!
!  Note: the outfld hashing logic will fail if any of the following are true:
!
!         1) The lower bound on the dimension of 'masterlist' is less than 1.
!
!         2) 'outfld' is called with field names that are not defined on
!            masterlist.  This applies to both initial/branch and restart
!            runs.
!
!         3) An inconsistency between a field's tape active flag
!            'masterlist(ff)%actflag(t)' and active fields read from
!            restart files.
!
!         4) Invoking function 'gen_hash_key' before the primary and secondary
!            hash tables have been created (routine bld_outfld_hash_tbls).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  User definable constants for hash and overflow tables.
!  Define size of primary hash table (specified as 2**size).
!
   integer, parameter :: tbl_hash_pri_sz_lg2 = 16
!
!  Define size of overflow hash table % of primary hash table.
!
   integer, parameter :: tbl_hash_oflow_percent = 20
!
!  Do *not* modify the parameters below.
!
   integer, parameter :: tbl_hash_pri_sz = 2**tbl_hash_pri_sz_lg2
   integer, parameter :: tbl_hash_oflow_sz = tbl_hash_pri_sz * (tbl_hash_oflow_percent/100.0_r8) 
!
!  The primary and overflow tables are organized to mimimize space (read:
!  try to maximimze cache line usage).
!
!  gen_hash_key(fieldname) will return an index on the interval
!  [0 ... tbl_hash_pri_sz-1].
!
!
!  Primary:
!  gen_hash_key(fieldname)-------+     +----------+
!                                |     |   -ii    | 1 ------>tbl_hash_oflow(ii)
!                                |     +----------+
!                                +-->  |    ff    | 2 ------>masterlist(ff)
!                                      +----------+
!                                      |          | ...
!                                      +----------+
!                                      |          | tbl_hash_pri_sz
!                                      +----------+
!
!  Overflow (if tbl_hash_pri() < 0):
!  tbl_hash_pri(gen_hash_key(fieldname))
!                         |
!                         |            +----------+
!                         |            |     1    | 1  (one entry on O.F. chain)
!                         |            +----------+
!                         |            |    ff_m  | 2
!                         |            +----------+
!                         +--------->  |     3    | 3  (three entries on chain)
!                                      +----------+
!                                      |    ff_x  | 4
!                                      +----------+
!                                      |    ff_y  | 5
!                                      +----------+
!                                      |    ff_z  | 6
!                                      +----------+
!                                      |          | ...
!                                      +----------+
!                                      |          | tbl_hash_oflow_sz
!                                      +----------+
!
!
   integer, dimension(0:tbl_hash_pri_sz-1) :: tbl_hash_pri ! Primary hash table
   integer, dimension(tbl_hash_oflow_sz) :: tbl_hash_oflow ! Overflow hash table
!
!  Constants used in hashing function gen_hash_key.
!  Note: if the constants in table 'tbl_gen_hash_key' below are modified,
!        changes are required to routine 'gen_hash_key' because of specific
!        logic in the routine that optimizes character strings of length 8.
!

   integer, parameter :: gen_hash_key_offset = z'000053db'

   integer, parameter :: tbl_max_idx = 15  ! 2**N - 1
   integer, dimension(0:tbl_max_idx) :: tbl_gen_hash_key = &
   (/61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1/)

   ! Used in output interpolation of unstructured grids
   logical :: interpolate, ncoloutput


!
!
! Generic procedures
!
!
! Public entities
!

!
! Filename specifiers for history, initial files and restart history files
! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = tape number)
!
   character(len=max_string_len) :: rhfilename_spec = '%c.cam.rh%t.%y-%m-%d-%s.nc' ! history restart
   character(len=max_string_len), public :: hfilename_spec(ptapes) = (/ (' ', i=1, ptapes) /) ! filename specifyer

! needed by history_default
   public :: init_masterlinkedlist

! Needed by initext

   public :: nhtfrq, mfilt, inithist, ctitle

! Needed by read_namelist

   public :: fincl, fincllonlat, fexcl, fwrtpr
   public :: pflds, ptapes, empty_htapes, ndens
   public :: avgflag_pertape, collect_column_output
   integer, public :: lcltod_start(ptapes) ! start time of day for local time averaging (sec)
   integer, public :: lcltod_stop(ptapes)  ! stop time of day for local time averaging, stop > start is wrap around (sec)

! Needed by stepon

   public :: hstwr
   public :: nfils

! Functions
   public :: init_restart_history      ! Write restart history data
   public :: write_restart_history     ! Write restart history data
   public :: read_restart_history      ! Read restart history data
   public :: wshist                    ! Write files out
   public :: outfld                    ! Output a field
   public :: intht                     ! Initialization
   public :: wrapup                    ! process history files at end of run
   public :: write_inithist            ! logical flag to allow dump of IC history buffer to IC file
   public :: addfld                    ! Add a field to history file
   public :: add_default               ! Add the default fields
   public :: get_hfilepath             ! Return history filename
   public :: get_ptapes                ! Return the number of tapes being used
   public :: get_hist_restart_filepath ! Return the full filepath to the history restart file
   public :: hist_fld_active           ! Determine if a field is active on any history file
   public :: hist_fld_col_active       ! Determine if a field is active on any history file at
                                       ! each column in a chunk


CONTAINS

  subroutine init_masterlinkedlist()
    use dyn_grid, only : get_dyn_grid_parm
    integer :: plev, plevp

    nullify(masterlinkedlist)
    nullify(tape)
    plev = get_dyn_grid_parm('plev')
    plevp = get_dyn_grid_parm('plevp')

    call register_hist_mdim('lev',plev)
    call register_hist_mdim('ilev',plevp)

  end subroutine init_masterlinkedlist


   subroutine intht ()
!
!----------------------------------------------------------------------- 
! 
! Purpose: Initialize history file handler for initial or continuation run.
!          For example, on an initial run, this routine initializes "ptapes"
!          history files.  On a restart or regeneration  run, this routine 
!          only initializes history files declared beyond what existed on the 
!          previous run.  Files which already existed on the previous run have 
!          already been initialized (i.e. named and opened) in routine RESTRT.
! 
! Method: Loop over tapes and fields per tape setting appropriate variables and
!         calling appropriate routines
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use ioFileMod
      use shr_sys_mod, only: shr_sys_getenv
      use time_manager, only: get_prev_time, get_curr_time
      use cam_control_mod, only : nsrest
      use dycore,           only: dycore_is
      use sat_hist,        only: sat_hist_init
#if (defined SPMD)
      use spmd_utils, only : mpichar, mpicom
#endif
!
!-----------------------------------------------------------------------
!
! Local workspace
!
      integer :: t, f              ! tape, field indices
      integer :: begdim1           ! on-node dim1 start index
      integer :: enddim1           ! on-node dim1 end index
      integer :: begdim2           ! on-node dim2 start index
      integer :: enddim2           ! on-node dim2 end index
      integer :: begdim3           ! on-node chunk or lat start index
      integer :: enddim3           ! on-node chunk or lat end index
      type (dim_index_3d) :: dimind  ! 3-D dimension index
      integer :: day, sec          ! day and seconds from base date
      integer :: i, ii, ff         ! index
      integer :: rcode             ! shr_sys_getenv return code
      type(master_entry), pointer :: listentry
      character(len=32) :: fldname ! temp variable used to produce a left justified field name
                                   ! in the formatted logfile output

!
! Print master field list
!

      if (masterproc) then
         write(iulog,*) ' '
         write(iulog,*)' ******* MASTER FIELD LIST *******'
      end if
      listentry=>masterlinkedlist
      f=0
      do while(associated(listentry))
         f=f+1
         if(masterproc) then
            fldname = listentry%field%name
            write(iulog,9000) f, fldname, listentry%field%units, listentry%field%numlev, &
                             listentry%avgflag(1), trim(listentry%field%long_name)
9000        format(i5, 1x, a32, 1x, a16, 1x, i4, 1x, a1, 2x, a)
         end if
         listentry=>listentry%next_entry
      end do
      nfmaster = f
      if(masterproc) write(iulog,*)'BLDFLD:nfmaster=',nfmaster

!
!  Now that masterlinkedlist is defined and we are performing a restart run
!  (after all addfld calls), construct primary and secondary hashing tables.
!
      if ( nsrest==1 ) then
         call bld_outfld_hash_tbls()
         call bld_htapefld_indices()
         return
      end if
!
! Get users logname and machine hostname
!
      if ( masterproc )then
         logname = ' '
         call shr_sys_getenv ('LOGNAME',logname,rcode)
         host = ' '
         call shr_sys_getenv ('HOST',host,rcode)
      end if
#ifdef SPMD
! PIO requires netcdf attributes have consistant values on all tasks
      call mpibcast(logname, len(logname), mpichar, 0, mpicom)
      call mpibcast(host, len(host), mpichar, 0, mpicom)
#endif
!
! Override averaging flag for all fields on a particular tape if namelist input so specifies
!
      do t=1,ptapes
         if (avgflag_pertape(t) /= ' ') then
            call h_override (t)
         end if
      end do
!
! Define field list information for all history files.  
!
      call fldlst ()
!
! Loop over max. no. of history files permitted  
!
      if ( nsrest .eq. 3 ) then
         call get_prev_time(day, sec)  ! elapased time since reference date
      else
         call get_curr_time(day, sec)  ! elapased time since reference date
      end if
      do t=1,ptapes
         nfils(t) = 0            ! no. of time samples in hist. file no. t

! Time at beginning of current averaging interval.
         
         beg_time(t) = day + sec/86400._r8
      end do

!
! Initialize history variables
!
      do t=1,ptapes
         do f=1,nflds(t)
            begdim1  = tape(t)%hlist(f)%field%begdim1
            enddim1  = tape(t)%hlist(f)%field%enddim1
            begdim2  = tape(t)%hlist(f)%field%begdim2
            enddim2  = tape(t)%hlist(f)%field%enddim2
            begdim3  = tape(t)%hlist(f)%field%begdim3
            enddim3  = tape(t)%hlist(f)%field%enddim3
            allocate(tape(t)%hlist(f)%hbuf(begdim1:enddim1,begdim2:enddim2,begdim3:enddim3))
            tape(t)%hlist(f)%hbuf = 0._r8
            if(tape(t)%hlist(f)%field%flag_xyfill .or. (avgflag_pertape(t) .eq. 'L')) then
               allocate (tape(t)%hlist(f)%nacs(begdim1:enddim1,begdim3:enddim3))
            else
               allocate (tape(t)%hlist(f)%nacs(1,begdim3:enddim3))
            end if
            tape(t)%hlist(f)%nacs(:,:)=0           
            tape(t)%hlist(f)%field%vector_compliment=0            
         end do
      end do
      ! Setup vector pairs for unstructured grid interpolation
      call setup_interpolation_and_define_vector_compliments()
      !  Initialize the sat following history subsystem
      call sat_hist_init()
      
      return
   end subroutine intht


   subroutine setup_interpolation_and_define_vector_compliments()
     use dycore, only : dycore_is
     use interp_mod,       only: setup_history_interpolation, var_is_vector_uvar, var_is_vector_vvar
     use sat_hist, only : is_satfile
     integer :: t, f, ii, ff

     if(dycore_is('UNSTRUCTURED')) then
        do t=1,ptapes
           do f=1,nflds(t)     
              tape(t)%hlist(f)%field%vector_compliment=0
           end do
        end do
! Do not interpolate IC history and sat hist files
        call setup_history_interpolation(ptapes-2)
        do t=1,ptapes
           if(.not. is_satfile(t) .and. .not. is_initfile(t)) then
              do f=1,nflds(t)
                 ii=var_is_vector_uvar(tape(t)%hlist(f)%field%name)
                 if(ii>0) then
                    do ff=1,nflds(t)
                       if(ii==var_is_vector_vvar(tape(t)%hlist(ff)%field%name)) then
                          tape(t)%hlist(f)%field%vector_compliment=ff
                          tape(t)%hlist(ff)%field%vector_compliment=-f
                       end if
                    end do
                 end if
              end do
           end if
        end do
     end if
   end subroutine setup_interpolation_and_define_vector_compliments





   subroutine restart_vars_setnames()

     restartvars(1)%name = 'rgnht'
     restartvars(1)%type = pio_int
     restartvars(1)%ndims = 1     
     restartvars(1)%dims(1) = 1  ! ptapes

     restartvars(2)%name = 'nhtfrq'
     restartvars(2)%type = pio_int
     restartvars(2)%ndims = 1
     restartvars(2)%dims(1) = 1  ! ptapes


     restartvars(3)%name = 'nflds'
     restartvars(3)%type = pio_int
     restartvars(3)%ndims = 1
     restartvars(3)%dims(1) = 1  ! ptapes

     restartvars(4)%name = 'nfils'
     restartvars(4)%type = pio_int
     restartvars(4)%ndims = 1
     restartvars(4)%dims(1) = 1  ! ptapes

     restartvars(5)%name = 'mfilt'
     restartvars(5)%type = pio_int
     restartvars(5)%ndims = 1
     restartvars(5)%dims(1) = 1  ! ptapes

     restartvars(6)%name = 'nfpath'
     restartvars(6)%type = pio_char
     restartvars(6)%ndims = 2
     restartvars(6)%dims(1) = 2  ! max_string_len
     restartvars(6)%dims(2) = 1  ! ptapes

     restartvars(7)%name = 'cpath'
     restartvars(7)%type = pio_char
     restartvars(7)%ndims = 2
     restartvars(7)%dims(1) = 2  ! max_string_len
     restartvars(7)%dims(2) = 1  ! ptapes

     restartvars(8)%name = 'nhfil'
     restartvars(8)%type = pio_char
     restartvars(8)%ndims = 2
     restartvars(8)%dims(1) = 2  ! max_string_len
     restartvars(8)%dims(2) = 1  ! ptapes

     restartvars(9)%name = 'ndens'
     restartvars(9)%type = pio_int
     restartvars(9)%ndims = 1
     restartvars(9)%dims(1) = 1  ! ptapes

     restartvars(10)%name = 'fincllonlat'
     restartvars(10)%type = pio_char
     restartvars(10)%ndims = 3
     restartvars(10)%dims(1) = 5  ! max_chars
     restartvars(10)%dims(2) = 4  ! pflds
     restartvars(10)%dims(3) = 1  ! ptapes

     restartvars(11)%name = 'ncprec'
     restartvars(11)%type = pio_int
     restartvars(11)%ndims = 1
     restartvars(11)%dims(1) = 1  ! ptapes

     restartvars(12)%name = 'ngroup'
     restartvars(12)%type = pio_int
     restartvars(12)%ndims = 1
     restartvars(12)%dims(1) = 1  ! ptapes

     restartvars(13)%name = 'beg_time'
     restartvars(13)%type = pio_double
     restartvars(13)%ndims = 1
     restartvars(13)%dims(1) = 1  ! ptapes

     restartvars(14)%name = 'fincl'
     restartvars(14)%type = pio_char
     restartvars(14)%ndims = 3
     restartvars(14)%dims(1) = 3  ! fieldname_lenp2
     restartvars(14)%dims(2) = 4  ! pflds
     restartvars(14)%dims(3) = 1  ! ptapes

     restartvars(15)%name = 'fexcl'
     restartvars(15)%type = pio_char
     restartvars(15)%ndims = 3
     restartvars(15)%dims(1) = 3  ! fieldname_lenp2
     restartvars(15)%dims(2) = 4  ! pflds
     restartvars(15)%dims(3) = 1  ! ptapes

     restartvars(16)%name = 'field_name'
     restartvars(16)%type = pio_char
     restartvars(16)%ndims = 3
     restartvars(16)%dims(1) = 6  ! max_fieldname_len
     restartvars(16)%dims(2) = 7  ! maxnflds
     restartvars(16)%dims(3) = 1  ! ptapes


     restartvars(17)%name = 'decomp_type'
     restartvars(17)%type = pio_int
     restartvars(17)%ndims = 2
     restartvars(17)%dims(1) = 7  ! maxnflds
     restartvars(17)%dims(2) = 1  ! ptapes


     restartvars(18)%name = 'numlev'
     restartvars(18)%type = pio_int
     restartvars(18)%ndims = 2
     restartvars(18)%dims(1) = 7  ! maxnflds
     restartvars(18)%dims(2) = 1  ! ptapes

     restartvars(19)%name = 'hrestpath'
     restartvars(19)%type = pio_char
     restartvars(19)%ndims = 2
     restartvars(19)%dims(1) = 2  ! max_string_len
     restartvars(19)%dims(2) = 1  ! ptapes


     restartvars(20)%name = 'hwrt_prec'
     restartvars(20)%type = pio_int
     restartvars(20)%ndims = 2
     restartvars(20)%dims(1) = 7  ! maxnflds
     restartvars(20)%dims(2) = 1  ! ptapes

     restartvars(21)%name = 'avgflag'
     restartvars(21)%type = pio_char
     restartvars(21)%ndims = 2
     restartvars(21)%dims(1) = 7  ! maxnflds
     restartvars(21)%dims(2) = 1  ! ptapes

     restartvars(22)%name = 'sampling_seq'
     restartvars(22)%type = pio_char
     restartvars(22)%ndims = 3
     restartvars(22)%dims(1) = 5  ! max_chars
     restartvars(22)%dims(2) = 7  ! maxnflds
     restartvars(22)%dims(3) = 1  ! ptapes

     restartvars(23)%name = 'long_name'
     restartvars(23)%type = pio_char
     restartvars(23)%ndims = 3
     restartvars(23)%dims(1) = 5  ! max_chars
     restartvars(23)%dims(2) = 7  ! maxnflds
     restartvars(23)%dims(3) = 1  ! ptapes

     restartvars(24)%name = 'units'
     restartvars(24)%type = pio_char
     restartvars(24)%ndims = 3
     restartvars(24)%dims(1) = 5  ! max_chars
     restartvars(24)%dims(2) = 7  ! maxnflds
     restartvars(24)%dims(3) = 1  ! ptapes

     restartvars(25)%name = 'xyfill'
     restartvars(25)%type = pio_int
     restartvars(25)%ndims = 2
     restartvars(25)%dims(1) = 7  ! maxnflds
     restartvars(25)%dims(2) = 1  ! ptapes     

     restartvars(26)%name = 'lcltod_start'
     restartvars(26)%type = pio_int
     restartvars(26)%ndims = 1
     restartvars(26)%dims(1) = 1  ! ptapes     

     restartvars(27)%name = 'lcltod_stop'
     restartvars(27)%type = pio_int
     restartvars(27)%ndims = 1
     restartvars(27)%dims(1) = 1  ! ptapes     

     restartvars(28)%name = 'fillvalue'
     restartvars(28)%type = pio_double
     restartvars(28)%ndims = 2
     restartvars(28)%dims(1) = 7  ! maxnflds
     restartvars(28)%dims(2) = 1  ! ptapes     

     restartvars(29)%name = 'mdims'
     restartvars(29)%type = pio_int
     restartvars(29)%ndims = 3
     restartvars(29)%dims(1) = 8  ! maxvarmdims
     restartvars(29)%dims(2) = 7  ! maxnflds
     restartvars(29)%dims(3) = 1  ! ptapes     


   end subroutine restart_vars_setnames

   subroutine restart_dims_setnames()
     restartdims(1)%name = 'ptapes'
     restartdims(1)%len  = ptapes

     restartdims(2)%name = 'max_string_len'
     restartdims(2)%len  = max_string_len

     restartdims(3)%name = 'fieldname_lenp2'
     restartdims(3)%len  = fieldname_lenp2

     restartdims(4)%name = 'pflds'
     restartdims(4)%len  = pflds

     restartdims(5)%name = 'max_chars'
     restartdims(5)%len  = max_chars

     restartdims(6)%name = 'max_fieldname_len'
     restartdims(6)%len  = max_fieldname_len

     restartdims(7)%name = 'maxnflds'
     restartdims(7)%len  = maxval(nflds)

     restartdims(8)%name = 'maxvarmdims'
     restartdims(8)%len  = maxvarmdims

   end subroutine restart_dims_setnames


   subroutine init_restart_history (File)

!--------------------------------------------------------------------------------------------------
!
! Arguments
!
      type(file_desc_t), intent(inout) :: File                      ! Pio file Handle
!
! Local 
!
      integer :: dimids(4), ndims
      integer :: ierr, i, k

      ! Don't need to write restart data if we have written the file this step
      where (hstwr(:))
         rgnht(:) = .false.
      elsewhere
         rgnht(:) = .true.
      end where

      if(maxval(nflds)>0) then
         call restart_vars_setnames()
         call restart_dims_setnames()

	 call pio_seterrorhandling(File, PIO_BCAST_ERROR)
         do i=1,restartdimcnt
	    ! it's possible that one or more of these have been defined elsewhere
	    ierr = pio_inq_dimid(File,restartdims(i)%name, restartdims(i)%dimid)
	    if(ierr/=PIO_NOERR) then
               ierr = pio_def_dim(File,restartdims(i)%name, restartdims(i)%len, &
                    restartdims(i)%dimid)
            end if
         end do
	 call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

         do i=1,restartvarcnt
            ndims= restartvars(i)%ndims
            do k=1,ndims
               dimids(k)=restartdims(restartvars(i)%dims(k))%dimid
            end do
            allocate(restartvars(i)%vdesc)
            ierr = pio_def_var(File, restartvars(i)%name, restartvars(i)%type, dimids(1:ndims), restartvars(i)%vdesc)
         end do
      end if
    end subroutine init_restart_history

    function restartvar_getdesc(name) result(vdesc)
      character(len=*), intent(in) :: name
      type(var_desc_t), pointer :: vdesc
      character(len=max_chars) :: errmsg
      integer :: i
      
      nullify(vdesc)
      do i=1,restartvarcnt
         if(name .eq. restartvars(i)%name) then
            vdesc=>restartvars(i)%vdesc
            exit
         end if
      end do
      if(.not.associated(vdesc)) then
         errmsg = 'Could not find restart variable '//name
         call endrun(errmsg)
      end if
    end function restartvar_getdesc


!#######################################################################

    subroutine write_restart_history ( File, & 
         yr_spec, mon_spec, day_spec, sec_spec )

      implicit none
      !--------------------------------------------------------------------------------------------------
      !
      ! Arguments
      !
      type(file_desc_t), intent(inout) :: file         ! PIO restart file pointer
      integer, intent(in), optional :: yr_spec         ! Simulation year
      integer, intent(in), optional :: mon_spec        ! Simulation month
      integer, intent(in), optional :: day_spec        ! Simulation day
      integer, intent(in), optional :: sec_spec        ! Seconds into current simulation day
      !
      ! Local workspace
      !
      integer :: ierr, t, f
      integer :: rgnht_int(ptapes), start(2), count(2), startc(3)
      type(var_desc_t), pointer :: vdesc

      type(var_desc_t), pointer ::  field_name_desc, decomp_type_desc, numlev_desc, &
           hrestpath_desc, avgflag_desc, sseq_desc, longname_desc, units_desc, hwrt_prec_desc, &
           xyfill_desc, mdims_desc, fillval_desc

      integer, allocatable :: xyfill(:,:), allmdims(:,:,:)
      integer :: maxnflds


      maxnflds = maxval(nflds)
      allocate(xyfill(maxnflds, ptapes))
      xyfill = 0
      !
      !-----------------------------------------------------------------------
      ! Write the history restart data if necessary
      !-----------------------------------------------------------------------

      rgnht_int(:) = 0

      if(.not.allocated(restarthistory_tape)) allocate(restarthistory_tape(ptapes))

      do t=1,ptapes
         ! No need to write history IC restart because it is always instantaneous
         if (is_initfile(file_index=t)) rgnht(t) = .false.
         ! No need to write restart data for empty files
         if (nflds(t) == 0) rgnht(t) = .false.
         if(rgnht(t)) then
            rgnht_int(t) = 1
            restarthistory_tape(t)%hlist => history_tape(t)%hlist

            if(associated(history_tape(t)%column)) then
               restarthistory_tape(t)%column => history_tape(t)%column
            end if
            if(associated(history_tape(t)%column_st)) then
               restarthistory_tape(t)%column_st => history_tape(t)%column_st
            end if
         end if
      end do

      if(maxval(nflds)<=0) return
      
      call wshist(rgnht)

      vdesc => restartvar_getdesc('fincl')
      ierr= pio_put_var(File, vdesc, fincl(:,1:ptapes))

      vdesc => restartvar_getdesc('fincllonlat')
      ierr= pio_put_var(File, vdesc, fincllonlat(:,1:ptapes))

      vdesc => restartvar_getdesc('fexcl')
      ierr= pio_put_var(File, vdesc, fexcl(:,1:ptapes))

      vdesc => restartvar_getdesc('rgnht')
      ierr= pio_put_var(File, vdesc, rgnht_int(1:ptapes))

      vdesc => restartvar_getdesc('nhtfrq')
      ierr= pio_put_var(File, vdesc, nhtfrq(1:ptapes))

      vdesc => restartvar_getdesc('nflds')
      ierr= pio_put_var(File, vdesc, nflds(1:ptapes))

      vdesc => restartvar_getdesc('nfils')
      ierr= pio_put_var(File, vdesc, nfils(1:ptapes))

      vdesc => restartvar_getdesc('mfilt')
      ierr= pio_put_var(File, vdesc, mfilt(1:ptapes))

      vdesc => restartvar_getdesc('nfpath')
      ierr= pio_put_var(File, vdesc, nfpath(1:ptapes))

      vdesc => restartvar_getdesc('cpath')
      ierr= pio_put_var(File, vdesc,  cpath(1:ptapes))

      vdesc => restartvar_getdesc('nhfil')
      ierr= pio_put_var(File, vdesc, nhfil(1:ptapes))

      vdesc => restartvar_getdesc('ndens')
      ierr= pio_put_var(File, vdesc, ndens(1:ptapes))
      vdesc => restartvar_getdesc('ncprec')
      ierr= pio_put_var(File, vdesc, ncprec(1:ptapes))
      vdesc => restartvar_getdesc('ngroup')
      ierr= pio_put_var(File, vdesc, ngroup(1:ptapes))
      vdesc => restartvar_getdesc('beg_time')
      ierr= pio_put_var(File, vdesc, beg_time(1:ptapes))

      vdesc => restartvar_getdesc('hrestpath')
      ierr = pio_put_var(File, vdesc, hrestpath(1:ptapes))

      vdesc => restartvar_getdesc('lcltod_start')
      ierr = pio_put_var(File, vdesc, lcltod_start(1:ptapes))

      vdesc => restartvar_getdesc('lcltod_stop')
      ierr = pio_put_var(File, vdesc, lcltod_stop(1:ptapes))

      field_name_desc => restartvar_getdesc('field_name')
      decomp_type_desc => restartvar_getdesc('decomp_type')
      numlev_desc => restartvar_getdesc('numlev')
      hwrt_prec_desc => restartvar_getdesc('hwrt_prec')

      sseq_desc => restartvar_getdesc('sampling_seq')
      longname_desc => restartvar_getdesc('long_name')
      units_desc => restartvar_getdesc('units')
      avgflag_desc => restartvar_getdesc('avgflag')
      xyfill_desc => restartvar_getdesc('xyfill')

      mdims_desc => restartvar_getdesc('mdims')
      fillval_desc => restartvar_getdesc('fillvalue')

      tape=>history_tape

      allocate(allmdims(maxvarmdims,maxval(nflds),ptapes))
      allmdims=-1

      startc(1)=1
      do t = 1,ptapes
         start(2)=t
	 startc(3)=t
         do f=1,nflds(t)
            start(1)=f
	    startc(2)=f
            ierr = pio_put_var(File, field_name_desc,startc,tape(t)%hlist(f)%field%name)
            ierr = pio_put_var(File, decomp_type_desc,start,tape(t)%hlist(f)%field%decomp_type)
            ierr = pio_put_var(File, numlev_desc,start,tape(t)%hlist(f)%field%numlev)

            ierr = pio_put_var(File, hwrt_prec_desc,start,tape(t)%hlist(f)%hwrt_prec)
            ierr = pio_put_var(File, sseq_desc,startc,tape(t)%hlist(f)%field%sampling_seq)
            ierr = pio_put_var(File, longname_desc,startc,tape(t)%hlist(f)%field%long_name)
            ierr = pio_put_var(File, units_desc,startc,tape(t)%hlist(f)%field%units)
            ierr = pio_put_var(File, avgflag_desc,start, tape(t)%hlist(f)%avgflag)

            ierr = pio_put_var(File, fillval_desc,start, tape(t)%hlist(f)%field%fillvalue)
            if(associated(tape(t)%hlist(f)%field%mdims)) then
               allmdims(1:size(tape(t)%hlist(f)%field%mdims),f,t) = tape(t)%hlist(f)%field%mdims
            end if
            if(tape(t)%hlist(f)%field%flag_xyfill)             xyfill(f,t) = 1
         end do
      end do
      ierr = pio_put_var(File, xyfill_desc, xyfill)
      ierr = pio_put_var(File, mdims_desc, allmdims)

      deallocate(xyfill, allmdims)
      return

    end subroutine write_restart_history


!#######################################################################

   subroutine read_restart_history (File)
     use pio, only : pio_inquire_variable, pio_inq_dimid, pio_inquire, pio_inq_varid, &
          pio_inquire_dimension
     
      use ppgrid,        only: begchunk, endchunk
      use phys_grid,     only: get_ncols_p
      use rgrid,         only: nlon
      use dycore,        only: dycore_is
      use cam_pio_utils, only: cam_pio_openfile, phys_decomp, dyn_decomp, dyn_stagger_decomp, get_decomp
      use ioFileMod,     only: getfil
      use sat_hist,     only: sat_hist_define, sat_hist_init
    
      use shr_sys_mod,   only: shr_sys_getenv
#if (defined SPMD)
      use spmd_utils,    only : iam, mpicom, mpichar
#endif
!
!-----------------------------------------------------------------------
!
! Arguments
!
      type(file_desc_t), intent(inout) :: File            ! unit number
!
! Local workspace
!
      integer t, f                     ! tape, field indices
      integer c                        ! chunk or lat index
      integer lenc                     ! length of useful character data
      integer numlev                   ! number of vertical levels (dimension and loop)
      integer begdim2                  ! on-node vert start index
      integer enddim2                  ! on-node vert end index
      integer ioerr                    ! error code from read()
      integer begdim1                  ! on-node dim1 start index
      integer enddim1                  ! on-node dim1 end index
      integer begdim3                  ! on-node chunk or lat start index
      integer enddim3                  ! on-node chunk or lat end index
      integer ncol                     ! number of active columns per chunk
      integer lenarr                   ! global size of array to be read
      integer flag_xyfill_int


      integer rgnht_int(ptapes)
      integer :: ierr

      character(len=max_string_len)  :: locfn       ! Local filename
      character(len=max_fieldname_len), allocatable :: tmpname(:,:)
      integer, allocatable :: decomp(:,:), tmpnumlev(:,:)
      integer, pointer :: nacs(:,:)    ! accumulation counter
      character(len=max_fieldname_len) :: fname_tmp ! local copy of field name
      character(len=max_chars), pointer :: tmpstr

      real(r8) :: sum_r8
      real(r4) :: sum_r4

      type (dim_index_3d) :: dimind, dimind_xyz, dimind_xzy    ! 3-D dimension index
      integer :: ii, i, j, k, ptapes_dimid, vsize
      integer :: beglatxy, endlatxy, beglat, endlat, beglonxy, endlonxy, beglon, endlon

      type(var_desc_t) :: vdesc, longname_desc, units_desc, avgflag_desc, sseq_desc, fillval_desc
      type(io_desc_t), pointer :: iodesc
      real(r8), allocatable :: tmpfldr8(:)
      real(r4), allocatable :: tmpfldr4(:)
      integer, allocatable :: tmpfldi(:), tmpprec(:,:)
      character(len=1), allocatable:: tmpavgflag(:, :)
      integer, allocatable :: xyfill(:,:), allmdims(:,:,:)
      integer :: nacsdimcnt, nacsval
      integer :: maxnflds, dimid
      type(master_entry), pointer :: listentry


      integer levdim             ! level dimension id
      character(len=16) :: name
      integer :: ndims, dimids(8)
      integer :: tmpdims(8), dimcnt
      integer :: value
      integer :: mtapes, mdimcnt
!
! Get users logname and machine hostname
!
      if ( masterproc )then
         logname = ' '
         call shr_sys_getenv ('LOGNAME',logname,ierr)
         host = ' '
         call shr_sys_getenv ('HOST',host,ierr)
      end if
#ifdef SPMD
      ! PIO requires netcdf attributes have consistant values on all tasks
      call mpibcast(logname, len(logname), mpichar, 0, mpicom)
      call mpibcast(host, len(host), mpichar, 0, mpicom)
#endif


      beglonxy = get_dyn_grid_parm('beglonxy')
      endlonxy = get_dyn_grid_parm('endlonxy')
      beglatxy = get_dyn_grid_parm('beglatxy')
      endlatxy = get_dyn_grid_parm('endlatxy')
      beglat = get_dyn_grid_parm('beglat')
      endlat = get_dyn_grid_parm('endlat')
      beglon = get_dyn_grid_parm('beglon')
      endlon = get_dyn_grid_parm('endlon')


      call pio_seterrorhandling(File, PIO_BCAST_ERROR)

      ierr = pio_inq_dimid(File, 'ptapes', ptapes_dimid)
      if(ierr/= PIO_NOERR) then
         if(masterproc) write(iulog,*) 'Not reading history info from restart file', ierr
         return   ! no history info in restart file
      end if
      call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

      ierr = pio_inq_dimlen(File, ptapes_dimid, mtapes)

      ierr = pio_inq_dimid(File, 'maxnflds', dimid)
      ierr = pio_inq_dimlen(File, dimid, maxnflds)

      ierr = pio_inq_dimid(File, 'maxvarmdims', dimid)
      ierr = pio_inq_dimlen(File, dimid, maxvarmdims)
      

      ierr = pio_inq_varid(File, 'rgnht', vdesc)
      ierr = pio_get_var(File, vdesc, rgnht_int(1:mtapes))      

      ierr = pio_inq_varid(File, 'nhtfrq', vdesc)
      ierr = pio_get_var(File, vdesc, nhtfrq(1:mtapes))

      ierr = pio_inq_varid(File, 'nflds', vdesc)
      ierr = pio_get_var(File, vdesc, nflds(1:mtapes))
      ierr = pio_inq_varid(File, 'nfils', vdesc)
      ierr = pio_get_var(File, vdesc, nfils(1:mtapes))
      ierr = pio_inq_varid(File, 'mfilt', vdesc)
      ierr = pio_get_var(File, vdesc, mfilt(1:mtapes))

      ierr = pio_inq_varid(File, 'nfpath', vdesc)
      ierr = pio_get_var(File, vdesc, nfpath(1:mtapes))
      ierr = pio_inq_varid(File, 'cpath', vdesc)
      ierr = pio_get_var(File, vdesc, cpath(1:mtapes))
      ierr = pio_inq_varid(File, 'nhfil', vdesc)
      ierr = pio_get_var(File, vdesc, nhfil(1:mtapes))
      ierr = pio_inq_varid(File, 'hrestpath', vdesc)
      ierr = pio_get_var(File, vdesc, hrestpath(1:mtapes))


      ierr = pio_inq_varid(File, 'ndens', vdesc)
      ierr = pio_get_var(File, vdesc, ndens(1:mtapes))
      ierr = pio_inq_varid(File, 'ncprec', vdesc)
      ierr = pio_get_var(File, vdesc, ncprec(1:mtapes))
      ierr = pio_inq_varid(File, 'ngroup', vdesc)
      ierr = pio_get_var(File, vdesc, ngroup(1:mtapes))
      ierr = pio_inq_varid(File, 'beg_time', vdesc)
      ierr = pio_get_var(File, vdesc, beg_time(1:mtapes))


      ierr = pio_inq_varid(File, 'fincl', vdesc)
      ierr = pio_get_var(File, vdesc, fincl(:,1:mtapes))

      ierr = pio_inq_varid(File, 'fincllonlat', vdesc)
      ierr = pio_get_var(File, vdesc, fincllonlat(:,1:mtapes))

      ierr = pio_inq_varid(File, 'fexcl', vdesc)
      ierr = pio_get_var(File, vdesc, fexcl(:,1:mtapes))

      ierr = pio_inq_varid(File, 'lcltod_start', vdesc)
      ierr = pio_get_var(File, vdesc, lcltod_start(1:mtapes))

      ierr = pio_inq_varid(File, 'lcltod_stop', vdesc)
      ierr = pio_get_var(File, vdesc, lcltod_stop(1:mtapes))




      allocate(tmpname(maxnflds, mtapes), decomp(maxnflds, mtapes), tmpnumlev(maxnflds,mtapes))
      ierr = pio_inq_varid(File, 'field_name', vdesc)
      ierr = pio_get_var(File, vdesc, tmpname)

      ierr = pio_inq_varid(File, 'decomp_type', vdesc)
      ierr = pio_get_var(File, vdesc, decomp)
      ierr = pio_inq_varid(File, 'numlev', vdesc)
      ierr = pio_get_var(File, vdesc, tmpnumlev)

      allocate(tmpprec(maxnflds,mtapes))
      ierr = pio_inq_varid(File, 'hwrt_prec',vdesc)
      ierr = pio_get_var(File, vdesc, tmpprec(:,:))

      allocate(xyfill(maxnflds,mtapes))
      ierr = pio_inq_varid(File, 'xyfill', vdesc)
      ierr = pio_get_var(File, vdesc, xyfill)

      allocate(allmdims(maxvarmdims,maxnflds,mtapes))
      ierr = pio_inq_varid(File, 'mdims', vdesc)
      ierr = pio_get_var(File, vdesc, allmdims)


      ierr = pio_inq_varid(File, 'avgflag', avgflag_desc)      

      ierr = pio_inq_varid(File, 'long_name', longname_desc)      
      ierr = pio_inq_varid(File, 'units', units_desc)      
      ierr = pio_inq_varid(File, 'sampling_seq', sseq_desc)      

      ierr = pio_inq_varid(File, 'fillvalue', fillval_desc)      


      rgnht(:)=.false.

      allocate(history_tape(mtapes))

      tape => history_tape

      do t=1,mtapes
         nullify(tape(t)%column)
         nullify(tape(t)%column_st)

         if(rgnht_int(t)==1) rgnht(t)=.true.
         

         call strip_null(nfpath(t))
         call strip_null(cpath(t))
         call strip_null(hrestpath(t))
         allocate(tape(t)%hlist(nflds(t)))
         
         do f=1,nflds(t)
            nullify(tape(t)%hlist(f)%field%mdims)
            ierr = pio_get_var(File,fillval_desc, (/f,t/), tape(t)%hlist(f)%field%fillvalue)
            ierr = pio_get_var(File,avgflag_desc, (/f,t/), tape(t)%hlist(f)%avgflag)
            ierr = pio_get_var(File,longname_desc, (/1,f,t/), tape(t)%hlist(f)%field%long_name)
            ierr = pio_get_var(File,units_desc, (/1,f,t/), tape(t)%hlist(f)%field%units)
            tape(t)%hlist(f)%field%sampling_seq(1:max_chars) = ' '
            ierr = pio_get_var(File,sseq_desc, (/1,f,t/), tape(t)%hlist(f)%field%sampling_seq)
            call strip_null(tape(t)%hlist(f)%field%sampling_seq)
            if(xyfill(f,t) ==1) then
               tape(t)%hlist(f)%field%flag_xyfill=.true.
            else
               tape(t)%hlist(f)%field%flag_xyfill=.false.
            end if
            call strip_null(tmpname(f,t))
            tape(t)%hlist(f)%field%name = tmpname(f,t)
            tape(t)%hlist(f)%field%decomp_type = decomp(f,t)
            tape(t)%hlist(f)%field%numlev = tmpnumlev(f,t)
            tape(t)%hlist(f)%hwrt_prec = tmpprec(f,t)

            mdimcnt = count(allmdims(:,f,t)>0)
            if(mdimcnt>0) then
               allocate(tape(t)%hlist(f)%field%mdims(mdimcnt))
               tape(t)%hlist(f)%field%mdims = allmdims(1:mdimcnt,f,t)
            end if

         end do
      end do
      deallocate(tmpname, tmpnumlev, tmpprec, decomp, xyfill)

      do t=1,ptapes
         do f=1,nflds(t)
            nullify(tape(t)%hlist(f)%field%colperdim3)
            numlev   = tape(t)%hlist(f)%field%numlev
            select case (tape(t)%hlist(f)%field%decomp_type)
            case (phys_decomp)
               tape(t)%hlist(f)%field%begdim2= 1
               tape(t)%hlist(f)%field%enddim2= numlev
               tape(t)%hlist(f)%field%begdim3 = begchunk
               tape(t)%hlist(f)%field%enddim3 = endchunk
               allocate (tape(t)%hlist(f)%field%colperdim3(begchunk:endchunk))
               do c=begchunk,endchunk
                  ncol = get_ncols_p(c)
                  tape(t)%hlist(f)%field%colperdim3(c) = ncol
               end do
               tape(t)%hlist(f)%field%begdim1  = 1
               tape(t)%hlist(f)%field%enddim1  = pcols
            case (dyn_decomp)               
               tape(t)%hlist(f)%field%begdim1  = beglonxy
               tape(t)%hlist(f)%field%enddim1  = endlonxy
               tape(t)%hlist(f)%field%begdim3  = beglatxy
               tape(t)%hlist(f)%field%enddim3  = endlatxy
               tape(t)%hlist(f)%field%begdim2  = 1
               tape(t)%hlist(f)%field%enddim2  = numlev
               
               if (.not. (dycore_is('LR') .or. dycore_is('UNSTRUCTURED'))) then
                  allocate (tape(t)%hlist(f)%field%colperdim3(beglat:endlat))
                  do c=beglat,endlat
                     tape(t)%hlist(f)%field%colperdim3(c) = nlon(c)
                  end do
               end if
            case (dyn_stagger_decomp)
                  tape(t)%hlist(f)%field%begdim1 = beglonxy
                  tape(t)%hlist(f)%field%enddim1 = endlonxy
                  tape(t)%hlist(f)%field%begdim2 = 1
                  tape(t)%hlist(f)%field%enddim2 = numlev
                  tape(t)%hlist(f)%field%begdim3 = beglatxy
                  tape(t)%hlist(f)%field%enddim3 = endlatxy
            case default
               write(iulog,*)'READ_RESTART_HISTORY: bad decomp_type=',tape(t)%hlist(f)%field%decomp_type
               call endrun ()
            end select

            begdim1  = tape(t)%hlist(f)%field%begdim1
            enddim1  = tape(t)%hlist(f)%field%enddim1
            begdim2  = tape(t)%hlist(f)%field%begdim2
            enddim2  = tape(t)%hlist(f)%field%enddim2
            begdim3  = tape(t)%hlist(f)%field%begdim3
            enddim3  = tape(t)%hlist(f)%field%enddim3

            allocate(tape(t)%hlist(f)%hbuf(begdim1:enddim1,begdim2:enddim2,begdim3:enddim3))

            nullify(tape(t)%hlist(f)%varid)
            nullify(tape(t)%hlist(f)%nacs)
            if(tape(t)%hlist(f)%field%flag_xyfill .or. (avgflag_pertape(t)=='L')) then
               allocate (tape(t)%hlist(f)%nacs(begdim1:enddim1,begdim3:enddim3))
            else
               allocate(tape(t)%hlist(f)%nacs(1,begdim3:enddim3))
            end if
            ! initialize all buffers to zero - this will be overwritten later by the
            ! data in the history restart file if it exists.
            call h_zero(f,t)


         end do
      end do
!
!-----------------------------------------------------------------------
! Read history restart files
!-----------------------------------------------------------------------
!
! Loop over the total number of history files declared and
! read the pathname for any history restart files
! that are present (if any). Test to see if the run is a restart run
! AND if any history buffer regen files exist (rgnht=.T.). Note, rgnht 
! is preset to false, reset to true in routine WSDS if hbuf restart files
! are written and saved in the master restart file. Each history buffer
! restart file is then obtained.
! Note: some f90 compilers (e.g. SGI) complain about I/O of 
! derived types which have pointer components, so explicitly read each one.
! 
      do t=1,mtapes
         if (rgnht(t)) then
!
! Open history restart file
!
            call getfil (hrestpath(t), locfn)
            call cam_pio_openfile(tape(t)%File, locfn, 0)
!
! Read history restart file
!
            !     Read the mdims on file - they are the last dims defined and 'lev' is always the first.
            ierr=pio_inq_dimid (tape(t)%File, 'lev', levdim)
            ierr = pio_inquire(tape(t)%File, ndims)
            
            do i=levdim,ndims
               ierr = pio_inquire_dimension(tape(t)%File, i, name, value)
               call register_hist_mdim(name,value)
            enddo

            do f=1,nflds(t)  
               begdim1    =  tape(t)%hlist(f)%field%begdim1
               enddim1    =  tape(t)%hlist(f)%field%enddim1
               begdim2    =  tape(t)%hlist(f)%field%begdim2
               enddim2    =  tape(t)%hlist(f)%field%enddim2
               begdim3    =  tape(t)%hlist(f)%field%begdim3
               enddim3    =  tape(t)%hlist(f)%field%enddim3
               numlev     =  tape(t)%hlist(f)%field%numlev

               fname_tmp = strip_suffix(tape(t)%hlist(f)%field%name)
               if(masterproc) write(iulog, *) 'Reading history variable ',fname_tmp
               ierr = pio_inq_varid(tape(t)%File, fname_tmp, vdesc)


               if(.not. associated(tape(t)%hlist(f)%field%mdims)) then
                  ierr = pio_inquire_variable(tape(t)%File, vdesc, ndims=ndims,dimids=dimids)
                  dimcnt=0
                  do i=1,ndims
                     if(dimids(i)>=levdim) then
                        dimcnt=dimcnt+1
                        tmpdims(dimcnt)=dimids(i)-levdim+1
                     end if
                  end do
                  if(dimcnt>0) then
                     allocate(tape(t)%hlist(f)%field%mdims(dimcnt))
                     tape(t)%hlist(f)%field%mdims(:) = tmpdims(1:dimcnt)
                     if(dimcnt > maxvarmdims) maxvarmdims=dimcnt
                  end if
               end if
               call get_decomp(iodesc, tape(t)%hlist(f)%field, pio_double)
               call pio_setframe(vdesc, int(1,kind=PIO_OFFSET))
               call pio_read_darray(tape(t)%File, vdesc, iodesc, tape(t)%hlist(f)%hbuf, ierr)

               ierr = pio_inq_varid(tape(t)%File, trim(fname_tmp)//'_nacs', vdesc)
	       ierr = pio_inq_varndims(tape(t)%File, vdesc, nacsdimcnt)

               if(nacsdimcnt>0) then
                  call get_decomp(iodesc, tape(t)%hlist(f)%field, pio_int, numlev_in=1)
                  allocate(tape(t)%hlist(f)%nacs(begdim1:enddim1,begdim3:enddim3))

                  nacs       => tape(t)%hlist(f)%nacs(:,:)
                  call pio_read_darray(tape(t)%File, vdesc, iodesc, nacs, ierr)
               else
                  allocate(tape(t)%hlist(f)%nacs(1,begdim3:enddim3))
                  ierr = pio_get_var(tape(t)%File, vdesc, nacsval)
                  tape(t)%hlist(f)%nacs(1,:)= nacsval
               end if



            end do
!          
! Done reading this history restart file
!
            call pio_closefile (tape(t)%File)

        end if  ! rgnht(t)
        call column_init(t, nflds(t))
      end do     ! end of do mtapes loop

      


!
! If the history files are partially complete (contain less than
! mfilt(t) time samples, then get the files and open them.)
!
! NOTE:  No need to perform this operation for IC history files or empty files
!

      do t=1,mtapes
         if (is_initfile(file_index=t)) then
            ! Initialize filename specifier for IC file
            hfilename_spec(t) = '%c.cam' // trim(inst_suffix) // '.i.%y-%m-%d-%s.nc'
            nfils(t) = 0
         else if (nflds(t) == 0) then
            nfils(t) = 0
         else
            if (nfils(t) > 0) then
               call getfil (cpath(t), locfn)
               call cam_pio_openfile(tape(t)%File, locfn, PIO_WRITE)
               call h_inquire (t)
               if(is_satfile(t)) then
                  !  Initialize the sat following history subsystem
                  call sat_hist_init()
                  call sat_hist_define(tape(t)%File)
               end if
            end if
            !
            ! If the history file is full, close the current unit
            !
            if (nfils(t) >= mfilt(t)) then
               if (masterproc) then
                  write(iulog,*)'READ_RESTART_HISTORY: nf_close(',t,')=',nhfil(t), mfilt(t)
               end if
               do f=1,nflds(t)
                  deallocate(tape(t)%hlist(f)%varid)
                  nullify(tape(t)%hlist(f)%varid)
               end do
               call pio_closefile(tape(t)%File)
               nfils(t) = 0
            end if
         end if
      end do

      ! Setup vector pairs for unstructured grid interpolation
      call setup_interpolation_and_define_vector_compliments()

      if(mtapes/=ptapes .and. masterproc) then
         write(iulog,*) ' WARNING: Restart file ptapes setting ',mtapes,' not equal to model setting ',ptapes
      end if

      return
   end subroutine read_restart_history

!#######################################################################

   character(len=max_string_len) function get_hfilepath( tape )
!
!----------------------------------------------------------------------- 
! 
! Purpose: Return full filepath of history file for given tape number
! This allows public read access to the filenames without making
! the filenames public data.
!
!----------------------------------------------------------------------- 
!
  integer, intent(in) :: tape  ! Tape number

  get_hfilepath = cpath( tape )
  end function get_hfilepath

!#######################################################################

   character(len=max_string_len) function get_hist_restart_filepath( tape )
!
!----------------------------------------------------------------------- 
! 
! Purpose: Return full filepath of restart file for given tape number
! This allows public read access to the filenames without making
! the filenames public data.
!
!----------------------------------------------------------------------- 
!
  integer, intent(in) :: tape  ! Tape number

  get_hist_restart_filepath = hrestpath( tape )
  end function get_hist_restart_filepath

!#######################################################################

  integer function get_ptapes( )
!
!----------------------------------------------------------------------- 
! 
! Purpose: Return the number of tapes being used.
! This allows public read access to the number of tapes without making
! ptapes public data.
!
!----------------------------------------------------------------------- 
!
  get_ptapes = ptapes
  end function get_ptapes

  recursive function get_entry_by_name(listentry, name) result(entry)
    type(master_entry),  pointer :: listentry
    character(len=*), intent(in) :: name ! variable name
    type(master_entry), pointer :: entry
    
    if(associated(listentry)) then
       if(listentry%field%name .eq. name) then
          entry => listentry
       else
          entry=>get_entry_by_name(listentry%next_entry, name)
       end if
    else
       nullify(entry)
    end if
  end function get_entry_by_name


!#######################################################################

   subroutine fldlst ()

     use dycore, only : dycore_is
!
!----------------------------------------------------------------------- 
! 
! Purpose: Define the contents of each history file based on namelist input for initial or branch
! run, and restart data if a restart run.
!          
! Method: Use arrays fincl and fexcl to modify default history tape contents.
!         Then sort the result alphanumerically for later use by OUTFLD to
!         allow an n log n search time.
!
!---------------------------Local variables-----------------------------
!
      integer t, f                   ! tape, field indices
      integer ff                     ! index into include, exclude and fprec list
      character(len=fieldname_len) :: name ! field name portion of fincl (i.e. no avgflag separator)
      character(len=max_fieldname_len) :: mastername ! name from masterlist field
      character(len=max_chars) :: tmpcolumn_name,latlonname,latlonnamep1 ! tmp char fields
      character(len=1) :: avgflag    ! averaging flag
      character(len=1) :: prec_wrt   ! history buffer write precision flag
      character(len=6) :: prec_str
      character(len=32) :: fldname

      type (hentry) :: tmp           ! temporary used for swapping
      type (column_info) :: tmpcolumn ! temporary used for swapping

      type(master_entry), pointer :: listentry

      
!
! First ensure contents of fincl, fexcl, and fwrtpr are all valid names
!
      do t=1,ptapes
         f = 1
         do while (f < pflds .and. fincl(f,t) /= ' ')
            name = getname (fincl(f,t))
            mastername=''
            listentry => get_entry_by_name(masterlinkedlist, name)
            if(associated(listentry)) mastername = listentry%field%name
            if (name /= mastername) then
               write(iulog,*)'FLDLST: ', trim(name), ' in fincl(', f, ') not found'
               call endrun
            end if
            f = f + 1
         end do

         f = 1
         do while (f < pflds .and. fexcl(f,t) /= ' ')
            mastername=''
            listentry => get_entry_by_name(masterlinkedlist, fexcl(f,t))
            if(associated(listentry)) mastername = listentry%field%name

            if (fexcl(f,t) /= mastername) then
               write(iulog,*)'FLDLST: ', fexcl(f,t), ' in fexcl(', f, ') not found'
               call endrun
            end if
            f = f + 1
         end do

         f = 1
         do while (f < pflds .and. fwrtpr(f,t) /= ' ')
            name = getname (fwrtpr(f,t))
            mastername=''
            listentry => get_entry_by_name(masterlinkedlist, name)
            if(associated(listentry)) mastername = listentry%field%name
            if (name /= mastername) then
               write(iulog,*)'FLDLST: ', trim(name), ' in fwrtpr(', f, ') not found'
               call endrun
            end if
            do ff=1,f-1                 ! If duplicate entry is found, stop
               if (trim(name) == trim(getname(fwrtpr(ff,t)))) then
                  write(iulog,*)'FLDLST: Duplicate field ', name, ' in fwrtpr'
                  call endrun
               end if
            end do
            f = f + 1
         end do
      end do

      nflds(:) = 0
      ! IC history file is to be created, set properties
      if(is_initfile()) then
         hfilename_spec(ptapes) = '%c.cam' // trim(inst_suffix) // '.i.%y-%m-%d-%s.nc'

         ncprec(ptapes) = pio_double
         ndens (ptapes) = 1
         mfilt (ptapes) = 1
      end if




      do t=1,ptapes
!
! Add the field to the tape if specified via namelist (FINCL[1-ptapes]), or if
! it is on by default and was not excluded via namelist (FEXCL[1-ptapes]).
! Also set history buffer accumulation and output precision values according
! to the values specified via namelist (FWRTPR[1-ptapes])
! or, if not on the list, to the default values given by ndens(t).
!
         listentry => masterlinkedlist
         do while(associated(listentry))
            mastername = listentry%field%name
            call list_index (fincl(1,t), mastername, ff)

            if (ff > 0) then
               nflds(t) = nflds(t)+1
            else if ((.not. empty_htapes) .or. (is_initfile(file_index=t))) then
               call list_index (fexcl(1,t), mastername, ff)
               if (ff == 0 .and. listentry%actflag(t)) then
                  nflds(t) = nflds(t)+1
               end if
            end if
            listentry=>listentry%next_entry
         end do
      enddo
!
! Determine total number of active history tapes
!
      if (masterproc) then
         do t=1,ptapes
            if (nflds(t)  ==  0) then
               write(iulog,*)'FLDLST: Tape ',t,' is empty'
            end if
         end do
      endif
      allocate(history_tape(ptapes))
      tape=>history_tape


      do t=1,ptapes
         nullify(tape(t)%hlist)
         nullify(tape(t)%column)
         nullify(tape(t)%column_st)
! Now we have a field count and can allocate
         if(nflds(t)> 0) allocate(tape(t)%hlist(nflds(t)))

         do ff=1,nflds(t)
            nullify(tape(t)%hlist(ff)%hbuf)
            nullify(tape(t)%hlist(ff)%nacs)
            nullify(tape(t)%hlist(ff)%varid)
         end do


         nflds(t) = 0 ! recount to support array based method
         listentry => masterlinkedlist
         do while(associated(listentry))
            mastername = listentry%field%name

            call list_index (fwrtpr(1,t), mastername, ff)
            if (ff > 0) then
               prec_wrt = getflag(fwrtpr(ff,t))
            else
               prec_wrt = ' '
            end if

            call list_index (fincl(1,t), mastername, ff)

            if (ff > 0) then
               avgflag = getflag (fincl(ff,t))
               call inifld (t, listentry, avgflag,  prec_wrt)
            else if ((.not. empty_htapes) .or. (is_initfile(file_index=t))) then
               call list_index (fexcl(1,t), mastername, ff)
               if (ff == 0 .and. listentry%actflag(t)) then
                  call inifld (t, listentry, ' ', prec_wrt)
               else
                  listentry%actflag(t) = .false.
               end if
            else
               listentry%actflag(t) = .false.
            end if
            listentry=>listentry%next_entry

         end do
!
! If column output is specified make sure there are some fields defined
! for that tape
!
         if (nflds(t) .eq. 0 .and. fincllonlat(1,t) .ne. ' ') then
            write(iulog,*) 'FLDLST: Column output is specified for tape ',t,' but no fields defined for that tape.'
            call endrun()
         else
            call column_init(t, nflds(t))
         end if
!
! Specification of tape contents now complete.  Sort each list of active 
! entries for efficiency in OUTFLD.  Simple bubble sort.
!
         do f=nflds(t)-1,1,-1
            do ff=1,f

               if (tape(t)%hlist(ff)%field%name > tape(t)%hlist(ff+1)%field%name) then
            
                  tmp = tape(t)%hlist(ff)
                  tape(t)%hlist(ff  ) = tape(t)%hlist(ff+1)
                  tape(t)%hlist(ff+1) = tmp

               else if (tape(t)%hlist(ff  )%field%name == tape(t)%hlist(ff+1)%field%name) then
                  
                  write(iulog,*)'FLDLST: Duplicate field ', tape(t)%hlist(ff  )%field%name
                  write(iulog,*)'t,ff,name=',t,ff,tape(t)%hlist(ff  )%field%name
                  call endrun
            
               end if

            end do
         end do
!
! Bubble sort columns check for duplicates and rebuild field_column_names in newly sorted order
!
         if (ngroup(t) .gt. 0) then
            do i=ngroup(t)-1,1,-1
               do ff=1,i
                  latlonname=trim(tape(t)%column(ff)%lon_name) // "_" // trim(tape(t)%column(ff)%lat_name)
                  latlonnamep1=trim(tape(t)%column(ff+1)%lon_name) // "_" // trim(tape(t)%column(ff+1)%lat_name)
                  if (trim(latlonname) > trim(latlonnamep1)) then
                     tmpcolumn = tape(t)%column(ff)
                     tape(t)%column(ff) = tape(t)%column(ff+1)
                     tape(t)%column(ff+1) = tmpcolumn
                     if(dycore_is('LR')) then
                        tmpcolumn = tape(t)%column_st(ff)
                        tape(t)%column_st(ff) = tape(t)%column_st(ff+1)
                        tape(t)%column_st(ff+1) = tmpcolumn
                     end if
                  else if (trim(latlonname) == trim(latlonnamep1)) then
                     write(iulog,*)'FLDLST: Duplicate column entry for tape ',t,'duplicate column name is ',trim(latlonname)
                     call endrun
                  end if
               end do
            end do
            do f=1,nflds(t)
               do i=1,ngroup(t)
                  if (ngroup(t) .gt. 0) then
                     tape(t)%hlist(f)%field_column_name(i) = trim(tape(t)%hlist(f)%field%name) // "_" // &
                          trim(tape(t)%column(i)%lon_name) // "_" // trim(tape(t)%column(i)%lat_name)
                  else 
                     tape(t)%hlist(f)%field_column_name(1) = ' '
                  end if
               end do
            end do
         end if

         if (masterproc) then
            if (nflds(t) > 0) then
               write(iulog,*) ' '
               write(iulog,*)'FLDLST: History file ', t, ' contains ', nflds(t), ' fields'

               if (is_initfile(file_index=t)) then
                  write(iulog,*) ' Write frequency:                 ',inithist,' (INITIAL CONDITIONS)'
               else if (nflds(t) > 0) then
                  if (nhtfrq(t) == 0 .and. t == 1) then
                     write(iulog,*) ' Write frequency:                  MONTHLY'
                  else
                     write(iulog,*) ' Write frequency:                 ',nhtfrq(t)
                  end if
               end if

               write(iulog,*) ' Filename specifier:              ', trim(hfilename_spec(t))

               prec_str = 'double'
               if (ndens(t) == 2) prec_str = 'single'
               write(iulog,*) ' Output precision:                ', prec_str
               write(iulog,*) ' Number of time samples per file: ', mfilt(t)

               write(iulog,*)' Included fields are:'
            end if
   
            do f=1,nflds(t)
               if (ngroup(t) .gt. 0) then
                  if (f .eq. 1) write(iulog,*) '   Fields on this tape will be output as column data (FIELD_LON_LAT)'
                  do i = 1, ngroup(t)
                     ff = (f-1)*ngroup(t) + i
                     fldname = tape(t)%hlist(f)%field_column_name(i)
                     write(iulog,9000) ff, fldname, tape(t)%hlist(f)%field%units, tape(t)%hlist(f)%field%numlev, &
                             tape(t)%hlist(f)%avgflag, tape(t)%hlist(f)%field%long_name

                  end do
               else
                  fldname = tape(t)%hlist(f)%field%name
                  write(iulog,9000) f, fldname, tape(t)%hlist(f)%field%units, tape(t)%hlist(f)%field%numlev, &
                             tape(t)%hlist(f)%avgflag, tape(t)%hlist(f)%field%long_name

               end if

            end do

         end if

      end do    ! do t=1,ptapes

9000  format(i5, 1x, a32, 1x, a16, 1x, i4, 1x, a1, 2x, 256a)
!
! Packing density, ndens: With netcdf, only 1 (nf_double) and 2 (pio_real)
! are allowed
!
      do t=1,ptapes
         if (ndens(t) == 1) then
            ncprec(t) = pio_double
         else if (ndens(t) == 2) then
            ncprec(t) = pio_real
         else
            call endrun ('FLDLST: ndens must be 1 or 2')
         end if

      end do
!
!  Now that masterlinkedlist is defined, construct primary and secondary hashing
!  tables.
!
      call bld_outfld_hash_tbls()
      call bld_htapefld_indices()

      return
   end subroutine fldlst

!#######################################################################
   subroutine inifld (t, listentry, avgflag, prec_wrt)
!
!----------------------------------------------------------------------- 
! 
! Purpose: Add a field to the active list for a history tape
! 
! Method: Copy the data from the master field list to the active list for the tape
!         Also: define mapping arrays from (col,chunk) -> (lon,lat)
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------


!
! Arguments
!
      integer, intent(in) :: t            ! history tape index

      type(master_entry), pointer :: listentry

      character*1, intent(in) :: avgflag  ! averaging flag
      character*1, intent(in) :: prec_wrt ! history output precision flag
!
! Local workspace
!
      integer :: i                  ! index
      integer :: n                  ! field index on defined tape


!
! Ensure that it is not to late to add a field to the history tape
!
      if (htapes_defined) then
         call endrun ('INIFLD: Attempt to add field '//listentry%field%name//' after history files set')
      end if


      nflds(t) = nflds(t) + 1
      n = nflds(t)
!
! Copy field info.
!
      if(n > size(tape(t)%hlist)) then
         write(iulog,*) 'tape field miscount error ', n, size(tape(t)%hlist)
         call endrun()
      end if

      tape(t)%hlist(n)%field = listentry%field

      select case (prec_wrt)
      case (' ')
         if (ndens(t) == 1) then
            tape(t)%hlist(n)%hwrt_prec = 8
         else
            tape(t)%hlist(n)%hwrt_prec = 4
         end if
      case ('4')
         tape(t)%hlist(n)%hwrt_prec = 4
         if (masterproc) then
            write(iulog,*) 'INIFLD: Output data type for ', tape(t)%hlist(n)%field%name, &
                       ' is real*4'
         end if
      case ('8')
         tape(t)%hlist(n)%hwrt_prec = 8
         if (masterproc) then
            write(iulog,*) 'INIFLD: Output data type for ', tape(t)%hlist(n)%field%name, &
                       ' is real*8'
         end if
      case default
         call endrun ('INIFLD: unknown prec_wrt='//prec_wrt)
      end select
!
! Override the default averaging (masterlist) averaging flag if non-blank
!
      if (avgflag == ' ') then
         tape(t)%hlist(n)%avgflag = listentry%avgflag(t)
         tape(t)%hlist(n)%time_op = listentry%time_op(t)
      else
         tape(t)%hlist(n)%avgflag = avgflag
         select case (avgflag)
         case ('A')
            tape(t)%hlist(n)%time_op = 'mean'
         case ('B')
            tape(t)%hlist(n)%time_op = 'mean00z'
         case ('I')
            tape(t)%hlist(n)%time_op = ' '
         case ('X')
            tape(t)%hlist(n)%time_op = 'maximum'
         case ('M')
            tape(t)%hlist(n)%time_op = 'minimum'
         case('L')
            tape(t)%hlist(n)%time_op = LT_DESC
         case default
            call endrun ('INIFLD: unknown avgflag='//avgflag)
         end select
      end if

#ifdef HDEBUG
      write(iulog,*)'HDEBUG: ',__LINE__,' field ', tape(t)%hlist(n)%field%name, ' added as ', 'field number ', n,' on tape ', t
      write(iulog,*)'units=',tape(t)%hlist(n)%field%units
      write(iulog,*)'numlev=',tape(t)%hlist(n)%field%numlev
      write(iulog,*)'avgflag=',tape(t)%hlist(n)%avgflag
      write(iulog,*)'time_op=',tape(t)%hlist(n)%time_op
      write(iulog,*)'hwrt_prec=',tape(t)%hlist(n)%hwrt_prec
#endif

      return
   end subroutine inifld


   subroutine column_init(t, nflds)
     use dyn_grid, only: get_dyn_grid_parm_real2d, get_dyn_grid_parm_real1d
     use dycore, only : dycore_is
     use spmd_utils, only : mpi_sum, mpi_integer, iam, mpicom
     use phys_grid, only : get_ncols_p, get_lat_all_p, get_lon_all_p
     use ppgrid,        only: begchunk, endchunk

     integer, intent(in) :: t, nflds
     integer :: ff, i

     real(r8), allocatable :: lontmp   (:)    ! Temp array holding longitude values 
     real(r8), allocatable :: lontmp_st(:)   ! Temp array holding longitude values 
     character(len=max_chars) :: lonlatname(pflds), lonname(pflds), latname(pflds) ! variable name
     integer :: lonind   (pflds,2) !   beginning and ending longitude range
     integer :: latind   (pflds,2) !   beginning and ending latitude range

     integer :: lonind_st(pflds,2) !   beginning and ending staggered longitude range
     integer :: latind_st(pflds,2) !   beginning and ending staggeredlatitude range
     real(r8), pointer :: londeg(:,:), londeg_st(:,:), latdeg(:), latdeg_st(:)
     character(len=max_chars) :: tmpname(pflds) ! variable name (will not be used)

     integer :: group              !   number of group
     integer :: splon, plat, plon
     integer, pointer :: hmask(:,:)
     integer :: lats(pcols), lons(pcols)
     integer :: cnt, cnttot, ierr, lchnk, ncol, j
     integer :: beglatxy,endlatxy,beglonxy,endlonxy
     

     plon = get_dyn_grid_parm('plon')
     splon = get_dyn_grid_parm('splon')
     plat = get_dyn_grid_parm('plat')

     beglonxy = get_dyn_grid_parm('beglonxy')
     beglatxy = get_dyn_grid_parm('beglatxy')
     endlonxy = get_dyn_grid_parm('endlonxy')
     endlatxy = get_dyn_grid_parm('endlatxy')


     londeg => get_dyn_grid_parm_real2d('londeg')
     latdeg => get_dyn_grid_parm_real1d('latdeg')
     if(dycore_is('LR')) then
        londeg_st => get_dyn_grid_parm_real2d('londeg_st')
        latdeg_st => get_dyn_grid_parm_real1d('latdeg_st')
     end if
     !
     ! Setup column information if this field will be written as group
     ! First verify the column information in the namelist
     !

     ff = 1

     allocate(lontmp(plon))
     if(dycore_is('LR')) then
        allocate(lontmp_st(splon))
     end if
     do while (fincllonlat(ff,t) /= ' ')

        lonlatname(ff) = trim(fincllonlat(ff,t))
        call getlatind(lonlatname(ff),latind(ff,1),latind(ff,2),latname(ff),latdeg     ,plat)
	
        do i = 1,plon
           lontmp(i) = londeg(i,1)
           if (lontmp(i) .lt. 0._r8) lontmp(i) = lontmp(i) + 360._r8
        end do
        call getlonind(lonlatname(ff),lonind(ff,1),lonind(ff,2),lonname(ff),lontmp,plon)
#ifdef HDEBUG
        write(iulog,*)'HDEBUG: ',__LINE__,' closest lon lat range for group ',lonlatname(ff),' is ',lonind(ff,:),latind(ff,:)
        write(iulog,*) ff,t,londeg(lonind(ff,1),1),londeg(lonind(ff,2),1),latdeg(latind(ff,1)),latdeg(latind(ff,2))
#endif
        if(dycore_is('LR')) then
           call getlatind(lonlatname(ff),latind_st(ff,1),latind_st(ff,2),tmpname(ff),latdeg_st     ,plat-1)
           do i = 1,splon
              lontmp_st(i) = londeg_st(i,1)
              if (lontmp_st(i) .lt. 0._r8) lontmp_st(i) = lontmp_st(i) + 360._r8
           end do
           call getlonind(lonlatname(ff),lonind_st(ff,1),lonind_st(ff,2),tmpname(ff),lontmp_st,splon )
#ifdef HDEBUG
           write(iulog,*)'HDEBUG: ',__LINE__,' closest staggered lon lat range for group ',lonlatname(ff), &
                ' is ',lonind_st(ff,:),latind_st(ff,:)
           write(iulog,*) londeg_st(lonind_st(ff,1),1),londeg_st(lonind_st(ff,2),1), &
                latdeg_st(latind_st(ff,1)),latdeg_st(latind_st(ff,2))
#endif
        end if
        ff = ff + 1
     end do
     deallocate(lontmp)

     if(dycore_is('LR')) then      
        deallocate(lontmp_st)
     endif

     group = ff-1
     ngroup(t)=0

     if (group .gt. 0) then 
        if(collect_column_output(t)) then
           allocate(tape(t)%column(1))
           tape(t)%column(1)%num_lats = 0
           tape(t)%column(1)%num_lons = 0

           allocate(hmask(plon, plat))
         
           hmask = 0
           cnt = 0
           do j=1,plat
              do i=1,plon
                 do ff=1,group
                    if(hmask(i,j)==0 .and. i>=lonind(ff,1) .and. i <= lonind(ff,2) .and. &
                       j>=latind(ff,1) .and. j <= latind(ff,2)) then
                       cnt = cnt+1
                       hmask(i,j)=cnt
                    end if
                 end do                 
              end do
           end do

           tape(t)%column(1)%column_cnt = cnt
           allocate(tape(t)%column(1)%hmask(pcols,begchunk:endchunk))
           tape(t)%column(1)%hmask=0
           allocate(tape(t)%column(1)%hmask_dyn(endlonxy-beglonxy+1,endlatxy-beglatxy+1))
           tape(t)%column(1)%hmask_dyn=hmask(beglonxy:endlonxy,beglatxy:endlatxy)

           do lchnk=begchunk,endchunk
              ncol = get_ncols_p(lchnk)
              call get_lat_all_p(lchnk, pcols, lats)
              call get_lon_all_p(lchnk, pcols, lons)
              do i=1,ncol
                 tape(t)%column(1)%hmask(i,lchnk)=hmask(lons(i),lats(i))
              end do
           end do
           tape(t)%column(1)%htape=t
           deallocate(hmask)


        else
           ngroup(t)=group
           allocate (tape(t)%column   (group))
           if(dycore_is('LR')) then
              allocate (tape(t)%column_st(group))
           endif
           do i=1,nflds
              allocate (tape(t)%hlist(i)%field_column_name(group))
           end do
           do ff = 1, group
              nullify(tape(t)%column(ff)%hmask)
              tape(t)%column(ff)%lat_name=trim(latname(ff))
              tape(t)%column(ff)%lon_name=trim(lonname(ff))
              tape(t)%column(ff)%columnlon(:)=lonind(ff,:)
              tape(t)%column(ff)%columnlat(:)=latind(ff,:)
              tape(t)%column(ff)%num_lats = &
                   tape(t)%column(ff)%columnlat(2)-tape(t)%column(ff)%columnlat(1)+1
              tape(t)%column(ff)%num_lons = &
                   tape(t)%column(ff)%columnlon(2)-tape(t)%column(ff)%columnlon(1)+1
              do i=1,nflds
                 tape(t)%hlist(i)%field_column_name(ff) = &
                      trim(tape(t)%hlist(i)%field%name) // "_" // trim(lonname(ff)) // "_" // trim(latname(ff))
              end do
              if(dycore_is('LR')) then
                 tape(t)%column_st(ff)%lat_name=trim(latname(ff)) // "_st"
                 tape(t)%column_st(ff)%lon_name=trim(lonname(ff)) // "_st"
                 tape(t)%column_st(ff)%columnlon(:)=lonind_st(ff,:)
                 tape(t)%column_st(ff)%columnlat(:)=latind_st(ff,:)
                 tape(t)%column_st(ff)%num_lats = &
                      tape(t)%column_st(ff)%columnlat(2)-tape(t)%column_st(ff)%columnlat(1)+1
                 tape(t)%column_st(ff)%num_lons = &
                      tape(t)%column_st(ff)%columnlon(2)-tape(t)%column_st(ff)%columnlon(1)+1
              endif
           
           end do
        end if

     else

        allocate (tape(t)%column(1))
        nullify(tape(t)%column(1)%hmask)
        do i=1,nflds
           allocate (tape(t)%hlist(i)%field_column_name(1))
           tape(t)%hlist(i)%field_column_name(1) = ' '
        end do
        tape(t)%column(1)%lat_name=' '
        tape(t)%column(1)%lon_name=' '
        tape(t)%column(1)%columnlon(:)=0
        tape(t)%column(1)%columnlat(:)=0
        tape(t)%column(1)%num_lats=0
        tape(t)%column(1)%num_lons=0
     end if


   end subroutine column_init








!#######################################################################

   subroutine strip_null(str)
     character(len=*), intent(inout) :: str
     do i=1,len(str)
        if(ichar(str(i:i))==0) str(i:i)=' '
     end do
   end subroutine strip_null

   character(len=max_fieldname_len) function strip_suffix (name)
!
!---------------------------------------------------------- 
! 
! Purpose:  Strip "&IC" suffix from fieldnames if it exists
!          
!----------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name
!
! Local workspace
!
      integer :: n
!
!-----------------------------------------------------------------------
!
      strip_suffix = ' '

      do n = 1,fieldname_len
         strip_suffix(n:n) = name(n:n)
         if(name(n+1:n+1         ) == ' '                       ) return
         if(name(n+1:n+fieldname_suffix_len) == fieldname_suffix) return
      end do

      strip_suffix(fieldname_len+1:max_fieldname_len) = name(fieldname_len+1:max_fieldname_len)

      return

   end function strip_suffix

!#######################################################################

   character(len=fieldname_len) function getname (inname)
!
!----------------------------------------------------------------------- 
! 
! Purpose: retrieve name portion of inname
!          
! Method:  If an averaging flag separater character is present (":") in inname, 
!          lop it off
! 
!-------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: inname
!
! Local workspace
!
      integer :: length
      integer :: i
   
      length = len (inname)
      
      if (length < fieldname_len .or. length > fieldname_lenp2) then
         write(iulog,*) 'GETNAME: bad length=',length
         call endrun
      end if
   
      getname = ' '
      do i=1,fieldname_len
         if (inname(i:i) == ':') exit
         getname(i:i) = inname(i:i)
      end do
      
      return
   end function getname

!#######################################################################

   subroutine getlatind(inname,beglatind,endlatind,latname,latdeg,mlat)
!
!----------------------------------------------------------------------- 
! 
! Purpose: retrieve closest model lat index given north south latitude string
!          
! Author: John Truesdale
! 
!-------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in)  :: inname
      character(len=*), intent(out) :: latname
      integer, intent(out) :: beglatind,endlatind
      integer, intent(in) :: mlat
      real(r8), intent(in) :: latdeg(mlat)
!
! Local workspace
!
      character(len=max_chars) str1,str2,tmpstr(4)
      integer :: i,j,marker,latind(2),tmplen(4)
      real(r8) :: latdegree(2)

!
!-------------------------------------------------------------------------------
!
      str1=' '
      str2=' '
      tmpstr(:)= ' '

!
! make sure _ separator is present
!      
      if (scan(inname,'_').eq.0) then
         write(iulog,*)'GETLATIND: Improperly formatted column string.  Missing underscore character (xxxE_yyyS) ', &
                   inname,scan(inname,'_')
         call endrun
      end if
!
! split inname string into lat lon strings
!      
      marker=scan(inname,'_')
      str1=inname(:marker-1)
      str2=inname(marker+1:)

!
! split ranges of lats( or lons) into seperate substrings. Substrings 1,2 will contain lats(lons) 
! from portion of string before underscore
! if a single column and not a range is specified substrings 1 and 2 will be the same
! and substrings 3,4 will contain lats(lons) from portion of main string after underscore character
! if a single column and not a range is specified substrings 3 and 4 will be the same
!
      if (scan(str1,':') .ne. 0) then
         marker=scan(str1,':')
         tmpstr(1)=str1(:marker-1)
         tmpstr(2)=str1(marker+1:)
      else
         tmpstr(1)=str1
         tmpstr(2)=str1
      end if
      if (scan(str2,':') .ne. 0) then
         marker=scan(str2,':')
         tmpstr(3)=str2(:marker-1)
         tmpstr(4)=str2(marker+1:)
      else
         tmpstr(3)=str2
         tmpstr(4)=str2
      end if

! check format of substrings - (number followed by single character north/south east/west designation)

      do i = 1,4
         tmplen(i)=len_trim(tmpstr(i))
         if (verify(tmpstr(i),"0123456789.").ne.tmplen(i) .or. verify(tmpstr(i)(tmplen(i):tmplen(i)),"eEwWnNsS") .ne. 0) then
            write(iulog,*)'GETLATIND (2): Improperly formatted column string. ',&
                 inname,verify(tmpstr(i),"0123456789."),tmplen(i),&
                 verify(tmpstr(i)(tmplen(i):tmplen(i)),"eEwWnNsS")
            call endrun
         end if
      end do

! find latitude substrings and put into temporary work space tmpstr(1:2)

      if (verify(tmpstr(1)(tmplen(1):tmplen(1)),"nNsSs").eq.0  .and. verify(tmpstr(2)(tmplen(2):tmplen(2)),"nNsSs").eq.0 ) then

      else if (verify(tmpstr(3)(tmplen(3):tmplen(3)),"nNsSs").eq.0 .and. verify(tmpstr(3)(tmplen(3):tmplen(3)),"nNsSs").eq.0) then
         tmpstr(1) = tmpstr(3)
         tmplen(1)=tmplen(3)
         tmpstr(2) = tmpstr(4)
         tmplen(2)=tmplen(4)
      else
         call endrun ('GETLATIND (3): Improperly formatted column string. '//inname)
      end if
!
! convert lat substrings to real and if in southern hemisphere make negative
!
      do i = 1,2
         read(tmpstr(i)(1:tmplen(i)-1),*) latdegree(i)
         if (verify(tmpstr(i)(tmplen(i):tmplen(i)),'sS').eq.0) then
            latdegree(i)=(-1._r8)*latdegree(i)
         end if
!
! Make sure specified latitudes is in bounds
!
         if (latdegree(i) .lt. -90._r8 .or. latdegree(i) .gt. 90._r8) then
            write(iulog,*)'GETLATIND: latitude for column namelist is out of range (-90 .. 90) value=',latdegree(i)
            call endrun
         endif
!
! Find closest lat index for each substring
!
         latind(i)=mlat
         do j = 1, mlat-1
            if ( abs(latdeg(j)-latdegree(i)) .lt.  &
                 abs(latdeg(j+1)-latdegree(i))) then
               latind(i)=j
               exit
            endif
         enddo
      end do
!
! output begining and ending latitude indices. If just a column is specified then beginning latitude index will be the same as the
! ending latitude index  

!      
      if (latind(1) .le. latind(2) ) then
         beglatind =latind(1)
         endlatind =latind(2)
      else 
         beglatind = latind(2)
         endlatind = latind(1)
      end if
      if (beglatind .eq. endlatind) then
         latname = 'LAT_'//trim(tmpstr(1))
      else
         if (latind(1) .le. latind(2) ) then
            latname = 'LAT_'//trim(tmpstr(1)) // "_to_" // trim(tmpstr(2))
         else 
            latname = 'LAT_'//trim(tmpstr(2)) // "_to_" // trim(tmpstr(1))
         end if
      end if
      ! one last sanity check
      if(scan(latname,' ') > 0) then
         latname(scan(latname,' '):) = ' '
      end if

      return
   end subroutine getlatind
!#######################################################################

   subroutine getlonind (inname,beglonind,endlonind,lonname,londeg,mlon)
!
!----------------------------------------------------------------------- 
! 
! Purpose: retrieve closes model lat index given north south latitude string
!          
! Method:  If an averaging flag separater character is present (":") in inname, 
!          lop it off
! 
!-------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*) inname,lonname
      integer :: beglonind,endlonind
      integer :: mlon
      real(r8) :: londeg(mlon)
!
! Local workspace
!
      character(len=max_chars) str1,str2,tmpstr(4)
      integer :: i,j,marker,lonind(2),tmplen(4)
      real(r8) :: londegree(2)
      real(r8) :: min, tmp
!
!-------------------------------------------------------------------------------
!
      str1=' '
      str2=' '
      tmpstr(:)= ' '
!
! make sure _ separator is present
!      
      if (scan(inname,'_').eq.0) then
         write(iulog,*)'GETLONIND: Improperly formatted column string.  Missing underscore character (xxxE_yyyS) ',&
              inname,scan(inname,'_')
         call endrun
      end if
!
! split string in to lat lon strings
!      
      marker=scan(inname,'_')
      str1=inname(:marker-1)
      str2=inname(marker+1:)

!
! split ranges of lats( or lons) into seperate substrings. Substrings 1,2 will contain lats(lons)
! from portion of string before underscore
! if a single column and not a range is specified substrings 1 and 2 will be the same
! and substrings 3,4 will contain lats(lons) from portion of main string after underscore character
! if a single column and not a range is specified substrings 3 and 4 will be the same
!
      if (scan(str1,':') .ne. 0) then
         marker=scan(str1,':')
         tmpstr(1)=str1(:marker-1)
         tmpstr(2)=str1(marker+1:)
      else
         tmpstr(1)=str1
         tmpstr(2)=str1
      end if

      if (scan(str2,':') .ne. 0) then
         marker=scan(str2,':')
         tmpstr(3)=str2(:marker-1)
         tmpstr(4)=str2(marker+1:)
      else
         tmpstr(3)=str2
         tmpstr(4)=str2
      end if

! check format of substrings - (number followed by single character north/south east/west designation)

      do i = 1,4
         tmplen(i)=len_trim(tmpstr(i))
         if (verify(tmpstr(i),"0123456789.").ne.tmplen(i) .or. verify(tmpstr(i)(tmplen(i):tmplen(i)),"eEwWnNsS") .ne. 0) then
            write(iulog,*)'GETLONIND (2): Improperly formatted column string. ', &
                 inname,verify(tmpstr(i),"0123456789."),tmplen(i), &
                 verify(tmpstr(i)(tmplen(i):tmplen(i)),"eEwWnNsS")
            call endrun
         end if
      end do

! find longitude substrings and put into temporary work space tmpstr(1:2)

      if (verify(tmpstr(1)(tmplen(1):tmplen(1)),"eEwW").eq.0  .and. verify(tmpstr(2)(tmplen(2):tmplen(2)),"eEwW").eq.0 ) then

      else if (verify(tmpstr(3)(tmplen(3):tmplen(3)),"eEwW").eq.0 .and. verify(tmpstr(3)(tmplen(3):tmplen(3)),"eEwW").eq.0) then
         tmpstr(1) = tmpstr(3)
         tmplen(1)=tmplen(3)
         tmpstr(2) = tmpstr(4)
         tmplen(2)=tmplen(4)
      else
         call endrun ('GETLONIND (3): Improperly formatted column string. '//inname)
      end if
!
! convert lon substrings to real and make sure its degrees east
!
      do i = 1,2
         read(tmpstr(i)(1:tmplen(i)-1),*) londegree(i)
         if (verify(tmpstr(i)(tmplen(i):tmplen(i)),'wW').eq.0) then
            londegree(i) = 360._r8 - londegree(i)
         end if
!
! Make sure specified longitudes are in bounds
!
         if (londegree(i) .lt. 0._r8 .or. londegree(i) .gt. 360._r8) then
            write(iulog,*)'GETLONIND: longitude for column namelist is out of range (0 .. 360) value=',londegree(i)
            call endrun
         endif
!
! Find closest lon index for each substring.  If just a column is specified then beginning longitude index will be the same as the
! ending longitude index  
!
         min  = 1.e+36_r8
         do j = 1, mlon
            tmp = abs(londeg(j)-londegree(i))
            if ( tmp .le. min) then
               min       = tmp
               lonind(i) = j
            endif
         end do
      end do
!
! output begining and ending longitude indicies. If just a column is specified then beginning longitude index
! will be the same as the ending longitude index  
!      
      if (lonind(1) .le. lonind(2) ) then
         beglonind =lonind(1)
         endlonind =lonind(2)
      else 
         beglonind = lonind(2)
         endlonind = lonind(1)
      end if


      if (beglonind .eq. endlonind) then
         lonname = 'LON_'//trim(tmpstr(1))
      else
         if (lonind(1) .le. lonind(2) ) then
            lonname = 'LON_'//trim(tmpstr(1)) // "_to_" // trim(tmpstr(2))
         else 
            lonname = 'LON_'//trim(tmpstr(2)) // "_to_" // trim(tmpstr(1))
         end if
      end if


      return
    end subroutine getlonind

!#######################################################################

   character(len=1) function getflag (inname)
!
!----------------------------------------------------------------------- 
! 
! Purpose: retrieve flag portion of inname
!          
! Method:  If an averaging flag separater character is present (":") in inname, 
!          return the character after it as the flag
! 
!-------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: inname   ! character string
!
! Local workspace
!
      integer :: length         ! length of inname
      integer :: i              ! loop index

      length = len (inname)

      if (length /= fieldname_lenp2) then
         write(iulog,*) 'GETFLAG: bad length=',length
         call endrun
      end if

      getflag = ' '
      do i=1,fieldname_lenp2-1
         if (inname(i:i) == ':') then
            getflag = inname(i+1:i+1)
            exit
         end if
      end do

      return
   end function getflag

!#######################################################################

   subroutine list_index (list, name, index)
!
! Input arguments
!
      character(len=*) , intent(in) :: list(pflds) ! input list of names, possibly ":" delimited
      character(len=max_fieldname_len), intent(in) :: name ! name to be searched for
!
! Output arguments
!
      integer, intent(out) :: index               ! index of "name" in "list"
!
! Local workspace
!
      character(len=fieldname_len) :: listname    ! input name with ":" stripped off.
      integer f                       ! field index

      index = 0
      do f=1,pflds
!
! Only list items
!
         listname = getname (list(f))
         if (listname == ' ') exit
         if (listname == name) then
            index = f
            exit
         end if
      end do
      
      return
   end subroutine list_index

!#######################################################################

   subroutine outfld (fname, field, idim, c)
     use cam_history_buffers, only: hbuf_accum_inst, hbuf_accum_add, &
          hbuf_accum_add00z, hbuf_accum_max, hbuf_accum_min, &
          hbuf_accum_addlcltime, dim_index_2d
     
!
!----------------------------------------------------------------------- 
! 
! Purpose: Accumulate (or take min, max, etc. as appropriate) input field
!          into its history buffer for appropriate tapes
! 
! Method: Check 'masterlist' whether the requested field 'fname' is active
!         on one or more history tapes, and if so do the accumulation.
!         If not found, return silently.
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: fname ! Field name--should be 8 chars long

      integer, intent(in) :: idim           ! Longitude dimension of field array
      integer, intent(in) :: c              ! chunk (physics) or latitude (dynamics) index

      real(r8), intent(in) :: field(idim,*) ! Array containing field values
!
! Local variables
!
      integer :: t, f                ! tape, field indices
      integer :: fl, fu              ! upper, lower indices used in binary search thru sorted list
      integer :: begdim1             ! on-node dim1 start index
      integer :: enddim1             ! on-node dim1 end index
      integer :: begdim2             ! on-node dim2 start index
      integer :: enddim2             ! on-node dim2 end index
      integer :: begdim3
      integer :: enddim3
      integer :: endi                ! ending longitude index (reduced grid)

      character*(max_fieldname_len) :: fname_loc  ! max-char equivalent of fname
      character*1 :: avgflag         ! averaging flag
      
      real(r8),pointer :: hbuf(:,:)      ! history buffer
      integer, pointer :: nacs(:)    ! accumulation counter
      type (dim_index_2d) :: dimind  ! 2-D dimension index
      logical :: flag_xyfill         ! non-applicable xy points flagged with fillvalue
      integer :: ff                  ! masterlist index pointer
      integer :: ierr, ncol
!-----------------------------------------------------------------------

      tape => history_tape

      fname_loc = fname

      ff = get_masterlist_indx(fname_loc)

!
!  If ( ff < 0 ), the field is not defined on the masterlist. This check
!  is necessary because of coding errors calling outfld without first defining
!  the field on masterlist.
!
      if ( ff < 0 ) then
         return
      end if
!
!  Next, check to see whether this field is active on one or more history
!  tapes.
!
      if ( .not. masterlist(ff)%thisentry%act_sometape )  then
         return
      end if
!
! Note, the field may be on any or all of the history files (primary
! and auxiliary).
!
!      write(iulog,*)'fname_loc=',fname_loc
      do 40 t=1,ptapes
         if ( .not. masterlist(ff)%thisentry%actflag(t)) cycle
         f = masterlist(ff)%thisentry%htapeindx(t)
!
! Update history buffer
!
         begdim1 = tape(t)%hlist(f)%field%begdim1
         enddim1 = tape(t)%hlist(f)%field%enddim1
         begdim2 = tape(t)%hlist(f)%field%begdim2
         enddim2 = tape(t)%hlist(f)%field%enddim2
         begdim3 = tape(t)%hlist(f)%field%begdim3
         enddim3 = tape(t)%hlist(f)%field%enddim3

         avgflag = tape(t)%hlist(f)%avgflag
         flag_xyfill = tape(t)%hlist(f)%field%flag_xyfill

         nacs   => tape(t)%hlist(f)%nacs(:,c)

         
         hbuf => tape(t)%hlist(f)%hbuf(:,:,c)

         if(associated(tape(t)%hlist(f)%field%colperdim3)) then         
            endi    = tape(t)%hlist(f)%field%colperdim3(c)
         else
            endi = tape(t)%hlist(f)%field%enddim1 - tape(t)%hlist(f)%field%begdim1 + 1
         end if
         dimind = dim_index_2d (1,endi,begdim2,enddim2)

         select case (avgflag)

         case ('I') ! Instantaneous
            call hbuf_accum_inst (hbuf, field, nacs, dimind, idim, flag_xyfill,tape(t)%hlist(f)%field%fillvalue)

         case ('A') ! Time average

            call hbuf_accum_add (hbuf, field, nacs, dimind, idim, flag_xyfill,tape(t)%hlist(f)%field%fillvalue)

         case ('B') ! Time average only 00z values

            call hbuf_accum_add00z (hbuf, field, nacs, dimind, idim, flag_xyfill,tape(t)%hlist(f)%field%fillvalue)

         case ('X') ! Maximum over time

            call hbuf_accum_max (hbuf, field, nacs, dimind, idim, flag_xyfill,tape(t)%hlist(f)%field%fillvalue)

         case ('M') ! Minimum over time

            call hbuf_accum_min (hbuf, field, nacs, dimind, idim, flag_xyfill,tape(t)%hlist(f)%field%fillvalue)

         case ('L')

            call hbuf_accum_addlcltime(hbuf, field, nacs, dimind, idim, flag_xyfill, &
                                       tape(t)%hlist(f)%field%fillvalue, &
                                       c, tape(t)%hlist(f)%field%decomp_type, &
                                       lcltod_start(t),lcltod_stop(t))

        case default

            call endrun ('OUTFLD: invalid avgflag='//avgflag)

         end select
40    continue

      return
   end subroutine outfld

!#######################################################################

   logical function is_initfile (file_index)
!
!------------------------------------------------------------------------ 
! 
! Purpose: to determine:
!
!   a) if an IC file is active in this model run at all
!       OR,
!   b) if it is active, is the current file index referencing the IC file
!      (IC file is always at ptapes)
! 
!------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in), optional :: file_index ! index of file in question

      is_initfile = .false.

      if (present(file_index)) then
         if (inithist /= 'NONE' .and. file_index == ptapes) is_initfile = .true.
      else
         if (inithist /= 'NONE'                           ) is_initfile = .true.
      end if

      return

   end function is_initfile

!#######################################################################

   integer function strcmpf (name1, name2)
!
!----------------------------------------------------------------------- 
! 
! Purpose: Return the lexical difference between two strings
! 
! Method: Use ichar() intrinsic as we loop through the names
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      character(len=max_fieldname_len), intent(in) :: name1, name2 ! strings to compare
      integer n                                     ! loop index
   
      do n=1,max_fieldname_len
         strcmpf = ichar(name1(n:n)) - ichar(name2(n:n))
         if (strcmpf /= 0) exit
      end do

      return
   end function strcmpf

!#######################################################################

   subroutine h_inquire (t)
     use pio, only : pio_inquire, pio_inq_dimid, pio_inquire_dimension, pio_inq_varid, &
          pio_inquire_variable, pio_inq_attlen
     use dycore, only : dycore_is
!
!----------------------------------------------------------------------- 
! 
! Purpose: Ensure that the proper variables are on a history file
! 
! Method: Issue the appropriate netcdf wrapper calls
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: t   ! tape index
!
! Local workspace
!
      integer f,ff               ! field index
      integer ret                ! return value from function call


      integer :: ierr
      integer :: len
      integer :: mdimsize

!
!
! Dimension id's
!
      tape => history_tape


      
!
! Create variables for model timing and header information 
!
      if(.not. is_satfile(t)) then
         ierr=pio_inq_varid (tape(t)%File,'ndcur   ',    tape(t)%ndcurid)
         ierr=pio_inq_varid (tape(t)%File,'nscur   ',    tape(t)%nscurid)
         ierr=pio_inq_varid (tape(t)%File,'nsteph  ',    tape(t)%nstephid)

         ierr=pio_inq_varid (tape(t)%File,'time_bnds',   tape(t)%tbndid)
         ierr=pio_inq_varid (tape(t)%File,'date_written',tape(t)%date_writtenid)
         ierr=pio_inq_varid (tape(t)%File,'time_written',tape(t)%time_writtenid)
#if ( defined BFB_CAM_SCAM_IOP )
         ierr=pio_inq_varid (tape(t)%File,'tsec    ',tape(t)%tsecid)
         ierr=pio_inq_varid (tape(t)%File,'bdate   ',tape(t)%bdateid)
#endif
         if (.not. is_initfile(file_index=t) ) then
            ! Don't write the GHG/Solar forcing data to the IC file.  It is never
            ! read from that file so it's confusing to have it there.
            ierr=pio_inq_varid (tape(t)%File,'co2vmr  ',    tape(t)%co2vmrid)
            ierr=pio_inq_varid (tape(t)%File,'ch4vmr  ',    tape(t)%ch4vmrid)
            ierr=pio_inq_varid (tape(t)%File,'n2ovmr  ',    tape(t)%n2ovmrid)
            ierr=pio_inq_varid (tape(t)%File,'f11vmr  ',    tape(t)%f11vmrid)
            ierr=pio_inq_varid (tape(t)%File,'f12vmr  ',    tape(t)%f12vmrid)
            ierr=pio_inq_varid (tape(t)%File,'sol_tsi ',    tape(t)%sol_tsiid)
            if (cam_chempkg_is('waccm_mozart').or.cam_chempkg_is('waccm_ghg')) then
               ierr=pio_inq_varid (tape(t)%File,'f107    ',    tape(t)%f107id)
               ierr=pio_inq_varid (tape(t)%File,'f107a   ',    tape(t)%f107aid)
               ierr=pio_inq_varid (tape(t)%File,'kp      ',    tape(t)%kpid)
               ierr=pio_inq_varid (tape(t)%File,'ap      ',    tape(t)%apid)
            endif
         end if
      end if
      ierr=pio_inq_varid (tape(t)%File,'date    ',    tape(t)%dateid)
      ierr=pio_inq_varid (tape(t)%File,'datesec ',    tape(t)%datesecid)
      ierr=pio_inq_varid (tape(t)%File,'time    ',    tape(t)%timeid)


!
! Obtain variable name from ID which was read from restart file
!
      do f=1,nflds(t)
         if(.not.associated(tape(t)%hlist(f)%varid)) then
            allocate(tape(t)%hlist(f)%varid(max(1,ngroup(t))))
         end if
!
! If this field will be put out as columns then get column names for field
!
         if (ngroup(t) .gt. 0) then
            do i=1,ngroup(t)
               ierr=pio_inq_varid (tape(t)%File, tape(t)%hlist(f)%field_column_name(i), tape(t)%hlist(f)%varid(i))
               ierr=pio_get_att(tape(t)%File, tape(t)%hlist(f)%varid(i), 'basename',tape(t)%hlist(f)%field%name)
            end do
         else   
            ierr=pio_inq_varid (tape(t)%File,tape(t)%hlist(f)%field%name, tape(t)%hlist(f)%varid(1))
         end if
         if(tape(t)%hlist(f)%field%numlev>1) then
            ierr = pio_inq_attlen(tape(t)%File,tape(t)%hlist(f)%varid(1),'mdims', mdimsize)
            if(.not. associated(tape(t)%hlist(f)%field%mdims)) then
               allocate(tape(t)%hlist(f)%field%mdims(mdimsize))
            end if
            ierr=pio_get_att(tape(t)%File,tape(t)%hlist(f)%varid(1),'mdims', &
                 tape(t)%hlist(f)%field%mdims(1:mdimsize))
            if(mdimsize>maxvarmdims) maxvarmdims=mdimsize
         end if

      end do

      if(masterproc) then
         write(iulog,*)'H_INQUIRE: Successfully opened netcdf file '
      end if

      return
   end subroutine h_inquire

!#######################################################################

   subroutine add_default (name, tindex, flag)
!
!----------------------------------------------------------------------- 
! 
! Purpose: Add a field to the default "on" list for a given history file
! 
! Method: 
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name  ! field name
      character(len=1), intent(in) :: flag  ! averaging flag

      integer, intent(in) :: tindex         ! history tape index
!
! Local workspace
!
      integer :: t            ! file index
      integer :: f            ! field index
      type(master_entry), pointer :: listentry

      if (htapes_defined) then
         call endrun ('ADD_DEFAULT: Attempt to add hist default '//trim(name)//' after history files set')
      end if
!
! Check validity of input arguments
!
      if (tindex > ptapes) then
         write(iulog,*)'ADD_DEFAULT: tape index=', tindex, ' is too big'
         call endrun
      end if

! Add to IC file if tindex = 0, reset to ptapes
      if (tindex == 0) then
         t = ptapes
         if ( .not. is_initfile(file_index=t) ) return
      else
         t = tindex
      end if

      if (flag /= ' ' .and. flag /= 'A' .and. flag /= 'B' .and. flag /= 'I' .and. &
          flag /= 'X' .and. flag /= 'M' .and. flag /= 'L') then

         call endrun ('ADD_DEFAULT: unknown averaging flag='//flag)
      end if
!
! Look through master list for input field name.  When found, set active
! flag for that tape to true.  Also set averaging flag if told to use other
! than default.
!
      listentry => get_entry_by_name(masterlinkedlist, trim(name))
      if(.not.associated(listentry)) then
         call endrun ('ADD_DEFAULT: field='//name//' not found')
      end if
      listentry%actflag(t) = .true.
      if (flag /= ' ') then
         listentry%avgflag(t) = flag
         select case (flag)
         case ('A')
            listentry%time_op(t) = 'mean'
         case ('B')
            listentry%time_op(t) = 'mean00z'
         case ('I')
            listentry%time_op(t) = ' '
         case ('X')
            listentry%time_op(t) = 'maximum'
         case ('M')
            listentry%time_op(t) = 'minimum'
         case ('L')
            listentry%time_op(t) = LT_DESC
         case default
            call endrun ('ADD_DEFAULT: unknown avgflag='//flag)
         end select
      end if
      
      return
   end subroutine add_default

!#######################################################################

   subroutine h_override (t)
!
!----------------------------------------------------------------------- 
! 
! Purpose: Override default history tape contents for a specific tape
!
! Method: Copy the flag into the master field list
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: t         ! history tape index
!
! Local workspace
!
      character(len=1) :: avgflg       ! lcl equiv of avgflag_pertape(t) (to address xlf90 compiler bug)

      type(master_entry), pointer :: listentry


      avgflg = avgflag_pertape(t)

      
      listentry=>masterlinkedlist
      do while(associated(listentry))
         listentry%avgflag(t) = avgflag_pertape(t)
         select case (avgflg)
         case ('A')
            listentry%time_op(t) = 'mean'
         case ('B')
            listentry%time_op(t) = 'mean00z'
         case ('I')
            listentry%time_op(t) = ' '
         case ('X')
            listentry%time_op(t) = 'maximum'
         case ('M')
            listentry%time_op(t) = 'minimum'
         case ('L')
            listentry%time_op(t) = LT_DESC
            
         case default
            call endrun ('H_OVERRIDE: unknown avgflag='//avgflag_pertape(t))
         end select
         listentry=>listentry%next_entry
      end do

   end subroutine h_override
         
!#######################################################################

   subroutine h_define (t, restart)
     !
     !----------------------------------------------------------------------- 
     ! 
     ! Purpose: Define contents of history file t
     ! 
     ! Method: Issue the required netcdf wrapper calls to define the history file contents
     ! 
     !-----------------------------------------------------------------------
     use pspect, only : ptrn, ptrk, ptrm
     use rgrid, only : nlon, wnummax

     use dyn_grid, only : get_dyn_grid_parm_real2d, get_dyn_grid_parm_real1d
     use dycore, only : dycore_is
     use time_manager, only: get_step_size, get_ref_date, timemgr_get_calendar_cf
     use filenames,    only: caseid
     use string_utils, only: to_upper
     use abortutils,   only: endrun
     use dycore,       only: dycore_is
     use physconst,     only: pi
     use spmd_utils,    only: iam, nsmps, npes, mpicom
     use cam_pio_utils, only: cam_pio_createfile, get_phys_decomp
     use interp_mod,    only: get_interp_lat, get_interp_lon, latlon_interpolation, add_interp_attributes
     use cosp_share,   only: srmid_cosp, dbzemid_cosp, scol_cosp, sza_cosp, htmisrmid_cosp, prsmid_cosp, &
          taumid_cosp, htmid_cosp, taumid_cosp_modis, prslim_cosp, taulim_cosp, taulim_cosp_modis, &
          dbzelim_cosp, srlim_cosp, htmisrlim_cosp, htlim_cosp
     use sat_hist,     only: sat_hist_define
     
     use ppgrid, only : begchunk,endchunk, pver, pverp
     use phys_grid, only : get_rlat_all_p, get_rlon_all_p     

     !-----------------------------------------------------------------------

     !
     ! Input arguments
     !
     integer, intent(in) :: t   ! tape index
     logical, intent(in) :: restart
     !
     ! Local workspace
     !
     integer :: i, j            ! longitude, latitude indices
     integer :: f               ! field index
     integer :: ff              ! varid index for fields output by column
     integer :: numlev          ! number of vertical levels (dimension and loop)
     integer :: ncreal          ! real data type for output
     integer :: dtime           ! timestep size
     integer :: ndbase = 0      ! days component of base time
     integer :: nsbase = 0      ! seconds component of base time
     integer :: nbdate          ! base date in yyyymmdd format
     integer :: bdate           ! base date in yyyymmdd format
     integer :: nbsec           ! time of day component of base date [seconds]
     integer :: yr, mon, day    ! year, month, day components of a date


     integer :: slatdim         ! staggered latitude dimension
     integer :: slondim         ! staggered longitude dimension
     type(var_desc_t) :: slatvar         ! variable id for staggered lat
     type(var_desc_t) :: slonvar         ! variable id for staggered lon
     type(var_desc_t) :: wsid            ! Staggered latitude weight ID

     real(r8), pointer:: londeg_st(:,:)    ! Staggered grid point array (lon)
     real(r8), pointer:: latdeg_st(: )   ! Staggered grid point array (lat)
     real(r8), pointer:: w_staggered(:)  ! Staggered location weights

     real(r8) :: ailev(pverp)      ! interface level values
     real(r8), pointer :: latdeg(:)    ! degrees gaussian latitudes 
     real(r8), pointer :: latdeg2d(:,:)    ! degrees latitudes for ncols
     real(r8), pointer :: londeg(:,:)    ! degrees longitude
     real(r8), pointer :: w(:)    ! grid weights
     real(r8), parameter :: radtodeg = 180.0_r8/pi
     real(r8), pointer :: harea(:)

     character(len=max_chars) str ! character temporary 
     character(len=max_chars) :: fname_tmp ! local copy of field name
     character(len=max_chars) calendar             ! Calendar type
     !
     ! netcdf variables
     !
     integer ret                ! function return value
     integer :: timdim             ! unlimited dimension id
     type(var_desc_t) :: latvar             ! latitude variable id
     type(var_desc_t) :: lonvar             ! longitude variable id
     type(var_desc_t) :: areavar             ! cell area variable id
     type(var_desc_t), pointer :: glatvar(:)     ! column latitude variable id
     type(var_desc_t), pointer :: glonvar(:)     ! column longitude variable id
     type(var_desc_t) :: ps0var             ! variable id for PS0
     integer :: chardim            ! character dimension id

     real(r8) :: alev(pver)        ! level values (pascals)
     real(r8), allocatable :: alon(:)        ! longitude values (degrees)
     integer  :: dimenchar(2)       ! character dimension ids


     integer :: nacsdims(2)        ! dimension ids for nacs (used in restart file)
     integer :: hdims, tdims
     integer londim             ! longitude dimension id
     integer latdim             ! latitude dimension id
     integer, allocatable:: grouplondim(:) ! longitude dimension id
     integer, allocatable:: grouplatdim(:) ! latitude dimension id
     integer levdim             ! level dimension id
     integer ilevdim            ! interface dimension id

     integer bnddim            ! bounds dimension id
     type(var_desc_t):: levvar             ! level variable id
     type(var_desc_t):: ilevvar            ! intfc variable id

     type(var_desc_t) :: hyaiid             ! hybrid A coef. intfc var id
     type(var_desc_t) :: hybiid             ! hybrid B coef. intfc var id
     type(var_desc_t) :: hyamid             ! hybrid A coef. level var id
     type(var_desc_t) :: hybmid             ! hybrid B coef. level var id
     type(var_desc_t) :: ntrmid             ! M truncation parameter var id
     type(var_desc_t) :: ntrnid             ! N truncation parameter var id
     type(var_desc_t) :: ntrkid             ! K truncation parameter var id



     integer :: dimindex(8)         ! dimension ids for variable declaration

     integer, pointer :: mdims(:)
     integer :: mdimsize


     type(var_desc_t), pointer :: glatvar_st(:)     ! column staggered latitude variable id
     type(var_desc_t), pointer :: glonvar_st(:)     ! column staggered longitude variable id
     integer, allocatable :: grouplondim_st(:) ! staggered longitude dimension id
     integer, allocatable :: grouplatdim_st(:) ! staggered latitude dimension id
     integer :: splon, plon

     integer ndims
     integer dim1s,dim2s        ! global size of the first and second horizontal dim.
     integer ncol
     integer ncoldim            ! horizontal dim for non rectangular grids
     integer :: ierr
     type(var_desc_t), pointer :: varid

     integer, allocatable :: mdimids(:)
     !
     integer :: amode

     type(var_desc_t) :: cosp_prs_var      	! COSP ISCCP mean pressure variable id
     type(var_desc_t) :: cosp_prs_bnds      	! COSP ISCCP mean pressure variable bounds
     type(var_desc_t) :: cosp_tau_var      	! COSP ISCCP mean optical depth variable id
     type(var_desc_t) :: cosp_tau_bnds      	! COSP ISCCP mean optical depth varible bounds
     type(var_desc_t) :: cosp_tau_modis_var    ! COSP MODIS mean optical depth variable id
     type(var_desc_t) :: cosp_tau_modis_bnds    ! COSP MODIS mean optical depth variable bounds
     type(var_desc_t) :: cosp_ht_var      	! COSP mean ht variable id
     type(var_desc_t) :: cosp_ht_bnds      	! COSP mean ht variable bounds
     type(var_desc_t) :: cosp_dbze_var         ! COSP mean dbze variable id
     type(var_desc_t) :: cosp_dbze_bnds         ! COSP mean dbze variable bounds
     type(var_desc_t) :: cosp_sr_var   	! COSP mean sr variable id
     type(var_desc_t) :: cosp_sr_bnds   	! COSP mean sr variable bounds
     type(var_desc_t) :: cosp_scol_var         ! COSP mean scol variable id
     type(var_desc_t) :: cosp_htmisr_var   	! COSP mean htmisr variable id
     type(var_desc_t) :: cosp_htmisr_bnds   	! COSP mean htmisr variable bounds
     type(var_desc_t) :: cosp_sza_var   	! COSP mean sza variable id

     type(io_desc_t), pointer :: iodesc

     integer :: lchnk, plat


     if(restart) then
        tape => restarthistory_tape
        if(masterproc) write(iulog,*)'Opening netcdf history restart file ', trim(hrestpath(t))
     else
        tape => history_tape
        if(masterproc) write(iulog,*)'Opening netcdf history file ', trim(nhfil(t))
     end if

     amode = PIO_CLOBBER


     if(restart) then
        call cam_pio_createfile (tape(t)%File, hrestpath(t), amode)
     else
        call cam_pio_createfile (tape(t)%File, nhfil(t), amode)
     end if
     if(is_satfile(t)) then
        ierr = pio_def_dim (tape(t)%File, 'ncol', pio_unlimited, timdim)  

        ierr=pio_def_var (tape(t)%File,'lat',pio_double,(/timdim/),latvar)
        ierr=pio_put_att (tape(t)%File, latvar, 'long_name', 'latitude')
        ierr=pio_put_att (tape(t)%File, latvar, 'units', 'degrees_north')

        ierr=pio_def_var (tape(t)%File,'lon',pio_double,(/timdim/),lonvar)
        ierr=pio_put_att (tape(t)%File, lonvar,'long_name','longitude')
        ierr=pio_put_att (tape(t)%File, lonvar,'units','degrees_east')

     else
        !
        ! Setup netcdf file - create the dimensions of lat,lon,time,level
        !     
        interpolate = (latlon_interpolation(t) .and. .not. restart )
        ncoloutput = ((dycore_is('UNSTRUCTURED') .and. .not. interpolate) &
             .or. associated(tape(t)%column(1)%hmask) &
             .and. .not. restart)

        call get_horiz_grid_dim_d(dim1s,dim2s)

        ! First define the horizontal grid dims
        
        latdim = -1
        londim = -1
        ncoldim = -1

        plat = get_dyn_grid_parm('plat')

        if(dycore_is('UNSTRUCTURED')) then
           if(interpolate) then
              dim1s=get_dyn_grid_parm('nlon')
              dim2s=get_dyn_grid_parm('nlat')

              ret = pio_def_dim (tape(t)%File, 'lat', dim2s, latdim)
              ret = pio_def_dim (tape(t)%File, 'lon', dim1s, londim)            
              
              call add_interp_attributes(tape(t)%File)

           else
              ncol = dim1s*dim2s
              ret = pio_def_dim (tape(t)%File, 'ncol', ncol, ncoldim)
           end if
           ret = pio_put_att(tape(t)%File, PIO_GLOBAL, 'np', get_dyn_grid_parm('np'))
           ret = pio_put_att(tape(t)%File, PIO_GLOBAL, 'ne', get_dyn_grid_parm('ne'))
        else if(ncoloutput) then
           ncol = tape(t)%column(1)%column_cnt
           ret =  pio_def_dim (tape(t)%File, 'ncol', ncol, ncoldim)
        else if(ngroup(t)==0 .or. restart) then
           ret = pio_def_dim (tape(t)%File, 'lat', dim2s, latdim)
           ret = pio_def_dim (tape(t)%File, 'lon', dim1s, londim)
           if(dycore_is('LR')) then
              ret = pio_def_dim (tape(t)%File, 'slat', plat-1, slatdim)
              splon = get_dyn_grid_parm('splon')
              ret = pio_def_dim (tape(t)%File, 'slon', splon, slondim)
           end if
        endif
        

        ! Define the unlimited time dim
        ret = pio_def_dim (tape(t)%File, 'time', pio_unlimited, timdim)
        ret = pio_def_dim (tape(t)%File, 'nbnd', 2, bnddim)
        ret = pio_def_dim (tape(t)%File, 'chars', 8, chardim)

        !
        ! create dimensions for groups of contiguous columns
        ! variables with for single columns are created later in this routine
        !

        if (.not. restart .and. ngroup(t).ne.0) then
           allocate(grouplatdim(ngroup(t)), grouplondim(ngroup(t)))
           allocate(glatvar(ngroup(t)), glonvar(ngroup(t)))
           do i = 1, ngroup(t)

              if(scan(trim(tape(t)%column(i)%lat_name),' ')> 0) then
                 tape(t)%column(i)%lat_name(scan(tape(t)%column(i)%lat_name, ' '):)=' '
              end if
              if(scan(trim(tape(t)%column(i)%lon_name),' ')> 0) then
                 tape(t)%column(i)%lon_name(scan(tape(t)%column(i)%lon_name, ' '):)=' '
              end if

              ret = pio_def_dim (tape(t)%File, trim(tape(t)%column(i)%lat_name), tape(t)%column(i)%num_lats, grouplatdim(i))
              ierr=pio_def_var (tape(t)%File,tape(t)%column(i)%lat_name,pio_double,grouplatdim(i:i),glatvar(i))
              ierr=pio_put_att (tape(t)%File, glatvar(i),'long_name','latitude')
              ierr=pio_put_att (tape(t)%File, glatvar(i),'units','degrees_north')

              ret = pio_def_dim (tape(t)%File, trim(tape(t)%column(i)%lon_name), tape(t)%column(i)%num_lons, grouplondim(i))
              ierr=pio_def_var (tape(t)%File,tape(t)%column(i)%lon_name,pio_double,grouplondim(i:i),glonvar(i))
              ierr=pio_put_att (tape(t)%File, glonvar(i),'long_name','longitude')
              ierr=pio_put_att (tape(t)%File, glonvar(i),'units','degrees_east')
           end do
           if(dycore_is('LR')) then
              !
              ! For staggered (FV) fields, include staggered dimension on h-file
              !
              do f=1,nflds(t)
                 fname_tmp = tape(t)%hlist(f)%field%name
                 if(fname_tmp .eq. 'US' .or. fname_tmp .eq. 'FU_S') then
                    if(.not.allocated(grouplatdim_st)) then
                       allocate(grouplatdim_st(ngroup(t)))
                       allocate(glatvar_st(ngroup(t)))

                       do i = 1, ngroup(t)

                          if(scan(trim(tape(t)%column_st(i)%lat_name),' ')> 0) then
                             tape(t)%column_st(i)%lat_name(scan(tape(t)%column_st(i)%lat_name, ' '):)=' '
                          end if
                          if(scan(trim(tape(t)%column_st(i)%lon_name),' ')> 0) then
                             tape(t)%column_st(i)%lon_name(scan(tape(t)%column_st(i)%lon_name, ' '):)=' '
                          end if

                          ret = pio_def_dim (tape(t)%File, tape(t)%column_st(i)%lat_name, &
                               tape(t)%column_st(i)%num_lats, grouplatdim_st(i))
                          ierr=pio_def_var (tape(t)%File,tape(t)%column_st(i)%lat_name,pio_double, &
                               grouplatdim_st(i:i),glatvar_st(i))
                          ierr=pio_put_att (tape(t)%File, glatvar_st(i),'long_name','staggered latitude')
                          ierr=pio_put_att (tape(t)%File, glatvar_st(i),'units','degrees_north')
                       end do
                    end if
                 else if(fname_tmp .eq. 'VS' .or. fname_tmp .eq. 'FV_S') then
                    if(.not.allocated(grouplondim_st)) then
                       allocate(grouplondim_st(ngroup(t)))
                       allocate(glonvar_st(ngroup(t)))
                       do i = 1, ngroup(t)
                          if(scan(trim(tape(t)%column_st(i)%lat_name),' ')> 0) then
                             tape(t)%column_st(i)%lat_name(scan(tape(t)%column_st(i)%lat_name, ' '):)=' '
                          end if
                          if(scan(trim(tape(t)%column_st(i)%lon_name),' ')> 0) then
                             tape(t)%column_st(i)%lon_name(scan(tape(t)%column_st(i)%lon_name, ' '):)=' '
                          end if
                          ret = pio_def_dim (tape(t)%File, tape(t)%column_st(i)%lon_name, &
                               tape(t)%column_st(i)%num_lons, grouplondim_st(i))
                          ierr=pio_def_var (tape(t)%File,tape(t)%column_st(i)%lon_name,pio_double, &
                               grouplondim_st(i:i),glonvar_st(i))
                          ierr=pio_put_att (tape(t)%File, glonvar_st(i),'long_name','staggered longitude')
                          ierr=pio_put_att (tape(t)%File, glonvar_st(i),'units','degrees_east')
                       end do
                    end if
                 end if
              end do
           endif
        end if
     endif   ! is satfile
     ! Restart requires that the def_dim call for registeredmdims be the last in the history file.



     ! Now define the registered mdims
     allocate(mdimids(registeredmdims))
     do i=1,registeredmdims
        ret = pio_def_dim (tape(t)%File, trim(hist_mdims(i)%name), hist_mdims(i)%value, mdimids(i))        
        ! define variables to label the dimensions, use the same names as the
        ! dimensions

        select case(trim(hist_mdims(i)%name))
        case ('lev')
           levdim = mdimids(i)
           ierr=pio_def_var (tape(t)%File,'lev',pio_double,(/levdim/),levvar)
           str = 'hybrid level at midpoints (1000*(A+B))'
           ierr=pio_put_att (tape(t)%File, levvar, 'long_name', trim(str))
           str = 'level'
           ierr=pio_put_att (tape(t)%File, levvar, 'units', trim(str))
           ierr=pio_put_att (tape(t)%File, levvar, 'positive', 'down')
           ierr=pio_put_att (tape(t)%File, levvar, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
           ierr=pio_put_att (tape(t)%File, levvar, 'formula_terms', 'a: hyam b: hybm p0: P0 ps: PS')

           ierr=pio_def_var (tape(t)%File,'hyam',PIO_DOUBLE,(/levdim/),hyamid)
           str = 'hybrid A coefficient at layer midpoints'
           ierr=pio_put_att (tape(t)%File, hyamid, 'long_name', trim(str))

           ierr=pio_def_var (tape(t)%File,'hybm',PIO_DOUBLE,(/levdim/),hybmid)
           str = 'hybrid B coefficient at layer midpoints'
           ierr=pio_put_att (tape(t)%File, hybmid, 'long_name', trim(str))

        case ('ilev')
           ilevdim = mdimids(i)
           ierr=pio_def_var (tape(t)%File,'ilev',pio_double,(/ilevdim/),ilevvar)
           str = 'hybrid level at interfaces (1000*(A+B))'
           ierr=pio_put_att (tape(t)%File, ilevvar, 'long_name', trim(str))
           str = 'level'
           ierr=pio_put_att (tape(t)%File, ilevvar, 'units', trim(str))
           ierr=pio_put_att (tape(t)%File, ilevvar, 'positive', 'down')
           ierr=pio_put_att (tape(t)%File, ilevvar, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
           ierr=pio_put_att (tape(t)%File, ilevvar, 'formula_terms', 'a: hyai b: hybi p0: P0 ps: PS')

           ierr=pio_def_var (tape(t)%File,'hyai',PIO_DOUBLE,(/ilevdim/),hyaiid)
           str = 'hybrid A coefficient at layer interfaces'
           ierr=pio_put_att (tape(t)%File, hyaiid, 'long_name', trim(str))

           ierr=pio_def_var (tape(t)%File,'hybi',PIO_DOUBLE,(/ilevdim/),hybiid)
           str = 'hybrid B coefficient at layer interfaces'
           ierr=pio_put_att (tape(t)%File, hybiid, 'long_name', trim(str))

	end select
        if(.not. is_satfile(t) .or. is_initfile(t)) then
           select case(trim(hist_mdims(i)%name))
           case('cosp_prs')
              ierr=pio_def_var (tape(t)%File, 'cosp_prs', PIO_DOUBLE,  mdimids(i:i), cosp_prs_var)
              str = 'COSP Mean ISCCP pressure'
              ierr=pio_put_att (tape(t)%File, cosp_prs_var, 'long_name', trim(str))
              str = 'mb'
              ierr=pio_put_att (tape(t)%File, cosp_prs_var, 'units', trim(str))
              str = 'cosp_prs_bnds'
              ierr=pio_put_att (tape(t)%File, cosp_prs_var, 'bounds', trim(str))
              ierr=pio_def_var (tape(t)%File, str, PIO_DOUBLE, (/bnddim,mdimids(i:i)/), cosp_prs_bnds)

           case('cosp_tau')
              ierr=pio_def_var (tape(t)%File, 'cosp_tau', PIO_DOUBLE,  mdimids(i:i), cosp_tau_var)
              str = 'COSP Mean ISCCP optical depth'
              ierr=pio_put_att (tape(t)%File, cosp_tau_var, 'long_name', trim(str))
              str = 'unitless'
              ierr=pio_put_att (tape(t)%File, cosp_tau_var, 'units', trim(str))
              str = 'cosp_tau_bnds'
              ierr=pio_put_att (tape(t)%File, cosp_tau_var, 'bounds', trim(str))
              ierr=pio_def_var (tape(t)%File, str, PIO_DOUBLE, (/bnddim,mdimids(i:i)/), cosp_tau_bnds)
           case('cosp_tau_modis')
              ierr=pio_def_var (tape(t)%File, 'cosp_tau_modis', PIO_DOUBLE,  mdimids(i:i), cosp_tau_modis_var)
              str = 'COSP Mean MODIS optical depth'
              ierr=pio_put_att (tape(t)%File, cosp_tau_modis_var, 'long_name', trim(str))
              str = 'unitless'
              ierr=pio_put_att (tape(t)%File, cosp_tau_modis_var, 'units', trim(str))
              str = 'cosp_tau_modis_bnds'
              ierr=pio_put_att (tape(t)%File, cosp_tau_modis_var, 'bounds', trim(str))
              ierr=pio_def_var (tape(t)%File, str, PIO_DOUBLE, (/bnddim,mdimids(i:i)/), cosp_tau_modis_bnds)
           case('cosp_ht')
              ierr=pio_def_var (tape(t)%File, 'cosp_ht', PIO_DOUBLE, mdimids(i:i), cosp_ht_var)
              str = 'COSP Mean Height for lidar and radar simulator outputs'
              ierr=pio_put_att (tape(t)%File, cosp_ht_var, 'long_name', trim(str))
              str = 'm'
              ierr=pio_put_att (tape(t)%File, cosp_ht_var, 'units', trim(str))
              str = 'cosp_ht_bnds'
              ierr=pio_put_att (tape(t)%File, cosp_ht_var, 'bounds', trim(str))
              ierr=pio_def_var (tape(t)%File, str, PIO_DOUBLE, (/bnddim,mdimids(i:i)/), cosp_ht_bnds)
           case('cosp_dbze')
              ierr=pio_def_var (tape(t)%File, 'cosp_dbze', PIO_DOUBLE,  mdimids(i:i),cosp_dbze_var)
              str = 'COSP Mean dBZe for radar simulator CFAD output'
              ierr=pio_put_att (tape(t)%File, cosp_dbze_var, 'long_name', trim(str))
              str = 'dbze'
              ierr=pio_put_att (tape(t)%File, cosp_dbze_var, 'units', trim(str))
              str = 'cosp_dbze_bnds'
              ierr=pio_put_att (tape(t)%File, cosp_dbze_var, 'bounds', trim(str))
              ierr=pio_def_var (tape(t)%File, str, PIO_DOUBLE, (/bnddim,mdimids(i:i)/), cosp_dbze_bnds)
           case('cosp_sr')
              ierr=pio_def_var (tape(t)%File, 'cosp_sr', PIO_DOUBLE,  mdimids(i:i),cosp_sr_var)
              str = 'COSP Mean Scattering Ratio for lidar simulator CFAD output'
              ierr=pio_put_att (tape(t)%File, cosp_sr_var, 'long_name', trim(str))
              str = '1'
              ierr=pio_put_att (tape(t)%File, cosp_sr_var, 'units', trim(str))
              str = 'cosp_sr_bnds'
              ierr=pio_put_att (tape(t)%File, cosp_sr_var, 'bounds', trim(str))
              ierr=pio_def_var (tape(t)%File, str, PIO_DOUBLE, (/bnddim,mdimids(i:i)/), cosp_sr_bnds)           
           case('cosp_scol')
              ierr=pio_def_var (tape(t)%File, 'cosp_scol', PIO_DOUBLE,  mdimids(i:i),cosp_scol_var)
              str = 'COSP subcolumn'
              ierr=pio_put_att (tape(t)%File, cosp_scol_var, 'long_name', trim(str))
              str = 'number'
              ierr=pio_put_att (tape(t)%File, cosp_scol_var, 'units', trim(str))

           case('cosp_htmisr')
              ierr=pio_def_var (tape(t)%File, 'cosp_htmisr', PIO_DOUBLE,  mdimids(i:i),cosp_htmisr_var)
              str = 'COSP MISR height'
              ierr=pio_put_att (tape(t)%File, cosp_htmisr_var, 'long_name', trim(str))
              str = 'km'
              ierr=pio_put_att (tape(t)%File, cosp_htmisr_var, 'units', trim(str))
              str = 'cosp_htmisr_bnds'
              ierr=pio_put_att (tape(t)%File, cosp_htmisr_var, 'bounds', trim(str))
              ierr=pio_def_var (tape(t)%File, str, PIO_DOUBLE, (/bnddim,mdimids(i:i)/), cosp_htmisr_bnds)
           case('cosp_sza')
              ierr=pio_def_var (tape(t)%File, 'cosp_sza', PIO_DOUBLE,  mdimids(i:i),cosp_sza_var)
              str = 'COSP Parasol SZA'
              ierr=pio_put_att (tape(t)%File, cosp_sza_var, 'long_name', trim(str))
              str = 'degrees'
              ierr=pio_put_att (tape(t)%File, cosp_sza_var, 'units', trim(str))
              !           case('cosp_nbnds')
              !          case('cosp_nhtml')
           case('lev','ilev')   ! just prevent the message for these cases.
           case default
              if(masterproc) write(iulog,*) 'Warning: no dimvar found for registered mdim named: ',hist_mdims(i)%name
           end select
        end if
     end do

     ! There should be no further calls to pio_def_dim after the registeredmdims loop

     ierr=pio_def_var (tape(t)%File,'P0',pio_double,ps0var)
     str = 'reference pressure'
     ierr=pio_put_att (tape(t)%File, ps0var, 'long_name', trim(str))
     ierr=pio_put_att (tape(t)%File, ps0var, 'units', 'Pa')

     call get_ref_date(yr, mon, day, nbsec)
     nbdate = yr*10000 + mon*100 + day
     ierr=pio_def_var (tape(t)%File,'time',pio_double,(/TIMDIM/),tape(t)%timeid)
     ierr=pio_put_att (tape(t)%File, tape(t)%timeid, 'long_name', 'time')
     str = 'days since ' // date2yyyymmdd(nbdate) // ' ' // sec2hms(nbsec)
     ierr=pio_put_att (tape(t)%File, tape(t)%timeid, 'units', trim(str))

     calendar = timemgr_get_calendar_cf()
     ierr=pio_put_att (tape(t)%File, tape(t)%timeid, 'calendar', calendar)

     
     ierr=pio_def_var (tape(t)%File,'date    ',pio_int,(/timdim/),tape(t)%dateid)
     str = 'current date (YYYYMMDD)'
     ierr=pio_put_att (tape(t)%File, tape(t)%dateid, 'long_name', trim(str))


     ierr=pio_def_var (tape(t)%File,'datesec ',pio_int,(/timdim/), tape(t)%datesecid)
     str = 'current seconds of current date'
     ierr=pio_put_att (tape(t)%File, tape(t)%datesecid, 'long_name', trim(str))


     !     
     ! Character header information 
     !
     str = 'CF-1.0'
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'Conventions', trim(str))
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'source', 'CAM')
#if ( defined BFB_CAM_SCAM_IOP )
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'CAM_GENERATED_FORCING','create SCAM IOP dataset')
#endif
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'case',caseid)
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'title',ctitle)
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'logname',logname)
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'host', host)
     ierr= pio_put_att (tape(t)%File, PIO_GLOBAL, 'Version', &
          '$Name$')
     ierr= pio_put_att (tape(t)%File, PIO_GLOBAL, 'revision_Id', &
          '$Id$')
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'initial_file', ncdata)
     ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'topography_file', bnd_topo)

     if(.not. is_satfile(t)) then
        if(ngroup(t)==0 .or. restart) then
           if(ncoloutput) then
              ierr=pio_def_var (tape(t)%File,'lat',pio_double,(/NCOLDIM/),latvar)
              if(dycore_is('UNSTRUCTURED')) then
                 ierr=pio_def_var (tape(t)%File,'area',pio_double,(/NCOLDIM/),areavar)
              end if
           else      
              ierr=pio_def_var (tape(t)%File,'lat',pio_double,(/LATDIM/),latvar)
           end if
           ierr=pio_put_att (tape(t)%File, latvar, 'long_name', 'latitude')
           ierr=pio_put_att (tape(t)%File, latvar, 'units', 'degrees_north')
           
           if(ncoloutput ) then
              ierr=pio_def_var (tape(t)%File,'lon',pio_double,(/NCOLDIM/),lonvar)
           else
              ierr=pio_def_var (tape(t)%File,'lon',pio_double,(/LONDIM/),lonvar)
           end if
           ierr=pio_put_att (tape(t)%File, lonvar,'long_name','longitude')
           ierr=pio_put_att (tape(t)%File, lonvar,'units','degrees_east')

           ! If a staggered grid is in use, output the lat and lon arrays variables.
           if(dycore_is('LR') .and. .not. ncoloutput) then
              
              ierr=pio_def_var (tape(t)%File, 'slat', pio_double,(/SLATDIM/),slatvar)
              ierr=pio_put_att (tape(t)%File, slatvar, 'long_name', 'staggered latitude')
              ierr=pio_put_att (tape(t)%File, slatvar, 'units', 'degrees_north')

              ierr=pio_def_var (tape(t)%File,'slon',pio_double,(/SLONDIM/),slonvar)
              ierr=pio_put_att (tape(t)%File, slonvar, 'long_name', 'staggered longitude')
              ierr=pio_put_att(tape(t)%File, slonvar, 'units', 'degrees_east')

              ierr=pio_def_var (tape(t)%File,'w_stag',pio_double,(/slatdim/),wsid)
              str = 'staggered latitude weights'
              ierr=pio_put_att (tape(t)%File, wsid, 'long_name', trim(str))
           endif
        end if

        ierr=pio_put_att (tape(t)%File, tape(t)%timeid, 'bounds', 'time_bnds')

        ierr=pio_def_var (tape(t)%File,'time_bnds',pio_double,(/bnddim,timdim/),tape(t)%tbndid)
        ierr=pio_put_att (tape(t)%File, tape(t)%tbndid, 'long_name', 'time interval endpoints')
        !
        ! Character
        !
        dimenchar(1) = chardim
        dimenchar(2) = timdim
        ierr=pio_def_var (tape(t)%File,'date_written',PIO_CHAR,dimenchar, tape(t)%date_writtenid)
        ierr=pio_def_var (tape(t)%File,'time_written',PIO_CHAR,dimenchar, tape(t)%time_writtenid)
        !
        ! Integer Header
        !
        ierr=pio_def_var (tape(t)%File,'ntrm',PIO_INT,ntrmid)
        str = 'spectral truncation parameter M'
        ierr=pio_put_att (tape(t)%File, ntrmid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'ntrn',PIO_INT, ntrnid)
        str = 'spectral truncation parameter N'
        ierr=pio_put_att (tape(t)%File, ntrnid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'ntrk',PIO_INT,ntrkid)
        str = 'spectral truncation parameter K'
        ierr=pio_put_att (tape(t)%File, ntrkid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'ndbase',PIO_INT,tape(t)%ndbaseid)
        str = 'base day'
        ierr=pio_put_att (tape(t)%File, tape(t)%ndbaseid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'nsbase',PIO_INT,tape(t)%nsbaseid)
        str = 'seconds of base day'
        ierr=pio_put_att (tape(t)%File, tape(t)%nsbaseid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'nbdate',PIO_INT,tape(t)%nbdateid)
        str = 'base date (YYYYMMDD)'
        ierr=pio_put_att (tape(t)%File, tape(t)%nbdateid, 'long_name', trim(str))

#if ( defined BFB_CAM_SCAM_IOP )
        ierr=pio_def_var (tape(t)%File,'bdate',PIO_INT,tape(t)%bdateid)
        str = 'base date (YYYYMMDD)'
        ierr=pio_put_att (tape(t)%File, tape(t)%bdateid, 'long_name', trim(str))
#endif
        ierr=pio_def_var (tape(t)%File,'nbsec',PIO_INT,tape(t)%nbsecid)
        str = 'seconds of base date'
        ierr=pio_put_att (tape(t)%File, tape(t)%nbsecid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'mdt',PIO_INT,tape(t)%mdtid)
        ierr=pio_put_att (tape(t)%File, tape(t)%mdtid, 'long_name', 'timestep')
        ierr=pio_put_att (tape(t)%File, tape(t)%mdtid, 'units', 's')


        if(latdim > -1) then
           if(.not. is_initfile(file_index=t) .and. .not. dycore_is('UNSTRUCTURED') ) then
              ierr=pio_def_var (tape(t)%File,'nlon',PIO_INT,(/latdim/),tape(t)%nlonid)
              str = 'number of longitudes'
              ierr=pio_put_att (tape(t)%File, tape(t)%nlonid, 'long_name', trim(str))

              ierr=pio_def_var (tape(t)%File,'wnummax',PIO_INT,(/latdim/),tape(t)%wnummaxid)
              str = 'cutoff Fourier wavenumber'
              ierr=pio_put_att (tape(t)%File, tape(t)%wnummaxid, 'long_name', trim(str))
           end if
        !
        ! Floating point time-invariant
        !
           ierr=pio_def_var (tape(t)%File,'gw',PIO_DOUBLE,(/latdim/),gwid)
           str = 'gauss weights'
           ierr=pio_put_att (tape(t)%File, gwid, 'long_name', trim(str))
        end if

        !
        ! Create variables for model timing and header information 
        !

        ierr=pio_def_var (tape(t)%File,'ndcur   ',pio_int,(/timdim/),tape(t)%ndcurid)
        str = 'current day (from base day)'
        ierr=pio_put_att (tape(t)%File, tape(t)%ndcurid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'nscur   ',pio_int,(/timdim/),tape(t)%nscurid)
        str = 'current seconds of current day'
        ierr=pio_put_att (tape(t)%File, tape(t)%nscurid, 'long_name', trim(str))


        if (.not. is_initfile(file_index=t)) then
           ! Don't write the GHG/Solar forcing data to the IC file.
           ierr=pio_def_var (tape(t)%File,'co2vmr  ',pio_double,(/timdim/),tape(t)%co2vmrid)
           str = 'co2 volume mixing ratio'
           ierr=pio_put_att (tape(t)%File, tape(t)%co2vmrid, 'long_name', trim(str))

           ierr=pio_def_var (tape(t)%File,'ch4vmr  ',pio_double,(/timdim/),tape(t)%ch4vmrid)
           str = 'ch4 volume mixing ratio'
           ierr=pio_put_att (tape(t)%File, tape(t)%ch4vmrid, 'long_name', trim(str))

           ierr=pio_def_var (tape(t)%File,'n2ovmr  ',pio_double,(/timdim/),tape(t)%n2ovmrid)
           str = 'n2o volume mixing ratio'
           ierr=pio_put_att (tape(t)%File, tape(t)%n2ovmrid, 'long_name', trim(str))

           ierr=pio_def_var (tape(t)%File,'f11vmr  ',pio_double,(/timdim/),tape(t)%f11vmrid)
           str = 'f11 volume mixing ratio'
           ierr=pio_put_att (tape(t)%File, tape(t)%f11vmrid, 'long_name', trim(str))

           ierr=pio_def_var (tape(t)%File,'f12vmr  ',pio_double,(/timdim/),tape(t)%f12vmrid)
           str = 'f12 volume mixing ratio'
           ierr=pio_put_att (tape(t)%File, tape(t)%f12vmrid, 'long_name', trim(str))

           ierr=pio_def_var (tape(t)%File,'sol_tsi ',pio_double,(/timdim/),tape(t)%sol_tsiid)
           str = 'total solar irradiance'
           ierr=pio_put_att (tape(t)%File, tape(t)%sol_tsiid, 'long_name', trim(str))
           str = 'W/m2'
           ierr=pio_put_att (tape(t)%File, tape(t)%sol_tsiid, 'units', trim(str))

           if (cam_chempkg_is('waccm_mozart').or.cam_chempkg_is('waccm_ghg')) then
              ! solar / geomagetic activity indices...
              ierr=pio_def_var (tape(t)%File,'f107',pio_double,(/timdim/),tape(t)%f107id)
              str = '10.7 cm solar radio flux (F10.7)'
              ierr=pio_put_att (tape(t)%File, tape(t)%f107id, 'long_name', trim(str))
              str = '10^-22 W m^-2 Hz^-1'
              ierr=pio_put_att (tape(t)%File, tape(t)%f107id, 'units', trim(str))

              ierr=pio_def_var (tape(t)%File,'f107a',pio_double,(/timdim/),tape(t)%f107aid)
              str = '81-day centered mean of 10.7 cm solar radio flux (F10.7)'
              ierr=pio_put_att (tape(t)%File, tape(t)%f107aid, 'long_name', trim(str))

              ierr=pio_def_var (tape(t)%File,'kp',pio_double,(/timdim/),tape(t)%kpid)
              str = 'Daily planetary K geomagnetic index'
              ierr=pio_put_att (tape(t)%File, tape(t)%kpid, 'long_name', trim(str))

              ierr=pio_def_var (tape(t)%File,'ap',pio_double,(/timdim/),tape(t)%apid)
              str = 'Daily planetary A geomagnetic index'
              ierr=pio_put_att (tape(t)%File, tape(t)%apid, 'long_name', trim(str))
           endif
        end if


#if ( defined BFB_CAM_SCAM_IOP )
        ierr=pio_def_var (tape(t)%File,'tsec ',pio_int,(/timdim/), tape(t)%tsecid)
        str = 'current seconds of current date needed for scam'
        ierr=pio_put_att (tape(t)%File, tape(t)%tsecid, 'long_name', trim(str))
#endif
        ierr=pio_def_var (tape(t)%File,'nsteph  ',pio_int,(/timdim/),tape(t)%nstephid)
        str = 'current timestep'
        ierr=pio_put_att (tape(t)%File, tape(t)%nstephid, 'long_name', trim(str))
     end if



     !
     ! Create variables and attributes for field list
     !

     do f=1,nflds(t)
        select case (tape(t)%hlist(f)%avgflag)
        case ('A')
           tape(t)%hlist(f)%time_op(:) = 'mean'
        case ('B')
           tape(t)%hlist(f)%time_op(:) = 'mean00z'
        case ('I')
           tape(t)%hlist(f)%time_op(:) = ' '
        case ('X')
           tape(t)%hlist(f)%time_op(:) = 'maximum'
        case ('M')
           tape(t)%hlist(f)%time_op(:) = 'minimum'
        case('L')
           tape(t)%hlist(f)%time_op(:) = LT_DESC
        case default
           call endrun ('h_define: unknown avgflag='//tape(t)%hlist(f)%avgflag)
        end select

        numlev = tape(t)%hlist(f)%field%numlev

        if (tape(t)%hlist(f)%hwrt_prec == 8 .or. restart) then
           ncreal = pio_double
        else
           ncreal = pio_real 
        end if
        if(.not.associated(tape(t)%hlist(f)%varid)) then
           allocate(tape(t)%hlist(f)%varid(1:max(1,ngroup(t))))
        end if
        fname_tmp = strip_suffix(tape(t)%hlist(f)%field%name)

        if(is_satfile(t)) then
           hdims=0
           nfils(t)=1
           call sat_hist_define(tape(t)%File)
        else
           if(ncoloutput) then
              dimindex(1)=ncoldim
              hdims=1
              nacsdims(1) = ncoldim
           else
              nacsdims=(/londim,latdim/)
              dimindex(1)=londim               
              dimindex(2)=latdim               
              hdims=2
              if(dycore_is('LR') .and. numlev>1) then
                 if(fname_tmp .eq. 'US' .or. fname_tmp .eq. 'FU_S') then
                    dimindex(1:2) = (/londim,slatdim/)
                    nacsdims=(/londim,slatdim/)
                 else if(fname_tmp .eq. 'VS' .or. fname_tmp .eq. 'FV_S') then
                    dimindex(1:2) = (/slondim,latdim/)
                    nacsdims=(/slondim,latdim/)
                 else
                    dimindex(1:2) = (/londim,latdim/)
                    nacsdims=(/londim,latdim/)
                 end if
              end if
           end if
        endif

        !
        !  Create variables and atributes for fields written out as columns
        !         
        if(associated(tape(t)%hlist(f)%field%mdims)) then
           mdims => tape(t)%hlist(f)%field%mdims
           mdimsize = size(mdims)
        else if(tape(t)%hlist(f)%field%numlev>1) then
           call endrun('mdims not defined for variable '//trim(tape(t)%hlist(f)%field%name))
        else
           mdimsize=0   
        endif
        do i=0,ngroup(t)
           if(i==0 .and. ngroup(t)>0 .and. .not. restart) cycle
           fname_tmp = strip_suffix(tape(t)%hlist(f)%field%name)
           if(i==0) then  ! i==0 full field output
              if(.not.associated(tape(t)%hlist(f)%varid)) then
                 allocate(tape(t)%hlist(f)%varid(1))
              end if
              varid => tape(t)%hlist(f)%varid(1)
              if(restart) then
                 if (tape(t)%hlist(f)%field%flag_xyfill .or. (avgflag_pertape(t)=='L')) then
                    ierr=pio_def_var(tape(t)%File,trim(fname_tmp)//'_nacs',pio_int,nacsdims(1:hdims),&
                         tape(t)%hlist(f)%nacs_varid)
                 else
                    ierr=pio_def_var(tape(t)%File,trim(fname_tmp)//'_nacs',pio_int,&
                         tape(t)%hlist(f)%nacs_varid)
                 end if
              endif
           else   ! i>0 column group output
              if(restart) exit
              ff=(f-1)*ngroup(t)+i
              dimindex(1)=grouplondim(i)               
              dimindex(2)=grouplatdim(i)               
              hdims=2
              if(dycore_is('LR')) then
                 if(fname_tmp .eq. 'US' .or. fname_tmp .eq. 'FU_S') then
                    dimindex(2)=grouplatdim_st(i)
                 else if(fname_tmp .eq. 'VS' .or. fname_tmp .eq. 'FV_S') then
                    dimindex(1)=grouplondim_st(i)               
                 end if
              endif
              fname_tmp = tape(t)%hlist(f)%field_column_name(i)
              varid => tape(t)%hlist(f)%varid(i)
           endif
           if(mdimsize==0) then
              if(restart) then
                 tdims=hdims
              else
                 dimindex(hdims+1) = timdim
                 tdims=hdims+1
              end if
              ierr=pio_def_var(tape(t)%File,trim(fname_tmp),ncreal,dimindex(1:tdims),varid)
           else 
              do j=1,mdimsize
                 dimindex(hdims+j) = mdimids(mdims(j))
              end do
              if(restart) then
                 tdims=hdims+mdimsize
              else
                 tdims=hdims+mdimsize+1
                 dimindex(tdims)=timdim
              end if
              ierr=pio_def_var(tape(t)%File,trim(fname_tmp),ncreal,dimindex(1:tdims),varid)
              ierr=pio_put_att(tape(t)%FIle,varid,'mdims', mdims(1:mdimsize))
           end if

           str = tape(t)%hlist(f)%field%sampling_seq
           if (len_trim(str)>0) then
              ierr=pio_put_att (tape(t)%File, varid,'Sampling_Sequence', trim(str))
           end if


           if (tape(t)%hlist(f)%field%flag_xyfill) then
              ! Add both _FillValue and missing_value to cover expectations of various applications.
              ! The attribute type must match the data type.
              if (tape(t)%hlist(f)%hwrt_prec == 8) then
                 ierr=pio_put_att (tape(t)%File, varid, '_FillValue', tape(t)%hlist(f)%field%fillvalue)
                 ierr=pio_put_att (tape(t)%File, varid, 'missing_value', tape(t)%hlist(f)%field%fillvalue)
              else
                 ierr=pio_put_att (tape(t)%File, varid, '_FillValue', real(tape(t)%hlist(f)%field%fillvalue,r4))
                 ierr=pio_put_att (tape(t)%File, varid, 'missing_value', real(tape(t)%hlist(f)%field%fillvalue,r4))
              end if
           end if

           str = tape(t)%hlist(f)%field%units
           ierr=pio_put_att (tape(t)%File, varid, 'units', trim(str))

           str = tape(t)%hlist(f)%field%long_name
           ierr=pio_put_att (tape(t)%File, varid, 'long_name', trim(str))
           !
           ! Assign field attributes defining valid levels and averaging info
           !
           str = tape(t)%hlist(f)%time_op
           select case (str)
           case ('mean', 'maximum', 'minimum' )
              ierr=pio_put_att (tape(t)%File, varid,'cell_methods', 'time: '//str)
           end select
           if(i>0) ierr= pio_put_att(tape(t)%File, varid, 'basename',tape(t)%hlist(f)%field%name)
        end do
     end do
     !
     deallocate(mdimids)
     ret = pio_enddef(tape(t)%File)

     if(masterproc) then
        write(iulog,*)'H_DEFINE: Successfully opened netcdf file '
     endif
     !
     ! Write time-invariant portion of history header
     !
     ierr = pio_put_var (tape(t)%File, ps0var, (/ps0/))
     if(.not. is_satfile(t)) then
        londeg => get_dyn_grid_parm_real2d('londeg')
        latdeg => get_dyn_grid_parm_real1d('latdeg')
        if(dycore_is('LR')) then
           londeg_st => get_dyn_grid_parm_real2d('londeg_st')
           latdeg_st => get_dyn_grid_parm_real1d('latdeg_st')
        end if
        if(ncoloutput) then
           if(dycore_is('UNSTRUCTURED')) then
              ! save memory for these global arrays by doing them 1 by 1:
              allocate(harea(ncol))
              call get_horiz_grid_d(ncol, clon_d_out=harea)
              do j=1,ncol
                 harea(j) = harea(j)*radtodeg
              end do
              ierr = pio_put_var (tape(t)%File, lonvar, harea)

              call get_horiz_grid_d(ncol, area_d_out=harea)
              ierr = pio_put_var (tape(t)%File, areavar, harea)

              call get_horiz_grid_d(ncol, clat_d_out=harea)
              do j=1,ncol
                 harea(j) = harea(j)*radtodeg
              end do
              ierr = pio_put_var (tape(t)%File, latvar, harea)
              deallocate(harea)
           else
              allocate(latdeg2d(pcols, begchunk:endchunk),londeg(pcols,begchunk:endchunk))
              latdeg2d=0_r8
              londeg=0_r8
              do lchnk=begchunk,endchunk
                 call get_rlat_all_p(lchnk,pcols,latdeg2d(:,lchnk))
                 call get_rlon_all_p(lchnk,pcols,londeg(:,lchnk))
              end do
              latdeg2d = latdeg2d*radtodeg
              londeg = londeg*radtodeg
              call get_phys_decomp(iodesc, 1,(/1/),1,pio_double,column=tape(t)%column(1))
              call pio_write_darray(tape(t)%File, lonvar, iodesc, londeg, ierr)
              call pio_write_darray(tape(t)%File, latvar, iodesc, latdeg2d, ierr)
              deallocate(latdeg2d,londeg)
           end if
        else if(interpolate) then
           ierr = pio_put_var (tape(t)%File, latvar, get_interp_lat())
           ierr = pio_put_var (tape(t)%File, lonvar, get_interp_lon())
        else if(ngroup(t)==0 .or. restart) then

           ierr=pio_put_var (tape(t)%File, latvar, latdeg)
           plon = get_dyn_grid_parm('plon')
           allocate(alon(plon))
           do i=1,plon
              alon(i) = (i-1) * 360.0_r8 / plon
           end do

           if (single_column) alon(1)=scmlon
           ierr = pio_put_var (tape(t)%File, lonvar, alon)
           deallocate(alon)
           if(dycore_is('LR')) then
              latdeg_st => get_dyn_grid_parm_real1d('latdeg_st')
              
              ierr = pio_put_var (tape(t)%File, slatvar, latdeg_st)
              
              londeg_st => get_dyn_grid_parm_real2d('londeg_st')
              ierr = pio_put_var (tape(t)%File, slonvar, londeg_st(:,1))
              
              w_staggered => get_dyn_grid_parm_real1d('w_staggered')
              ierr = pio_put_var (tape(t)%File, wsid, w_staggered)
           endif
        end if
        !
        ! write out coordinate variables for columns
        !
        if (.not. restart .and. ngroup(t)>0) then
           do i = 1, ngroup(t)
              ierr = pio_put_var (tape(t)%File, glonvar(i), &
                   londeg(tape(t)%column(i)%columnlon(1):tape(t)%column(i)%columnlon(2),1))
              ierr = pio_put_var (tape(t)%File, glatvar(i), &
                   latdeg(tape(t)%column(i)%columnlat(1):tape(t)%column(i)%columnlat(2)))
           end do
           deallocate(glonvar, glatvar,grouplatdim,grouplondim)
           !
           ! For staggered (FV) fields, write out staggered coordinate variables for columns
           !
           do f=1,nflds(t)
              select case (tape(t)%hlist(f)%field%name)
              case ('US','FU_S')
                 do i = 1, ngroup(t)
                    ierr = pio_put_var(tape(t)%File, glatvar_st(i), &
                         latdeg_st(tape(t)%column_st(i)%columnlat(1):tape(t)%column_st(i)%columnlat(2)))
                 end do
              case ('VS','FV_S')
                 do i = 1, ngroup(t)
                    ierr = pio_put_var(tape(t)%File, glonvar_st(i), &
                         londeg_st(tape(t)%column_st(i)%columnlon(1):tape(t)%column_st(i)%columnlon(2),1))
                 end do
              end select
           end do
        end if


        if(latdim>-1) then
           w => get_dyn_grid_parm_real1d('w')
           ierr = pio_put_var (tape(t)%File, gwid, w)
        end if
        ierr = pio_put_var (tape(t)%File, ntrmid, (/ptrm/))
        ierr = pio_put_var (tape(t)%File, ntrnid, (/ptrn/))
        ierr = pio_put_var (tape(t)%File, ntrkid, (/ptrk/))
        dtime = get_step_size()
        ierr = pio_put_var (tape(t)%File, tape(t)%mdtid, (/dtime/))
        !
        ! Model date info
        !
        ierr = pio_put_var (tape(t)%File, tape(t)%ndbaseid, (/ndbase/))
        ierr = pio_put_var (tape(t)%File, tape(t)%nsbaseid, (/nsbase/))
        
        ierr = pio_put_var (tape(t)%File, tape(t)%nbdateid, (/nbdate/))
#if ( defined BFB_CAM_SCAM_IOP )
        ierr = pio_put_var (tape(t)%File, tape(t)%bdateid, (/nbdate/))
#endif
        ierr = pio_put_var (tape(t)%File, tape(t)%nbsecid, (/nbsec/))
        !
        ! Reduced grid info
        !

        if(.not. is_initfile(file_index=t) .and. .not. dycore_is('UNSTRUCTURED') .and. latdim>-1) then
           ierr = pio_put_var (tape(t)%File, tape(t)%nlonid, nlon)
           ierr = pio_put_var (tape(t)%File, tape(t)%wnummaxid, wnummax)
        end if
     end if
     if(allocated(grouplondim_st)) then
        deallocate(glonvar_st, grouplondim_st)
     end if
     if(allocated(grouplatdim_st)) then
        deallocate(glatvar_st, grouplatdim_st)
     end if
     !
     ! 0.01 converts Pascals to millibars
     !
     alev(:pver) = 0.01_r8*ps0*(hyam(:pver) + hybm(:pver))
     ailev(:pverp) = 0.01_r8*ps0*(hyai(:pverp) + hybi(:pverp))


     do i=1,registeredmdims
        select case(hist_mdims(i)%name)
        case ('lev')
           ierr = pio_put_var (tape(t)%File, levvar, alev)
           ierr = pio_put_var (tape(t)%File, hyamid, hyam)
           ierr = pio_put_var (tape(t)%File, hybmid, hybm)
        case ('ilev')
           ierr = pio_put_var (tape(t)%File, ilevvar, ailev)
           ierr = pio_put_var (tape(t)%File, hyaiid, hyai)
           ierr = pio_put_var (tape(t)%File, hybiid, hybi)
        end select
        if(.not. (is_satfile(t) .or. is_initfile(t))) then
           select case(hist_mdims(i)%name)
           case('cosp_prs')
              ierr = pio_put_var (tape(t)%file, cosp_prs_var, prsmid_cosp)
              ierr = pio_put_var (tape(t)%file, cosp_prs_bnds, prslim_cosp)
           case('cosp_tau')
              ierr = pio_put_var (tape(t)%file, cosp_tau_var, taumid_cosp)
              ierr = pio_put_var (tape(t)%file, cosp_tau_bnds, taulim_cosp)
           case('cosp_tau_modis')
              ierr = pio_put_var (tape(t)%file, cosp_tau_modis_var, taumid_cosp_modis)
              ierr = pio_put_var (tape(t)%file, cosp_tau_modis_bnds, taulim_cosp_modis)
           case('cosp_ht')
              ierr = pio_put_var (tape(t)%file, cosp_ht_var, htmid_cosp)
              ierr = pio_put_var (tape(t)%file, cosp_ht_bnds, htlim_cosp)
           case('cosp_dbze')
              ierr = pio_put_var (tape(t)%file, cosp_dbze_var, dbzemid_cosp)
              ierr = pio_put_var (tape(t)%file, cosp_dbze_bnds, dbzelim_cosp)
           case('cosp_sr')
              ierr = pio_put_var (tape(t)%file, cosp_sr_var, srmid_cosp)
              ierr = pio_put_var (tape(t)%file, cosp_sr_bnds, srlim_cosp)
           case('cosp_scol')
              ierr = pio_put_var (tape(t)%file, cosp_scol_var, scol_cosp)
           case('cosp_htmisr')              
              ierr = pio_put_var (tape(t)%file, cosp_htmisr_var, htmisrmid_cosp)
              ierr = pio_put_var (tape(t)%file, cosp_htmisr_bnds, htmisrlim_cosp)
           case('cosp_sza')
              ierr = pio_put_var (tape(t)%file, cosp_sza_var, sza_cosp)
           case('lev','ilev')  ! avoids the warning message
           case default
              if(masterproc) write(iulog,*) 'Warning: no dimvar found for registered mdim named: ',hist_mdims(i)%name
           end select
        end if
     end do




   end subroutine h_define

!#######################################################################


   subroutine h_normalize (f, t)

   use dycore, only: dycore_is

!
!----------------------------------------------------------------------- 
! 
! Purpose: Normalize fields on a history file by the number of accumulations
! 
! Method: Loop over fields on the tape.  Need averaging flag and number of
!         accumulations to perform normalization.
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer, intent(in) :: f       ! field index
      integer, intent(in) :: t       ! tape index
!
! Local workspace
!
      integer begdim2                ! on-node vert start index
      integer enddim2                ! on-node vert end index
      integer c                      ! chunk (or lat) index
      integer endi                   ! terminating column index
      integer begdim1                ! on-node dim1 start index
      integer enddim1                ! on-node dim1 end index  
      integer begdim3                ! on-node chunk or lat start index
      integer enddim3                ! on-node chunk or lat end index  
      integer :: ib, ie    ! beginning and ending indices of first dimension
      integer :: jb, je    ! beginning and ending indices of second dimension
      integer k

      logical :: flag_xyfill         ! non-applicable xy points flagged with fillvalue
      character*1 avgflag            ! averaging flag

      call t_startf ('h_normalize')

      begdim2 = tape(t)%hlist(f)%field%begdim2
      enddim2 = tape(t)%hlist(f)%field%enddim2
      avgflag = tape(t)%hlist(f)%avgflag
!
! normalize by number of accumulations for averaged case
!
      begdim1 = tape(t)%hlist(f)%field%begdim1
      enddim1 = tape(t)%hlist(f)%field%enddim1
      begdim3 = tape(t)%hlist(f)%field%begdim3
      enddim3 = tape(t)%hlist(f)%field%enddim3
      flag_xyfill = tape(t)%hlist(f)%field%flag_xyfill
      do c=begdim3,enddim3
 ! Some slices might not be entirely full  -- use endi instead of enddim2
         if(associated(tape(t)%hlist(f)%field%colperdim3)) then         
            endi    = tape(t)%hlist(f)%field%colperdim3(c)
         else
            endi = enddim1-begdim1+1 
         end if

         ib = begdim1
         ie = begdim1+endi-1
         jb = begdim2
         je = enddim2
         
         if (flag_xyfill) then
               do k=jb,je
                  where (tape(t)%hlist(f)%nacs(ib:ie,c) == 0)
                     tape(t)%hlist(f)%hbuf(ib:ie,k,c) = tape(t)%hlist(f)%field%fillvalue
                  endwhere
               end do
         end if

         if (avgflag == 'A' .or. avgflag == 'B' .or. avgflag == 'L') then
            if(flag_xyfill .or. avgflag == 'L') then
               do k=jb,je
                  where (tape(t)%hlist(f)%nacs(ib:ie,c) /= 0)
                     tape(t)%hlist(f)%hbuf(ib:ie,k,c) = &
                          tape(t)%hlist(f)%hbuf(ib:ie,k,c) &
                          / tape(t)%hlist(f)%nacs(ib:ie,c)
                  endwhere
               end do
            else if(tape(t)%hlist(f)%nacs(1,c) > 0) then
               do k=jb,je
                  tape(t)%hlist(f)%hbuf(ib:ie,k,c) = &
                       tape(t)%hlist(f)%hbuf(ib:ie,k,c) &
                       / tape(t)%hlist(f)%nacs(1,c)
               end do
            end if
         end if
      end do

      call t_stopf ('h_normalize')
      
      return
   end subroutine h_normalize

!#######################################################################

   subroutine h_zero (f, t)
   use dycore, only: dycore_is
!
!----------------------------------------------------------------------- 
! 
! Purpose: Zero out accumulation buffers for a tape
! 
! Method: Loop through fields on the tape
! 
!-----------------------------------------------------------------------
!
      integer, intent(in) :: f     ! field index
      integer, intent(in) :: t     ! tape index
!
! Local workspace
!
      integer begdim1              ! on-node 1st dim start index
      integer enddim1              ! on-node 1st dim end index
      integer begdim2              ! on-node vert start index
      integer enddim2              ! on-node vert end index
      integer c                    ! chunk index
      integer endi                 ! terminating column index
      integer begdim3              ! on-node chunk or lat start index
      integer enddim3              ! on-node chunk or lat end index  
      type (dim_index_3d) :: dimind ! 3-D dimension index
      
      call t_startf ('h_zero')

      begdim1 = tape(t)%hlist(f)%field%begdim1
      enddim1 = tape(t)%hlist(f)%field%enddim1
      begdim2 = tape(t)%hlist(f)%field%begdim2
      enddim2 = tape(t)%hlist(f)%field%enddim2
      begdim3 = tape(t)%hlist(f)%field%begdim3
      enddim3 = tape(t)%hlist(f)%field%enddim3

      if(associated(tape(t)%hlist(f)%field%colperdim3)) then         
         do c=begdim3,enddim3
            endi    = tape(t)%hlist(f)%field%colperdim3(c)
            tape(t)%hlist(f)%hbuf(begdim1:begdim1+endi-1,begdim2:enddim2,c)=0._r8
         end do
         tape(t)%hlist(f)%nacs(:,:) = 0
      else
         endi = tape(t)%hlist(f)%field%enddim1 - tape(t)%hlist(f)%field%begdim1 + 1
         tape(t)%hlist(f)%hbuf(begdim1:begdim1+endi-1,begdim2:enddim2,begdim3:enddim3)=0._r8
         tape(t)%hlist(f)%nacs(:,:) = 0
      end if

      call t_stopf ('h_zero')

      return
   end subroutine h_zero

!#######################################################################

   subroutine dump_field (f, t, restart)
     use interp_mod,    only: latlon_interpolation
     use cam_pio_utils, only: get_decomp, pio_subsystem     
     use dycore,        only: dycore_is
     use interp_mod, only : write_interpolated
     use spmd_utils, only : iam
     integer, intent(in) :: f, t
     logical, intent(in) :: restart
!
!----------------------------------------------------------------------- 
! 
! Purpose: Write a variable to a history tape
! 
! Method: If SPMD, first gather the data to the master processor.  
!         Next, transpose the data to COORDS order (the default).
!         Finally, issue the netcdf call to write the variable
! 
!-----------------------------------------------------------------------
     integer :: bsize, ierr
     
     real(r4), allocatable :: localbuf4(:)
     type(io_desc_t), pointer :: int_iodesc, iodesc
     type(var_desc_t), pointer :: varid, varid_vecv
     character(len=max_fieldname_len) :: fname
     type(column_info) :: tmpcolumn
     integer :: ff

     interpolate = (latlon_interpolation(t) .and. .not. restart )
     do i=1,max(ngroup(t),1)
        if(restart .and. i>1) exit
        varid=>tape(t)%hlist(f)%varid(i)

        if(restart) then
           call pio_setframe(varid, int(1,kind=PIO_Offset))
           if(.not. interpolate) then
              call get_decomp(iodesc, tape(t)%hlist(f)%field, pio_double)
           end if
        else
           call pio_setframe(varid, int(max(1,nfils(t)),kind=PIO_Offset))
           if(ngroup(t)>0) then
              tmpcolumn = tape(t)%column(i)
              nullify(tmpcolumn%hmask)
              if(dycore_is('LR')) then
                 fname = tape(t)%hlist(f)%field%name
                 if(fname .eq. 'US' .or. fname .eq. 'FU_S') then
                    tmpcolumn%columnlat(1) = tape(t)%column_st(i)%columnlat(1)+1
                    ! Adding 1 because staggered lat indexing for "xzybuf" starts at 2
                    tmpcolumn%columnlat(2) = tape(t)%column_st(i)%columnlat(2)+1
                    tmpcolumn%num_lats =  tape(t)%column_st(i)%num_lats
                 else if(fname .eq. 'VS' .or. fname .eq. 'FV_S') then
                    tmpcolumn%columnlon = tape(t)%column_st(i)%columnlon
                    tmpcolumn%num_lons =  tape(t)%column_st(i)%num_lons
                 end if
              end if
              allocate(iodesc)  ! no lookup table for columns
              call get_decomp(iodesc, tape(t)%hlist(f)%field, pio_double,column=tmpcolumn)

#ifdef HDEBUG
              write(iulog,*)'DUMP_FIELD:',__LINE__,' writing column time indx ', &
                   nfils(t), ' field id', tape(t)%hlist(f)%varid(i)%varid, &
                   ' name ', trim(tape(t)%hlist(f)%field_column_name(i)), &
                   'numlons=',tmpcolumn%num_lons,'numlats=', tmpcolumn%num_lats, ' base name ', &
                   trim(tape(t)%hlist(f)%field%name),' columnlon=', &
                   tmpcolumn%columnlon, ' columnlat=', tmpcolumn%columnlat
#endif
           else if(associated(tape(t)%column(1)%hmask)) then
              call get_decomp(iodesc, tape(t)%hlist(f)%field, pio_double, column=tape(t)%column(1))
           else if(.not. interpolate) then
              call get_decomp(iodesc, tape(t)%hlist(f)%field, pio_double)
           endif
        end if


        bsize=size(tape(t)%hlist(f)%hbuf)
        if(interpolate) then
           ff = tape(t)%hlist(f)%field%vector_compliment
           ! ff < 0 is the vector compliment of ff*(-1)
           if(ff > 0) then
              varid_vecv=>tape(t)%hlist(ff)%varid(i)
              call pio_setframe(varid_vecv, int(max(1,nfils(t)),kind=PIO_Offset))
              call write_interpolated(tape(t)%File, varid, varid_vecv, tape(t)%hlist(f)%hbuf,&
                   tape(t)%hlist(ff)%hbuf, tape(t)%hlist(f)%field%numlev,  pio_double, tape(t)%hlist(f)%field%decomp_type)
           else if(ff==0) then
              call write_interpolated(tape(t)%File, varid, tape(t)%hlist(f)%hbuf,&
                   tape(t)%hlist(f)%field%numlev,  pio_double, tape(t)%hlist(f)%field%decomp_type)
           endif
        else
           if(tape(t)%hlist(f)%field%numlev>1 .and. .not. associated(tape(t)%hlist(f)%field%mdims)) then
              call endrun('mdims not defined for variable '//trim(tape(t)%hlist(f)%field%name))
           end if
           call pio_write_darray(tape(t)%File, varid, iodesc, tape(t)%hlist(f)%hbuf, ierr)
        endif

        if(restart) then
           if(tape(t)%hlist(f)%field%flag_xyfill .or. (avgflag_pertape(t) == 'L')) then
              call get_decomp(int_iodesc, tape(t)%hlist(f)%field, pio_int, 1)
              call pio_write_darray(tape(t)%File, tape(t)%hlist(f)%nacs_varid, &
                   int_iodesc, reshape(tape(t)%hlist(f)%nacs, (/size(tape(t)%hlist(f)%nacs)/)), ierr)

           else
              ierr = pio_put_var(tape(t)%File, tape(t)%hlist(f)%nacs_varid,&
                   tape(t)%hlist(f)%nacs(:, tape(t)%hlist(f)%field%begdim3:tape(t)%hlist(f)%field%enddim3))
           end if
        end if
        if(.not. restart .and. ngroup(t)>0) then
           call pio_freedecomp(pio_subsystem, iodesc)
           deallocate(iodesc) ! no lookup table for columns
        end if

     end do
     
     return
   end subroutine dump_field

!#######################################################################

   logical function write_inithist ()
!
!-----------------------------------------------------------------------
! 
! Purpose: Set flags that will initiate dump to IC file when OUTFLD and
! WSHIST are called
! 
!-----------------------------------------------------------------------
!
      use time_manager, only: get_nstep, get_curr_date, get_step_size, is_last_step
!
! Local workspace
!
      integer :: yr, mon, day      ! year, month, and day components of
                                   ! a date
      integer :: nstep             ! current timestep number
      integer :: ncsec             ! current time of day [seconds]
      integer :: dtime             ! timestep size

!-----------------------------------------------------------------------

      write_inithist  = .false.

      if(is_initfile()) then

         nstep = get_nstep()
         call get_curr_date(yr, mon, day, ncsec)

         if    (inithist == '6-HOURLY') then
            dtime  = get_step_size()
            write_inithist = nstep /= 0 .and. mod( nstep, nint((6._r8*3600._r8)/dtime) ) == 0
         elseif(inithist == 'DAILY'   ) then
            write_inithist = nstep /= 0 .and. ncsec == 0
         elseif(inithist == 'MONTHLY' ) then
            write_inithist = nstep /= 0 .and. ncsec == 0 .and. day == 1
         elseif(inithist == 'YEARLY'  ) then
            write_inithist = nstep /= 0 .and. ncsec == 0 .and. day == 1 .and. mon == 1
         elseif(inithist == 'CAMIOP'  ) then
            write_inithist = nstep == 0 
         elseif(inithist == 'ENDOFRUN'  ) then
            write_inithist = nstep /= 0 .and. is_last_step()
         end if
      end if

      return
   end function write_inithist

!#######################################################################

   subroutine wshist (rgnht_in)
     !
     !----------------------------------------------------------------------- 
     ! 
     ! Purpose: Driver routine to write fields on history tape t
     ! 
     ! Method: For variables which do not need to be gathered (SPMD) just issue the netcdf call
     !         For those that do need to be gathered, call "dump_field" to do the operation.
     !         Finally, zero the history buffers for each field written.
     ! 
     ! Author: CCM Core Group
     ! 
     !-----------------------------------------------------------------------
     use dycore, only: dycore_is

     use time_manager,  only: get_nstep, get_curr_date, get_curr_time,get_step_size
     use chem_surfvals, only: chem_surfvals_get, chem_surfvals_co2_rad
     use solar_data,    only: sol_tsi
     use sat_hist,      only: sat_hist_write
     use mo_solar_parms,only: get_solar_parms

     logical, intent(in), optional :: rgnht_in(ptapes)
     !
     ! Local workspace
     !
     character(len=8) :: cdate  ! system date
     character(len=8) :: ctime  ! system time

     logical  :: rgnht(ptapes), restart
     integer t, f               ! tape, field indices
     integer ret                ! return value from netcdf call
     integer start              ! starting index required by nf_put_vara
     integer count1             ! count values required by nf_put_vara
     integer startc(2)          ! start values required by nf_put_vara (character)
     integer countc(2)          ! count values required by nf_put_vara (character)
#ifdef HDEBUG
     !      integer begdim3
     !      integer enddim3
#endif

     integer :: yr, mon, day      ! year, month, and day components of a date
     integer :: nstep             ! current timestep number
     integer :: ncdate            ! current date in integer format [yyyymmdd]
     integer :: ncsec             ! current time of day [seconds]
     integer :: ndcur             ! day component of current time
     integer :: nscur             ! seconds component of current time
     real(r8) :: time             ! current time
     real(r8) :: tdata(2)         ! time interval boundaries
     character(len=max_string_len) :: fname ! Filename
     logical :: prev              ! Label file with previous date rather than current
     integer :: ierr
     integer :: bsize
     integer :: dim1s,dim2s,ncol
#if ( defined BFB_CAM_SCAM_IOP )
     integer :: tsec             ! day component of current time
     integer :: dtime            ! seconds component of current time
#endif
     real(r8) :: kp, ap, f107, f107a

     !-----------------------------------------------------------------------

     if (cam_chempkg_is('waccm_mozart').or.cam_chempkg_is('waccm_ghg')) then
        ! get solar/geomagnetic activity data...
        call get_solar_parms( f107_s=f107, f107a_s=f107a, ap_s=ap, kp_s=kp )
     endif

     if(present(rgnht_in)) then
        rgnht=rgnht_in
        restart=.true.
        tape => restarthistory_tape
     else
        rgnht=.false.
        restart=.false.
        tape => history_tape
     end if

     nstep = get_nstep()
     call get_curr_date(yr, mon, day, ncsec)
     ncdate = yr*10000 + mon*100 + day
     call get_curr_time(ndcur, nscur)
     !
     ! Write time-varying portion of history file header
     !
     do t=1,ptapes
        if (nflds(t) == 0 .or. (restart .and.(.not.rgnht(t)))) cycle
        !
        ! Check if this is the IC file and if it's time to write.
        ! Else, use "nhtfrq" to determine if it's time to write
        ! the other history files.
        !
        if((.not. restart) .or. rgnht(t)) then
           if( is_initfile(file_index=t) ) then
              hstwr(t) =  write_inithist()
              prev     = .false.
           else 
              if (nhtfrq(t) == 0) then
                 hstwr(t) = nstep /= 0 .and. day == 1 .and. ncsec == 0
                 prev     = .true.
              else
                 hstwr(t) = mod(nstep,nhtfrq(t)) == 0
                 prev     = .false.
              end if
           end if
        end if
        if (hstwr(t) .or. (restart .and. rgnht(t))) then
           if(masterproc) then
              if(is_initfile(file_index=t)) then
                 write(iulog,100) yr,mon,day,ncsec
100              format('WSHIST: writing time sample to Initial Conditions h-file', &
                      ' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
              else if(is_satfile(t)) then
                 write(iulog,150) nfils(t),t,yr,mon,day,ncsec
150              format('WSHIST: writing sat columns ',i6,' to h-file ', &
                      i1,' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
              else if(hstwr(t)) then
                 write(iulog,200) nfils(t),t,yr,mon,day,ncsec
200              format('WSHIST: writing time sample ',i3,' to h-file ', &
                      i1,' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
              else if(restart .and. rgnht(t)) then
                 write(iulog,300) nfils(t),t,yr,mon,day,ncsec
300              format('WSHIST: writing history restart ',i3,' to hr-file ', &
                      i1,' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
              end if
              write(iulog,*)
           end if
           !
           ! Starting a new volume => define the metadata
           !
           if (nfils(t)==0 .or. (restart.and.rgnht(t))) then
              if(restart) then
                 rhfilename_spec = '%c.cam' // trim(inst_suffix) // '.rh%t.%y-%m-%d-%s.nc' 
                 fname = interpret_filename_spec( rhfilename_spec, number=(t-1))
                 hrestpath(t)=fname
              else if(is_initfile(file_index=t)) then
                 fname = interpret_filename_spec( hfilename_spec(t) )
              else
                 fname = interpret_filename_spec( hfilename_spec(t), number=(t-1), &
                      prev=prev )
              end if
              !
              ! Check that this new filename isn't the same as a previous or current filename
              !
              do f = 1, ptapes
                 if (masterproc.and. trim(fname) == trim(nhfil(f)) )then
                    write(iulog,*)'WSHIST: New filename same as old file = ', trim(fname)
                    write(iulog,*)'Is there an error in your filename specifiers?'
                    write(iulog,*)'hfilename_spec(', t, ') = ', hfilename_spec(t)
                    if ( t /= f )then
                       write(iulog,*)'hfilename_spec(', f, ') = ', hfilename_spec(f)
                    end if
                    call endrun
                 end if
              end do
              if(.not. restart) then
                 nhfil(t) = fname
                 if(masterproc) write(iulog,*)'WSHIST: nhfil(',t,')=',trim(nhfil(t))
                 cpath(t) = nhfil(t)
                 if ( len_trim(nfpath(t)) == 0 ) nfpath(t) = cpath(t)
              end if
              call h_define (t, restart)
           end if

           if(is_satfile(t)) then
              call sat_hist_write( tape(t), nflds(t), nfils(t))
           else
              if(restart) then
                 start=1
              else
                 nfils(t) = nfils(t) + 1
                 start = nfils(t)
              end if
              count1 = 1


              ierr = pio_put_var (tape(t)%File, tape(t)%ndcurid,(/start/), (/count1/),(/ndcur/))
              ierr = pio_put_var (tape(t)%File, tape(t)%nscurid,(/start/), (/count1/),(/nscur/))
              ierr = pio_put_var (tape(t)%File, tape(t)%dateid,(/start/), (/count1/),(/ncdate/))

              if (.not. is_initfile(file_index=t)) then
                 ! Don't write the GHG/Solar forcing data to the IC file.
                 ierr=pio_put_var (tape(t)%File, tape(t)%co2vmrid,(/start/), (/count1/),(/chem_surfvals_co2_rad(vmr_in=.true.)/))
                 ierr=pio_put_var (tape(t)%File, tape(t)%ch4vmrid,(/start/), (/count1/),(/chem_surfvals_get('CH4VMR')/))
                 ierr=pio_put_var (tape(t)%File, tape(t)%n2ovmrid,(/start/), (/count1/),(/chem_surfvals_get('N2OVMR')/))
                 ierr=pio_put_var (tape(t)%File, tape(t)%f11vmrid,(/start/), (/count1/),(/chem_surfvals_get('F11VMR')/))
                 ierr=pio_put_var (tape(t)%File, tape(t)%f12vmrid,(/start/), (/count1/),(/chem_surfvals_get('F12VMR')/))
                 ierr=pio_put_var (tape(t)%File, tape(t)%sol_tsiid,(/start/), (/count1/),(/sol_tsi/))

                 if (cam_chempkg_is('waccm_mozart').or.cam_chempkg_is('waccm_ghg')) then
                    ierr=pio_put_var (tape(t)%File, tape(t)%f107id, (/start/), (/count1/),(/ f107 /) )
                    ierr=pio_put_var (tape(t)%File, tape(t)%f107aid,(/start/), (/count1/),(/ f107a /) )
                    ierr=pio_put_var (tape(t)%File, tape(t)%kpid,   (/start/), (/count1/),(/ kp /) )
                    ierr=pio_put_var (tape(t)%File, tape(t)%kpid,   (/start/), (/count1/),(/ ap /) )
                 endif

              end if



              ierr = pio_put_var (tape(t)%File, tape(t)%datesecid,(/start/),(/count1/),(/ncsec/))
#if ( defined BFB_CAM_SCAM_IOP )
              dtime = get_step_size()
              tsec=dtime*nstep
              ierr = pio_put_var (tape(t)%File, tape(t)%tsecid,(/start/),(/count1/),(/tsec/))
#endif
              ierr = pio_put_var (tape(t)%File, tape(t)%nstephid,(/start/),(/count1/),(/nstep/))
              time = ndcur + nscur/86400._r8
              ierr=pio_put_var (tape(t)%File, tape(t)%timeid, (/start/),(/count1/),(/time/))

              startc(1) = 1
              startc(2) = start
              countc(1) = 2
              countc(2) = 1
              if (is_initfile(file_index=t)) then
                 tdata = time   ! Inithist file is always instantanious data
              else
                 tdata(1) = beg_time(t)
                 tdata(2) = time
              end if
              ierr=pio_put_var (tape(t)%File, tape(t)%tbndid, startc, countc, tdata)
              if(.not.restart) beg_time(t) = time  ! update beginning time of next interval
              startc(1) = 1
              startc(2) = start
              countc(1) = 8
              countc(2) = 1
              call datetime (cdate, ctime)
              ierr = pio_put_var (tape(t)%File, tape(t)%date_writtenid, startc, countc, (/cdate/))
              ierr = pio_put_var (tape(t)%File, tape(t)%time_writtenid, startc, countc, (/ctime/))

              if(.not. restart) then
                 !$OMP PARALLEL DO PRIVATE (F)
                 do f=1,nflds(t)
                    ! Normalized averaged fields
                    if (tape(t)%hlist(f)%avgflag /= 'I') then
                       call h_normalize (f, t)
                    end if
                 end do
              end if
              !
              ! Write field to history tape.  Note that this is NOT threaded due to netcdf limitations
              !
              call t_startf ('dump_field')
              do f=1,nflds(t)
                 call dump_field(f, t, restart)
              end do
              call t_stopf ('dump_field')
              !
              ! Zero history buffers and accumulators now that the fields have been written.
              !



              if(restart) then
                 do f=1,nflds(t)
                    if(associated(tape(t)%hlist(f)%varid)) then
                       deallocate(tape(t)%hlist(f)%varid)
                       nullify(tape(t)%hlist(f)%varid)
                    end if
                 end do
                 call pio_closefile(tape(t)%File)
              else
                 !$OMP PARALLEL DO PRIVATE (F)
                 do f=1,nflds(t)
                    call h_zero (f, t)
                 end do
              end if
           end if
        end if

     end do

     return
   end subroutine wshist

!#######################################################################

   subroutine addvar (ncid, name, xtype , ndims , dimids, vid)

!
!----------------------------------------------------------------------- 
! 
! Purpose: Issue the netcdf call to add a variable to the dataset
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
      type(file_desc_t), intent(inout) :: ncid       ! netcdf file id
      integer, intent(in) :: xtype      ! netcdf type flag
      integer, intent(in) :: ndims      ! number of dimensions
      integer, intent(in) :: dimids(:)  ! dimension ids
      type(var_desc_t), intent(out) :: vid        ! variable ids
      
      character(len=*), intent(in) :: name ! variable name
      integer :: ierr
!
! Local workspace
!
      type(master_entry), pointer :: listentry

      ierr=pio_def_var (ncid, name, xtype , dimids(1:ndims), vid)

      listentry=> get_entry_by_name(masterlinkedlist, name)
      if(associated(listentry)) then
         ierr=pio_put_att (ncid, vid, 'long_name', listentry%field%long_name)
         ierr=pio_put_att (ncid, vid, 'units', listentry%field%units)
      end if
         
!
! Field not found in masterlist: long_name and units attributes unknown so just
! return
!
      return
   end subroutine addvar

!#######################################################################

   subroutine addfld (fname, units, numlev, avgflag, long_name, &
                      decomp_type, flag_xyfill,sampling_seq, &
                      begdim1, enddim1, &
                      begdim2, enddim2,&
                      begdim3, enddim3, mdimnames, fill_value )
                 
!
!----------------------------------------------------------------------- 
! 
! Purpose: Add a field to the master field list
! 
! Method: Put input arguments of field name, units, number of levels, averaging flag, and 
!         long name into a type entry in the global master field list (masterlist).
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use ppgrid,        only: begchunk, endchunk, pver, pverp
      use rgrid,         only: nlon
      use phys_grid,     only: get_ncols_p, physgrid_set
      use dycore,        only: dycore_is
      use cam_pio_utils, only: phys_decomp, dyn_decomp, dyn_stagger_decomp
      use cam_history_support, only : fillvalue

!
! Arguments
!
      character(len=*), intent(in) :: fname      ! field name--should be "max_fieldname_len" characters long
                                                 ! or less
      character(len=*), intent(in) :: units      ! units of fname--should be 8 chars
      character(len=1), intent(in) :: avgflag    ! averaging flag
      character(len=*), intent(in) :: long_name  ! long name of field
      
      integer, intent(in) :: numlev              ! number of vertical levels (dimension and loop)
      integer, intent(in) :: decomp_type         ! decomposition type

      logical, intent(in), optional :: flag_xyfill ! non-applicable xy points flagged with fillvalue

      character(len=*), intent(in), optional :: sampling_seq ! sampling sequence - if not every timestep, 
                                                             ! how often field is sampled:  
                                                             ! every other; only during LW/SW radiation calcs, etc.
      integer, intent(in), optional :: begdim1, enddim1
      integer, intent(in), optional :: begdim2, enddim2
      integer, intent(in), optional :: begdim3, enddim3

      character(len=*), intent(in), optional :: mdimnames(:)
      real(r8), intent(in), optional :: fill_value


      integer dim1s,dim2s                        ! global size of the first and second horizontal dim.

!
! Local workspace
!
      character(len=max_fieldname_len) :: fname_tmp ! local copy of fname
      type(master_entry), pointer :: listentry

      integer :: c             ! chunk (physics ) or latitude (dynamics) index
      integer :: ncol          ! number of columns per chunk
      integer :: beglat, endlat
      integer :: beglatxy, endlatxy, beglonxy, endlonxy, beglon, endlon
      integer :: dimcnt

      if (htapes_defined) then
         call endrun ('ADDFLD: Attempt to add field '//trim(fname)//' after history files set')
      end if

      beglat = get_dyn_grid_parm('beglat')
      endlat = get_dyn_grid_parm('endlat')
      beglon = get_dyn_grid_parm('beglon')
      endlon = get_dyn_grid_parm('endlon')
      beglatxy = get_dyn_grid_parm('beglatxy')
      endlatxy = get_dyn_grid_parm('endlatxy')
      beglonxy = get_dyn_grid_parm('beglonxy')
      endlonxy = get_dyn_grid_parm('endlonxy')
!
! Ensure that required grid info is available now
!
      select case (decomp_type)
      case (phys_decomp)
         if (.not. physgrid_set) then
            call endrun ('ADDFLD: Attempt to add field '//fname//' to masterlist before physics grid set')
         end if
      case (dyn_decomp)
         if (.not. dyndecomp_set) then
            call endrun ('ADDFLD: Attempt to add field '//fname//' to masterlist before dynamics grid set')
         end if
      end select
      call get_horiz_grid_dim_d(dim1s,dim2s)

!
! Ensure that new field name is not all blanks
!
      if (len_trim(fname)==0) then
         call endrun('ADDFLD: blank field name not allowed')
      end if
!
! Ensure that new field name is not longer than allowed
! (strip "&IC" suffix if it exists)
!
      fname_tmp  = fname
      fname_tmp  = strip_suffix(fname_tmp)

      if (len_trim(fname_tmp) > fieldname_len) then
         write(iulog,*)'ADDFLD: field name cannot be longer than ',fieldname_len,' characters long'
         write(iulog,*)'Field name:  ',fname
         call endrun()
      end if
!
! Ensure that new field doesn't already exist
!
      listentry => get_entry_by_name(masterlinkedlist, fname)
      if(associated(listentry)) then
         call endrun ('ADDFLD:  '//fname//' already on list')
      end if

!
! Add field to Master Field List arrays fieldn and iflds
!
      allocate(listentry)
      listentry%field%name        = fname
      listentry%field%long_name   = long_name
      listentry%field%units       = units      
      listentry%field%numlev      = numlev
      listentry%field%decomp_type = decomp_type      

      nullify(listentry%field%mdims)
      nullify(listentry%field%colperdim3)

      listentry%htapeindx(:)=-1
      listentry%act_sometape = .false.
      listentry%actflag(:) = .false.

!
! Indicate sampling sequence of field (i.e., how often "outfld" is called)
! If not every timestep (default), then give a descriptor indicating the
! sampling pattern.  Currently, the only valid value is "rad_lwsw" for sampling
! during LW/SW radiation timesteps only
!
      if (present(sampling_seq)) then
         listentry%field%sampling_seq = sampling_seq
      else
         listentry%field%sampling_seq = ' '
      end if
!
! Whether to apply xy fillvalue: default is false
!
      if (present(flag_xyfill)) then
         listentry%field%flag_xyfill = flag_xyfill
      else
         listentry%field%flag_xyfill = .false.
      end if

!
!    Allow external packages to have fillvalues different than default
!

      if(present(fill_value)) then
         listentry%field%fillvalue = fill_value
      else
         listentry%field%fillvalue = fillvalue
      endif

!
! ID and size of extra dimensions
!

      if(present(mdimnames)) then
         dimcnt = size(mdimnames)
         allocate(listentry%field%mdims(dimcnt))
         call lookup_hist_mdim_indicies(mdimnames, listentry%field%mdims)
         if(dimcnt > maxvarmdims) maxvarmdims = dimcnt
      else if(numlev>1) then
         
         allocate(listentry%field%mdims(1))
         select case(numlev)
         case(pver)
            call lookup_hist_mdim_indicies((/'lev'/),listentry%field%mdims)
         case(pverp)
            call lookup_hist_mdim_indicies((/'ilev'/),listentry%field%mdims)
         case default
            write(iulog,*) 'No Vertical dimension = ',numlev,' defined for var ',trim(fname)
            call endrun()
         end select
      endif
!
! Dimension history info based on decomposition type (dynamics or physics)
!
      nullify (listentry%field%colperdim3)
      select case (decomp_type)
      case (phys_decomp)
         listentry%field%begdim1  = 1
         listentry%field%enddim1  = pcols
         listentry%field%begdim2  = 1
         listentry%field%enddim2  = numlev
         listentry%field%begdim3  = begchunk
         listentry%field%enddim3  = endchunk
         allocate (listentry%field%colperdim3(begchunk:endchunk))
         do c=begchunk,endchunk
            ncol = get_ncols_p(c)
            listentry%field%colperdim3(c) = ncol
         end do
      case (dyn_stagger_decomp)
            listentry%field%begdim1  = beglonxy
            listentry%field%enddim1  = endlonxy
            listentry%field%begdim2  = 1
            listentry%field%enddim2  = numlev
            listentry%field%begdim3  = beglatxy
            listentry%field%enddim3  = endlatxy
      case (dyn_decomp)
         if ( dycore_is('LR') )then
            listentry%field%begdim1  = beglonxy
            listentry%field%enddim1  = endlonxy
            listentry%field%begdim2  = 1
            listentry%field%enddim2  = numlev
            listentry%field%begdim3  = beglatxy
            listentry%field%enddim3  = endlatxy
         else
            listentry%field%begdim3  = beglat
            listentry%field%enddim3  = endlat
            listentry%field%begdim2  = 1
            listentry%field%enddim2  = numlev
            listentry%field%begdim1  = 1
            listentry%field%enddim1  = dim1s
            if(.not. dycore_is('UNSTRUCTURED')) then
               allocate (listentry%field%colperdim3(beglat:endlat))
               do c=beglat,endlat
                  listentry%field%colperdim3(c) = nlon(c)
               end do
            endif
         endif

      case default
         write(iulog,*)'ADDFLD: unknown decomp_type=', decomp_type
         call endrun ()
      end select

      if(present(begdim1)) then
         listentry%field%begdim1  = begdim1
      end if
      if(present(enddim1)) then
         listentry%field%enddim1  = enddim1
      end if
      if(present(begdim2)) then
         listentry%field%begdim2  = begdim2
      end if
      if(present(enddim2)) then
         listentry%field%enddim2  = enddim2
      end if
      if(present(begdim3)) then
         listentry%field%begdim3  = begdim3
      end if
      if(present(enddim3)) then
         listentry%field%enddim3  = enddim3
      end if



!
! These 2 fields are used only in master field list, not runtime field list
!
      listentry%avgflag(:) = avgflag
      listentry%actflag(:) = .false.

      select case (avgflag)
      case ('A')
         listentry%time_op(:) = 'mean'
      case ('B')
         listentry%time_op(:) = 'mean00z'
      case ('I')
         listentry%time_op(:) = ' '
      case ('X')
         listentry%time_op(:) = 'maximum'
      case ('M')
         listentry%time_op(:) = 'minimum'
      case ('L')
         listentry%time_op(:) = LT_DESC
      case default
         call endrun ('ADDFLD: unknown avgflag='//avgflag)
      end select
      nullify(listentry%next_entry)

      call add_entry_to_master(listentry)
      return
   end subroutine addfld

   subroutine add_entry_to_master( newentry)
     type(master_entry), target, intent(in) :: newentry
     type(master_entry), pointer :: listentry
 
     if(associated(masterlinkedlist)) then
        listentry => masterlinkedlist
        do while(associated(listentry%next_entry))
           listentry=>listentry%next_entry
        end do
        listentry%next_entry=>newentry
     else
        masterlinkedlist=>newentry
     end if

   end subroutine add_entry_to_master

!#######################################################################

   subroutine wrapup (rstwr, nlend)
!
!-----------------------------------------------------------------------
!
! Purpose: 
! Close history files.
! 
! Method: 
! This routine will close any full hist. files
! or any hist. file that has data on it when restart files are being 
! written.
! If a partially full history file was disposed (for restart 
! purposes), then wrapup will open that unit back up and position 
! it for appending new data. 
!
! Original version: CCM2
!
!-----------------------------------------------------------------------
!
      use pio, only : pio_file_is_open
      use shr_kind_mod,  only: r8 => shr_kind_r8
      use pspect
      use ioFileMod
      use time_manager,  only: get_nstep, get_curr_date, get_curr_time
      use cam_pio_utils, only: cam_pio_openfile

!
! Input arguments
!
      logical, intent(in) :: rstwr   ! true => restart files are written this timestep
      logical, intent(in) :: nlend   ! Flag if time to end

!
! Local workspace
!
      integer :: nstep          ! current timestep number
      integer :: ncdate         ! current date in integer format [yyyymmdd]
      integer :: ncsec          ! time of day relative to current date [seconds]
      integer :: ndcur          ! days component of current time
      integer :: nscur          ! seconds component of current time
      integer :: yr, mon, day   ! year, month, day components of a date

      logical lfill   (ptapes)  ! Is history file ready to dispose?
      logical hdispose(ptapes)  ! Primary file disposed
      logical lhdisp            ! true => history file is disposed
      logical lhfill            ! true => history file is full

      integer t                 ! History file number
      integer f
      integer :: ierr
      real(r8) tday             ! Model day number for printout
!-----------------------------------------------------------------------

      tape => history_tape

      nstep = get_nstep()
      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day
      call get_curr_time(ndcur, nscur)
!
!-----------------------------------------------------------------------
! Dispose history files.
!-----------------------------------------------------------------------
!
! Begin loop over ptapes (the no. of declared history files - primary
! and auxiliary).  This loop disposes a history file to Mass Store
! when appropriate.
!
      do t=1,ptapes
         if (nflds(t) == 0) cycle

         hdispose(t) = .false.
         lfill(t) = .false.
!
! Find out if file is full
!
         if (hstwr(t) .and. nfils(t) >= mfilt(t)) then
            lfill(t) = .true.
         endif
!
! Dispose history file if 
!    1) file is filled or 
!    2) this is the end of run and file has data on it or
!    3) restarts are being put out and history file has data on it
!
         if (lfill(t) .or. (nlend .and. nfils(t) >= 1) .or. (rstwr .and. nfils(t) >= 1)) then
!
! Dispose history file
!
            hdispose(t) = .true.
!
! Is this the 0 timestep data of a monthly run?
! If so, just close primary unit do not dispose. 
!
            if (masterproc) write(iulog,*)'WRAPUP: nf_close(',t,')=',trim(nhfil(t))
            if(pio_file_is_open(tape(t)%File)) then
               if (nlend .or. lfill(t)) then
                  do f=1,nflds(t)
                     if (associated(tape(t)%hlist(f)%varid)) then
                        deallocate(tape(t)%hlist(f)%varid)
                        nullify(tape(t)%hlist(f)%varid)
                     end if
                  end do
               end if
               call PIO_CloseFile (tape(t)%File)
            end if
            if (nhtfrq(t) /= 0 .or. nstep > 0) then

! 
! Print information concerning model output.
! Model day number = iteration number of history file data * delta-t / (seconds per day)
! 
               tday = ndcur + nscur/86400._r8
               if(masterproc) then
                  if (t==1) then
                     write(iulog,*)'   Primary history file'
                  else
                     write(iulog,*)'   Auxiliary history file number ', t-1
                  end if
                  write(iulog,9003)nstep,nfils(t),tday
                  write(iulog,9004)
               end if
                  !                      
! Auxilary files may have been closed and saved off without being full. 
! We must reopen the files and position them for more data.
! Must position auxiliary files if not full
!              
               if (.not.nlend .and. .not.lfill(t)) then
                  call cam_PIO_openfile (tape(t)%File, nhfil(t), PIO_WRITE)
                  call h_inquire(t)
               end if
            endif                 ! if 0 timestep of montly run****
         end if                      ! if time dispose history fiels***   
      end do                         ! do ptapes
!
! Reset number of files on each history tape
!
      do t=1,ptapes
         if (nflds(t) == 0) cycle
         lhfill = hstwr(t) .and. nfils(t) >= mfilt(t)
         lhdisp = lhfill .or. (nlend .and. nfils(t) >= 1) .or. &
              (rstwr .and. nfils(t) >= 1)
         if (lhfill.and.lhdisp) then
            nfils(t) = 0
         endif
      end do
      return
9003  format('    Output at NSTEP     = ',i10,/, &
             '    Number of time samples on this file = ',i10,/, &
             '    Model Day           = ',f10.2)
9004  format('---------------------------------------')
   end subroutine wrapup


   integer function gen_hash_key(string)
!dir$ INLINENEVER gen_hash_key
!
!-----------------------------------------------------------------------
!
! Purpose: Generate a hash key on the interval [0 .. tbl_hash_pri_sz-1]
!          given a character string.
!
! Algorithm is a variant of perl's internal hashing function.
!
!-----------------------------------------------------------------------
!
   implicit none
!
!  Arguments:
!
   character(len=*), intent(in) :: string
!
!  Local.
!
   integer :: hash
   integer :: i

   hash = gen_hash_key_offset

   if ( len(string) /= 19 ) then
!
!     Process arbitrary string length.
!
      do i = 1, len(string)
         hash = ieor(hash , (ichar(string(i:i)) * tbl_gen_hash_key(iand(i-1,tbl_max_idx))))
      end do
   else
!
!     Special case string length = 19
!
      hash = ieor(hash , ichar(string(1:1))   * 61)
      hash = ieor(hash , ichar(string(2:2))   * 59)
      hash = ieor(hash , ichar(string(3:3))   * 53)
      hash = ieor(hash , ichar(string(4:4))   * 47)
      hash = ieor(hash , ichar(string(5:5))   * 43)
      hash = ieor(hash , ichar(string(6:6))   * 41)
      hash = ieor(hash , ichar(string(7:7))   * 37)
      hash = ieor(hash , ichar(string(8:8))   * 31)
      hash = ieor(hash , ichar(string(9:9))   * 29)
      hash = ieor(hash , ichar(string(10:10)) * 23)
      hash = ieor(hash , ichar(string(11:11)) * 17)
      hash = ieor(hash , ichar(string(12:12)) * 13)
      hash = ieor(hash , ichar(string(13:13)) * 11)
      hash = ieor(hash , ichar(string(14:14)) * 7)
      hash = ieor(hash , ichar(string(15:15)) * 3)
      hash = ieor(hash , ichar(string(16:16)) * 1)
      hash = ieor(hash , ichar(string(17:17)) * 61)
      hash = ieor(hash , ichar(string(18:18)) * 59)
      hash = ieor(hash , ichar(string(19:19)) * 53)
   end if

   gen_hash_key = iand(hash, tbl_hash_pri_sz-1)

   return

   end function gen_hash_key

!#######################################################################

   integer function get_masterlist_indx(fldname)
!
!-----------------------------------------------------------------------
!
! Purpose: Return the the index of the field's name on the master file list.
!
!          If the field is not found on the masterlist, return -1.
!
!-----------------------------------------------------------------------
!
!  Arguments:
!
   character(len=*), intent(in) :: fldname
!
!  Local.
!
   integer :: hash_key
   integer :: ff
   integer :: ii
   integer :: io   ! Index of overflow chain in overflow table
   integer :: in   ! Number of entries on overflow chain

   hash_key = gen_hash_key(fldname)
   ff = tbl_hash_pri(hash_key)
   if ( ff < 0 ) then
      io = abs(ff)
      in = tbl_hash_oflow(io)
      do ii = 1, in
         ff = tbl_hash_oflow(io+ii)
         if ( masterlist(ff)%thisentry%field%name == fldname ) exit
      end do
   end if

   if (ff == 0) then
      ! fldname generated a hash key that doesn't have an entry in tbl_hash_pri.
      ! This means that fldname isn't in the masterlist
      call endrun ('GET_MASTERLIST_INDX: attemping to output field '//fldname//' not on master list')
   end if
   
   if (associated(masterlist(ff)%thisentry) .and. masterlist(ff)%thisentry%field%name /= fldname ) then
      call endrun ('GET_MASTERLIST_INDX: error finding field '//fldname//' on master list')
   end if

   get_masterlist_indx = ff
   return
   end function get_masterlist_indx
!#######################################################################

   subroutine bld_outfld_hash_tbls()
!
!-----------------------------------------------------------------------
!
! Purpose: Build primary and overflow hash tables for outfld processing.
!
! Steps:
!  1) Foreach field on masterlist, find all collisions.
!  2) Given the number of collisions, verify overflow table has sufficient
!     space.
!  3) Build primary and overflow indices.
!
!-----------------------------------------------------------------------
!
!  Local.
!
   integer :: ff
   integer :: ii
   integer :: itemp
   integer :: ncollisions
   integer :: hash_key
   type(master_entry), pointer :: listentry
!
!  1) Find all collisions.
!
   tbl_hash_pri = 0

   ff=0
   allocate(masterlist(nfmaster))
   listentry=>masterlinkedlist
   do while(associated(listentry))
      ff=ff+1
      masterlist(ff)%thisentry=>listentry
      listentry=>listentry%next_entry
   end do
   if(ff /= nfmaster) then
      write(iulog,*) 'nfmaster = ',nfmaster, ' ff=',ff
      call endrun('mismatch in expected size of nfmaster')
   end if

   
   do ff = 1, nfmaster
      hash_key = gen_hash_key(masterlist(ff)%thisentry%field%name)
      tbl_hash_pri(hash_key) = tbl_hash_pri(hash_key) + 1
   end do

!
!  2) Count number of collisions and define start of a individual
!     collision's chain in overflow table. A collision is defined to be any
!     location in tbl_hash_pri that has a value > 1.
!
   ncollisions = 0
   do ii = 0, tbl_hash_pri_sz-1
      if ( tbl_hash_pri(ii) > 1 ) then  ! Define start of chain in O.F. table
         itemp = tbl_hash_pri(ii)
         tbl_hash_pri(ii) = -(ncollisions + 1)
         ncollisions = ncollisions + itemp + 1
      end if
   end do

   if ( ncollisions > tbl_hash_oflow_sz ) then
      write(iulog,*) 'BLD_OUTFLD_HASH_TBLS: ncollisions > tbl_hash_oflow_sz', &
              ncollisions, tbl_hash_oflow_sz
      call endrun()
   end if

!
!  3) Build primary and overflow tables.
!     i - set collisions in tbl_hash_pri to point to their respective
!         chain in the overflow table.
!
   tbl_hash_oflow = 0

   do ff = 1, nfmaster
      hash_key = gen_hash_key(masterlist(ff)%thisentry%field%name)
      if ( tbl_hash_pri(hash_key) < 0 ) then
         ii = abs(tbl_hash_pri(hash_key))
         tbl_hash_oflow(ii) = tbl_hash_oflow(ii) + 1
         tbl_hash_oflow(ii+tbl_hash_oflow(ii)) = ff
      else
         tbl_hash_pri(hash_key) = ff
      end if
   end do

!
!  Dump out primary and overflow hashing tables.
!
!   if ( masterproc ) then
!      do ii = 0, tbl_hash_pri_sz-1
!         if ( tbl_hash_pri(ii) /= 0 ) write(iulog,666) 'tbl_hash_pri', ii, tbl_hash_pri(ii)
!      end do
!
!      do ii = 1, tbl_hash_oflow_sz
!         if ( tbl_hash_oflow(ii) /= 0 ) write(iulog,666) 'tbl_hash_oflow', ii, tbl_hash_oflow(ii)
!      end do
!
!      itemp = 0
!      ii = 1
!      do 
!         if ( tbl_hash_oflow(ii) == 0 ) exit
!         itemp = itemp + 1
!         write(iulog,*) 'Overflow chain ', itemp, ' has ', tbl_hash_oflow(ii), ' entries:'
!         do ff = 1, tbl_hash_oflow(ii)  ! dump out colliding names on this chain
!            write(iulog,*) '     ', ff, ' = ', tbl_hash_oflow(ii+ff), &
!                       ' ', masterlist(tbl_hash_oflow(ii+ff))%thisentry%field%name
!         end do
!         ii = ii + tbl_hash_oflow(ii) +1 !advance pointer to start of next chain
!      end do
!   end if

   return
666 format(1x, a, '(', i4, ')', 1x, i6)

   end subroutine bld_outfld_hash_tbls

!#######################################################################

   subroutine bld_htapefld_indices
!
!-----------------------------------------------------------------------
!
! Purpose: Set history tape field indicies in masterlist for each
!          field defined on every tape.
!
! Note: because of restart processing, the actflag field is cleared and
!       then set only for active output fields on the different history
!       tapes.
!
!-----------------------------------------------------------------------
!
!  Arguments:
!

!
!  Local.
!
   integer :: f
   integer :: t

!
!  Initialize htapeindx to an invalid value.
!
   type(master_entry), pointer :: listentry

   ! reset all the active flags to false 
   ! this is needed so that restarts work properly -- fvitt
   listentry=>masterlinkedlist
   do while(associated(listentry))
      listentry%actflag(:) = .false.
      listentry%act_sometape = .false.
      listentry=>listentry%next_entry
   end do

   do t = 1, ptapes
      do f = 1, nflds(t)
         listentry => get_entry_by_name(masterlinkedlist, tape(t)%hlist(f)%field%name)
         if(.not.associated(listentry)) then
            write(iulog,*) 'BLD_HTAPEFLD_INDICES: something wrong, field not found on masterlist'
            write(iulog,*) 'BLD_HTAPEFLD_INDICES: t, f, ff = ', t, f
            write(iulog,*) 'BLD_HTAPEFLD_INDICES: tape%name = ', tape(t)%hlist(f)%field%name
            call endrun
         end if
         listentry%act_sometape = .true.
         listentry%actflag(t) = .true.
         listentry%htapeindx(t) = f
      end do
   end do

!
! set flag indicating h-tape contents are now defined (needed by addfld)
!
   htapes_defined = .true.      

   return
   end subroutine bld_htapefld_indices

!#######################################################################

   logical function hist_fld_active (fname)
!
!------------------------------------------------------------------------
!
! Purpose: determine if a field is active on any history file
!
!------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: fname ! Field name
!
! Local variables
!
      character*(max_fieldname_len) :: fname_loc  ! max-char equivalent of fname
      integer :: ff                  ! masterlist index pointer
!-----------------------------------------------------------------------

      fname_loc = fname
      ff = get_masterlist_indx(fname_loc)
      if ( ff < 0 ) then
         hist_fld_active = .false.
      else 
         hist_fld_active = masterlist(ff)%thisentry%act_sometape
      end if

   end function hist_fld_active

!#######################################################################

function hist_fld_col_active(fname, lchnk)

! Determine whether each column in a field is active on any history file.
! The purpose of this routine is to provide information which would allow
! a diagnostic physics parameterization to only be run on a subset of
! columns in the case when only column or regional output is requested.
!
! **N.B.** The field is assumed to be using the physics decomposition.

   use phys_grid, only: get_ncols_p, get_lat_all_p, get_lon_all_p

! Arguments
   character(len=*), intent(in) :: fname ! Field name
   integer,          intent(in) :: lchnk ! chunk ID

! Return value
   logical :: hist_fld_col_active(pcols)

! Local variables
   character*(max_fieldname_len) :: fname_loc  ! max-char equivalent of fname

   integer :: ff         ! masterlist index pointer
   integer :: i
   integer :: t          ! history file (tape) index
   integer :: igp, iub    ! column group indices
   integer :: ilats(pcols)
   integer :: ilons(pcols)
   integer :: ncol, lon1, lon2, lat1, lat2

   type (active_entry), pointer :: tape(:)

!-----------------------------------------------------------------------

! Initialize to false.  Then look to see if and where active.
   hist_fld_col_active = .false.

! Check for name in the master list.
   fname_loc = fname
   ff = get_masterlist_indx(fname_loc)

! If not in master list then return.
   if ( ff < 0 ) return

! If in master list, but not active on any file then return
   if (.not. masterlist(ff)%thisentry%act_sometape) return

! Loop over history files and check for the field/column in each one

   tape => history_tape
   do t = 1, ptapes

! Is the field active in this file?  If not the cycle to next file.
      if (.not. masterlist(ff)%thisentry%actflag(t)) cycle

! Check whether this file has column output.
      if (tape(t)%column(1)%lat_name == ' ') then

! No column output has been requested.  In that case the field has
! global output which implies all columns are active.  No need to 
! check any other history files.
         hist_fld_col_active = .true.
         exit

      else if (associated(tape(t)%column(1)%hmask)) then

! All column output is collected into a single group for output.
         ncol = get_ncols_p(lchnk)
         do i = 1, ncol
            if (tape(t)%column(1)%hmask(i,lchnk) > 0) then
               hist_fld_col_active(i) = .true.
            end if
         end do

      else

! Column output has been requested, and multiple columns/regions 
! are possible.  Need to loop over the groups

! The column group info is stored as indices into the global array.
! Retrieve the global indices for the columns in the chunk.

         ncol = get_ncols_p(lchnk)
         call get_lat_all_p(lchnk, pcols, ilats)
         call get_lon_all_p(lchnk, pcols, ilons)

         do igp = 1, ngroup(t)
            lon1 = tape(t)%column(igp)%columnlon(1)
            lon2 = tape(t)%column(igp)%columnlon(2)
            lat1 = tape(t)%column(igp)%columnlat(1)
            lat2 = tape(t)%column(igp)%columnlat(2)

! Is the column in the group
            do i = 1, ncol
! only need to check if the column has not already been
! identified as active
               if (.not. hist_fld_col_active(i)) then
                  if ( ilons(i) >= lon1 .and. ilons(i) <= lon2 .and. &
                       ilats(i) >= lat1 .and. ilats(i) <= lat2 ) then
                     hist_fld_col_active(i) = .true.
                  end if
               end if
            end do ! columns in chunk

         end do ! groups in history file

      end if

   end do ! history files

end function hist_fld_col_active



end module cam_history
