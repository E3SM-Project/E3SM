module conditional_diag
!--------------------------------------------------------
! This is one of the modules that support conditional
! sampling and budget analysis in EAM.
!
! This module contains 
!  - derived types definitions
!  - subroutine for namelist handling
!  - subroutine for memory allocation
!
! History:
!  First version by Hui Wan, PNNL, March - May 2021
!--------------------------------------------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  implicit none

  private

  ! Derived types

  public cnd_diag_info_t   ! for metadata
  public cnd_diag_t        ! for physical quantities

  ! Variable(s) of derived type

  public cnd_diag_info  ! metadata

  ! Subroutines

  public cnd_diag_readnl
  public cnd_diag_alloc

  ! module parameters

  integer, parameter :: ncnd_max         = 10 ! max # of conditions allowed in a single simulation
  integer, parameter :: mname_maxlen     = 8  ! string length for metric name

  integer, parameter :: nqoi_max         = 20 ! max # of conditionally sampled QoIs in a single simulation
  integer, parameter :: qoiname_maxlen   = 8  ! string length for QoI name

  integer, parameter :: nchkpt_max       = 250 ! max # of active checkpoints in a single simulation
  integer, parameter :: chkptname_maxlen = 10  ! string length for checkpoint name

  ! what kind of dp (pressure layer thickness) to multiply QoI by

  integer, parameter         :: UNSET   = -1
  integer, parameter, public :: NODP    = 0
  integer, parameter, public :: PDEL    = 1
  integer, parameter, public :: PDELDRY = 2

  ! fillvalue for cells that are masked out 

  real(r8),parameter, public :: FILLVALUE = 0._r8

  !-------------------------------------------------------------------------------
  ! Derived type for metadata
  !-------------------------------------------------------------------------------
  type cnd_diag_info_t

    ! Do we want to write out the QoIs to history tapes?
    logical :: l_output_state = .false.

    ! Do we want to write out increments of the QoIs to history tapes? 
    logical :: l_output_incrm = .false.

    ! Which history tapes should contain a complete set of the conditionally sampled diagnostics?
    ! Note: we allow the user to specify multiple history tapes that will contain all output variables 
    ! associated with the conditional diagnostics. It is envisioned that these multiple history 
    ! tapes might have different output frequencies and/or averaging flags. Apart from that,
    ! like with all other output variables on the master list, the user can manually add or 
    ! remove individual variables using fincl or fexcl. 

    integer             :: ntape                         ! number of tapes
    integer,allocatable :: hist_tape_with_all_output(:)  ! tape indices

    ! Sampling conditions. 
    ! The current implementation allows the user to define multiple conditions,
    ! each of which will have its own metric, QoIs, and output.
    ! But to keep it simple (at least as a start), we assume that
    ! the QoIs and checkpoints to monitor are the same for different sampling conditions

    integer                      :: ncnd = 0             ! total # of sampling conditions used in this simulation
    character(len=mname_maxlen),&
                     allocatable :: metric_name(:)       ! shape = (ncnd); metric names
    integer,allocatable          :: metric_nver(:)       ! shape = (ncnd); # of vertical levels of the metrics
    integer,allocatable          :: metric_cmpr_type(:)  ! shape = (ncnd); see module parameters in conditional_diag_main.F90
    real(r8),allocatable         :: metric_threshold(:)  ! shape = (ncnd); threshold values for conditional sampling 
    real(r8),allocatable         :: metric_tolerance(:)  ! shape = (ncnd); tolerance for the "equal to" comparison type
    character(len=chkptname_maxlen),&
                     allocatable :: cnd_eval_chkpt(:)    ! shape = (ncnd); checkpoints at which the evaluation of 
                                                         ! the corresponding sampling condition will happen
    character(len=chkptname_maxlen),&
                     allocatable :: cnd_end_chkpt(:)     ! shape = (ncnd); checkpoints marking end-of-time-step
                                                         ! of the corresponding sampling conditions. 
                                                         ! At this checkpoint, the masking of QoIs will be done
                                                         ! and the masked values will be send to history buffer.

    ! QoIs to be monitored. Each QoI can have 1, pver, or pver+1 vertical levels
    integer                                   :: nqoi = 0
    character(len=qoiname_maxlen),allocatable :: qoi_name(:)       ! shape = (nqoi)
    integer,allocatable                       :: qoi_nver(:)       ! shape = (nqoi); # of vertical levels of the QoI
    integer,allocatable                       :: qoi_nver_save(:)  ! shape = (nqoi); # of vertical levels of the QoI
                                                                   ! if no vertical integral is requested; 1 otherwise.

    ! Active checkpoints at which the QoI will be monitored 
    integer                                      :: nchkpt = 0     ! total # of active checkpoints
    character(len=chkptname_maxlen), allocatable :: qoi_chkpt(:)   ! checkpoints at which QoIs will be monitored

    ! "Multiply by dp" options
    integer,allocatable :: x_dp(:,:)  ! shape = (nqoi,nchkpt); whether QoIs at each checkpoint and checkpoint 
                                      ! should be multiplied by dp, and if so, wet or dry.


  end type cnd_diag_info_t

  !-------------------------------------------------------------------------------
  ! Derived types for the sampling conditions and monitored QoIs
  !-------------------------------------------------------------------------------
  ! Values of a single QoI at different checkpoints and the 
  ! inrements relative to the previous active checkpoint

  type snapshots_and_increments_t

    real(r8), allocatable :: val(:,:,:) ! shape = (pcols,info%qoi_nver_save(iqoi),info%nchkpt) QoI (or vert. integral) values at active checkpoints
    real(r8), allocatable :: inc(:,:,:) ! shape = (pcols,info%qoi_nver_save(iqoi),info%nchkpt) QoI (or vert. integral) increments between adjacent active checkpoints
    real(r8), allocatable :: old(:,:)   ! shape = (pcols,info%qoi_nver_save(iqoi)) QoI (or vert. integral) values at the previous active checkpoint

  end type snapshots_and_increments_t

  !----------------------------------------------------------------------
  ! The collection of all QoIs monitored under the same condition, and
  ! the values of the metric and flag used for sampling

  type metric_and_qois_t

    real(r8),                        allocatable :: metric (:,:)     ! shape = (pcols, info%metric_nver(icnd))
    real(r8),                        allocatable :: flag   (:,:)     ! shape = (pcols, info%metric_nver(icnd))
    type(snapshots_and_increments_t),allocatable :: qoi    (:)       ! shape = (info%nqoi)

  end type metric_and_qois_t

  !----------------------------------------------------------------------
  ! A collection of multiple conditions (including the corresponding
  ! metrics and monitored QoIs)

  type cnd_diag_t

    type(metric_and_qois_t), allocatable :: cnd(:) ! shape = (info%ncnd)

  end type cnd_diag_t

!===============================================================================
! Module variables
!===============================================================================
  type(cnd_diag_info_t) :: cnd_diag_info

contains
!===============================================================================
! Procedures
!===============================================================================
subroutine cnd_diag_readnl(nlfile)

   use cam_history_support,only: ptapes
   use ppgrid,             only: pver
   use infnan,             only: nan, assignment(=), isnan
   use namelist_utils,     only: find_group_name
   use units,              only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   character(len=*), parameter :: subname = 'cnd_diag_readnl'

   integer :: ncnd, nchkpt, nqoi, ntape, ii, icnd, ichkpt

   ! Local variables for reading namelist

   integer :: unitn, ierr

   character(len=mname_maxlen)     :: metric_name     (ncnd_max)
   integer                         :: metric_nver     (ncnd_max)
   integer                         :: metric_cmpr_type(ncnd_max)
   real(r8)                        :: metric_threshold(ncnd_max)
   real(r8)                        :: metric_tolerance(ncnd_max)

   character(len=chkptname_maxlen) :: cnd_eval_chkpt  (ncnd_max)
   character(len=chkptname_maxlen) :: cnd_end_chkpt   (ncnd_max)

   character(len=chkptname_maxlen) :: qoi_chkpt(nchkpt_max)

   character(len=qoiname_maxlen)   :: qoi_name (nqoi_max)
   integer                         :: qoi_nver (nqoi_max)

   integer ::   qoi_x_dp(nqoi_max)
   integer :: chkpt_x_dp(nchkpt_max)

   logical :: l_output_state, l_output_incrm
   integer :: hist_tape_with_all_output(ptapes)  ! tape indices

   !-------
   namelist /conditional_diag_nl/     &
            metric_name, metric_nver, &
            metric_cmpr_type, metric_threshold, metric_tolerance, &
            cnd_eval_chkpt, cnd_end_chkpt,  & 
            qoi_chkpt, qoi_name, qoi_nver,  &
            qoi_x_dp, chkpt_x_dp, &
            l_output_state, l_output_incrm, &
            hist_tape_with_all_output

   !----------------------------------------
   !  Default values
   !----------------------------------------
   metric_name(:)      = ' '
   metric_nver(:)      = 0
   metric_cmpr_type(:) = -99
   metric_threshold(:) = nan
   metric_tolerance(:) = 0._r8

   cnd_eval_chkpt(:) = ' '
    cnd_end_chkpt(:) = ' '

   qoi_name (:) = ' '
   qoi_nver (:) = 0 
   qoi_chkpt(:) = ' '

     qoi_x_dp(:) = NODP
   chkpt_x_dp(:) = UNSET

   l_output_state = .false.
   l_output_incrm = .false.

   hist_tape_with_all_output(:) = -1

   !----------------------------------------
   ! Read namelist and check validity
   !----------------------------------------
   if (masterproc) then

      ! Read namelist

      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'conditional_diag_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, conditional_diag_nl, iostat=ierr)
         if (ierr /= 0) then
            write(iulog,*) 'read error ',ierr
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)

      ! Count user-specified sampling conditions

      ii = 0
      do while ( (ii+1) <= ncnd_max .and. metric_name(ii+1) /= ' ')
         ii = ii + 1
      end do
      ncnd = ii

      !----------------------------------------------------------------------
      ! If no condition has been specified, set the other counts to zero, too
      !----------------------------------------------------------------------
      if (ncnd==0) then

         nqoi   = 0
         nchkpt = 0
         ntape  = 0

      !----------------------------------------------------------------------
      ! If at least one condition has been sepecified, do some sanity check, 
      ! then parse additional namelist settings
      !----------------------------------------------------------------------
      else

         do ii = 1,ncnd
            if (trim(adjustl(metric_name(ii)))=='ALL') then

            ! Metric name 'ALL' is interpreted as selecting all grid cells.
            ! In this case, the metric will be set to a constant field of 1 and
            ! the flags will be set to ON.

            ! - set metric_nver(ii) to 1 to save memory

                metric_nver(ii) = 1

            ! - metric_cmpr_type and metric_threshold will no longer
            !   be needed; set them to some values to avoid program abort.

                metric_cmpr_type(ii) = 0
                metric_threshold(ii) = 1._r8
                metric_tolerance(ii) = 1._r8

            ! - if user did not specify cnd_end_chkpt(ii), set it to the non-empty cnd_eval_chkpt(ii).
            ! - if neither cnd_end_chkpt(ii) or cnd_eval_chkpt(ii) is specified, set 
            !   cnd_end_chkpt(ii) to 'PBCDIAG' as this is the checkpoint at which 
            !   where most of the standard model output variables are sent to history buffer.

                if ( cnd_end_chkpt(ii) == ' ' ) then 
                   if ( cnd_eval_chkpt(ii) /= ' ' ) then 
                      cnd_end_chkpt(ii) = cnd_eval_chkpt(ii) 
                   else
                      cnd_end_chkpt(ii) = 'PBCDIAG'
                   end if
                end if

            ! - cnd_eval_chkpt is not needed as we will sample all time steps; 
            !   set to cnd_end_chkpt so that it has a value, but the value really doesn't matter.

                cnd_eval_chkpt(ii) = cnd_end_chkpt(ii)

            end if
         end do
         ! metric_nver = -1 read in from namelist should be interpreted as pver for flexible configuration
         where ( metric_nver(1:ncnd) == -1 ) metric_nver(1:ncnd) = pver

         if (any( metric_nver     (1:ncnd) <= 0   )) call endrun(subname//' error: need positive metric_nver for each metric_name')
         if (any( metric_cmpr_type(1:ncnd) == -99 )) call endrun(subname//' error: need valid metric_cmpr_type for each metric_name')
         if (any( isnan(metric_threshold(1:ncnd)) )) call endrun(subname//' error: need valid metric_threshold for each metric_name')
         if (any( cnd_eval_chkpt  (1:ncnd) == ' ' )) call endrun(subname//' error: be sure to specify cnd_eval_chkpt for each metric_name')
         
         do ii = 1,ncnd
            if (cnd_end_chkpt(ii)==' ')  cnd_end_chkpt(ii) = cnd_eval_chkpt(ii)
         end do

         !-------------------------------------------------------
         ! Count QoIs to be monitored, then do some sanity check

         ii = 0
         do while ( (ii+1) <= nqoi_max .and. qoi_name(ii+1) /= ' ')
            ii = ii + 1
         end do
         nqoi = ii

         ! qoi_nver = -1 read in from namelist should be interpreted as pver for flexible configuration
         where ( qoi_nver(1:nqoi) == -1 ) qoi_nver(1:nqoi) = pver
         if (any(qoi_nver(1:nqoi)<=0)) call endrun(subname//'error: need positive qoi_nver for each qoi_name')

         !---------------------------------------------
         ! Count active checkpoints for QoI monitoring

         ii = 0
         do while ( (ii+1) <= nchkpt_max .and. qoi_chkpt(ii+1) /= ' ')
            ii = ii + 1
         end do
         nchkpt = ii

         !---------------------
         ! Enforce consistency

         if (nqoi==0) nchkpt = 0 ! If user did not specify any QoI, set nchkpt to 0 for consistency
         if (nchkpt==0) nqoi = 0 ! If user did not specify any checkpoint for QoI monitoring, set nqoi to 0 for consistency

         !--------------------------
         ! "multiply by dp" options

         ! ... for QoIs

         do ii = 1,nqoi

            ! Set user-provided unexpected values to NODP

            if ( qoi_x_dp(ii)/=PDEL     .and. &
                 qoi_x_dp(ii)/=PDEL+100 .and. &
                 qoi_x_dp(ii)/=PDELDRY  .and. &
                 qoi_x_dp(ii)/=PDELDRY+100    ) then

               qoi_x_dp(ii) = NODP
            end if

            ! Multiplication by dp can only be applied to QoIs defined at layer midpoints
            ! or as layer averages. Force selection to be NODP for other types of QoIs.

            if (qoi_nver(ii)/=pver) qoi_x_dp(ii) = NODP

         end do

         ! ... for checkpoints

         do ii = 1,nchkpt
            if (chkpt_x_dp(ii)/=PDEL .and. chkpt_x_dp(ii)/=PDELDRY) chkpt_x_dp(ii) = UNSET
         end do

         !---------------------------------------------------------------------------------
         ! Count history tapes that will each contain a full suite of the output variables

         ii = 0
         do while ( (ii+1) <= ptapes .and. hist_tape_with_all_output(ii+1) >= 0)
            ii = ii + 1
         end do
         ntape = ii

         ! If the user did not specify any tape, then at least add the output variables to h0

         if (ntape==0) then
            ntape = 1
            hist_tape_with_all_output(ntape) = 1
         end if

      end if ! ncnd = 0 or > 0

   end if ! masterproc
   !--------------------------------------

#ifdef SPMD
   call mpibcast(ncnd,   1, mpiint, 0, mpicom)
   call mpibcast(nqoi,   1, mpiint, 0, mpicom)
   call mpibcast(nchkpt, 1, mpiint, 0, mpicom)
   call mpibcast(ntape,  1, mpiint, 0, mpicom)
#endif

   cnd_diag_info% ncnd   = ncnd
   cnd_diag_info% nqoi   = nqoi
   cnd_diag_info% nchkpt = nchkpt
   cnd_diag_info% ntape  = ntape

   if (ncnd==0) then

      if (masterproc) then
         write(iulog,*)'==========================================================='
         write(iulog,*)'       *** Conditional diagnostics NOT requested ***'
         write(iulog,*)'==========================================================='
      end if

      return
   end if

#ifdef SPMD
   !--------------------------------------
   ! Broadcast namelist variables
   !--------------------------------------
   call mpibcast(metric_name,      ncnd_max*len(metric_name(1)), mpichar, 0, mpicom)
   call mpibcast(metric_nver,      ncnd_max,                     mpiint,  0, mpicom)
   call mpibcast(metric_cmpr_type, ncnd_max,                     mpiint,  0, mpicom)
   call mpibcast(metric_threshold, ncnd_max,                     mpir8,   0, mpicom)
   call mpibcast(metric_tolerance, ncnd_max,                     mpir8,   0, mpicom)

   call mpibcast(cnd_eval_chkpt,   ncnd_max*len(cnd_eval_chkpt(1)), mpichar, 0, mpicom)
   call mpibcast(cnd_end_chkpt,    ncnd_max*len(cnd_end_chkpt(1)),  mpichar, 0, mpicom)
   call mpibcast(qoi_chkpt,      nchkpt_max*len(qoi_chkpt(1)),      mpichar, 0, mpicom)

   call mpibcast(qoi_name,  nqoi_max*len(qoi_name(1)),  mpichar, 0, mpicom)
   call mpibcast(qoi_nver,  nqoi_max,                   mpiint,  0, mpicom)

   call mpibcast(  qoi_x_dp,   nqoi_max, mpiint,  0, mpicom)
   call mpibcast(chkpt_x_dp, nchkpt_max, mpiint,  0, mpicom)

   call mpibcast(l_output_state, 1, mpilog, 0, mpicom)
   call mpibcast(l_output_incrm, 1, mpilog, 0, mpicom)

   call mpibcast(hist_tape_with_all_output, ptapes, mpiint, 0, mpicom)
#endif

   !-------------------------------------------
   ! Pack information into to cnd_diag_info
   !-------------------------------------------
   ! Conditions for conditional sampling

   allocate( cnd_diag_info% metric_name(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_name')
   do ii = 1,ncnd
      cnd_diag_info% metric_name(ii) = trim(adjustl(metric_name(ii)))
   end do

   allocate( cnd_diag_info% metric_nver(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_nver')
   cnd_diag_info% metric_nver(1:ncnd) = metric_nver(1:ncnd)

   allocate( cnd_diag_info% metric_cmpr_type(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_cmpr_type')
   cnd_diag_info% metric_cmpr_type(1:ncnd) = metric_cmpr_type(1:ncnd)

   allocate( cnd_diag_info% metric_threshold(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_threshold')
   cnd_diag_info% metric_threshold(1:ncnd) = metric_threshold(1:ncnd)

   allocate( cnd_diag_info% metric_tolerance(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_tolerance')
   cnd_diag_info% metric_tolerance(1:ncnd) = metric_tolerance(1:ncnd)

   allocate( cnd_diag_info% cnd_eval_chkpt(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% cnd_eval_chkpt')
   cnd_diag_info% cnd_eval_chkpt(1:ncnd) = cnd_eval_chkpt(1:ncnd)

   allocate( cnd_diag_info% cnd_end_chkpt(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% cnd_end_chkpt')
   cnd_diag_info% cnd_end_chkpt(1:ncnd) = cnd_end_chkpt(1:ncnd)

   ! QoIs 

   allocate( cnd_diag_info% qoi_name(nqoi), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% qoi_name')
   do ii = 1,nqoi
      cnd_diag_info% qoi_name(ii) = trim(adjustl(qoi_name(ii)))
   end do

   allocate( cnd_diag_info% qoi_nver(nqoi), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% qoi_nver')
   cnd_diag_info% qoi_nver(1:nqoi) = qoi_nver(1:nqoi)

   ! Active checkpoints at which the QoIs will be monitored

   allocate( cnd_diag_info% qoi_chkpt(nchkpt), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% qoi_chkpt')
   do ii = 1,nchkpt
      cnd_diag_info% qoi_chkpt(ii) = trim(adjustl(qoi_chkpt(ii)))
   end do

   ! "multiply by dp" selections
   
   allocate( cnd_diag_info% x_dp(nqoi,nchkpt), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% x_dp')

   call set_x_dp( nqoi, qoi_x_dp, nchkpt, chkpt_x_dp, cnd_diag_info%x_dp )! in, ..., in, out

   ! qoi_nver_save(ii) is the # of vertical levels on which the ii-th QoI and/or its increment
   ! will be saved in the data structure for output. qoi_nver_save(ii) will be 1 if 
   ! vertical integral is requested by the user by choosing an appropriate value for qoi_x_dp(ii).
   ! Otherwise  qoi_nver_save(ii) ==  qoi_nver(ii).

   allocate( cnd_diag_info% qoi_nver_save(nqoi), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% qoi_nver_save')

   where (qoi_x_dp(1:nqoi)==NODP)
      cnd_diag_info%qoi_nver_save = cnd_diag_info%qoi_nver
   elsewhere
      cnd_diag_info%qoi_nver_save = 1
   end where

   ! output to history tape(s)

   cnd_diag_info%l_output_state = l_output_state
   cnd_diag_info%l_output_incrm = l_output_incrm

   allocate( cnd_diag_info% hist_tape_with_all_output(ntape), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% hist_tape_with_all_output')
   cnd_diag_info% hist_tape_with_all_output(1:ntape) = hist_tape_with_all_output(1:ntape)

   !-----------------------------------------------
   ! Send information to log file
   !-----------------------------------------------
   if (masterproc) then

      write(iulog,*)
      write(iulog,*)'==============================================================================='
      write(iulog,*)'               *** Conditional diagnostics requested ***'
      write(iulog,*)'==============================================================================='

      write(iulog,*)
      write(iulog,'(4x,a12,a6,a12,2a12,2a17)') 'metric','nlev','cmpr_type','threshold','tolerance', &
                                               'cnd_eval_chkpt','cnd_end_chkpt' 
      do ii = 1,cnd_diag_info%ncnd
       if (trim(cnd_diag_info% metric_name(ii))=='ALL') then

         write(iulog,'(i4.3,a12,a6,a12,2a12,  2a17)') ii,                                 &
                                             adjustr(cnd_diag_info% metric_name(ii)),     &
                                                     '-','-','-','-','-',                 &
                                             adjustr(cnd_diag_info% cnd_end_chkpt(ii))
       else
         write(iulog,'(i4.3,a12,i6,i12,2e12.2,2a17)') ii,                                 &
                                             adjustr(cnd_diag_info% metric_name(ii)),     &
                                                     cnd_diag_info% metric_nver(ii),      &
                                                     cnd_diag_info% metric_cmpr_type(ii), &
                                                     cnd_diag_info% metric_threshold(ii), &
                                                     cnd_diag_info% metric_tolerance(ii), &
                                             adjustr(cnd_diag_info% cnd_eval_chkpt(ii)),  &
                                             adjustr(cnd_diag_info% cnd_end_chkpt(ii))
       end if
      end do

      write(iulog,*)
      write(iulog,*)'--------------------------------------------------'
      write(iulog,'(4x,a12,a15)') 'qoi_chkpt','mult_by_dp'
      do ii = 1,cnd_diag_info%nchkpt
         write(iulog,'(i4.3,a12,i15)') ii, adjustr(cnd_diag_info%qoi_chkpt(ii)), chkpt_x_dp(ii)
      end do

      write(iulog,*)
      write(iulog,*)'--------------------------------------------------'
      write(iulog,'(4x,a12,a8,a15,a12)')'QoI_name','nlev', 'mult_by_dp', 'nlev_save'
      do ii = 1,cnd_diag_info%nqoi
         write(iulog,'(i4.3,a12,i8,i15,i12)') ii, adjustr(cnd_diag_info%qoi_name(ii)),  &
                                              cnd_diag_info%qoi_nver(ii), qoi_x_dp(ii), &
                                              cnd_diag_info%qoi_nver_save(ii)
      end do
      write(iulog,*)

      write(iulog,*)'--------------------------------------------------'
      write(iulog,*)' l_output_state = ',l_output_state
      write(iulog,*)' l_output_incrm = ',l_output_incrm
      write(iulog,*)
      write(iulog,'(a,12i5)')' hist_tape_with_all_output = ',hist_tape_with_all_output(1:ntape)
      write(iulog,*)'--------------------------------------------------'
      write(iulog,*)
      write(iulog,*)'      "multiply by dp" selections, final'
      write(iulog,*)
      write(iulog,'(10x,20a10)') ( adjustr(cnd_diag_info%qoi_name(ii)), ii=1,nqoi )
      do ichkpt = 1,nchkpt
         write(iulog,'(a10,20i10)') adjustr(cnd_diag_info%qoi_chkpt(ichkpt)), ( cnd_diag_info%x_dp(ii,ichkpt), ii=1,nqoi )
      end do
      write(iulog,*)'==============================================================================='
      write(iulog,*)

  end if  ! masterproc

end subroutine cnd_diag_readnl

!===============================================================
subroutine set_x_dp( nqoi, qoi_x_dp, nchkpt, chkpt_x_dp, x_dp_out )

   integer,intent(in) :: nqoi, nchkpt
   integer,intent(in) ::   qoi_x_dp(:)
   integer,intent(in) :: chkpt_x_dp(:)
   integer,intent(out) ::  x_dp_out(:,:)

   integer :: icnd, iqoi, ichkpt

   character(len=256) :: msg

        do iqoi = 1,nqoi

           !---------------------------------------------------------------------------
           ! If user chose NODP, PDEL, or PDELDRY for this QoI, then take that value
           !---------------------------------------------------------------------------
           select case ( qoi_x_dp(iqoi) )
           case (NODP,PDEL,PDELDRY)

              x_dp_out(iqoi,:) = qoi_x_dp(iqoi) 

           !---------------------------------------------------------------------------------
           ! If user set qoi_x_dp(iqoi) to PDEL+100 or PDELDRY+100, then use PDEL or PDELDRY
           ! but give precedence to chkpt_x_dp.
           !---------------------------------------------------------------------------------
           case ( PDEL+100, PDELDRY+100 )

              !--------------------------------------------------------------
              ! Loop through all checkpoints. If chkpt_x_dp(ichkpt) has been 
              ! set, then take that value; otherwise, use qoi_x_dp(iqoi)
              !--------------------------------------------------------------
              do ichkpt = 1,nchkpt
                 if (chkpt_x_dp(ichkpt)/=(-1)) then
                    x_dp_out(iqoi,ichkpt) = chkpt_x_dp(ichkpt)
                 else
                    x_dp_out(iqoi,ichkpt) = qoi_x_dp(iqoi) - 100
                 end if
              end do ! ichkpt

           case default 
           !-----------------------------------------------------------------------
           ! Subroutine cnd_diag_readnl is supposed to have set the default value
           ! of qoi_x_dp to NODP, so we do not expect values other than
           ! NODP, PDEL(+100), or PDELDRY(+100).
           !-----------------------------------------------------------------------
              write(msg,'(a,i2,a,i2,a)') "qoi_x_dp(",iqoi,") =",qoi_x_dp(iqoi),' is unexpected'
              call endrun(trim(msg))

           end select ! qoi_x_dp(iqoi)
        end do    ! iqoi

end subroutine set_x_dp

!===============================================================================
subroutine cnd_diag_alloc( phys_diag, begchunk, endchunk, pcols, cnd_diag_info )

  type(cnd_diag_t), pointer :: phys_diag(:)

  integer, intent(in) :: begchunk, endchunk
  integer, intent(in) :: pcols
  type(cnd_diag_info_t),intent(in) :: cnd_diag_info

  integer :: ierr, lchnk

  character(len=*), parameter :: subname = 'cnd_diag_alloc'
  !-----------------------------------------------------------------

  allocate(phys_diag(begchunk:endchunk), stat=ierr)
  if( ierr /= 0 ) then
     write(iulog,*) subname//': phys_diag allocation error = ',ierr
     call endrun(subname//': failed to allocate phys_diag array')
  end if

  if (cnd_diag_info%ncnd <= 0) return

  do lchnk=begchunk,endchunk
     call single_chunk_cnd_diag_alloc( phys_diag(lchnk), lchnk, pcols, cnd_diag_info )
  end do

  if (masterproc) then
     write(iulog,*)
     write(iulog,*) "====================================================================="
     write(iulog,*) " Finished memory allocation for conditional diagnostics including:"
     write(iulog,*) cnd_diag_info%ncnd,   " conditions"
     write(iulog,*) cnd_diag_info%nqoi,   " quantities of interest (QoIs)"
     write(iulog,*) cnd_diag_info%nchkpt, " active checkpoints to monitor QoIs"
     write(iulog,*) "====================================================================="
     write(iulog,*)
  end if

end subroutine cnd_diag_alloc


!-----------------------------------------------------------------------------
! Allocate memory for conditional diagnostics in a single grid chunk
!-----------------------------------------------------------------------------
subroutine single_chunk_cnd_diag_alloc( diag, lchnk, psetcols, cnd_diag_info )

  type(cnd_diag_t), intent(inout)   :: diag
  integer,intent(in)                :: lchnk

  integer, intent(in)               :: psetcols
  type(cnd_diag_info_t), intent(in) :: cnd_diag_info

  integer :: ierr = 0
  integer :: ncnd, icnd

  character(len=*), parameter :: subname = 'single_chunk_cnd_diag_alloc'

  !----------------------------------------------
  ! Allocate an array for all sampling conditions
  !----------------------------------------------
  ncnd = cnd_diag_info%ncnd

  allocate( diag%cnd(ncnd), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of diag%cnd')

  !----------------------------------------------------------------
  ! Allocate memory for metrics and QoIs of each condition
  !----------------------------------------------------------------
  do icnd = 1,ncnd

     call metrics_and_qois_alloc( diag%cnd(icnd),                  &
                                  cnd_diag_info%metric_nver(icnd), &
                                  cnd_diag_info%nchkpt,            &
                                  cnd_diag_info%nqoi,              &
                                  cnd_diag_info%qoi_nver_save,     &
                                  cnd_diag_info%l_output_state,    &
                                  cnd_diag_info%l_output_incrm,    &
                                  psetcols                         )
  end do

end subroutine single_chunk_cnd_diag_alloc


!-----------------------------------------------------------------------------
! Allocate memory for metrics and diagnostics for a single sampling condition
!-----------------------------------------------------------------------------
subroutine metrics_and_qois_alloc( cnd, metric_nver, nchkpt, nqoi, qoi_nver_save, &
                                   l_output_state, l_output_incrm, psetcols )

  type(metric_and_qois_t), intent(inout) :: cnd

  integer, intent(in) :: metric_nver, nchkpt
  integer, intent(in) :: nqoi
  integer, intent(in) :: qoi_nver_save(nqoi)
  logical, intent(in) :: l_output_state
  logical, intent(in) :: l_output_incrm
  integer, intent(in) :: psetcols

  integer :: iqoi
  integer :: ierr

  character(len=*), parameter :: subname = 'metrics_and_qois_alloc'

  ! metric and flag, which might have 1, pver, or pver+1 vertical levels

  allocate( cnd% metric(psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of cnd% metric')

  allocate( cnd% flag  (psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of cnd% flag')

  ! QoIs 

  if (nqoi > 0) then

     allocate( cnd% qoi( nqoi ), stat=ierr)
     if ( ierr /= 0 ) call endrun(subname//': allocation of cnd% qoi')

     ! snapshots of each QoI 
     if (l_output_state) then
      do iqoi = 1, nqoi

        allocate( cnd%qoi(iqoi)% val(psetcols,qoi_nver_save(iqoi),nchkpt), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi% val')

        cnd%qoi(iqoi)% val(:,:,:) = 0._r8

      end do !iqoi
     end if

     ! increments of each QoI
     if (l_output_incrm) then
      do iqoi = 1, nqoi

        allocate( cnd%qoi(iqoi)% old(psetcols,qoi_nver_save(iqoi)), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi% old')

        allocate( cnd%qoi(iqoi)% inc(psetcols,qoi_nver_save(iqoi),nchkpt), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi% inc')

        cnd%qoi(iqoi)% old(:,:)   = 0._r8
        cnd%qoi(iqoi)% inc(:,:,:) = 0._r8

      end do !iqoi
     end if

   end if

end subroutine metrics_and_qois_alloc

end module conditional_diag

