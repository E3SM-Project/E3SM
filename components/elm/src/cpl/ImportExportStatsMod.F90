module ImportExportStatsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Statistics (min, max, mean, standard deviation, NaN count, Inf count) of every
  ! field ELM imports from the coupler (x2l) and exports to the coupler (l2x), plus the
  ! gridcell at which the minimum and the maximum occurred.
  !
  ! The fields are discovered by walking the coupler attribute vectors, so a field added
  ! to seq_flds_x2l_fields / seq_flds_l2x_fields is picked up automatically and no code
  ! in this module or in elm_cpl_indices has to be extended.
  !
  ! The output goes to the land log file on the inst/daily/monthly/annual/all-time
  ! cadence of the budget diagnostics, and this module deliberately mirrors the shape of
  ! WaterBudgetMod / HeatBudgetMod so that the two read alike.
  !
  ! Accumulators are task-local and are not reduced until print time, which is both
  ! cheaper than the budget scheme and necessary for min/max to mean anything.
  !
  ! The accumulators are not written to the restart file, so a restart resets any
  ! partially accumulated period.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : shr_infnan_isnan, shr_infnan_isinf
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use elm_varctl     , only : iulog, do_import_export_stats
  use spmdMod        , only : masterproc
  use mct_mod        , only : mct_aVect, mct_avect_nRattr, mct_aVect_getRList2c
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ImportExportStats_Init
  public :: ImportExportStats_Accum
  public :: ImportExportStats_Print
  !
  ! !PUBLIC DATA:
  integer, parameter, public :: dir_x2l  = 1   ! coupler -> land (import)
  integer, parameter, public :: dir_l2x  = 2   ! land -> coupler (export)
  integer, parameter, public :: dir_size = dir_l2x

  !--- P for period ---
  integer, parameter :: p_inst = 1
  integer, parameter :: p_day  = 2
  integer, parameter :: p_mon  = 3
  integer, parameter :: p_ann  = 4
  integer, parameter :: p_inf  = 5

  integer, parameter :: p_size = p_inf

  character(len=8),parameter :: pname(p_size) = &
       (/&
       '    inst', &
       '   daily', &
       ' monthly', &
       '  annual', &
       'all_time'  &
       /)

  character(len=10),parameter :: dname(dir_size) = &
       (/&
       'IMPORT x2l', &
       'EXPORT l2x'  &
       /)

  integer, parameter :: name_len = 64   ! storage length for a coupler field name

  !--- accumulators for one direction; the L arrays are task-local, the G arrays hold
  !--- the reduced result and are meaningful on all tasks after a reduction ---
  type, private :: stats_type
     integer                          :: nfld = 0
     character(len=name_len), pointer :: fldname(:) => null()

     real(r8), pointer :: nsampL (:,:) => null()  ! count of finite samples
     real(r8), pointer :: fsumL  (:,:) => null()  ! sum of finite samples
     real(r8), pointer :: fsumsqL(:,:) => null()  ! sum of squares of finite samples
     real(r8), pointer :: fminL  (:,:) => null()  ! minimum over finite samples
     real(r8), pointer :: fmaxL  (:,:) => null()  ! maximum over finite samples
     real(r8), pointer :: gminL  (:,:) => null()  ! gridcell global index at fminL
     real(r8), pointer :: gmaxL  (:,:) => null()  ! gridcell global index at fmaxL
     real(r8), pointer :: nnanL  (:,:) => null()  ! count of NaN samples
     real(r8), pointer :: ninfL  (:,:) => null()  ! count of +/-Inf samples

     real(r8), pointer :: nsampG (:,:) => null()
     real(r8), pointer :: fsumG  (:,:) => null()
     real(r8), pointer :: fsumsqG(:,:) => null()
     real(r8), pointer :: fminG  (:,:) => null()
     real(r8), pointer :: fmaxG  (:,:) => null()
     real(r8), pointer :: gminG  (:,:) => null()
     real(r8), pointer :: gmaxG  (:,:) => null()
     real(r8), pointer :: nnanG  (:,:) => null()
     real(r8), pointer :: ninfG  (:,:) => null()
     real(r8), pointer :: latminG(:,:) => null()
     real(r8), pointer :: lonminG(:,:) => null()
     real(r8), pointer :: latmaxG(:,:) => null()
     real(r8), pointer :: lonmaxG(:,:) => null()
  end type stats_type

  type(stats_type) :: stats(dir_size)

  logical :: stats_active = .false.   ! .true. once Init has allocated everything

  ! The gridcell range this task owns, captured at Init. FillLoc must scan exactly the
  ! range that Accum indexed, so this is recorded once rather than re-derived.
  integer :: begg_stats = 0
  integer :: endg_stats = -1

  !----- formats -----
  character(*),parameter :: FH = "('    ',a30,4(a14),a13,2(a7),2(a32))"
  character(*),parameter :: FR = "('    ',a30,4(es14.5),f13.0,2(f7.0)," // &
                                 "2(1x,i12,' (',f7.2,',',f8.2,')'))"
  character(*),parameter :: FE = "('    ',a30,a)"

  ! Width of a table row: 4 + 30 + 4*14 + 13 + 2*7 + 2*32
  integer, parameter :: row_width = 181

  ! Left-justified so that it lines up with the field names, which are written with
  ! adjustl. A bare 'field' in an a30 edit descriptor would be right-justified.
  character(len=30),parameter :: hdr_field = 'field'

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ImportExportStats_Init(bounds, x2l, l2x)
    !
    ! !DESCRIPTION:
    ! Allocate the accumulators and record the name of every field in the two coupler
    ! attribute vectors. Must be called after the namelist has been read, so that
    ! do_import_export_stats is known, and after the attribute vectors exist.
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds ! gridcell bounds this task owns
    type(mct_aVect)  , intent(in) :: x2l    ! coupler -> land attribute vector
    type(mct_aVect)  , intent(in) :: l2x    ! land -> coupler attribute vector
    !
    ! !LOCAL VARIABLES:
    integer :: idir, k, ip, nfld
    character(len=*), parameter :: subname = 'ImportExportStats_Init'
    !-----------------------------------------------------------------------

#ifdef CPL_BYPASS
    if (do_import_export_stats) then
       call endrun(msg=subname//' ERROR: do_import_export_stats is not supported in '// &
            'CPL_BYPASS builds, where the export is compiled out and the import does '// &
            'not come from the coupler. '//errMsg(__FILE__, __LINE__))
    end if
#endif

    if (.not. do_import_export_stats) return

    begg_stats = bounds%begg
    endg_stats = bounds%endg

    do idir = 1, dir_size

       if (idir == dir_x2l) then
          nfld = mct_avect_nRattr(x2l)
       else
          nfld = mct_avect_nRattr(l2x)
       end if

       if (nfld <= 0) then
          call endrun(msg=subname//' ERROR: coupler attribute vector has no real '// &
               'attributes '//errMsg(__FILE__, __LINE__))
       end if

       stats(idir)%nfld = nfld

       allocate(stats(idir)%fldname(nfld))

       allocate(stats(idir)%nsampL (nfld,p_size))
       allocate(stats(idir)%fsumL  (nfld,p_size))
       allocate(stats(idir)%fsumsqL(nfld,p_size))
       allocate(stats(idir)%fminL  (nfld,p_size))
       allocate(stats(idir)%fmaxL  (nfld,p_size))
       allocate(stats(idir)%gminL  (nfld,p_size))
       allocate(stats(idir)%gmaxL  (nfld,p_size))
       allocate(stats(idir)%nnanL  (nfld,p_size))
       allocate(stats(idir)%ninfL  (nfld,p_size))

       allocate(stats(idir)%nsampG (nfld,p_size))
       allocate(stats(idir)%fsumG  (nfld,p_size))
       allocate(stats(idir)%fsumsqG(nfld,p_size))
       allocate(stats(idir)%fminG  (nfld,p_size))
       allocate(stats(idir)%fmaxG  (nfld,p_size))
       allocate(stats(idir)%gminG  (nfld,p_size))
       allocate(stats(idir)%gmaxG  (nfld,p_size))
       allocate(stats(idir)%nnanG  (nfld,p_size))
       allocate(stats(idir)%ninfG  (nfld,p_size))
       allocate(stats(idir)%latminG(nfld,p_size))
       allocate(stats(idir)%lonminG(nfld,p_size))
       allocate(stats(idir)%latmaxG(nfld,p_size))
       allocate(stats(idir)%lonmaxG(nfld,p_size))

       do k = 1, nfld
          if (idir == dir_x2l) then
             stats(idir)%fldname(k) = mct_aVect_getRList2c(k, x2l)
          else
             stats(idir)%fldname(k) = mct_aVect_getRList2c(k, l2x)
          end if
       end do

       do ip = 1, p_size
          call ImportExportStats_Zero(idir, ip)
       end do

    end do

    stats_active = .true.

    if (masterproc) then
       write(iulog,*)''
       write(iulog,'(a,i5,a,i5,a)') ' ImportExportStats: reporting statistics for ', &
            stats(dir_x2l)%nfld,' imported (x2l) and ', &
            stats(dir_l2x)%nfld,' exported (l2x) coupler fields'
       write(iulog,*)''
    end if

  end subroutine ImportExportStats_Init

  !-----------------------------------------------------------------------
  subroutine ImportExportStats_Zero(idir, ip)
    !
    ! !DESCRIPTION:
    ! Zero the task-local accumulators of one direction for one period.
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: idir
    integer, intent(in) :: ip
    !-----------------------------------------------------------------------

    stats(idir)%nsampL (:,ip) =  0._r8
    stats(idir)%fsumL  (:,ip) =  0._r8
    stats(idir)%fsumsqL(:,ip) =  0._r8
    stats(idir)%fminL  (:,ip) =  huge(1._r8)
    stats(idir)%fmaxL  (:,ip) = -huge(1._r8)
    stats(idir)%gminL  (:,ip) =  0._r8
    stats(idir)%gmaxL  (:,ip) =  0._r8
    stats(idir)%nnanL  (:,ip) =  0._r8
    stats(idir)%ninfL  (:,ip) =  0._r8

  end subroutine ImportExportStats_Zero

  !-----------------------------------------------------------------------
  subroutine ImportExportStats_Accum(idir, bounds, arr)
    !
    ! !DESCRIPTION:
    ! Accumulate one sample per (field, gridcell) into every period. The array is the raw
    ! coupler buffer, dimensioned (nfld, ngridcell), which is the layout of both the MCT
    ! (%rattr) and the MOAB (transposed temporary) paths.
    !
    ! Each field is reduced over gridcells once and the resulting partial statistics are
    ! then merged into each period, rather than touching all p_size columns per gridcell.
    !
    ! !USES:
    use elm_time_manager , only : get_nstep
    use GridcellType     , only : grc_pp
    !
    ! !ARGUMENTS:
    implicit none
    integer          , intent(in) :: idir       ! dir_x2l or dir_l2x
    type(bounds_type), intent(in) :: bounds     ! bounds
    real(r8)         , intent(in) :: arr(:,:)   ! (nfld, ngridcell) coupler buffer
    !
    ! !LOCAL VARIABLES:
    integer  :: k, i, g, ip, ngrid
    real(r8) :: v
    real(r8) :: cnsamp, cfsum, cfsumsq, cfmin, cfmax, cgmin, cgmax, cnnan, cninf
    character(len=*), parameter :: subname = 'ImportExportStats_Accum'
    !-----------------------------------------------------------------------

    if (.not. stats_active) return

    ! lnd_export is also called during initialization; that call must not be a sample.
    if (get_nstep() <= 0) return

    ngrid = bounds%endg - bounds%begg + 1

    if (bounds%begg /= begg_stats .or. bounds%endg /= endg_stats) then
       call endrun(msg=subname//' ERROR: gridcell bounds differ from the ones '// &
            'recorded at initialization '//errMsg(__FILE__, __LINE__))
    end if
    if (size(arr,1) /= stats(idir)%nfld) then
       call endrun(msg=subname//' ERROR: field count does not match the coupler '// &
            'attribute vector '//errMsg(__FILE__, __LINE__))
    end if
    if (size(arr,2) /= ngrid) then
       call endrun(msg=subname//' ERROR: gridcell count does not match bounds '// &
            errMsg(__FILE__, __LINE__))
    end if

    do k = 1, stats(idir)%nfld

       ! --- partial statistics for this field over this task's gridcells ---
       cnsamp  =  0._r8
       cfsum   =  0._r8
       cfsumsq =  0._r8
       cfmin   =  huge(1._r8)
       cfmax   = -huge(1._r8)
       cgmin   =  0._r8
       cgmax   =  0._r8
       cnnan   =  0._r8
       cninf   =  0._r8

       do i = 1, ngrid
          v = arr(k,i)
          if (shr_infnan_isnan(v)) then
             cnnan = cnnan + 1._r8
          else if (shr_infnan_isinf(v)) then
             cninf = cninf + 1._r8
          else
             g       = bounds%begg + i - 1
             cnsamp  = cnsamp  + 1._r8
             cfsum   = cfsum   + v
             cfsumsq = cfsumsq + v*v
             if (v < cfmin) then
                cfmin = v
                cgmin = real(grc_pp%gindex(g), r8)
             end if
             if (v > cfmax) then
                cfmax = v
                cgmax = real(grc_pp%gindex(g), r8)
             end if
          end if
       end do

       ! --- merge into every period ---
       do ip = 1, p_size
          stats(idir)%nsampL (k,ip) = stats(idir)%nsampL (k,ip) + cnsamp
          stats(idir)%fsumL  (k,ip) = stats(idir)%fsumL  (k,ip) + cfsum
          stats(idir)%fsumsqL(k,ip) = stats(idir)%fsumsqL(k,ip) + cfsumsq
          stats(idir)%nnanL  (k,ip) = stats(idir)%nnanL  (k,ip) + cnnan
          stats(idir)%ninfL  (k,ip) = stats(idir)%ninfL  (k,ip) + cninf
          if (cfmin < stats(idir)%fminL(k,ip)) then
             stats(idir)%fminL(k,ip) = cfmin
             stats(idir)%gminL(k,ip) = cgmin
          end if
          if (cfmax > stats(idir)%fmaxL(k,ip)) then
             stats(idir)%fmaxL(k,ip) = cfmax
             stats(idir)%gmaxL(k,ip) = cgmax
          end if
       end do

    end do

  end subroutine ImportExportStats_Accum

  !-----------------------------------------------------------------------
  subroutine ImportExportStats_Print(stats_print_inst, stats_print_daily, &
       stats_print_month, stats_print_ann, stats_print_ltann)
    !
    ! !DESCRIPTION:
    ! Reduce the task-local accumulators across all tasks, write the tables that are due
    ! this time step, then zero the periods whose boundary has just been crossed.
    !
    ! Every task must reach ImportExportStats_Reduce below, including tasks on which
    ! nothing is printed, because its reductions are collectives. The print levels come
    ! from the namelist and the dates from the model clock, both of which are identical
    ! on every task, so the control flow here is too.
    !
    ! !USES:
    use elm_time_manager , only : get_curr_date, get_prev_date, get_nstep
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: stats_print_inst
    integer, intent(in) :: stats_print_daily
    integer, intent(in) :: stats_print_month
    integer, intent(in) :: stats_print_ann
    integer, intent(in) :: stats_print_ltann
    !
    ! !LOCAL VARIABLES:
    integer :: idir, ip
    integer :: plev                     ! print level for this period
    integer :: year, mon, day, sec
    integer :: cdate
    logical :: at_bound(p_size)         ! this period's boundary was crossed this step
    logical :: do_print(p_size)         ! a table is due for this period
    logical :: any_print
    !-----------------------------------------------------------------------

    if (.not. stats_active) return

    if (get_nstep() <= 1) then
       call get_prev_date(year, mon, day, sec)
    else
       call get_curr_date(year, mon, day, sec)
    end if

    cdate = year*10000 + mon*100 + day

    ! A boundary has been crossed when the clock has just rolled over into the next day,
    ! month or year. p_inst turns over every step; p_inf never does.
    at_bound(:)      = .false.
    at_bound(p_inst) = .true.
    at_bound(p_day)  = (sec == 0)
    at_bound(p_mon)  = (sec == 0 .and. day == 1)
    at_bound(p_ann)  = (sec == 0 .and. day == 1 .and. mon == 1)
    at_bound(p_inf)  = .false.

    do ip = 1, p_size
       plev = 0
       if (ip == p_inst)                       plev = max(plev, stats_print_inst)
       if (ip == p_day  .and. at_bound(p_day)) plev = max(plev, stats_print_daily)
       if (ip == p_mon  .and. at_bound(p_mon)) plev = max(plev, stats_print_month)
       if (ip == p_ann  .and. at_bound(p_ann)) plev = max(plev, stats_print_ann)
       if (ip == p_inf  .and. at_bound(p_ann)) plev = max(plev, stats_print_ltann)

       do_print(ip) = (plev > 0)

       ! Nothing has accumulated into the longer periods yet on the first step.
       if (ip /= p_inst .and. get_nstep() == 1) do_print(ip) = .false.
    end do

    any_print = any(do_print(:))

    do idir = 1, dir_size

       ! Collective: entered by every task, and only when something will be written, so
       ! that the reduction count stays the same on every task.
       if (any_print) call ImportExportStats_Reduce(idir)

       if (masterproc) then
          do ip = 1, p_size
             if (do_print(ip)) call ImportExportStats_Write(idir, ip, cdate, sec)
          end do
       end if

       ! Zero the periods that have just turned over, whether or not they were printed,
       ! so that a period left at print level 0 does not accumulate for the whole run.
       do ip = 1, p_size
          if (at_bound(ip)) call ImportExportStats_Zero(idir, ip)
       end do

    end do

  end subroutine ImportExportStats_Print

  !-----------------------------------------------------------------------
  subroutine ImportExportStats_Reduce(idir)
    !
    ! !DESCRIPTION:
    ! Reduce one direction's task-local accumulators across all tasks into its G arrays.
    ! Collective: must be called by every task.
    !
    ! shr_mpi_min and shr_mpi_max have rank-0 and rank-1 specifics only, so the
    ! (nfld, p_size) arrays are packed into rank-1 buffers. Element (k,ip) of a packed
    ! array lives at (ip-1)*nfld + k, which ImportExportStats_FillLoc relies on.
    !
    ! !USES:
    use spmdMod     , only : mpicom
    use shr_mpi_mod , only : shr_mpi_sum, shr_mpi_min, shr_mpi_max
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: idir
    !
    ! !LOCAL VARIABLES:
    integer                :: nfld, nd
    real(r8), allocatable  :: lbuf(:), gbuf(:)
    real(r8), allocatable  :: cand(:,:), loc(:,:)
    character(*),parameter :: subName = '(ImportExportStats_Reduce)'
    !-----------------------------------------------------------------------

    nfld = stats(idir)%nfld
    nd   = nfld*p_size

    ! --- (1) counts and moments, in one reduction ---
    allocate(lbuf(5*nd), gbuf(5*nd))
    lbuf(     1:  nd) = reshape(stats(idir)%nsampL , (/nd/))
    lbuf(  nd+1:2*nd) = reshape(stats(idir)%fsumL  , (/nd/))
    lbuf(2*nd+1:3*nd) = reshape(stats(idir)%fsumsqL, (/nd/))
    lbuf(3*nd+1:4*nd) = reshape(stats(idir)%nnanL  , (/nd/))
    lbuf(4*nd+1:5*nd) = reshape(stats(idir)%ninfL  , (/nd/))
    call shr_mpi_sum(lbuf, gbuf, mpicom, subName, all=.true.)
    stats(idir)%nsampG  = reshape(gbuf(     1:  nd), (/nfld,p_size/))
    stats(idir)%fsumG   = reshape(gbuf(  nd+1:2*nd), (/nfld,p_size/))
    stats(idir)%fsumsqG = reshape(gbuf(2*nd+1:3*nd), (/nfld,p_size/))
    stats(idir)%nnanG   = reshape(gbuf(3*nd+1:4*nd), (/nfld,p_size/))
    stats(idir)%ninfG   = reshape(gbuf(4*nd+1:5*nd), (/nfld,p_size/))
    deallocate(lbuf, gbuf)

    allocate(lbuf(nd), gbuf(nd), cand(nfld,p_size))

    ! --- (2) global minimum ---
    lbuf = reshape(stats(idir)%fminL, (/nd/))
    call shr_mpi_min(lbuf, gbuf, mpicom, subName, all=.true.)
    stats(idir)%fminG = reshape(gbuf, (/nfld,p_size/))

    ! --- (3) global maximum ---
    lbuf = reshape(stats(idir)%fmaxL, (/nd/))
    call shr_mpi_max(lbuf, gbuf, mpicom, subName, all=.true.)
    stats(idir)%fmaxG = reshape(gbuf, (/nfld,p_size/))

    ! --- (4) gridcell global index of the minimum. Tasks that do not hold the global
    !         minimum contribute huge(), so the reduction picks the lowest gindex among
    !         the tasks that do, which makes the winner unique.
    !         The equality test on reals is deliberate and exact: MPI_MIN returns one of
    !         the contributed values bit for bit, so the owning task compares equal.
    where (stats(idir)%fminL == stats(idir)%fminG)
       cand = stats(idir)%gminL
    elsewhere
       cand = huge(1._r8)
    end where
    lbuf = reshape(cand, (/nd/))
    call shr_mpi_min(lbuf, gbuf, mpicom, subName, all=.true.)
    stats(idir)%gminG = reshape(gbuf, (/nfld,p_size/))
    where (stats(idir)%gminG >= huge(1._r8)) stats(idir)%gminG = 0._r8

    ! --- (5) and of the maximum ---
    where (stats(idir)%fmaxL == stats(idir)%fmaxG)
       cand = stats(idir)%gmaxL
    elsewhere
       cand = huge(1._r8)
    end where
    lbuf = reshape(cand, (/nd/))
    call shr_mpi_min(lbuf, gbuf, mpicom, subName, all=.true.)
    stats(idir)%gmaxG = reshape(gbuf, (/nfld,p_size/))
    where (stats(idir)%gmaxG >= huge(1._r8)) stats(idir)%gmaxG = 0._r8

    deallocate(lbuf, gbuf, cand)

    ! --- (6) latitude and longitude of those two gridcells. A global index belongs to
    !         exactly one task, so exactly one task contributes and the sum recovers it.
    allocate(loc(nd,4))
    call ImportExportStats_FillLoc(idir, nfld, begg_stats, endg_stats, loc)

    allocate(lbuf(4*nd), gbuf(4*nd))
    lbuf(     1:  nd) = loc(:,1)
    lbuf(  nd+1:2*nd) = loc(:,2)
    lbuf(2*nd+1:3*nd) = loc(:,3)
    lbuf(3*nd+1:4*nd) = loc(:,4)
    call shr_mpi_sum(lbuf, gbuf, mpicom, subName, all=.true.)
    stats(idir)%latminG = reshape(gbuf(     1:  nd), (/nfld,p_size/))
    stats(idir)%lonminG = reshape(gbuf(  nd+1:2*nd), (/nfld,p_size/))
    stats(idir)%latmaxG = reshape(gbuf(2*nd+1:3*nd), (/nfld,p_size/))
    stats(idir)%lonmaxG = reshape(gbuf(3*nd+1:4*nd), (/nfld,p_size/))
    deallocate(loc, lbuf, gbuf)

  end subroutine ImportExportStats_Reduce

  !-----------------------------------------------------------------------
  subroutine ImportExportStats_FillLoc(idir, nfld, begg, endg, loc)
    !
    ! !DESCRIPTION:
    ! For each (field, period), write this task's latitude and longitude into loc if this
    ! task owns the gridcell holding the global extremum, and zero otherwise.
    !
    ! !USES:
    use GridcellType , only : grc_pp
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: idir
    integer , intent(in)  :: nfld
    integer , intent(in)  :: begg, endg
    real(r8), intent(out) :: loc(:,:)   ! (nfld*p_size, 4): latmin, lonmin, latmax, lonmax
    !
    ! !LOCAL VARIABLES:
    integer :: k, ip, n, g, gwin
    !-----------------------------------------------------------------------

    loc(:,:) = 0._r8

    do ip = 1, p_size
       do k = 1, nfld
          n = (ip-1)*nfld + k

          gwin = nint(stats(idir)%gminG(k,ip))
          if (gwin > 0 .and. nint(stats(idir)%gminL(k,ip)) == gwin) then
             do g = begg, endg
                if (grc_pp%gindex(g) == gwin) then
                   loc(n,1) = grc_pp%latdeg(g)
                   loc(n,2) = grc_pp%londeg(g)
                   exit
                end if
             end do
          end if

          gwin = nint(stats(idir)%gmaxG(k,ip))
          if (gwin > 0 .and. nint(stats(idir)%gmaxL(k,ip)) == gwin) then
             do g = begg, endg
                if (grc_pp%gindex(g) == gwin) then
                   loc(n,3) = grc_pp%latdeg(g)
                   loc(n,4) = grc_pp%londeg(g)
                   exit
                end if
             end do
          end if

       end do
    end do

  end subroutine ImportExportStats_FillLoc

  !-----------------------------------------------------------------------
  subroutine ImportExportStats_Write(idir, ip, cdate, sec)
    !
    ! !DESCRIPTION:
    ! Write the table for one direction and one period. Called on masterproc only; the
    ! reduced G arrays must already be current.
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: idir
    integer, intent(in) :: ip
    integer, intent(in) :: cdate
    integer, intent(in) :: sec
    !
    ! !LOCAL VARIABLES:
    integer  :: k
    real(r8) :: mean, var, std
    !-----------------------------------------------------------------------

    write(iulog,*)''
    write(iulog,'(a,a,a,a,a,i8.8,i6)') ' ELM ',trim(dname(idir)), &
         ' FIELD STATS : period ',trim(pname(ip)),': date = ',cdate,sec
    write(iulog,FH) hdr_field, '  min', '  max', '  mean', '  stddev', &
         '  nsamp', '  nan', '  inf', &
         '   argmin (g, lat, lon)', '   argmax (g, lat, lon)'
    write(iulog,'(4x,a)') repeat('-',row_width-4)

    do k = 1, stats(idir)%nfld
       if (stats(idir)%nsampG(k,ip) > 0._r8) then
          mean = stats(idir)%fsumG(k,ip) / stats(idir)%nsampG(k,ip)
          ! One-pass variance. The max() guards the cancellation that occurs when a
          ! field is constant; it is adequate in double precision for a diagnostic.
          var  = max(0._r8, stats(idir)%fsumsqG(k,ip)/stats(idir)%nsampG(k,ip) - mean*mean)
          std  = sqrt(var)
          write(iulog,FR) adjustl(stats(idir)%fldname(k)), &
               stats(idir)%fminG(k,ip), stats(idir)%fmaxG(k,ip), mean, std, &
               stats(idir)%nsampG(k,ip), stats(idir)%nnanG(k,ip), stats(idir)%ninfG(k,ip), &
               nint(stats(idir)%gminG(k,ip)), &
               stats(idir)%latminG(k,ip), stats(idir)%lonminG(k,ip), &
               nint(stats(idir)%gmaxG(k,ip)), &
               stats(idir)%latmaxG(k,ip), stats(idir)%lonmaxG(k,ip)
       else
          write(iulog,FE) adjustl(stats(idir)%fldname(k)), &
               '    no finite samples'
       end if
    end do

    write(iulog,'(4x,a)') repeat('-',row_width-4)
    write(iulog,*)''

  end subroutine ImportExportStats_Write

end module ImportExportStatsMod
