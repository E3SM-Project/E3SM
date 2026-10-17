program test_hist_compute

  !-----------------------------------------------------------------------
  ! Unit test of histComputeMod against hist_reference_mod (the history
  ! algorithm before histComputeMod) on a synthetic subgrid.
  !
  ! Model arrays at every subgrid level, 1d and 2d, some with their lower
  ! bound reset to 1 and some multi-layer snow fields, are registered on three
  ! tapes (two averaged to gridcells, one on the native subgrid) with random
  ! averaging flags and scale types; one model array may appear on several
  ! tapes. Each step gets new values (10% missing), and every few steps new
  ! weights, active flags and snow layers. After every update, and every
  ! normalize at the end of a tape's interval, all buffers and counters must
  ! equal the reference bit for bit.
  !
  !   test_hist_compute          run the test (exit status 1 on a mismatch)
  !   test_hist_compute bench    time the reference and histComputeMod on a
  !                              larger subgrid
  !-----------------------------------------------------------------------

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use histComputeMod
  use fake_subgrid_mod
  use hist_reference_mod
#ifdef _OPENMP
  use omp_lib           , only : omp_get_max_threads
#endif

  implicit none

  integer, parameter :: ntapes = 3

  type model_array
     real(r8), allocatable :: a1(:), a2(:,:)
     integer :: level, rank, num2d, no_snow
  end type model_array

  type ref_buffer
     real(r8), allocatable :: hbuf(:,:)
     integer , allocatable :: nacs(:,:)
  end type ref_buffer

  type(model_array), allocatable, target :: arrays(:)
  type(hist_src1d) , allocatable :: src1(:)
  type(hist_src2d) , allocatable :: src2(:)
  type(hist_field_desc), allocatable :: fields(:)
  type(ref_buffer) , allocatable :: ref(:)
  integer :: tape_first(ntapes+1)
  integer :: seed_size, nfail, ncase
  integer, allocatable :: seed(:)
  character(len=16) :: arg
  real(r8) :: t_ref, t_new
  logical  :: bench = .false.   ! timing run: fixed missing values, fixed subgrid

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed(:) = 20261017
  call random_seed(put=seed)

  arg = ' '
  if (command_argument_count() > 0) call get_command_argument(1, arg)

  if (trim(arg) == 'bench') then
     bench = .true.
     call build_subgrid(ng=4000, offset=1000)
     call setup_fields(nsrc=300, nfld=600)
     call run(nsteps=24, extrasnow=.false., tile=256, check=.false.)
     write(*,'(a,i0,a,i0,a,i0,a)') 'bench: ', endg-begg+1, ' gridcells, ', endp-begp+1, &
          ' patches, ', size(fields), ' fields, 24 steps'
#ifdef _OPENMP
     write(*,'(a,i0)') '  OpenMP threads:        ', omp_get_max_threads()
#endif
     write(*,'(a,f9.3,a)') '  reference (old code): ', t_ref, ' s'
     write(*,'(a,f9.3,a)') '  histComputeMod:       ', t_new, ' s'
     write(*,'(a,f9.1)')   '  speedup:              ', t_ref / t_new
     stop
  end if

  call build_subgrid(ng=37, offset=10)
  call setup_fields(nsrc=90, nfld=180)
  write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'subgrid: ', endg-begg+1, ' gridcells, ', &
       endt-begt+1, ' topounits, ', endl-begl+1, ' landunits, ', endc-begc+1, ' columns, ', &
       endp-begp+1, ' patches'
  write(*,'(a,i0,a,i0,a)') 'fields:  ', size(fields), ' on 3 tapes, from ', size(arrays), ' model arrays'

  nfail = 0
  ncase = 0
  call run(nsteps=60, extrasnow=.false., tile=256, check=.true.)
  call run(nsteps=60, extrasnow=.true.,  tile=3,   check=.true.)
  if (nfail > 0) then
     write(*,'(a,i0,a)') 'FAIL: ', nfail, ' mismatches'
     error stop 1
  end if
  write(*,'(a,i0,a)') 'PASS: ', ncase, ' comparisons bit-for-bit with the reference'

contains

  !-----------------------------------------------------------------------
  subroutine setup_fields(nsrc, nfld)
    ! Random model arrays, and fields on three tapes that use them.
    integer, intent(in) :: nsrc, nfld
    integer :: i, f, t, beg, end, n, lb

    allocate(arrays(nsrc), src1(nsrc), src2(nsrc))
    do i = 1, nsrc
       associate(m => arrays(i))
       m%level   = rand_int(lev_grid, lev_patch)
       m%rank    = rand_int(1, 2)
       m%num2d   = 1
       m%no_snow = no_snow_unset
       if (m%rank == 2) m%num2d = rand_int(1, 7)
       if (m%rank == 2 .and. m%level >= lev_col .and. urand() < 0.4_r8) then
          m%num2d   = nlevsno
          m%no_snow = rand_int(no_snow_normal, no_snow_zero)
       end if
       beg = level_beg(m%level); end = level_end(m%level); n = end - beg + 1
       lb = beg
       if (urand() < 0.25_r8) lb = 1      ! pointer to an array section
       if (m%rank == 1) then
          allocate(m%a1(lb:lb+n-1))
          src1(i)%ptr => m%a1
       else
          if (m%no_snow /= no_snow_unset) then
             allocate(m%a2(lb:lb+n-1, -nlevsno+1:0))   ! as ELM's snow arrays
          else
             allocate(m%a2(lb:lb+n-1, m%num2d))
          end if
          src2(i)%ptr => m%a2
       end if
       end associate
    end do

    ! tapes 1 and 2 average to gridcells (hist_dov2xy), tape 3 doesn't
    allocate(fields(nfld), ref(nfld))
    f = 0
    do t = 1, ntapes
       tape_first(t) = f + 1
       do i = 1, nfld / ntapes
          f = f + 1
          fields(f)%hpindex = rand_int(1, nsrc)
          associate(m => arrays(fields(f)%hpindex), fd => fields(f))
          fd%level   = m%level
          fd%rank    = m%rank
          fd%num2d   = m%num2d
          fd%no_snow = m%no_snow
          fd%avg     = rand_int(avg_inst, avg_min)
          fd%to_gridcell = t < 3
          fd%c2l     = 1          ! mostly 'unity', as in ELM's field list
          fd%l2g     = 1
          if (urand() < 0.4_r8) fd%c2l = rand_int(2, n_c2l)
          if (urand() < 0.4_r8) fd%l2g = rand_int(2, n_l2g)
          if (fd%to_gridcell .and. fd%level /= lev_grid) then
             beg = begg; end = endg
          else
             beg = level_beg(fd%level); end = level_end(fd%level)
          end if
          allocate(ref(f)%hbuf(end-beg+1, fd%num2d), ref(f)%nacs(end-beg+1, fd%num2d))
          end associate
       end do
    end do
    tape_first(ntapes+1) = f + 1
  end subroutine setup_fields

  !-----------------------------------------------------------------------
  subroutine run(nsteps, extrasnow, tile, check)
    integer, intent(in) :: nsteps, tile
    logical, intent(in) :: extrasnow, check
    type(hist_subgrid_desc)  :: sg
    type(hist_subgrid_state) :: st
    integer :: step, f, t, ierr, bad
    character(len=128) :: msg
    integer(8) :: c0, c1, rate

    extrasnowlayers = extrasnow
    call subgrid_desc(sg)
    call subgrid_state(st)
    call hist_compute_set_subgrid(sg, ierr, msg)
    if (ierr /= 0) error stop 'set_subgrid failed'
    call hist_compute_build(fields, tape_first, src1, src2, spval, ierr, msg, tile=tile)
    if (ierr /= 0) then
       write(*,*) 'build failed for field ', ierr, ': ', trim(msg)
       error stop 1
    end if
    do f = 1, size(fields)
       ref(f)%hbuf(:,:) = 0._r8
       ref(f)%nacs(:,:) = 0
    end do

    t_ref = 0._r8; t_new = 0._r8
    do step = 1, nsteps
       call new_values()
       if (mod(step, 7) == 0 .and. .not. bench) call perturb_subgrid(0.2_r8)

       call system_clock(c0, rate)
       do f = 1, size(fields)
          call ref_update_field(fields(f), src1, src2, ref(f)%hbuf, ref(f)%nacs)
       end do
       call system_clock(c1)
       t_ref = t_ref + real(c1 - c0, r8) / rate

       call system_clock(c0)
       call hist_compute_update(src1, src2, st, bad)
       call system_clock(c1)
       t_new = t_new + real(c1 - c0, r8) / rate
       if (bad /= 0) error stop 'gridcell weights sum to more than 1'
       if (check) call compare('update', step)

       ! end of each tape's interval: normalize, compare, zero
       do t = 1, ntapes
          if (mod(step, 4 + 3*t) /= 0) cycle
          call hist_compute_normalize(t)
          do f = tape_first(t), tape_first(t+1) - 1
             call ref_normalize(fields(f), ref(f)%hbuf, ref(f)%nacs)
          end do
          if (check) call compare('normalize', step)
          call hist_compute_zero(t)
          do f = tape_first(t), tape_first(t+1) - 1
             ref(f)%hbuf(:,:) = 0._r8
             ref(f)%nacs(:,:) = 0
          end do
       end do
    end do
  end subroutine run

  !-----------------------------------------------------------------------
  subroutine new_values()
    ! New values of all model arrays, 10% missing.
    integer :: i, k, j

    do i = 1, size(arrays)
       associate(m => arrays(i))
       if (m%rank == 1) then
          do k = lbound(m%a1, 1), ubound(m%a1, 1)
             m%a1(k) = value(i + k)
          end do
       else
          do j = lbound(m%a2, 2), ubound(m%a2, 2)
             do k = lbound(m%a2, 1), ubound(m%a2, 1)
                m%a2(k, j) = value(i + k + j)
             end do
          end do
       end if
       end associate
    end do
  end subroutine new_values

  real(r8) function value(key)
    ! 10% missing: random each step (test), or at fixed points (bench, as in
    ! a model, where missing values mark points without the quantity)
    integer, intent(in) :: key
    logical :: missing
    if (bench) then
       missing = mod(key * 7919, 10) == 0
    else
       missing = urand() < 0.1_r8
    end if
    if (missing) then
       value = spval
    else
       value = 600._r8 * (urand() - 0.3_r8)
    end if
  end function value

  !-----------------------------------------------------------------------
  subroutine compare(what, step)
    character(len=*), intent(in) :: what
    integer         , intent(in) :: step
    real(r8), pointer :: hbuf(:,:)
    integer , pointer :: nacs(:,:)
    integer :: f

    do f = 1, size(fields)
       call hist_compute_views(f, hbuf, nacs)
       ncase = ncase + 1
       if (any(shape(hbuf) /= shape(ref(f)%hbuf))) then
          write(*,*) 'shape mismatch: field ', f
          nfail = nfail + 1
       else if (any(hbuf /= ref(f)%hbuf) .or. any(nacs /= ref(f)%nacs)) then
          nfail = nfail + 1
          if (nfail <= 10) then
             write(*,'(a,a,a,i0,a,i0,a,5i3,l2,a,i0)') 'mismatch after ', what, ' step ', step, &
                  ' field ', f, ' (level rank num2d avg snow gcell):', fields(f)%level, &
                  fields(f)%rank, fields(f)%num2d, fields(f)%avg, fields(f)%no_snow, &
                  fields(f)%to_gridcell, ' values differing: ', &
                  count(hbuf /= ref(f)%hbuf) + count(nacs /= ref(f)%nacs)
          end if
       end if
    end do
  end subroutine compare

end program test_hist_compute
