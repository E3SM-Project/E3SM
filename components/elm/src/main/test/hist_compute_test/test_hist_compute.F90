program test_hist_compute

  !-----------------------------------------------------------------------
  ! Unit test of histComputeMod and subgridAveViewMod against
  ! hist_reference_mod (the original serial ELM algorithm) on a synthetic
  ! subgrid.
  !
  ! - Averages to gridcells: p2g_range, c2g_range, l2g_range and t2g_range for every
  !   gridcell and scale type must equal the serial scatter loops.
  ! - History: model arrays at every subgrid level, 1d and 2d, some with their
  !   lower bound reset to 1 and some multi-layer snow fields, are registered
  !   on three tapes (two averaged to gridcells, one on the native subgrid)
  !   with random time operations (I, A, X, M, U, L) and scale types; one
  !   model array may appear on several tapes. Each step gets new values (10%
  !   missing), and every few steps new weights, active flags and snow layers.
  !   After every update, and every normalize at the end of a tape's
  !   interval, all buffers and counters must equal the reference bit for bit.
  !
  ! With OpenACC the model arrays and the subgrid live on the device, as in
  ! ELM: the test copies new values there and the buffers back.
  !
  !   test_hist_compute          run the test (exit status 1 on a mismatch)
  !   test_hist_compute bench    time the reference and histComputeMod on a
  !                              larger subgrid
  !-----------------------------------------------------------------------

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use subgridAveViewMod , only : subgrid_ave_view, p2g_range, c2g_range, l2g_range, t2g_range
  use subgridAveViewMod , only : p2g_range_unity, c2g_range_unity, scale_unity
  use histComputeMod
  use fake_subgrid_mod
  use hist_reference_mod

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
  type(hist_active_field), allocatable :: fields(:)
  type(ref_buffer) , allocatable :: ref(:)
  type(subgrid_ave_view) :: view
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
     call subgrid_view(view)
     call setup_fields(nsrc=300, nfld=600)
     call run(nsteps=24, extrasnow=.false., check=.false.)
     write(*,'(a,i0,a,i0,a,i0,a)') 'bench: ', endg-begg+1, ' gridcells, ', endp-begp+1, &
          ' patches, ', size(fields), ' fields, 24 steps'
     write(*,'(a,f9.3,a)') '  reference (old code): ', t_ref, ' s'
     write(*,'(a,f9.3,a)') '  histComputeMod:       ', t_new, ' s'
     write(*,'(a,f9.1)')   '  speedup:              ', t_ref / t_new
     stop
  end if

  call build_subgrid(ng=37, offset=10)
  call subgrid_view(view)
  call setup_fields(nsrc=90, nfld=180)
  write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'subgrid: ', endg-begg+1, ' gridcells, ', &
       endt-begt+1, ' topounits, ', endl-begl+1, ' landunits, ', endc-begc+1, ' columns, ', &
       endp-begp+1, ' patches'
  write(*,'(a,i0,a,i0,a)') 'fields:  ', size(fields), ' on 3 tapes, from ', size(arrays), ' model arrays'

  nfail = 0
  ncase = 0
  call run(nsteps=60, extrasnow=.false., check=.true.)
  call run(nsteps=60, extrasnow=.true.,  check=.true.)
  if (nfail > 0) then
     write(*,'(a,i0,a)') 'FAIL: ', nfail, ' mismatches'
     error stop 1
  end if
  write(*,'(a,i0,a)') 'PASS: ', ncase, ' comparisons bit-for-bit with the reference'

contains

  !-----------------------------------------------------------------------
  subroutine setup_fields(nsrc, nfld)
    ! Random model arrays (on the device with OpenACC), and fields on three
    ! tapes that use them.
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
          m%a1(:) = 0._r8
          !$acc enter data copyin(m%a1)
          src1(i)%ptr => m%a1
       else
          if (m%no_snow /= no_snow_unset) then
             allocate(m%a2(lb:lb+n-1, -nlevsno+1:0))   ! as ELM's snow arrays
          else
             allocate(m%a2(lb:lb+n-1, m%num2d))
          end if
          m%a2(:,:) = 0._r8
          !$acc enter data copyin(m%a2)
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
          fields(f)%source_index = rand_int(1, nsrc)
          associate(m => arrays(fields(f)%source_index), fld => fields(f))
          fld%source_level     = m%level
          fld%source_rank      = m%rank
          fld%num_levels       = m%num2d
          fld%no_snow_behavior = m%no_snow
          fld%time_op          = rand_int(time_inst, time_min_valid)
          fld%average_to_gridcell = t < 3 .and. m%level /= lev_grid
          fld%c2l_scale = 1          ! mostly 'unity', as in ELM's field list
          fld%l2g_scale = 1
          if (urand() < 0.4_r8) fld%c2l_scale = rand_int(2, n_c2l)
          if (urand() < 0.4_r8) fld%l2g_scale = rand_int(2, n_l2g)
          if (fld%average_to_gridcell) then
             beg = begg; end = endg
          else
             beg = level_beg(fld%source_level); end = level_end(fld%source_level)
          end if
          allocate(ref(f)%hbuf(end-beg+1, fld%num_levels), ref(f)%nacs(end-beg+1, fld%num_levels))
          end associate
       end do
    end do
    tape_first(ntapes+1) = f + 1
  end subroutine setup_fields

  !-----------------------------------------------------------------------
  subroutine run(nsteps, extrasnow, check)
    integer, intent(in) :: nsteps
    logical, intent(in) :: extrasnow, check
    integer :: step, f, t, ierr
    character(len=128) :: msg
    integer(8) :: c0, c1, rate

    extrasnowlayers = extrasnow
    call hist_compute_build(view, fields, tape_first, src1, src2, extrasnow, ierr, msg)
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
       if (check) call check_averages(step)

       call system_clock(c0, rate)
       do f = 1, size(fields)
          call ref_update_field(fields(f), src1, src2, ref(f)%hbuf, ref(f)%nacs)
       end do
       call system_clock(c1)
       t_ref = t_ref + real(c1 - c0, r8) / rate

       call system_clock(c0)
       call hist_compute_update(view, src1, src2)
       call system_clock(c1)
       t_new = t_new + real(c1 - c0, r8) / rate
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
          !$acc update device(m%a1)
       else
          do j = lbound(m%a2, 2), ubound(m%a2, 2)
             do k = lbound(m%a2, 1), ubound(m%a2, 1)
                m%a2(k, j) = value(i + k + j)
             end do
          end do
          !$acc update device(m%a2)
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
  subroutine check_averages(step)
    ! The averages to gridcells against the serial scatter loops, for every
    ! scale type, on random values of each level (10% missing): over all
    ! gridcells on the host, and in chunks of gridcells in a kernel (on the
    ! device with OpenACC), as subgridAveMod calls them.
    integer, intent(in) :: step
    integer, parameter :: chunk_len = 16
    real(r8) :: pv(begp:endp), cv(begc:endc), lv(begl:endl), tv(begt:endt)
    real(r8) :: ref_avg(begg:endg), avg(begg:endg)
    integer  :: i, ic2l, il2g, nbad, chunk, nchunk, g0, g1

    do i = begp, endp ; pv(i) = value(i) ; end do
    do i = begc, endc ; cv(i) = value(i) ; end do
    do i = begl, endl ; lv(i) = value(i) ; end do
    do i = begt, endt ; tv(i) = value(i) ; end do

    nbad = 0
    do ic2l = 1, n_c2l
       do il2g = 1, n_l2g
          call ref_p2g(pv, ref_avg, ic2l, il2g)
          call p2g_range(view, begg, endg, pv, begp, ic2l, il2g, avg)
          nbad = nbad + count(avg /= ref_avg)
          call ref_c2g(cv, ref_avg, ic2l, il2g)
          call c2g_range(view, begg, endg, cv, begc, ic2l, il2g, avg)
          nbad = nbad + count(avg /= ref_avg)
       end do
    end do
    ! the unity kernels: bit for bit the general ones with the unity ids
    call ref_p2g(pv, ref_avg, scale_unity, scale_unity)
    call p2g_range_unity(view, begg, endg, pv, begp, avg)
    nbad = nbad + count(avg /= ref_avg)
    call ref_c2g(cv, ref_avg, scale_unity, scale_unity)
    call c2g_range_unity(view, begg, endg, cv, begc, avg)
    nbad = nbad + count(avg /= ref_avg)
    do il2g = 1, n_l2g
       call ref_l2g(lv, ref_avg, il2g)
       call l2g_range(view, begg, endg, lv, begl, il2g, avg)
       nbad = nbad + count(avg /= ref_avg)
    end do
    call ref_t2g(tv, ref_avg)
    call t2g_range(view, begg, endg, tv, begt, avg)
    nbad = nbad + count(avg /= ref_avg)

    ! the same averages in chunks of gridcells, one gang each
    nchunk = (endg - begg) / chunk_len + 1
    call ref_p2g(pv, ref_avg, n_c2l, n_l2g)
    !$acc parallel loop gang default(present) copyin(pv) copyout(avg) private(g0, g1)
    do chunk = 1, nchunk
       g0 = begg + (chunk-1) * chunk_len ; g1 = min(g0 + chunk_len - 1, endg)
       call p2g_range(view, g0, g1, pv, begp, n_c2l, n_l2g, avg(g0:g1))
    end do
    nbad = nbad + count(avg /= ref_avg)
    call ref_c2g(cv, ref_avg, n_c2l, n_l2g)
    !$acc parallel loop gang default(present) copyin(cv) copyout(avg) private(g0, g1)
    do chunk = 1, nchunk
       g0 = begg + (chunk-1) * chunk_len ; g1 = min(g0 + chunk_len - 1, endg)
       call c2g_range(view, g0, g1, cv, begc, n_c2l, n_l2g, avg(g0:g1))
    end do
    nbad = nbad + count(avg /= ref_avg)
    call ref_p2g(pv, ref_avg, scale_unity, scale_unity)
    !$acc parallel loop gang default(present) copyin(pv) copyout(avg) private(g0, g1)
    do chunk = 1, nchunk
       g0 = begg + (chunk-1) * chunk_len ; g1 = min(g0 + chunk_len - 1, endg)
       call p2g_range_unity(view, g0, g1, pv, begp, avg(g0:g1))
    end do
    nbad = nbad + count(avg /= ref_avg)
    call ref_c2g(cv, ref_avg, scale_unity, scale_unity)
    !$acc parallel loop gang default(present) copyin(cv) copyout(avg) private(g0, g1)
    do chunk = 1, nchunk
       g0 = begg + (chunk-1) * chunk_len ; g1 = min(g0 + chunk_len - 1, endg)
       call c2g_range_unity(view, g0, g1, cv, begc, avg(g0:g1))
    end do
    nbad = nbad + count(avg /= ref_avg)
    call ref_l2g(lv, ref_avg, n_l2g)
    !$acc parallel loop gang default(present) copyin(lv) copyout(avg) private(g0, g1)
    do chunk = 1, nchunk
       g0 = begg + (chunk-1) * chunk_len ; g1 = min(g0 + chunk_len - 1, endg)
       call l2g_range(view, g0, g1, lv, begl, n_l2g, avg(g0:g1))
    end do
    nbad = nbad + count(avg /= ref_avg)
    call ref_t2g(tv, ref_avg)
    !$acc parallel loop gang default(present) copyin(tv) copyout(avg) private(g0, g1)
    do chunk = 1, nchunk
       g0 = begg + (chunk-1) * chunk_len ; g1 = min(g0 + chunk_len - 1, endg)
       call t2g_range(view, g0, g1, tv, begt, avg(g0:g1))
    end do
    nbad = nbad + count(avg /= ref_avg)

    ncase = ncase + 1
    if (nbad > 0) then
       nfail = nfail + 1
       write(*,'(a,i0,a,i0,a)') 'averages differ from the serial loops at step ', step, ': ', &
            nbad, ' gridcell values'
    end if
  end subroutine check_averages

  !-----------------------------------------------------------------------
  subroutine compare(what, step)
    character(len=*), intent(in) :: what
    integer         , intent(in) :: step
    real(r8), pointer :: hbuf(:,:)
    integer , pointer :: nacs(:,:)
    integer :: f, t

    do t = 1, ntapes
       call hist_compute_update_host(t)
    end do
    do f = 1, size(fields)
       call hist_compute_views(f, hbuf, nacs)
       ncase = ncase + 1
       if (any(shape(hbuf) /= shape(ref(f)%hbuf))) then
          write(*,*) 'shape mismatch: field ', f
          nfail = nfail + 1
       else if (any(hbuf /= ref(f)%hbuf) .or. any(nacs /= ref(f)%nacs)) then
          nfail = nfail + 1
          if (nfail <= 10) then
             write(*,'(a,a,a,i0,a,i0,a,5i3,l2,a,i0,a,i0,a,es9.2)') 'mismatch after ', what, ' step ', step, &
                  ' field ', f, ' (level rank levels time_op snow gcell):', fields(f)%source_level, &
                  fields(f)%source_rank, fields(f)%num_levels, fields(f)%time_op, &
                  fields(f)%no_snow_behavior, fields(f)%average_to_gridcell, ' values differing: ', &
                  count(hbuf /= ref(f)%hbuf), ', counters: ', count(nacs /= ref(f)%nacs), &
                  ', max relative difference: ', &
                  maxval(abs(hbuf - ref(f)%hbuf) / max(abs(ref(f)%hbuf), tiny(1._r8)))
          end if
       end if
    end do
  end subroutine compare

end program test_hist_compute
