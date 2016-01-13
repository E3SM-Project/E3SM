program parallel_config_eval
  implicit none
  integer :: ewn, nsn, lhalo, uhalo, l_npe, u_npe, stride_npe
  integer :: i, maxsize, maxdiff, config_cnt

  read (5,*) ewn, nsn, lhalo, uhalo, l_npe, u_npe, stride_npe
  if (stride_npe < 1) stride_npe = 1
  do i=l_npe,u_npe,stride_npe
    write(6,*) 
    call flush(6)
    call distributed_grid(ewn,nsn,lhalo,uhalo,i,maxsize,maxdiff,config_cnt)
    if (maxsize > 1) then
       write(6,1) ewn,nsn,lhalo,uhalo,i,maxsize,maxdiff,config_cnt
    else
       write(6,2) ewn,nsn,lhalo,uhalo,i
    endif
    call flush(6)
1 format("WORKED: Grid: (", I6, ",", I6, ") Halo: (", I2, ",", I2, ") Tasks:", I6, " Max Block Size:", I6, " Max Side Diff:", I2, " Config Cnt:", I2)
2 format("FAILED: Grid: (", I6, ",", I6, ") Halo: (", I2, ",", I2, ") Tasks:", I6)
  enddo
end

!pw  subroutine distributed_grid(ewn,nsn)
!pw++
  subroutine distributed_grid(ewn,nsn,lhalo,uhalo,tasks,maxsize,maxdiff,config_cnt)
!pw--
    implicit none
    integer :: ewn,nsn

!pw++
    integer, intent(in)  :: lhalo,uhalo,tasks
    integer, intent(out) :: maxsize,maxdiff,config_cnt
    integer :: global_ewn,global_nsn
    integer :: ProcsEW, this_rank
    integer :: global_col_offset, global_row_offset
    integer :: ewlb, ewub, nslb, nsub
    integer :: local_ewn, local_nsn
    integer :: own_ewn, own_nsn
    integer :: west, east, south, north
    integer :: l_ewn(0:tasks-1), l_nsn(0:tasks-1)
!pw--

    integer :: best,i,j,metric
    integer :: ewrank,ewtasks,nsrank,nstasks
    real(8) :: rewtasks,rnstasks

    ! begin

    maxsize = 1
    maxdiff = max(ewn,nsn)
    global_ewn = ewn
    global_nsn = nsn

    ewtasks = 0
    nstasks = 0
    best = huge(best)
    do i = 1,min(tasks,global_ewn)
       j = tasks/i
       if (j<=global_nsn.and.i*j==tasks) then ! try to use all tasks
          metric = abs(i*global_nsn-j*global_ewn) ! zero if ewn/nsn == i/j
          if (metric<best) then
             best = metric
             ewtasks = i
             nstasks = j
          end if
       end if
    end do
!pw    if (ewtasks*nstasks/=tasks) call parallel_stop(__FILE__,__LINE__)
!pw++
    if (ewtasks*nstasks/=tasks) then
      write(6,*) "FAILED: Decomposition algorithm failed"
      call flush(6)
      return
    endif
!pw--

    ! Store critical value for creating global IDs.  Defines grid distribution.
    ProcsEW = ewtasks

    ! For globalID calculations determine processor's global grid index offsets
    ! sum block sizes for row blocks preceding this_rank
    ! Do not include halo offsets in global calculations
    ! (There are ProcsEW processors per row.)
!pw++
    do this_rank=0,tasks-1
!pw--
    global_col_offset = 0
    do ewrank=0,mod(this_rank, ProcsEW)-1
      rewtasks = 1/real(ewtasks,8)
      ewlb = nint(ewrank*global_ewn*rewtasks)+1
      ewub = nint((ewrank+1)*global_ewn*rewtasks)
      own_ewn = ewub-ewlb+1
      global_col_offset = global_col_offset + own_ewn
    enddo

    ! sum block sizes for column blocks preceding this_rank
    ! (Integer division required for this_rank/ProcsEW)
    global_row_offset = 0
    do nsrank=0,(this_rank/ProcsEW)-1
      rnstasks = 1/real(nstasks,8)
      nslb = nint(nsrank*global_nsn*rnstasks)+1
      nsub = nint((nsrank+1)*global_nsn*rnstasks)
      own_nsn = nsub-nslb+1
      global_row_offset = global_row_offset + own_nsn
    enddo

    ! Set local processor's grid indices, including halo offsets
    ewrank = mod(this_rank,ewtasks)
    rewtasks = 1/real(ewtasks,8)
    ewlb = nint(ewrank*global_ewn*rewtasks)+1-lhalo
    ewub = nint((ewrank+1)*global_ewn*rewtasks)+uhalo
    local_ewn = ewub-ewlb+1
    own_ewn = local_ewn-lhalo-uhalo
!pw    ewn = local_ewn

    nsrank = this_rank/ewtasks
    rnstasks = 1/real(nstasks,8)
    nslb = nint(nsrank*global_nsn*rnstasks)+1-lhalo
    nsub = nint((nsrank+1)*global_nsn*rnstasks)+uhalo
    local_nsn = nsub-nslb+1
    own_nsn = local_nsn-lhalo-uhalo
!pw    nsn = local_nsn

    west = this_rank-1
    if ((west/ewtasks<this_rank/ewtasks).or.(west<0)) west = west+ewtasks
    east = this_rank+1
    if (east/ewtasks>this_rank/ewtasks) east = east-ewtasks
    south = this_rank-ewtasks
    if (south<0) south = south+tasks
    north = this_rank+ewtasks
    if (north>=tasks) north = north-tasks

    ! Check that haven't split up the problem too much.  Idea is that do not want halos overlapping in either dimension.
    ! local_* - lhalo - uhalo is the actual number of non-halo cells on a processor.
    if ((local_nsn - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
!pw        write(*,*) "NS halos overlap on processor ", this_rank
!pw        call parallel_stop(__FILE__, __LINE__)
!pw++
!pw
        write(*,*) "FAILED: NS halos overlap (pocessor,own_nsn):", this_rank, own_nsn
        return
!pw--
    endif

    if ((local_ewn  - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
!pw        write(*,*) "EW halos overlap on processor ", this_rank
!pw        call parallel_stop(__FILE__, __LINE__)
!pw++
        write(*,*) "FAILED: EW halos overlap (processor,own_ewn):", this_rank, own_ewn
        return
!pw--
    endif
!pw++
    l_ewn(this_rank) = own_ewn
    l_nsn(this_rank) = own_nsn
    enddo
!pw--

    ! Print grid geometry
    !write(*,*) "Process ", this_rank, " Total = ", tasks, " ewtasks = ", ewtasks, " nstasks = ", nstasks
    !write(*,*) "Process ", this_rank, " ewrank = ", ewrank, " nsrank = ", nsrank
    !write(*,*) "Process ", this_rank, " l_ewn = ", local_ewn, " o_ewn = ", own_ewn
    !write(*,*) "Process ", this_rank, " l_nsn = ", local_nsn, " o_nsn = ", own_nsn
    !write(*,*) "Process ", this_rank, " ewlb = ", ewlb, " ewub = ", ewub
    !write(*,*) "Process ", this_rank, " nslb = ", nslb, " nsub = ", nsub
    !write(*,*) "Process ", this_rank, " east = ", east, " west = ", west
    !write(*,*) "Process ", this_rank, " north = ", north, " south = ", south
    !write(*,*) "Process ", this_rank, " ew_vars = ", own_ewn, " ns_vars = ", own_nsn
!pw    call distributed_print_grid(own_ewn, own_nsn)
!pw++
    call distributed_print_grid(l_ewn, l_nsn, tasks, maxsize, maxdiff, config_cnt)
!pw--
  end subroutine distributed_grid

!pw  subroutine distributed_print_grid(l_ewn,l_nsn)
!pw++
  subroutine distributed_print_grid(l_ewn,l_nsn,tasks,maxsize,maxdiff,config_cnt)
!pw--
    ! Gathers and prints the overall grid layout by processor counts.
    implicit none

!pw    integer :: l_ewn, l_nsn
!pw++
    integer, intent(in) :: l_ewn(0:tasks-1), l_nsn(0:tasks-1)
    integer, intent(in) :: tasks
    integer, intent(out) :: maxsize, maxdiff, config_cnt
!pw--
    integer :: i,j,curr_count
!pw    integer,dimension(2) :: mybounds
    integer,dimension(:,:),allocatable :: bounds

    ! begin
!pw    mybounds(1) = l_ewn
!pw    mybounds(2) = l_nsn

!pw    if (main_task) then
       allocate(bounds(2,tasks))
!pw    else
!pw       allocate(bounds(1,1))
!pw    end if
!pw    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,mpi_integer,main_rank,comm)
!pw++
    do i=1,tasks
       bounds(1,i) = l_ewn(i-1)
       bounds(2,i) = l_nsn(i-1)
    enddo
!pw--
!pw    if (main_task) then
       do i = 1,tasks
          if (bounds(1,i) .ne. -1) then
             ! total up number of processors with matching distribution
             curr_count = 1
             do j = i+1,tasks
                if ((bounds(1,i) .eq. bounds(1,j)) .and. (bounds(2,i) .eq. bounds(2,j))) then
                   ! if matching current distribution, increment counter
                   curr_count = curr_count + 1
                   bounds(1,j) = -1  ! mark so not counted later
                   bounds(2,j) = -1
                endif
             enddo
             write(*,*) "Layout(EW,NS) = ", bounds(1,i), bounds(2,i), " total procs = ", curr_count
          endif
       end do
!pw++
       maxsize = 1
       maxdiff = 0
       config_cnt = 0
       do i = 1,tasks
          if (bounds(1,i) .ne. -1) then
             if (bounds(1,i)*bounds(2,i) > maxsize) then
                maxsize = bounds(1,i)*bounds(2,i)
             endif
             if (abs(bounds(1,i)-bounds(2,i)) > maxdiff) then
                maxdiff = abs(bounds(1,i)-bounds(2,i))
             endif
             config_cnt = config_cnt + 1
          endif
       enddo
       return
!pw--
!pw    end if
    ! automatic deallocation
  end subroutine distributed_print_grid

