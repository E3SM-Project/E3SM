! ===========================================
! Module to support hybrid programming model
! hybrid_t is assumed to be a private struct
! ===========================================
module hybrid_mod

use parallel_mod  , only : parallel_t, copy_par
use thread_mod    , only : omp_set_num_threads, omp_get_thread_num 
use thread_mod    , only : horz_num_threads, vert_num_threads, tracer_num_threads
use dimensions_mod, only : nlev, qsize, ntrac

implicit none
private

  type, private :: hybrid_p
     integer :: ibeg, iend
     integer :: kbeg, kend
     integer :: qbeg, qend
  end type

  type, public :: hybrid_t
     type (parallel_t) :: par
     integer           :: ithr
     integer           :: nthreads
     integer           :: ibeg, iend
     integer           :: kbeg, kend
     integer           :: qbeg, qend
     logical           :: masterthread
  end type 

  integer, allocatable :: work_pool_horz(:,:)
  integer, allocatable :: work_pool_vert(:,:)
  integer, allocatable :: work_pool_trac(:,:)

  integer :: nelemd_save
  logical :: init_ranges = .true.
  integer :: region_num_threads
  character(len=64) :: region_name

  public :: PrintHybrid
  public :: set_region_num_threads
  private :: set_loop_ranges
  public :: get_loop_ranges
  public :: init_loop_ranges
  public :: threadOwnsTracer, threadOwnsVertlevel
  public :: config_thread_region

  interface config_thread_region 
      module procedure config_thread_region_par
      module procedure config_thread_region_hybrid
  end interface
  interface PrintHybrid 
      module procedure PrintHybridnew
  end interface

contains

  subroutine PrintHybridnew(hybt,vname)
    type (hybrid_t) :: hybt
    character(len=*) :: vname

    write(*,21) vname, hybt%par%rank, hybt%ithr, hybt%nthreads, &
                hybt%ibeg, hybt%iend,hybt%kbeg,hybt%kend, &
                hybt%qbeg, hybt%qend
21  format('PrintHybrid: (',a, ', rank: ',i8, ', ithrd: ',i4,',  nthreads: ',i4, &
           ', i{beg,end}: ',2(i4),', k{beg,end}: ',2(i4),', q{beg,end}: ',2(i4),')')

  end subroutine PrintHybridnew

  
  function config_thread_region_hybrid(old,region_name) result(new)
     type (hybrid_t), intent(in) :: old
     character(len=*), intent(in) :: region_name
     type (hybrid_t) :: new

     integer :: ithr
     integer :: kbeg_range, kend_range, qbeg_range, qend_range
     

     ithr = omp_get_thread_num()

     if ( TRIM(region_name) == 'serial') then 
         region_num_threads = 1
         new%ibeg = old%ibeg;      new%iend = old%iend
         new%kbeg = old%kbeg;      new%kend = old%kend
         new%qbeg = old%qbeg;      new%qend = old%qend
     endif
     if ( TRIM(region_name) == 'vertical') then
         region_num_threads = vert_num_threads
         call set_thread_ranges_1D ( work_pool_vert, kbeg_range, kend_range, ithr )
         new%ibeg = old%ibeg;      new%iend = old%iend
         new%kbeg = kbeg_range;    new%kend = kend_range
         new%qbeg = old%qbeg;      new%qend = old%qend
      endif

      if ( TRIM(region_name) == 'tracer' ) then
         region_num_threads = tracer_num_threads
         call set_thread_ranges_1D ( work_pool_trac, qbeg_range, qend_range, ithr)
         new%ibeg = old%ibeg;      new%iend = old%iend
         new%kbeg = old%kbeg;      new%kend = old%kend
         new%qbeg = qbeg_range;    new%qend = qend_range
      endif


      if ( TRIM(region_name) == 'vertical_and_tracer' ) then
         region_num_threads = vert_num_threads*tracer_num_threads
         call set_thread_ranges_2D ( work_pool_vert, work_pool_trac, kbeg_range, kend_range, &
                                                       qbeg_range, qend_range, ithr )
         new%ibeg = old%ibeg;      new%iend = old%iend
         new%kbeg = kbeg_range;    new%kend = kend_range
         new%qbeg = qbeg_range;    new%qend = qend_range
      endif


      new%par          = old%par      ! relies on parallel_mod copy constructor
      new%nthreads     = old%nthreads * region_num_threads
      if( region_num_threads .ne. 1 ) then 
          new%ithr         = old%ithr * region_num_threads + ithr
      else 
          new%ithr         = old%ithr
      endif
      new%masterthread = old%masterthread
!  Do we want to make this following call?
!      call omp_set_num_threads(new%nthreads)

  end function config_thread_region_hybrid

  function config_thread_region_par(par,region_name) result(hybrid)
      type (parallel_t) , intent(in) :: par
      character(len=*), intent(in) :: region_name
      type (hybrid_t)                :: hybrid
      ! local 
      integer    :: ithr
      integer    :: ibeg_range, iend_range
      integer    :: kbeg_range, kend_range
      integer    :: qbeg_range, qend_range
      integer    :: nthreads

      ithr            = omp_get_thread_num()

      if ( TRIM(region_name) == 'serial') then
         region_num_threads = 1
         if ( .NOT. allocated(work_pool_horz) ) allocate(work_pool_horz(horz_num_threads,2))
         call set_thread_ranges_1D ( work_pool_horz, ibeg_range, iend_range, ithr )
         hybrid%ibeg = 1;          hybrid%iend = nelemd_save
         hybrid%kbeg = 1;          hybrid%kend = nlev
         hybrid%qbeg = 1;          hybrid%qend = qsize
      endif

      if ( TRIM(region_name) == 'horizontal') then
         region_num_threads = horz_num_threads 
         call set_thread_ranges_1D ( work_pool_horz, ibeg_range, iend_range, ithr )
         hybrid%ibeg = ibeg_range; hybrid%iend = iend_range
         hybrid%kbeg = 1;          hybrid%kend = nlev
         hybrid%qbeg = 1;          hybrid%qend = qsize
      endif

      if ( TRIM(region_name) == 'vertical') then
         region_num_threads = vert_num_threads 
         call set_thread_ranges_1D ( work_pool_vert, kbeg_range, kend_range, ithr )
         hybrid%ibeg = 1;          hybrid%iend = nelemd_save
         hybrid%kbeg = kbeg_range; hybrid%kend = kend_range
         hybrid%qbeg = 1;          hybrid%qend = qsize
      endif
  
      if ( TRIM(region_name) == 'tracer' ) then
         region_num_threads = tracer_num_threads
         call set_thread_ranges_1D ( work_pool_trac, qbeg_range, qend_range, ithr)
         hybrid%ibeg = 1;          hybrid%iend = nelemd_save
         hybrid%kbeg = 1;          hybrid%kend = nlev
         hybrid%qbeg = qbeg_range; hybrid%qend = qend_range
      endif

    
      if ( TRIM(region_name) == 'vertical_and_tracer' ) then
         region_num_threads = vert_num_threads*tracer_num_threads
         call set_thread_ranges_2D ( work_pool_vert, work_pool_trac, kbeg_range, kend_range, &
                                                       qbeg_range, qend_range, ithr )
         hybrid%ibeg = 1;          hybrid%iend = nelemd_save
         hybrid%kbeg = kbeg_range; hybrid%kend = kend_range
         hybrid%qbeg = qbeg_range; hybrid%qend = qend_range
      endif
      call omp_set_num_threads(region_num_threads)

!      hybrid%par          = par      ! relies on parallel_mod copy constructor
      call copy_par(hybrid%par,par)
      hybrid%nthreads     = region_num_threads
      hybrid%ithr         = ithr
      hybrid%masterthread = (par%masterproc .and. ithr==0)

  end function config_thread_region_par

  subroutine init_loop_ranges(nelemd)

      integer, intent(in) :: nelemd
      integer :: ith, beg_index, end_index

      
      if ( init_ranges ) then
        nelemd_save=nelemd
!JMD#ifdef _OPENMP
        if ( .NOT. allocated(work_pool_horz) ) allocate(work_pool_horz(horz_num_threads,2))
        do ith=0,horz_num_threads-1
          call create_work_pool( 1, nelemd, horz_num_threads, ith, beg_index, end_index )
          work_pool_horz(ith+1,1) = beg_index
          work_pool_horz(ith+1,2) = end_index
        end do

        if ( .NOT. allocated(work_pool_vert) ) allocate(work_pool_vert(vert_num_threads,2))
        do ith=0,vert_num_threads-1
          call create_work_pool( 1, nlev, vert_num_threads, ith, beg_index, end_index )
          work_pool_vert(ith+1,1) = beg_index
          work_pool_vert(ith+1,2) = end_index
        end do

        if ( .NOT. allocated(work_pool_trac) ) allocate(work_pool_trac(tracer_num_threads,2))
        do ith=0,tracer_num_threads-1
          call create_work_pool( 1, qsize, tracer_num_threads, ith, beg_index, end_index )
          work_pool_trac(ith+1,1) = beg_index
          work_pool_trac(ith+1,2) = end_index
        end do

!JMD#endif
        init_ranges = .false.
      endif

  end subroutine init_loop_ranges

 subroutine set_region_num_threads( local_name )

  character(len=*), intent(in) :: local_name

  region_name = local_name

#ifdef _OPENMP

  if ( TRIM(region_name) == 'horizontal') then
    region_num_threads = horz_num_threads 
    call omp_set_num_threads(region_num_threads)
    return
  endif

  if ( TRIM(region_name) == 'vertical') then
    region_num_threads = vert_num_threads 
    call omp_set_num_threads(region_num_threads)
    return
  endif
  
  if ( TRIM(region_name) == 'tracer' ) then
    region_num_threads = tracer_num_threads
    call omp_set_num_threads(region_num_threads)
    return
  endif
    
  if ( TRIM(region_name) == 'vertical_and_tracer' ) then
    region_num_threads = vert_num_threads*tracer_num_threads
    call omp_set_num_threads(region_num_threads)
    return
  endif
 
#endif
    
  end subroutine set_region_num_threads

  subroutine set_loop_ranges (pybrid)

  type (hybrid_p) :: pybrid

  integer :: ibeg_range, iend_range
  integer :: kbeg_range, kend_range
  integer :: qbeg_range, qend_range
  integer :: idthread

#ifdef _OPENMP
  idthread = omp_get_thread_num()

  if ( TRIM(region_name) == 'horizontal' ) then
    call set_thread_ranges_1D ( work_pool_horz, ibeg_range, iend_range, idthread )
    pybrid%ibeg = ibeg_range; pybrid%iend = iend_range
    pybrid%kbeg = 1;          pybrid%kend = nlev
    pybrid%qbeg = 1;          pybrid%qend = qsize
  endif

  if ( TRIM(region_name) == 'vertical' ) then
    call set_thread_ranges_1D ( work_pool_vert, kbeg_range, kend_range, idthread )
    !FIXME: need to set ibeg, iend as well
    pybrid%kbeg = kbeg_range; pybrid%kend = kend_range
    pybrid%qbeg = 1;          pybrid%qend = qsize
  endif

  if ( TRIM(region_name) == 'tracer' ) then
    call set_thread_ranges_1D ( work_pool_trac, qbeg_range, qend_range, idthread )
    !FIXME: need to set ibeg, iend as well
    pybrid%kbeg = 1;          pybrid%kend = nlev
    pybrid%qbeg = qbeg_range; pybrid%qend = qend_range
  endif

  if ( TRIM(region_name) == 'vertical_and_tracer' ) then
    call set_thread_ranges_2D ( work_pool_vert, work_pool_trac, kbeg_range, kend_range, &
                                                       qbeg_range, qend_range, idthread )
    !FIXME: need to set ibeg, iend as well
    pybrid%kbeg = kbeg_range; pybrid%kend = kend_range
    pybrid%qbeg = qbeg_range; pybrid%qend = qend_range
  endif

#else
  call reset_loop_ranges(pybrid, region_name)
#endif 

  end subroutine set_loop_ranges

  subroutine get_loop_ranges (pybrid, ibeg, iend, kbeg, kend, qbeg, qend)

  type (hybrid_t), intent(in) :: pybrid
  integer, optional, intent(out) :: ibeg, iend, kbeg, kend, qbeg, qend

  if ( present(ibeg) ) then
    ibeg = pybrid%ibeg
  endif
  if ( present(iend) ) then
    iend = pybrid%iend
  endif
  if ( present(kbeg) ) then
    kbeg = pybrid%kbeg
  endif
  if ( present(kend) ) then
    kend = pybrid%kend
  endif
  if ( present(qbeg) ) then
    qbeg = pybrid%qbeg
  endif
  if ( present(qend) ) then
    qend = pybrid%qend
  endif

  end subroutine get_loop_ranges

  function threadOwnsVertlevel(hybrid,value) result(found) 

   type (hybrid_t), intent(in) :: hybrid
   integer, intent(in) :: value
   logical :: found

   found = .false.
   if ((value >= hybrid%kbeg) .and. (value <= hybrid%kend)) then 
      found = .true. 
   endif
 
  end function threadOwnsVertlevel

  function threadOwnsTracer(hybrid,value) result(found) 

   type (hybrid_t), intent(in) :: hybrid
   integer, intent(in) :: value
   logical :: found

   found = .false.
   if ((value >= hybrid%qbeg) .and. (value <= hybrid%qend)) then 
      found = .true. 
   endif
 
  end function threadOwnsTracer

  subroutine reset_loop_ranges (pybrid, region_name)

  type (hybrid_p)              :: pybrid
  character(len=*), intent(in) :: region_name

  if ( TRIM(region_name) == 'vertical' ) then
    pybrid%kbeg = 1; pybrid%kend = nlev
  endif

  if ( TRIM(region_name) == 'tracer' ) then
    pybrid%qbeg = 1; pybrid%qend = qsize
  endif

  if ( TRIM(region_name) == 'vertical_and_tracer' ) then
    pybrid%kbeg = 1; pybrid%kend = nlev
    pybrid%qbeg = 1; pybrid%qend = qsize
  endif

  end subroutine reset_loop_ranges 

  subroutine set_thread_ranges_3D ( work_pool_x, work_pool_y, work_pool_z, &
                       beg_range_1, end_range_1, beg_range_2, end_range_2, &
                                    beg_range_3, end_range_3, idthread )

  integer, intent (in   ) :: work_pool_x(:,:)
  integer, intent (in   ) :: work_pool_y(:,:)
  integer, intent (in   ) :: work_pool_z(:,:)
  integer, intent (inout) :: beg_range_1 
  integer, intent (inout) :: end_range_1 
  integer, intent (inout) :: beg_range_2 
  integer, intent (inout) :: end_range_2 
  integer, intent (inout) :: beg_range_3 
  integer, intent (inout) :: end_range_3 
  integer, intent (inout) :: idthread

  integer :: index(3)
  integer :: i, j, k, ind, irange, jrange, krange

  ind = 0

  krange = SIZE(work_pool_z,1)
  jrange = SIZE(work_pool_y,1)
  irange = SIZE(work_pool_x,1)
  do k = 1, krange
    do j = 1, jrange
      do i = 1, irange
        if( ind == idthread ) then
          index(1) = i
          index(2) = j
          index(3) = k
        endif
        ind = ind + 1
      enddo
    enddo
  enddo
  beg_range_1 = work_pool_x(index(1),1)
  end_range_1 = work_pool_x(index(1),2)
  beg_range_2 = work_pool_y(index(2),1)
  end_range_2 = work_pool_y(index(2),2)
  beg_range_3 = work_pool_z(index(3),1)
  end_range_3 = work_pool_z(index(3),2)

!  write(6,1000) idthread, beg_range_1, end_range_1, &
!                          beg_range_2, end_range_2, &
!                          beg_range_3, end_range_3
!  call flush(6)
1000 format( 'set_thread_ranges_3D', 7(i4) )

  end subroutine set_thread_ranges_3D

  subroutine set_thread_ranges_2D( work_pool_x, work_pool_y, beg_range_1, end_range_1, &
                                                    beg_range_2, end_range_2, idthread )

  integer, intent (in   ) :: work_pool_x(:,:)
  integer, intent (in   ) :: work_pool_y(:,:)
  integer, intent (inout) :: beg_range_1 
  integer, intent (inout) :: end_range_1 
  integer, intent (inout) :: beg_range_2 
  integer, intent (inout) :: end_range_2 
  integer, intent (inout) :: idthread

  integer :: index(2)
  integer :: i, j, ind, irange, jrange

  ind = 0

  jrange = SIZE(work_pool_y,1)
  irange = SIZE(work_pool_x,1)
  do j = 1, jrange
    do i = 1, irange
      if( ind == idthread ) then
        index(1) = i
        index(2) = j
      endif
      ind = ind + 1
    enddo
  enddo
  beg_range_1 = work_pool_x(index(1),1)
  end_range_1 = work_pool_x(index(1),2)
  beg_range_2 = work_pool_y(index(2),1)
  end_range_2 = work_pool_y(index(2),2)

! write(6,1000) idthread, beg_range_1, end_range_1, &
!                         beg_range_2, end_range_2
! call flush(6)

1000 format( 'set_thread_ranges_2D', 7(i4) )

  end subroutine set_thread_ranges_2D

  subroutine set_thread_ranges_1D( work_pool, beg_range, end_range, idthread )

  integer, intent (in   ) :: work_pool(:,:)
  integer, intent (inout) :: beg_range
  integer, intent (inout) :: end_range
  integer, intent (inout) :: idthread

  integer :: index
  integer :: i, j, ind, irange

  ind = 0
  
  irange = SIZE(work_pool) 
  do i = 1, irange
    if( ind == idthread ) then
      index = i 
    endif
    ind = ind + 1
  enddo
  beg_range = work_pool(index,1)
  end_range = work_pool(index,2)

! write(6,1000) idthread, beg_range, end_range
! call flush(6)
1000 format( 'set_thread_ranges_1D', 7(i4) )

  end subroutine set_thread_ranges_1D

  subroutine create_work_pool( start_domain, end_domain, ndomains, ipe, beg_index, end_index )

  integer, intent(in) :: start_domain, end_domain
  integer, intent(in) :: ndomains, ipe
  integer, intent(out) ::beg_index, end_index

  integer :: beg(0:ndomains)
  integer :: length
  integer :: n

  length = end_domain - start_domain + 1
  beg(0) = start_domain

  do n=1,ndomains-1
     if (n.le.mod(length,ndomains)) then
        beg(n)=beg(n-1)+(length-1)/ndomains+1
     else
        beg(n)=beg(n-1)+length/ndomains
     end if
  end do
 
  beg(ndomains) = start_domain + length

  beg_index = beg(ipe)
  end_index = beg(ipe+1) - 1

  end subroutine create_work_pool

end module hybrid_mod
