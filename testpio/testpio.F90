#ifdef BGP
#define BGx
#endif
#ifdef BGL
#define BGx
#endif
!>
!! @file testpio.F90
!! An example of how PIO can be used
!<
program testpio

  ! Modules from PIO package that are used by this application

  use kinds_mod
  !>
  !! Use PIO methods and data structures
  !<
  use pio             ! _EXTERNAL
  use utils_mod
#ifdef TIMING
  use perf_mod        ! _EXTERNAL
#endif

  use pio_support, only : CheckMPIReturn ! _EXTERNAL
  use pio_support, only : piodie ! _EXTERNAL
  use pio_support, only : pio_writedof,pio_readdof ! _EXTERNAL

  ! Modules from testpio suite that are used by this application

  use gdecomp_mod, only: gdecomp_type, gdecomp_DOF, gdecomp_read_nml
  use alloc_mod       ! _EXTERNAL
  use check_mod
  use namelist_mod

  implicit none

  include 'mpif.h'    ! _EXTERNAL

  ! Code name, used in debug prints and passed to called routines for flow tracing
  character(len=*), parameter :: myname='testpio'

  integer(i4) :: my_task, ierr
  integer(i4) :: iostat
  integer(i4) :: indx
  integer(i4) :: mode

  integer(i4) :: ip,numPhases 
  character(len=*), parameter :: TestR8CaseName = 'r8_test'
  character(len=*), parameter :: TestR4CaseName = 'r4_test'
  character(len=*), parameter :: TestI4CaseName  = 'i4_test'
  character(len=*), parameter :: TestComboCaseName = 'combo_test'

  type (iosystem_desc_t) :: PIOSYS
  type (File_desc_t)  :: File, File_r8,File_r4,File_i4
  type (Var_desc_t)   :: vard_i4, &
       vard_r8c,vard_r4c,vard_i4c, &
       vard_i4i,vard_i4j,vard_i4k,vard_i4m,vard_i4dof
  type(var_desc_t), pointer :: vard_r8(:), vard_r4(:)

  type (IO_desc_t)    :: IOdesc_r8,IOdesc_r4,IOdesc_i4
  type (gdecomp_type) :: gdecomp

  real(r8), parameter :: MBYTES = 1.0e-6

  integer(i4) :: cbad, ivar
  integer(i4) :: i,j,is,ie,itmp,it,n,i1,j1,k1
  integer(i4) :: gDims3D(3)

  character(6) :: ew_type,ns_type
  character(len=10) :: varname

  integer(i4) :: varid,dimid_x,dimid_y,dimid_z

  integer(kind=PIO_OFFSET) :: one

  integer, parameter :: ntest = 5
  integer(i4), dimension(ntest) :: num_agg
  data num_agg/ 8,12,16,24,32/

  integer(i4),pointer :: test_i4wr(:),test_i4rd(:),diff_i4(:)
  integer(i4),pointer :: test_i4i(:),test_i4j(:),test_i4k(:),test_i4m(:),test_i4dof(:)
  real(r4),   pointer :: test_r4wr(:),test_r4rd(:),diff_r4(:)
  real(r8),   pointer :: test_r8wr(:),test_r8rd(:),diff_r8(:)

  logical :: TestR8    = .true.,   &
       TestR4    = .false.,  &
       TestInt   = .false.,  &
       TestCombo = .false.

  logical :: writePhase, readPhase
  logical, parameter :: splitPhase = .true.
  integer :: numPhase

  real(r8) :: lsum,lsum2,gsum
  real(r8) :: st,et  ! start/end times for timing
  real(r8) :: dt_write_r8, dt_write_r4, dt_write_i4 ! individual write times
  real(r8) :: dt_read_r8, dt_read_r4, dt_read_i4 ! individual read times
  ! Arrays to hold globally reduced read/write times--one element per time trial
  real(r8), dimension(:), pointer :: gdt_write_r8, gdt_write_r4, gdt_write_i4 
  real(r8), dimension(:), pointer :: gdt_read_r8, gdt_read_r4, gdt_read_i4

  integer(i4) :: nprocs
  integer(i4) :: lLength    ! local number of words in the computational decomposition 

  integer(i4), parameter ::  nml_in = 10
  character(len=*), parameter :: nml_filename = 'testpio_in'

  integer(i4)  :: master_task
  logical      :: log_master_task
  integer(i4)  :: nml_error
  integer(i4)  :: sdof,sdof_sum,sdof_min,sdof_max

  ! memory tracking stuff
  integer(i4)  :: msize,mrss,mrss0,mrss1,mrss2
  real(r8)     :: mb_blk
  real(r8),allocatable :: mem_tmp(:)
  integer(i4),allocatable :: lmem(:),gmem(:,:)

  integer(i4), pointer  :: compDOF(:), ioDOF(:)
  integer(i4), pointer  :: ioDOFR(:),ioDOFW(:)

  integer(i4) :: startIO(3),countIO(3), &
       startCOMP(3), countCOMP(3), &
       start(3), count(3)
  integer(i4) :: lenblocks, glenr8, glenr4, gleni4
  integer(kind=PIO_OFFSET) :: startpio(3), countpio(3)

  character(len=80) :: fname, fname_r8,fname_r4,fname_i4
  logical, parameter :: Debug = .FALSE.
  integer :: mpi_comm_compute, mpi_comm_io, mpi_icomm_cio

  character(len=3) :: citer

  ! Initialize MPI

  call MPI_INIT(ierr)
  call CheckMPIReturn('Call to MPI_INIT()',ierr,__FILE__,__LINE__)
  

  !    call enable_abort_on_exit

  call MPI_COMM_RANK(MPI_COMM_WORLD,my_task,ierr)
  call CheckMPIReturn('Call to MPI_COMM_RANK()',ierr,__FILE__,__LINE__)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  call CheckMPIReturn('Call to MPI_COMM_SIZE()',ierr,__FILE__,__LINE__)
  master_task = 0
  if (my_task == master_task) then
     log_master_task = .true.
  else
     log_master_task = .false.
  endif

  if(Debug)    print *,'testpio: before call to t_initf'
#ifdef TIMING
  !---------------------------------------------------------------
  ! timing library
  !---------------------------------------------------------------
  if(Debug)    print *,'testpio: point #1'
  call t_initf(nml_filename, logprint=.false., logunit=6, &
       mpicom=MPI_COMM_WORLD, MasterTask=log_master_task)
  if(Debug)    print *,'testpio: point #2'
  call t_startf('testpio_total')

  if(Debug)    print *,'testpio: after call to t_startf'
  !-----------------------------------------------------------------------------
  ! Memory test
  !-----------------------------------------------------------------------------
  call get_memusage(msize,mrss0)
  if(Debug)    print *,'testpio: after get_memusage #3'
  allocate(mem_tmp(1024*1024))    ! 1 MWord, 8 MB
  mem_tmp = -1.0
  call get_memusage(msize,mrss1)
  if(Debug)    print *,'testpio: after get_memusage #4'
  deallocate(mem_tmp)
  call get_memusage(msize,mrss2)
  if(Debug)    print *,'testpio: after get_memusage #5'
  mb_blk = 0.0_r8
  if(Debug)    print *,'testpio: after get_memusage #6'
  if (mrss1 - mrss0 > 0) then
     mb_blk = (8.0_r8)/((mrss1-mrss0)*1.0_r8)
  endif
  if (my_task == master_task) then
     write(*,*) myname,' 8 MB memory   alloc in MB is ',(mrss1-mrss0)*mb_blk
     write(*,*) myname,' 8 MB memory dealloc in MB is ',(mrss1-mrss2)*mb_blk
     write(*,*) myname,' Memory block size conversion in bytes is ',mb_blk*1024_r8*1024.0_r8
  endif

#endif

  !----------------------------------------------------------------
  ! Read in namelist and set File IO Type and Format accordingly...
  !----------------------------------------------------------------

  if(Debug)    print *,'testpio: before call to readTestPIO_Namelist'
  if(my_task == master_task) then 
     call ReadTestPIO_Namelist(nml_in, nprocs, nml_filename, myname, nml_error)
  endif
  if(Debug) print *,'testpio: before call to broadcast_namelist'
  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call Broadcast_Namelist(myname, my_task, master_task, MPI_COMM_WORLD, ierr)
  if(Debug) print *,'testpio: after call to broadcast_namelist'

  !-------------------------------------
  ! Checks
  !-------------------------------------

#if !defined(BGx) 
  if (num_iotasks <= 0) then
     write(*,*) trim(myname),' ERROR: ioprocs invalid num_iotasks=',num_iotasks
     call piodie(__FILE__,__LINE__)
  endif
#endif

  ! ----------------------------------------------------------------
  ! if stride is and num_iotasks is incompatible than reset stride
  ! ----------------------------------------------------------------
  if (base + num_iotasks * (stride-1) > nprocs-1) then
     write(*,*) trim(myname),' ERROR: num_iotasks, base and stride too large', &
          ' base=',base,' num_iotasks=',num_iotasks,' stride=',stride,' nprocs=',nprocs
     call piodie(__FILE__,__LINE__)
  endif

  !--------------------------------------
  ! Initalizes the parallel IO subsystem 
  !--------------------------------------
  call PIO_setDebugLevel(DebugLevel)

  if(Debug)    print *,'testpio: before call to PIO_init'

  if(async) then
!     print *,__FILE__,__LINE__,nprocs, num_iotasks, stride, base
     call split_comm(mpi_comm_world,nprocs, num_iotasks, stride, base, &
          mpi_comm_compute, mpi_comm_io, mpi_icomm_cio)
     call PIO_init(mpi_comm_compute, mpi_comm_io, mpi_icomm_cio, PIOSYS)
  else
     mpi_comm_compute = mpi_comm_world

     call PIO_init(my_task, MPI_COMM_COMPUTE, num_iotasks, num_aggregator, stride, &
          rearr_type, PIOSYS, base)
  end if
  if(Debug)    print *,'testpio: after call to PIO_init'

  gDims3D(1) = nx_global
  gDims3D(2) = ny_global
  gDims3D(3) = nz_global

  !-----------------------------------------
  ! Compute compDOF based on namelist input
  !-----------------------------------------

  startCOMP = 0
  if (trim(compdof_input) == 'namelist') then
     if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #1'
     call gdecomp_read_nml(gdecomp,nml_filename,'comp',PIOSYS%comp_rank,PIOSYS%num_tasks,gDims3D(1:3))
     if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #2'
     call gdecomp_DOF(gdecomp,PIOSYS%comp_rank,compDOF,start,count)
     if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #3'
  else
     call pio_readdof(trim(compdof_input),compDOF,MPI_COMM_COMPUTE,75)
     sdof = size(compDOF)
     start = gDims3D(1:3)   ! min tmp
     count = 0         ! max tmp
     do n = 1,sdof
        call c1dto3d(compdof(n),gDims3D(1),gDims3D(2),gDims3D(3),i1,j1,k1)
        start(1) = min(start(1),i1)
        start(2) = min(start(2),j1)
        start(3) = min(start(3),k1)
        count(1) = max(count(1),i1)
        count(2) = max(count(2),j1)
        count(3) = max(count(3),k1)
     enddo
     do n = 1,3
        count(n) = max(count(n)-start(n)+1,0)
     enddo
     if (count(1)*count(2)*count(3) == sdof) then
        ! start and count seem consistent with compDOF
     else
        ! start and count NOT consistent with compDOF, zero them out
        start = 0
        count = 0
     endif
  endif
  startCOMP(1:3) = start(1:3)
  countCOMP(1:3) = count(1:3)
  if (trim(compdof_output) /= 'none') then
     call pio_writedof(trim(compdof_output),compDOF,MPI_COMM_COMPUTE,75)
  endif

  sdof = size(compDOF)
  call MPI_REDUCE(sdof,sdof_sum,1,MPI_INTEGER,MPI_SUM,master_task,MPI_COMM_COMPUTE,ierr)
  call CheckMPIReturn('Call to MPI_REDUCE SUM',ierr,__FILE__,__LINE__)
  call MPI_REDUCE(sdof,sdof_min,1,MPI_INTEGER,MPI_MIN,master_task,MPI_COMM_COMPUTE,ierr)
  call CheckMPIReturn('Call to MPI_REDUCE MIN',ierr,__FILE__,__LINE__)
  call MPI_REDUCE(sdof,sdof_max,1,MPI_INTEGER,MPI_MAX,master_task,MPI_COMM_COMPUTE,ierr)
  call CheckMPIReturn('Call to MPI_REDUCE MAX',ierr,__FILE__,__LINE__)
  if (my_task == master_task) then
     write(6,*) trim(myname),' total nprocs = ',nprocs
     write(6,*) trim(myname),' compDOF sum/min/max = ',sdof_sum,sdof_min,sdof_max
  endif

  if (mod(my_task,((nprocs/32)+1)) == 0) then
     if(Debug)       write(6,*) trim(myname),' my_task,sdof,start,count = ',my_task,sdof,start,count
  endif

  !--------------------------------
  ! calculate ioDOF
  !--------------------------------

  startIO = 0
  countIO = 0
  if (trim(rearr) == 'none') then
     ioDOF   => compDOF
     startIO(1:3) = startCOMP(1:3)
     countIO(1:3) = countCOMP(1:3)
  elseif (trim(rearr) == 'box') then
     ! do nothing
     if (trim(iodof_input) == 'namelist') then
        if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #4'
        call gdecomp_read_nml(gdecomp,nml_filename,'io',PIOSYS%io_rank,PIOSYS%num_iotasks,gDims3D(1:3))
        if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #5'
        call gdecomp_DOF(gdecomp,PIOSYS%io_rank,ioDOF,start,count)
        if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #6'
        startIO(1:3) = start(1:3)
        countIO(1:3) = count(1:3)
     endif
  elseif (trim(rearr) == 'mct') then
     call gdecomp_read_nml(gdecomp,nml_filename,'io',PIOSYS%io_rank,PIOSYS%num_iotasks,gDims3D(1:3))
     call gdecomp_DOF(gdecomp,PIOSYS%io_rank,ioDOF,start,count)
     startIO(1:3) = start(1:3)
     countIO(1:3) = count(1:3)
  else
     call piodie(__FILE__,__LINE__,' rearr '//trim(rearr)//' not supported')
  endif

  ioDOFR => ioDOF
  ioDOFW => ioDOF
  startpio = startIO
  countpio = countIO
  lenblocks = countIO(1)

  if(Debug) print *,'comp_rank,io_rank: ',piosys%comp_rank,piosys%io_rank,' ioDOF ',ioDOF
  if(Debug) print *,'comp_rank: ',PIOSYS%comp_rank,SIZE(compDOF),SIZE(ioDOF)

  !-----------------------------------------------------
  ! number of words on each computational processor owns
  !-----------------------------------------------------

  lLength = size(compDOF)

  !----------------------
  ! allocate and set test arrays 
  !----------------------

  if(TestR8 .or. TestCombo) then 
     call alloc_check(test_r8wr,lLength,'testpio:test_r8wr')
  endif
  if(TestR4 .or. TestCombo) then 
     call alloc_check(test_r4wr,lLength,'testpio:test_r4wr' )
  endif
  if(TestInt .or. TestCombo) then 
     call alloc_check(test_i4wr,lLength,'testpio:test_i4wr')
  endif
  if(TestInt) then 
    call alloc_check(test_i4i ,lLength,'testpio:test_i4i ')
    call alloc_check(test_i4j ,lLength,'testpio:test_i4j ')
    call alloc_check(test_i4k ,lLength,'testpio:test_i4k ')
    call alloc_check(test_i4m ,lLength,'testpio:test_i4m ')
    call alloc_check(test_i4dof,lLength,'testpio:test_i4dof')
  endif

  do n = 1,lLength
     call c1dto3d(compdof(n),gDims3D(1),gDims3D(2),gDims3D(3),i1,j1,k1)
     if(TestInt) then 
	test_i4dof(n) = compdof(n)
        test_i4i(n) = i1
        test_i4j(n) = j1
        test_i4k(n) = k1
        test_i4m(n) = my_task
     endif
     if(TestR8 .or. TestCombo) then 
         test_r8wr(n) = 10.0_r8*cos(20.*real(i1,kind=r8)/real(gDims3D(1),kind=r8))* &
          cos(10.*real(j1,kind=r8)/real(gDims3D(2),kind=r8))* &
          (1.0+1.0*real(j1,kind=r8)/real(gDims3D(2),kind=r8))* &
          cos(25.*real(k1,kind=r8)/real(gDims3D(3),kind=r8))
     endif
     if(TestR4 .or. TestCombo) then 
         test_r4wr(n) = 10.0_r4*cos(20.*real(i1,kind=r4)/real(gDims3D(1),kind=r4))* &
          cos(10.*real(j1,kind=r4)/real(gDims3D(2),kind=r4))* &
          (1.0+1.0*real(j1,kind=r4)/real(gDims3D(2),kind=r4))* &
          cos(25.*real(k1,kind=r4)/real(gDims3D(3),kind=r4))
     endif
     if(TestInt .or. TestCombo) then 
         test_i4wr(n) = nint(10.0_r8*cos(20.*real(i1,kind=r8)/real(gDims3D(1),kind=r8))* &
          cos(10.*real(j1,kind=r8)/real(gDims3D(2),kind=r8))* &
          (1.0+1.0*real(j1,kind=r8)/real(gDims3D(2),kind=r8))* &
          cos(25.*real(k1,kind=r8)/real(gDims3D(3),kind=r8))*1000.0_r8)
     endif
  enddo
  if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #10'


  if(TestR8  .or. TestCombo)  call alloc_check(test_r8rd,lLength,'testpio:test_r8rd')
  if(TestInt .or. TestCombo) call alloc_check(test_i4rd,lLength,'testpio:test_i4rd')
  if(TestR4  .or. TestCombo)  call alloc_check(test_r4rd,lLength,'testpio:test_r4rd')

  if(TestR8  .or. TestCombo) test_r8rd(:) = 1000.00
  if(TestR4  .or. TestCombo) test_r4rd(:) = 1000.00
  if(TestInt .or. TestCombo) test_i4rd(:) = 1000

  if(Debug) then
     write(*,'(a,2(a,i8))') myname,':: Before call to OpenFile().  comp_rank=',piosys%comp_rank, &
          ' io_rank=',piosys%io_rank
  endif

  !--------------------------------
  ! allocate arrays for holding globally-reduced timing information
  !--------------------------------

  call alloc_check(gdt_write_r8, maxiter, ' testpio:gdt_write_r8 ')
  call alloc_check(gdt_read_r8, maxiter, ' testpio:gdt_read_r8 ')
  call alloc_check(gdt_write_r4, maxiter, ' testpio:gdt_write_r4 ')
  call alloc_check(gdt_read_r4, maxiter, ' testpio:gdt_read_r4 ')
  call alloc_check(gdt_write_i4, maxiter, ' testpio:gdt_write_i4 ')
  call alloc_check(gdt_read_i4, maxiter, ' testpio:gdt_read_i4 ')
  if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #11'

  if(splitPhase) then 
    numPhases = 2
  else
    numPhases = 1
  endif
  do ip=1,numPhases
     if(numPhases == 1) then 
        readPhase = .true.
        writePhase = .true.
     else
        if(ip == 1) then 
	   writePhase = .true.
	   readPhase = .false.
	else
	   writePhase = .false.
	   readPhase = .true.
        endif
     endif
     if(log_master_task) print *,'{write,read}Phase:  ',writePhase,readPhase
  do it=1,maxiter

     !-------------------------------------------------------
     ! Explain the distributed array decomposition to PIOlib
     !-------------------------------------------------------

     if (trim(rearr) == 'box') then
        !JMD print *,__FILE__,__LINE__,gdims3d,minval(compdof),maxval(compdof)

        if (trim(iodof_input) == 'namelist') then
           if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #7'
           if(TestR8 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_double,  gDims3D,compDOF,IOdesc_r8,startpio,countpio)
           if(TestR4 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_real,    gDims3D,compDOF,IOdesc_r4,startpio,countpio)
           if(TestInt .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_int,     gDims3D,compDOF,IOdesc_i4,startpio,countpio)
           if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #8'
        else
           if(TestR8 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_double,  gDims3D,compDOF,IOdesc_r8)
           if(TestR4 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_real,    gDims3D,compDOF,IOdesc_r4)
           if(TestInt .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_int,     gDims3D,compDOF,IOdesc_i4)
        endif
     else
        if(iofmtd.eq.'nc') then ! netCDF
           if (num_iodofs == 1) then
              if(TestR8 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_double, gDims3D,lenblocks,compDOF,ioDOF,startpio,countpio,IOdesc_r8)
              if(TestR4 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_real,   gDims3D,lenblocks,compDOF,ioDOF,startpio,countpio,IOdesc_r4)
              if(TestInt .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_int,    gDims3D,lenblocks,compDOF,ioDOF,startpio,countpio,IOdesc_i4)
           elseif (num_iodofs == 2) then
              if(TestR8 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_double, gDims3D,lenblocks,compDOF,ioDOFR,ioDOFW,startpio,countpio,IOdesc_r8)
              if(TestR4 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_real,   gDims3D,lenblocks,compDOF,ioDOFR,ioDOFW,startpio,countpio,IOdesc_r4)
              if(TestInt .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_int,    gDims3D,lenblocks,compDOF,ioDOFR,ioDOFW,startpio,countpio,IOdesc_i4)
           else
              call piodie(__FILE__,__LINE__,' num_iodofs not 1 or 2')
           endif
        else
           ! tcraig: there are cases where lenblocks is not valid here like different size IO blocks
           if (num_iodofs == 1) then
              if(TestR8 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_double, gDims3D,lenblocks,compDOF,ioDOF,IOdesc_r8)
              if(TestR4 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_real,   gDims3D,lenblocks,compDOF,ioDOF,IOdesc_r4)
              if(TestInt .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_int,    gDims3D,lenblocks,compDOF,ioDOF,IOdesc_i4)
           elseif (num_iodofs == 2) then
              if(TestR8 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_double, gDims3D,lenblocks,compDOF,ioDOFR,ioDOFW,IOdesc_r8)
              if(TestR4 .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_real,   gDims3D,lenblocks,compDOF,ioDOFR,ioDOFW,IOdesc_r4)
              if(TestInt .or. TestCombo) call PIO_initDecomp(PIOSYS,PIO_int,    gDims3D,lenblocks,compDOF,ioDOFR,ioDOFW,IOdesc_i4)
           else
              call piodie(__FILE__,__LINE__,' num_iodofs not 1 or 2')
           endif
        endif
     endif
     if(Debug)       print *,'iam: ',PIOSYS%comp_rank,'testpio: point #9'

     if(Debug) then
        write(*,'(a,2(a,i8))') myname,':: After call to initDecomp.  comp_rank=',piosys%comp_rank, &
             ' io_rank=',piosys%io_rank
     endif

     call PIO_getnumiotasks(PIOSYS,num_iotasks)
     !------------
     ! Open file{s} 
     !------------
     write(citer,'(i3.3)') it

     fname    = TRIM(dir)//'foo.'//citer//TRIM(Iofmtd)
     fname_r8 = TRIM(dir)//'foo.r8.'//citer//TRIM(Iofmtd)
     fname_r4 = TRIM(dir)//'foo.r4.'//citer//TRIM(Iofmtd)
     fname_i4 = TRIM(dir)//'foo.i4.'//citer//TRIM(Iofmtd)
     !   print *, __FILE__,__LINE__,'>',fname,'<'
     !   print *, __FILE__,__LINE__,'>',fname_r8,'<'
     !   print *, __FILE__,__LINE__,'>',fname_i4,'<'
     !   print *, __FILE__,__LINE__,'>',fname_r4,'<'
#if defined(_NETCDF) || defined(_PNETCDF)
     mode = pio_64bit_offset
#else
     mode = 0
#endif

if(writePhase) then 
     if(TestCombo) then
        if(Debug) write(*,'(2a,i8)') myname,':: Combination Test:  Creating File...it=',it
        ierr = PIO_CreateFile(PIOSYS,File,iotype,trim(fname), mode)
        call check_pioerr(ierr,__FILE__,__LINE__,' combo createfile')
     endif

     if(TestR8) then
        if (Debug) write(*,'(2a,i8)') myname,':: REAL*8 Test:  Creating File...it=',it
        ierr = PIO_CreateFile(PIOSYS,File_r8,iotype,trim(fname_r8), mode)
        call check_pioerr(ierr,__FILE__,__LINE__,' r8 createfile')
     endif

     if(TestR4) then
        if(Debug) write(*,'(2a,i8)') myname,':: REAL*4 Test:  Creating File...,it=',it
        ierr = PIO_CreateFile(PIOSYS,File_r4,iotype,trim(fname_r4), mode)
        call check_pioerr(ierr,__FILE__,__LINE__,' r4 createfile')
     endif

     if(TestInt) then
        if(Debug) write(*,'(2a,i8)') myname,':: INTEGER*4 Test:  Creating File...,it=',it
        ierr = PIO_CreateFile(PIOSYS,File_i4,iotype,trim(fname_i4), mode)
        call check_pioerr(ierr,__FILE__,__LINE__,' i4 createfile')
     endif

     ! Set Frame to '1' in the PIO descriptor file

     allocate(vard_r8(nvars), vard_r4(nvars))

     one = 1
     do ivar=1,nvars
        call PIO_SetFrame(vard_r8(ivar),one)
        call PIO_SetFrame(vard_r4(ivar),one)
     end do

     call PIO_SetFrame(vard_i4,one)
     call PIO_SetFrame(vard_r8c,one)
     call PIO_SetFrame(vard_r4c,one)
     call PIO_SetFrame(vard_i4c,one)
     call PIO_SetFrame(vard_i4i,one)
     call PIO_SetFrame(vard_i4j,one)
     call PIO_SetFrame(vard_i4k,one)
     call PIO_SetFrame(vard_i4m,one)
     call PIO_SetFrame(vard_i4dof,one)

     !---------------------------
     ! Code specifically for netCDF files 
     !---------------------------
     if(iotype == iotype_pnetcdf .or. & 
          iotype == iotype_netcdf .or. &
          iotype == PIO_iotype_netcdf4p .or. &
          iotype == PIO_iotype_netcdf4c) then

        if(TestR8) then 
           !-----------------------------------
           ! for the single record real*8 file 
           !-----------------------------------
           call WriteHeader(File_r8,nx_global,ny_global,nz_global,dimid_x,dimid_y,dimid_z)

           do ivar = 1, nvars
              write(varname,'(a,i5.5)') 'field',ivar
              iostat = PIO_def_var(File_r8,varname,PIO_double,(/dimid_x,dimid_y,dimid_z/),vard_r8(ivar))
              call check_pioerr(iostat,__FILE__,__LINE__,' r8 defvar')
           end do
           iostat = PIO_enddef(File_r8)
           call check_pioerr(iostat,__FILE__,__LINE__,' r8 enddef')
        endif

        if(TestR4) then 
           !-----------------------------------
           ! for the single record real*4 file 
           !-----------------------------------
           call WriteHeader(File_r4,nx_global,ny_global,nz_global,dimid_x,dimid_y,dimid_z)

           do ivar = 1, nvars
              write(varname,'(a,i5.5)') 'field',ivar
              iostat = PIO_def_var(File_r4,varname,PIO_real,(/dimid_x,dimid_y,dimid_z/),vard_r4(ivar))
              call check_pioerr(iostat,__FILE__,__LINE__,' r4 defvar')
           end do
           iostat = PIO_enddef(File_r4)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4 enddef')
        endif

        if(TestInt) then 
           !-----------------------------------
           ! for the single record integer file 
           !-----------------------------------
           call WriteHeader(File_i4,nx_global,ny_global,nz_global,dimid_x,dimid_y,dimid_z)

           iostat = PIO_def_var(File_i4,'field',PIO_int,(/dimid_x,dimid_y,dimid_z/),vard_i4)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4 defvar')
           iostat = PIO_def_var(File_i4,'fi',PIO_int,(/dimid_x,dimid_y,dimid_z/),vard_i4i)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4i defvar')
           iostat = PIO_def_var(File_i4,'fj',PIO_int,(/dimid_x,dimid_y,dimid_z/),vard_i4j)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4j defvar')
           iostat = PIO_def_var(File_i4,'fk',PIO_int,(/dimid_x,dimid_y,dimid_z/),vard_i4k)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4k defvar')
           iostat = PIO_def_var(File_i4,'my_task',PIO_int,(/dimid_x,dimid_y,dimid_z/),vard_i4m)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4m defvar')
           iostat = PIO_def_var(File_i4,'fdof',PIO_int,(/dimid_x,dimid_y,dimid_z/),vard_i4dof)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4dof defvar')
           iostat = PIO_enddef(File_i4)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4 enddef')
        endif

        if(TestCombo) then 
           !-----------------------------------
           ! for the multi record file 
           !-----------------------------------
           call WriteHeader(File,nx_global,ny_global,nz_global,dimid_x,dimid_y,dimid_z)
           iostat = PIO_def_var(File,'field_r8',PIO_double,(/dimid_x,dimid_y,dimid_z/),vard_r8c)
           call check_pioerr(iostat,__FILE__,__LINE__,' combo r8 defvar')
           iostat = PIO_def_var(File,'field_r4',PIO_real,(/dimid_x,dimid_y,dimid_z/),vard_r4c)
           call check_pioerr(iostat,__FILE__,__LINE__,' combo r4 defvar')
           iostat = PIO_def_var(File,'field_i4',PIO_int,(/dimid_x,dimid_y,dimid_z/),vard_i4c)
           call check_pioerr(iostat,__FILE__,__LINE__,' combo i4 defvar')
           iostat = PIO_enddef(File)
           call check_pioerr(iostat,__FILE__,__LINE__,' combo enddef')
        endif

     endif ! if(iotype == iotype_pnetcdf .or. iotype == iotype_netcdf ) then

     if(Debug) then
        write(*,'(a,2(a,i8))') myname,':: After call to OpenFile.  comp_rank=',piosys%comp_rank, &
             ' io_rank=',piosys%io_rank
     endif

     !-------------------------
     ! Time the parallel write 
     !-------------------------

     if(TestR8) then
        dt_write_r8 = 0.
        call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
        call CheckMPIReturn('Call to MPI_BARRIER()',ierr,__FILE__,__LINE__)
        st = MPI_Wtime()
#ifdef TIMING
        call t_startf('testpio_write')
#endif
        do ivar=1,nvars
           print *,__FILE__,__LINE__,ivar
           call PIO_write_darray(File_r8,vard_r8(ivar), iodesc_r8, test_r8wr, iostat)
           call check_pioerr(iostat,__FILE__,__LINE__,' r8 write_darray')
        end do
        print *,__FILE__,__LINE__
#ifdef TIMING
        call t_stopf('testpio_write')
#endif
        et = MPI_Wtime()
        dt_write_r8 = dt_write_r8 + (et - st)/nvars
        call PIO_CloseFile(File_r8)
     endif

     if(TestR4) then
        dt_write_r4 = 0.
        call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
        call CheckMPIReturn('Call to MPI_BARRIER()',ierr,__FILE__,__LINE__)
        st = MPI_Wtime()
#ifdef TIMING
        call t_startf('testpio_write')
#endif
        do ivar=1,nvars
           call PIO_write_darray(File_r4,vard_r4(ivar),iodesc_r4, test_r4wr,iostat)
           call check_pioerr(iostat,__FILE__,__LINE__,' r4 write_darray')
        end do
#ifdef TIMING
        call t_stopf('testpio_write')
#endif
        et = MPI_Wtime()
        dt_write_r4 = dt_write_r4 + (et - st)/nvars
        call PIO_CloseFile(File_r4)
     endif

     if(TestInt) then 
        dt_write_i4 = 0.
        call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
        call CheckMPIReturn('Call to MPI_BARRIER()',ierr,__FILE__,__LINE__)
        st = MPI_Wtime()
#ifdef TIMING
        call t_startf('testpio_write')
#endif
        call PIO_write_darray(File_i4,vard_i4,iodesc_i4,test_i4wr,iostat)
#ifdef TIMING
        call t_stopf('testpio_write')
#endif
        et = MPI_Wtime()
        dt_write_i4 = dt_write_i4 + et - st
        call check_pioerr(iostat,__FILE__,__LINE__,' i4 write_darray')
        call PIO_write_darray(File_i4,vard_i4i,iodesc_i4,test_i4i,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' i4i write_darray')
        call PIO_write_darray(File_i4,vard_i4j,iodesc_i4,test_i4j,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' i4j write_darray')
        call PIO_write_darray(File_i4,vard_i4k,iodesc_i4,test_i4k,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' i4k write_darray')
        call PIO_write_darray(File_i4,vard_i4m,iodesc_i4,test_i4m,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' i4m write_darray')
        call PIO_write_darray(File_i4,vard_i4dof,iodesc_i4,test_i4dof,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' i4dof write_darray')
        call PIO_CloseFile(File_i4) 
     endif

     if(TestCombo) then 
        call PIO_write_darray(File,vard_r8c,iodesc_r8,test_r8wr,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' combo r8 write_darray')
        call PIO_write_darray(File,vard_r4c,iodesc_r4, test_r4wr,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' combo r4 write_darray')
        call PIO_write_darray(File,vard_i4c,iodesc_i4, test_i4wr,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' combo i4 write_darray')
        call PIO_CloseFile(File)
     endif

     if(Debug) then
        write(*,'(a,2(a,i8))') myname,':: After calls to PIO_write_darray.  comp_rank=',piosys%comp_rank, & 
             ' io_rank=',piosys%io_rank
     endif

endif
     call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
     call CheckMPIReturn('Call to MPI_BARRIER()',ierr,__FILE__,__LINE__)

if (readPhase) then 
     !-------------------------------------
     ! Open the file back up and check data
     !-------------------------------------

     if(TestR8) then 
        ierr = PIO_OpenFile(PIOSYS, File_r8, iotype, fname_r8)
        call check_pioerr(ierr,__FILE__,__LINE__,' r8 openfile')
     endif

     if(TestR4) then 
        ierr = PIO_OpenFile(PIOSYS,File_r4,iotype, fname_r4)
        call check_pioerr(ierr,__FILE__,__LINE__,' r4 openfile')
     endif

     if(TestInt) then
        ierr = PIO_OpenFile(PIOSYS,File_i4,iotype, fname_i4)
        call check_pioerr(ierr,__FILE__,__LINE__,' int openfile')
     endif

     !   if(TestCombo) ierr = PIO_OpenFile(PIOSYS,File,iotype,fname)
     if(Debug) then
        write(*,'(2a,i8)') myname,':: After calls to PIO_OpenFile.  my_task=',my_task
     endif

     if(iotype == iotype_pnetcdf .or. &
          iotype == iotype_netcdf ) then
        do ivar=1,nvars
           if(TestR8) then 
              iostat = PIO_inq_varid(File_r8,'field00001',vard_r8(ivar))
              call check_pioerr(iostat,__FILE__,__LINE__,' r8 inq_varid')
           endif

           if(TestR4) then 
              iostat = PIO_inq_varid(File_r4,'field00001',vard_r4(ivar))
              call check_pioerr(iostat,__FILE__,__LINE__,' r4 inq_varid')
           endif
        end do
        if(TestInt) then 
           iostat = PIO_inq_varid(File_i4,'field',vard_i4)
           call check_pioerr(iostat,__FILE__,__LINE__,' i4 inq_varid')
        endif

     endif ! if((iotype == iotype_pnetcdf) .or (iotype == iotype_netcdf))...
     do ivar=1,nvars
        call PIO_SetFrame(vard_r8(ivar),one)
        call PIO_SetFrame(vard_r4(ivar),one)
     end do
     call PIO_SetFrame(vard_i4,one)
     call PIO_SetFrame(vard_r8c,one)
     call PIO_SetFrame(vard_r4c,one)
     call PIO_SetFrame(vard_i4c,one)

     !-------------------------
     ! Time the parallel  read
     !-------------------------
     if(TestR8) then 
        dt_read_r8 = 0.
        call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
        call CheckMPIReturn('Call to MPI_BARRIER()',ierr,__FILE__,__LINE__)
        st = MPI_Wtime()
#ifdef TIMING
        call t_startf('testpio_read')
#endif
	do ivar=1,nvars
          call PIO_read_darray(File_r8,vard_r8(ivar),iodesc_r8,test_r8rd,iostat)
          call check_pioerr(iostat,__FILE__,__LINE__,' r8 read_darray')
        enddo
#ifdef TIMING
        call t_stopf('testpio_read')
#endif
        et = MPI_Wtime()
        dt_read_r8 = dt_read_r8 + (et - st)/nvars        
        call check_pioerr(iostat,__FILE__,__LINE__,' r8 read_darray')
     endif

     if(TestR4) then
        dt_read_r4 = 0.
        call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
        call CheckMPIReturn('Call to MPI_BARRIER()',ierr,__FILE__,__LINE__)
        st = MPI_Wtime()
#ifdef TIMING
        call t_startf('testpio_read')
#endif
	do ivar=1,nvars
           call PIO_read_darray(File_r4,vard_r4(ivar),iodesc_r4,test_r4rd,iostat)
           call check_pioerr(iostat,__FILE__,__LINE__,' r4 read_darray')
        enddo
#ifdef TIMING
        call t_stopf('testpio_read')
#endif
        et = MPI_Wtime()
        dt_read_r4 = dt_read_r4 + (et - st)/nvars
        call check_pioerr(iostat,__FILE__,__LINE__,' r4 read_darray')
     endif

     if(TestInt) then
        dt_read_i4 = 0.
        call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
        call CheckMPIReturn('Call to MPI_BARRIER()',ierr,__FILE__,__LINE__)
        st = MPI_Wtime()
#ifdef TIMING
        call t_startf('testpio_read')
#endif
        call PIO_read_darray(File_i4,vard_i4,iodesc_i4, test_i4rd,iostat)
#ifdef TIMING
        call t_stopf('testpio_read')
#endif
        et = MPI_Wtime()
        dt_read_i4 = dt_read_i4 + et - st
        call check_pioerr(iostat,__FILE__,__LINE__,' i4 read_darray')
     endif

     !-------------------------------
     ! Print the maximum memory usage 
     !-------------------------------
     call alloc_print_usage(0,'testpio: after calls to PIO_read_darray')

#ifdef TESTMEM
     stop 
#endif

     if(Debug) then
        write(*,'(a,2(a,i8))') myname,':: After PIO_read_darray tests, my_task=', &
             my_task,', it=',it
     endif

     !-------------------
     ! close the file up 
     !-------------------
     if(TestR8) call PIO_CloseFile(File_r8)
     if(TestR4) call PIO_CloseFile(File_r4)
     if(TestInt) call PIO_CloseFile(File_i4)

     call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
     call CheckMPIReturn('Call to MPI_BARRIER()',ierr,__FILE__,__LINE__)

     if(Debug) then
        write(*,*) myname,':: my_task=',my_task,'test_r8wr= ',test_r8wr
        if(TestR8 .or. TestCombo) write(*,*) myname,':: my_task=',my_task,'test_r8rd= ',test_r8rd
     endif

     !-----------------------------
     ! Perform correctness testing 
     !-----------------------------
     if(TestR8) then
        call checkpattern(fname_r8,test_r8wr,test_r8rd,lLength,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' checkpattern r8 test')
     endif
     
     if( TestR4) then
        call checkpattern(fname_r4,test_r4wr,test_r4rd,lLength,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' checkpattern r4 test')
     endif

     if(TestInt) then
        call checkpattern(fname_i4, test_i4wr,test_i4rd,lLength,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' checkpattern i4 test')
     endif

     if(TestCombo) then 

        !-------------------------------------
        !  Open up and read the combined file 
        !-------------------------------------
           
        ierr = PIO_OpenFile(PIOSYS,File,iotype,fname)
        call check_pioerr(ierr,__FILE__,__LINE__,' combo test read openfile')
           
        if(iofmtd(1:2).eq.'nc') then
           iostat = PIO_inq_varid(File,'field_r8',vard_r8c)
           call check_pioerr(iostat,__FILE__,__LINE__,' combo test r8 inq_varid')
           iostat = PIO_inq_varid(File,'field_r4',vard_r4c)
           call check_pioerr(iostat,__FILE__,__LINE__,' combo test r4 inq_varid')
           iostat = PIO_inq_varid(File,'field_i4',vard_i4c)
           call check_pioerr(iostat,__FILE__,__LINE__,' combo test i4 inq_varid')
        endif
           
        call PIO_read_darray(File,vard_r8c,iodesc_r8,test_r8rd,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' combo test r4 pio_read_darray')
        call PIO_read_darray(File,vard_r4c,iodesc_r4,test_r4rd,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' combo test r4 pio_read_darray')
        call PIO_read_darray(File,vard_i4c,iodesc_i4,test_i4rd,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' combo test i4 pio_read_darray')

        call PIO_CloseFile(File)
           
        !-----------------------------
        ! Check the combined file 
        !-----------------------------
        call checkpattern(fname,test_r8wr,test_r8rd,lLength,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' checkpattern test_r8 ')
        
        call checkpattern(fname,test_r4wr,test_r4rd,lLength,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' checkpattern test_r4 ')
        
        call checkpattern(fname,test_i4wr,test_i4rd,lLength,iostat)
        call check_pioerr(iostat,__FILE__,__LINE__,' checkpattern test_i4 ')

     endif
     !---------------------------------------
     ! Print out the performance measurements 
     !---------------------------------------
     call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
endif

     if(TestR8) then
        ! Maximum read/write times
        if(readPhase)  call GetMaxTime(dt_read_r8, gdt_read_r8(it), MPI_COMM_COMPUTE, ierr)
        if(writePhase) call GetMaxTime(dt_write_r8, gdt_write_r8(it), MPI_COMM_COMPUTE, ierr)
     endif

     if(TestR4) then
        ! Maximum read/write times
        if(readPhase)  call GetMaxTime(dt_read_r4, gdt_read_r4(it), MPI_COMM_COMPUTE, ierr)
        if(writePhase) call GetMaxTime(dt_write_r4, gdt_write_r4(it), MPI_COMM_COMPUTE, ierr)
     endif

     if(TestInt) then
        ! Maximum read/write times
        if(readPhase)  call GetMaxTime(dt_read_i4, gdt_read_i4(it), MPI_COMM_COMPUTE, ierr)
        if(writePhase) call GetMaxTime(dt_write_i4, gdt_write_i4(it), MPI_COMM_COMPUTE, ierr)
     endif


     if(TestR8 .or. TestCombo) glenr8=iodesc_r8%glen
     if(TestR4 .or. TestCombo) glenr4=iodesc_r4%glen
     if(TestInt .or. TestCombo) gleni4=iodesc_i4%glen
     if(TestR8 .or. TestCombo) call pio_freedecomp(PIOSYS, iodesc_r8)
     if(TestInt .or. TestCombo) call pio_freedecomp(PIOSYS, iodesc_r4)
     if(TestR4 .or. TestCombo) call pio_freedecomp(PIOSYS, iodesc_i4)
  enddo ! do it=1,maxiter
  enddo ! do ip=1,numphase

  !--------------------------------
  ! Clean up initialization memory 
  !   note: make sure DOFs are not used later
  !--------------------------------
  if (PIOSYS%comp_rank >= 0) call dealloc_check(compDOF)
  if (trim(rearr) == 'mct') then
     if (PIOSYS%io_rank >= 0) call dealloc_check(ioDOF)
  endif

  !----------------------------------
  ! Print summary bandwidth statistics 
  !----------------------------------

  if(TestR8 .and. (piosys%io_rank == 0) ) then
     call WriteTimeTrialsStats(casename,TestR8CaseName, fname_r8, glenr8, gdt_read_r8, gdt_write_r8, maxiter) 
  endif

  if(TestR4 .and. (piosys%io_rank == 0) ) then
     call WriteTimeTrialsStats(casename,TestR4CaseName, fname_r4, glenr4, gdt_read_r4, gdt_write_r4, maxiter) 
  endif

  if(TestInt .and. (piosys%io_rank == 0) ) then
     call WriteTimeTrialsStats(casename,TestI4CaseName, fname_i4, gleni4, gdt_read_i4, gdt_write_i4, maxiter) 
  endif

  !-------------------------------
  ! Print timers and memory usage 
  !-------------------------------

#ifdef TIMING
  call t_stopf('testpio_total')
  call t_prf('timing.testpio',MPI_COMM_COMPUTE)
  call t_finalizef()
  call get_memusage(msize,mrss)
  allocate(lmem(2),gmem(2,0:nprocs-1))
  lmem(1) = msize
  lmem(2) = mrss
  call mpi_gather(lmem,2,MPI_INTEGER,gmem,2,MPI_INTEGER,0,MPI_COMM_COMPUTE,ierr)
  if (my_task == master_task) then
     do n = 0,nprocs-1
        write(*,'(2a,i8,a,2f10.2)') myname,' my_task=',n,' : (hw, usage) memory (MB) = ',gmem(1,n)*mb_blk,gmem(2,n)*mb_blk
     enddo
!     indx = MAXLOC(gmem(1,:),dim=1) - 1 ! offset the location of the maximum memory usage by one
     indx = MAXLOC(gmem(2,:),dim=1) - 1
     write(*,'(2a,i8,a,2f10.2)') myname,' my_task=',indx,' : (hw, usage) MAX memory (MB) = ',gmem(1,indx)*mb_blk,gmem(2,indx)*mb_blk
  endif
  deallocate(lmem,gmem)
#endif

  call MPI_Barrier(MPI_COMM_COMPUTE,ierr)
  if (my_task == master_task) then
     print *,' '
     print *,'testpio completed successfully'
     print *,' '
  endif

  call PIO_finalize(PIOSYS,ierr)
  call MPI_Finalize(ierr)
  call CheckMPIReturn('Call to MPI_FINALIZE()',ierr,__FILE__,__LINE__)

  !=============================================================================
contains
  !=============================================================================

  subroutine GetMaxTime(dtLocal, gdtMax, comm, ierror)

    implicit none

    real(r8),    intent(IN)  :: dtLocal
    real(r8),    intent(OUT) :: gdtMax
    integer(i4), intent(IN)  :: comm
    integer(i4), intent(OUT) :: ierror
    real(r8) :: local_temp

    local_temp = max(dtlocal, MPI_Wtick())

    call MPI_Allreduce(Local_temp, gdtMax, 1,MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierror)

    call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierror,__FILE__,__LINE__)

  end subroutine GetMaxTime

  !=============================================================================

  subroutine WriteStats(CaseName, FileName, glen, trialNo, dtRead, dtWrite)

    implicit none

    character(len=*),  intent(IN) :: CaseName
    character(len=80), intent(IN) :: FileName
    integer(i4),       intent(in) :: glen
    integer(i4),       intent(IN) :: trialNo
    real(r8),          intent(IN) :: dtRead
    real(r8),          intent(IN) :: dtWrite

    character(len=*), parameter :: myname_=myname//'::WriteStats'
    integer :: datumSize

    select case(CaseName)
    case(TestR8CaseName)
       datumSize = r8
    case(TestR4CaseName)
       datumSize = r4
    case(TestI4CaseName)
       datumSize = i4
    case(TestComboCaseName)
       write(*,'(4a)') myname_,':: Case ',CaseName,' not supported.  Returning without writing output.'
    case default
       write(*,'(4a)') myname_,':: Case ',CaseName,' not supported.  Returning without writing output.'
    end select

    print *,'-----------------------------------------'
    print *,myname,':: Timings for ',CaseName,' Trial Number=',trialNo
    print *,'Total Procs: ',nprocs,' IO Procs: ',num_iotasks, &
         ' Aggregators: ',num_aggregator,' Stride: ',stride
    print *,'Record bytes: ',INT(glen*datumSize,kind=i8), &
         INT(glen*datumSize,kind=i8)
    print *,'-----------------------------------------'
    print *,'Type of I/O performed: ',Iofmtd
    print *,' File name: ',FileName
    print *,'-----------------------------------------'
    print *, CaseName,' Trial No. ',trialNo,' Read time: ',dtRead, &
         'Read Mbytes/sec: ',MBYTES*glen*datumSize/dtRead
    print *, CaseName,' Trial No. ',trialNo,' Write time: ',dtWrite, &
         'Write Mbytes/sec: ',MBYTES*glen*datumSize/dtWrite
    print *,'-----------------------------------------'

  end subroutine WriteStats

  !=============================================================================

  subroutine WriteTimeTrialsStats(casename,TestName, FileName, glen, ReadTimes, WriteTimes, nTrials) 

    implicit none

    character(len=*),          intent(IN) :: casename
    character(len=*),          intent(IN) :: TestName
    character(len=80),         intent(IN) :: FileName
    integer(i4),          intent(IN) :: glen
    real(r8),    dimension(:), pointer ::    ReadTimes
    real(r8),    dimension(:), pointer ::    WriteTimes
    integer(i4),               intent(IN) :: nTrials

    character(len=*), parameter :: myname_=myname//'::WriteTimeTrialsStats'

    real(r8), parameter :: tiny = 1.e-10
    real(r8) :: ReadBWAvg, ReadBWStdDev, ReadBWStdErrMean, ReadBWMin, ReadBWMax
    real(r8) :: WriteBWAvg, WriteBWStdDev, WriteBWStdErrMean, WriteBWMin, WriteBWMax
    real(r8) :: WriteTimeAvg,ReadTimeAvg
    real(r8) :: TotalMBytes
    integer ::  datumSize, i, nDOF
    real(r8), dimension(:), pointer :: ReadBW, WriteBW

    if(nTrials .ne. size(ReadTimes)) then
       write(*,'(a,2(a,i8))') myname_,':: ERROR--nTrials = ',nTrials,' size(ReadTimes)=',size(ReadTimes)
       call piodie(__FILE__,__LINE__)
    endif

    if(nTrials .ne. size(WriteTimes)) then
       write(*,'(a,2(a,i8))') myname_,':: ERROR--nTrials = ',nTrials,' size(WriteTimes)=',size(WriteTimes)
       call piodie(__FILE__,__LINE__)
    endif

    select case(TestName)
    case(TestR8CaseName)
       datumSize = r8
    case(TestR4CaseName)
       datumSize = r4
    case(TestI4CaseName)
       datumSize = i4
    case(TestComboCaseName)
       datumSize = 0
    case default
       write(*,'(4a)') myname_,':: TestName ',TestName,' not supported.  Returning without writing output.'
       return
    end select

    TotalMBytes = MBYTES * glen * datumSize

    write(*,*)'-----------------------------------------'
    write(*,*) 'Timing Output :: case = ',trim(casename),'  test = ',trim(TestName)
    write(*,'(5a,g14.6)') &
         '  I/O format = ',trim(Iofmtd),'  File name = ',trim(FileName), &
         '  Record Size Mbytes = ',TotalMBytes
    write(*,*) '  Trials = ',nTrials,'  Total Procs = ',nprocs
    write(*,*) '  IO Procs=',num_iotasks,'  Aggregators=',num_aggregator,&
         '  Base=',base,'  Stride=',stride
    do i = 1,nTrials
       if(writetimes(i)>0) then
          write(*,101) 'n=',i,'  write (Mb/sec)=',TotalMBytes/WriteTimes(i), &
               '  write_time(sec)=',WriteTimes(i), &
               trim(casename),trim(TestName)
       else
          write(*,101) 'n=',i,'  write (Mb/sec)=',TotalMBytes, &
               '?  write_time(sec)=',WriteTimes(i), &
               trim(casename),trim(TestName)
       end if
    enddo
    if (nTrials > 1) write(*,*) '         -------------------'
    do i = 1,nTrials
	if(Readtimes(i)>0) then
           write(*,101) 'n=',i,'   read (Mb/sec)=',TotalMBytes/ReadTimes(i), &
                '   read_time(sec)=',ReadTimes(i), &
                trim(casename),trim(TestName)
        else
           write(*,101) 'n=',i,'   read (Mb/sec)=',TotalMBytes, &
                '?   read_time(sec)=',ReadTimes(i), &
                trim(casename),trim(TestName)
        end if
    enddo

101 format(3x,a,i5,a,f9.1,a,e12.4,2x,a,2x,a)

    if (nTrials > 1) then

       call alloc_check(ReadBW, nTrials, myname_//':: ReadBW')
       call alloc_check(WriteBW, nTrials, myname_//':: WriteBW')

       ! Compute mean read/write bandwidths
       ReadBWAvg = 0.
       ReadTimeAvg = 0.
       WriteBWAvg = 0.
       WriteTimeAvg = 0.
       do i=1,nTrials
          ReadBW(i) = TotalMBytes / (ReadTimes(i) + tiny)
          WriteBW(i) = TotalMBytes / (WriteTimes(i) + tiny)
          ReadBWAvg = ReadBWAvg + ReadBW(i)
          ReadTimeAvg = ReadTimeAvg + 1.0e3*ReadTimes(i)
          WriteBWAvg = WriteBWAvg + WriteBW(i)
          WriteTimeAvg = WriteTimeAvg + 1.0e3*WriteTimes(i)
       enddo

       ReadBWAvg = ReadBWAvg / float(nTrials)
       ReadTimeAvg = ReadTimeAvg / float(nTrials)
       WriteBWAvg = WriteBWAvg / float(nTrials)
       WriteTimeAvg = WriteTimeAvg / float(nTrials)

       ! Compute Standard Deviation and Std Error of the Mean
       ReadBWStdDev = 0.
       WriteBWStdDev = 0.
       do i=1,nTrials
          ReadBWStdDev = ReadBWStdDev + (ReadBWAvg - ReadBW(i))**2
          WriteBWStdDev = WriteBWStdDev + (WriteBWAvg - WriteBW(i))**2
       enddo

       ! Compute std. deviation and std. error of the mean
       nDOF = max(1,nTrials-1) ! sample number of degrees-of-freedom
       ReadBWStdDev = sqrt( ReadBWStdDev / float(nDOF) )
       WriteBWStdDev = sqrt( WriteBWStdDev / float(nDOF) )
       ReadBWStdErrMean = ReadBWStdDev / sqrt( float(nTrials) )
       WriteBWStdErrMean = WriteBWStdDev / sqrt( float(nTrials) )

       ! Determine minimum and maximum BW values
       ReadBWMin = minval(ReadBW)
       ReadBWMax = maxval(ReadBW)
       WriteBWMin = minval(WriteBW)
       WriteBWMax = maxval(WriteBW)

       write(*,*) '         -------------------'
       write(*,*) '  Summary BW Stats (MB/sec) for ',nTrials,' trials of ',trim(casename),' ',trim(TestName)
       write(*,102)  'write avg=',WriteBWAvg,' +/-',WriteBWStderrMean, &
            ' min=',WriteBWMin,' max=',WriteBWMax,' stddev=',WriteBWStdDev, &
            trim(casename),trim(TestName)
       write(*,102)  'read  avg=',ReadBWAvg ,' +/-',ReadBWStderrMean, &
            ' min=',ReadBWMin ,' max=',ReadBWMax ,' stddev=',ReadBWStdDev, &
            trim(casename),trim(TestName)
       write(*,103)   'Write Time Avg (usec) =',WriteTimeAvg
       write(*,103)   'Read Time Avg (usec) =',ReadTimeAvg

102    format(3x,5(a,f9.1),1x,a,1x,a)
103    format(3x,a,f9.1)

       call dealloc_check(ReadBW, myname_//':: ReadBW')
       call dealloc_check(WriteBW, myname_//':: WriteBW')

    endif    ! (nTrials > 1)
    write(*,*) '-----------------------------------------'

  end subroutine WriteTimeTrialsStats

  !=======================================================================
#ifdef TIMING
  subroutine get_memusage(msize,mrss)

    integer :: msize,mrss

    integer :: mshare,mtext,mdatastack
    integer :: ierr
    integer :: GPTLget_memusage

    ierr = GPTLget_memusage (msize, mrss, mshare, mtext, mdatastack)

  end subroutine get_memusage
#endif
  !=============================================================================

  subroutine c1dto3d(gindex,nx,ny,nz,i,j,k)

    implicit none
    integer,intent(in) :: gindex,nx,ny,nz
    integer,intent(out) :: i,j,k

    k = (gindex                          - 1) / (nx*ny) + 1
    j = (gindex - (k-1)*nx*ny            - 1) / (nx)    + 1
    i = (gindex - (k-1)*nx*ny - (j-1)*nx - 1)           + 1

  end subroutine c1dto3d

  !=============================================================================
  subroutine check_pioerr(ierr, file, line, str1, str2)

    implicit none
    integer(i4),intent(in) :: ierr
    character(len=*),intent(in) :: file
    integer(i4),intent(in) :: line
    character(len=*),optional,intent(in) :: str1
    character(len=*),optional,intent(in) :: str2

    character(len=256) lstr1
    character(len=256) lstr2
    character(len=*),parameter :: myname_='check_pioerr'

    lstr1 = ''
    if (present(str1)) then
       lstr1 = trim(str1)
    endif
    lstr2 = trim(lstr1)
    if (present(str2)) then
       lstr2 = trim(str2)
    endif

    if(ierr /= PIO_noerr) then
       write(*,*) trim(myname_),':: ERROR on my_task=',my_task,' ierr=',ierr,'  ',trim(lstr1)
       call piodie(file,line,trim(lstr2))
    endif

  end subroutine check_pioerr
  !=============================================================================

end program testpio

