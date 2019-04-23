module dead_mod

  use shr_kind_mod  , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_const_mod , only : shr_const_pi, shr_const_rearth
  use shr_file_mod  , only : shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod   , only : shr_mpi_bcast
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use dead_data_mod , only : dead_grid_lat, dead_grid_lon, dead_grid_area, dead_grid_mask, dead_grid_frac, dead_grid_index
  use dead_data_mod , only : dead_grid_total

  implicit none
  private

  public :: dead_setNewGrid
  public :: dead_read_inparms

!===============================================================================
contains
!===============================================================================

  subroutine dead_read_inparms(model, mpicom, my_task, master_task, &
       inst_index, inst_suffix, inst_name, logunit, &
       nxg, nyg, decomp_type, nproc_x, seg_len, flood)

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)       , intent(in)    :: model
    integer(IN)            , intent(in)    :: mpicom      ! mpi communicator
    integer(IN)            , intent(in)    :: my_task     ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task ! task number of master task
    integer                , intent(in)    :: inst_index  ! number of current instance (ie. 1)
    character(len=16)      , intent(in)    :: inst_suffix ! char string associated with instance
    character(len=16)      , intent(in)    :: inst_name   ! fullname of current instance (ie. "lnd_0001")
    integer(IN)            , intent(in)    :: logunit     ! logging unit number
    integer(IN)		   , intent(out)   :: nproc_x
    integer(IN)            , intent(out)   :: seg_len
    integer(IN)            , intent(out)   :: nxg         ! global dim i-direction
    integer(IN)            , intent(out)   :: nyg         ! global dim j-direction
    integer(IN)            , intent(out)   :: decomp_type ! decomposition type
    logical                , intent(out)   :: flood       ! rof flood flag

    !--- local variables ---
    character(CL) :: fileName    ! generic file name
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: unitn       ! Unit for namelist file

    !--- formats ---
    character(*), parameter :: F00   = "('(dead_read_inparms) ',8a)"
    character(*), parameter :: F01   = "('(dead_read_inparms) ',a,a,4i8)"
    character(*), parameter :: F02   = "('(dead_read_inparms) ',a,L2)"
    character(*), parameter :: F03   = "('(dead_read_inparms) ',a,a,i8,a)"
    character(*), parameter :: subName = "(dead_read_inpamrs) "
    !-------------------------------------------------------------------------------

    ! read the input parms (used to configure model)
    nxg            =  -9999
    nyg            =  -9999
    nproc_x        =  -9999
    seg_len        =  -9999
    decomp_type    =  -9999

    if (my_task == master_task) then
       unitn = shr_file_getUnit()
       open(unitn, file='x'//model//'_in'//trim(inst_suffix), status='old' )
       read(unitn,*) nxg
       read(unitn,*) nyg
       read(unitn,*) decomp_type
       read(unitn,*) nproc_x
       read(unitn,*) seg_len
       if (model.eq.'rof') then
          read(unitn,*) flood
       end if
       close (unitn)
       call shr_file_freeunit(unitn)
    endif

    call shr_mpi_bcast(nxg        , mpicom,'x'//model//' nxg')
    call shr_mpi_bcast(nyg        , mpicom,'x'//model//' nyg')
    call shr_mpi_bcast(decomp_type, mpicom,'x'//model//' decomp_type')
    call shr_mpi_bcast(nproc_x    , mpicom,'x'//model//' nproc_x')
    call shr_mpi_bcast(seg_len    , mpicom,'x'//model//' seg_len')
    if (model.eq.'rof') then
       call shr_mpi_bcast(flood    , mpicom,'xrof flood')
    end if

    if (my_task == master_task) then
       write(logunit,*)' Read in X'//model//' input from file= x'//model//'_in'
       write(logunit,F00) model
       write(logunit,F00) model,'         Model  :  ',model
       write(logunit,F01) model,'           NGX  :  ',nxg
       write(logunit,F01) model,'           NGY  :  ',nyg
       write(logunit,F01) model,' Decomposition  :  ',decomp_type
       write(logunit,F03) model,' Num pes in X   :  ',nproc_x,'  (type 3 only)'
       write(logunit,F03) model,' Segment Length :  ',seg_len,'  (type 11 only)'
       write(logunit,F01) model,'    inst_index  :  ',inst_index
       write(logunit,F00) model,'    inst_name   :  ',trim(inst_name)
       write(logunit,F00) model,'    inst_suffix :  ',trim(inst_suffix)
       if (model.eq.'rof') then
          write(logunit,F02) ' Flood mode     :  ',flood
       endif
       write(logunit,F00) model
       call shr_sys_flush(logunit)
    end if
  end subroutine dead_read_inparms

  !===============================================================================
  subroutine dead_setNewGrid(decomp_type, nxg, nyg, totpe, mype, logunit, lsize,  &
                             gbuf, seg_len, nproc_x)

    implicit none

    ! !DESCRIPTION:
    !     This sets up some defaults.  The user may want to overwrite some
    !     of these fields in the main program after initialization in complete.

    ! !INPUT/OUTPUT PARAMETERS:

    integer(IN) ,intent(in)          :: decomp_type !
    integer(IN) ,intent(in)          :: nxg,nyg     ! global grid sizes
    integer(IN) ,intent(in)          :: totpe       ! total number of pes
    integer(IN) ,intent(in)          :: mype        ! local pe number
    integer(IN) ,intent(in)          :: logunit     ! output logunit
    integer(IN) ,intent(out)         :: lsize       ! local grid sizes
    real(R8)    ,pointer             :: gbuf(:,:)   ! output data
    integer(IN) ,intent(in),optional :: seg_len     ! seg len decomp setting
    integer(IN) ,intent(in),optional :: nproc_x     ! 2d decomp setting

    !--- local ---
    integer(IN)             :: ierr            ! error code
    logical                 :: found
    integer(IN)             :: i,j,ig,jg
    integer(IN)             :: n,ng,is,ie,js,je,nx,ny      ! indices
    integer(IN)             :: npesx,npesy,mypex,mypey,nxp,nyp
    real   (R8)             :: hscore,bscore
    real   (R8)             :: dx,dy,deg2rad,ys,yc,yn,area,re
    integer(IN),allocatable :: gindex(:)

    !--- formats ---
    character(*), parameter :: F00   = "('(dead_setNewGrid) ',8a)"
    character(*), parameter :: F01   = "('(dead_setNewGrid) ',a,4i8)"
    character(*), parameter :: F02   = "('(dead_setNewGrid) ',a,4es13.6)"
    character(*), parameter :: subName = "(dead_setNewGrid) "
    !-------------------------------------------------------------------------------

    if (decomp_type == 1 .or. &
         decomp_type == 2 .or. &
         decomp_type == 3 .or. &
         decomp_type == 4 .or. &
         decomp_type == 11) then
       ! valid
    else
       !-------------------------------------------------------------------------
       ! invalid decomposition type
       !-------------------------------------------------------------------------
       if (mype == 0) then
          write(logunit,F01) 'ERROR: invalid decomp_type = ',decomp_type
       end if
       call shr_sys_abort(subName//'invalid decomp_type')
    endif

    if (nxg*nyg == 0) then
       lsize = 0
       allocate(gbuf(lsize,dead_grid_total))
       !      gbuf = -888.0_R8
       if (mype == 0) then
          write(logunit,*) subname,' grid size is zero, lsize = ',lsize
       end if
       return
    endif

    found = .false.

    if (decomp_type == 1) then  ! 1d decomp by lat
       npesx = 1
       npesy = totpe
       found = .true.
    elseif (decomp_type == 2) then  ! 1d decomp by lon
       npesx = totpe
       npesy = 1
       found = .true.
    elseif (decomp_type == 3) then  ! 2d decomp
       if (present(nproc_x)) then
          if ( nproc_x > 0 ) then
             npesx=nproc_x
             npesy=totpe/npesx
             if ( npesx*npesy /= totpe) then
                write(logunit,F00) 'ERROR: uneven decomposition'
                call shr_sys_abort(subName//'uneven decomp')
             end if
             found = .true.
          endif
       endif
       if (.not.found) then  ! narrow blocks
          do nx = 1,totpe
             ny = totpe/nx
             if (nx*ny == totpe) then
                npesx = nx
                npesy = ny
                found = .true.
             endif
          enddo
       endif
    elseif (decomp_type == 4) then  ! 2d evenly divisible square block decomp
       hscore = nxg*nyg
       do nx = 1,totpe
          ny = totpe/nx
          if (nx*ny == totpe .and. mod(nxg,nx) == 0 .and. mod(nyg,ny) == 0) then
             bscore = ((nxg*ny*1.0_r8) / (nyg*nx*1.0_r8)) - 1.0_r8
             bscore = bscore * bscore
             if (bscore < hscore .or. .not.found) then
                hscore = bscore
                npesx = nx
                npesy = ny
                found = .true.
             endif
          endif
       enddo
    endif

    if (found) then
       nx = nxg/npesx
       mypex = mod(mype,npesx)
       mypey = mype/npesx
       is = (mypex    ) * (nx) + 1
       ie = (mypex + 1) * (nx)

       ny = nyg/npesy
       js = (mypey    ) * (ny) + 1
       je = (mypey + 1) * (ny)

       nxp = nxg - (nx*npesx)       ! extra lons not accounted for yet
       nyp = nyg - (ny*npesy)       ! extra lats not accounted for yet

       is = is + min(mypex,nxp)      ! add them to first few pes and shift everything
       ie = ie + min(mypex+1,nxp)
       js = js + min(mypey,nyp)      ! add them to first few pes and shift everything
       je = je + min(mypey+1,nyp)

       lsize = (ie - is + 1) * (je - js + 1)

       allocate(gindex(lsize))
       n = 0
       do j = js,je
          do i = is,ie
             n = n + 1
             gindex(n) = (j-1)*nxg + i
          enddo
       enddo
    endif

    if (.not.found) then
       !-------------------------------------------------------------------------
       ! type 11 general segment decomp
       !-------------------------------------------------------------------------
       nx = nxg*nyg / (totpe*13) + 1   ! 13 segments per pe (arbitrary)
       ! nx override with seg_len
       if (present(seg_len)) then
          if (seg_len > 0) nx = seg_len
       endif

       n = 0
       i = 0
       lsize = 0
       do while (n < nxg*nyg)
          ny = min(nx,nxg*nyg-n)
          do j = 1,ny
             n = n + 1
             if (mype == mod(i,totpe)) then
                lsize = lsize + 1
             endif
          enddo
          i = i + 1
       enddo

       allocate(gindex(lsize))

       n = 0
       i = 0
       lsize = 0
       do while (n < nxg*nyg)
          ny = min(nx,nxg*nyg-n)
          do j = 1,ny
             n = n + 1
             if (mype == mod(i,totpe)) then
                lsize = lsize + 1
                gindex(lsize) = n
             endif
          enddo
          i = i + 1
       enddo

       if (mype == 0) then
          write(logunit,*) 'dead_setNewGrid decomp seg ',mype,lsize,nx
       end if

       found = .true.

    endif

    if ( .not.found ) then
       write(logunit,F01) 'ERROR: with decomp nxg,nyg,totpe=',nxg,nyg,totpe
       call shr_sys_abort(subName//'decomp')
    endif

    deg2rad = shr_const_pi / 180.0_R8
    re = shr_const_rearth

    allocate(gbuf(lsize,dead_grid_total))
    gbuf = -888.0_R8
    if (mype == 0) then
       write(logunit,*) subname,' Decomp is ',decomp_type,' lsize = ',lsize
    end if

    n=0
    dx = 360.0_R8/nxg * deg2rad
    do n = 1,lsize
       ig = mod((gindex(n)-1),nxg) + 1
       jg =     (gindex(n)-1)/nxg  + 1

       ys = -90.0_R8 + (jg-1.0_R8)*180.0_R8/(nyg)
       yc = -90.0_R8 + (jg-0.5_R8)*180.0_R8/(nyg)
       yn = -90.0_R8 + (jg-0.0_R8)*180.0_R8/(nyg)
       dy = sin(yn*deg2rad) - sin(ys*deg2rad)
       area = dx*dy*re*re

       gbuf(n,dead_grid_lon  ) = (ig-1.0_R8)*360.0_R8/(nxg)
       gbuf(n,dead_grid_lat  ) = yc
       gbuf(n,dead_grid_index) = gindex(n)
       gbuf(n,dead_grid_area ) = area
       gbuf(n,dead_grid_mask ) = 1
       gbuf(n,dead_grid_frac ) = 1.0_R8
    enddo

    deallocate(gindex)

  end subroutine dead_setNewGrid

end module dead_mod
