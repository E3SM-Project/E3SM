module restart_dynamics

!----------------------------------------------------------------------- 
! 
! Purpose: Read and write restart file
!
! !HISTORY:
!   2000.01.01   ??????     Creation
!   2006.04.13   Sawyer     Removed dependency on prognostics
!   2006.07.01   Sawyer     Transitioned q3 tracers to T_TRACERS
!   2008.10.02   Edwards    Added pio support
!
!-----------------------------------------------------------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8, r4=> shr_kind_r4
  use pmgrid,         only: plon, plat, plev
  use constituents,   only: pcnst
  use dyn_comp,       only: dyn_import_t, dyn_export_t
  use abortutils,     only: endrun
  use fv_control_mod, only: tmass0
  use dynamics_vars,  only: T_FVDYCORE_STATE
  use cam_logfile,    only: iulog
#if ( defined OFFLINE_DYN )
  use metdata,        only: write_met_restart, read_met_restart
#endif
  use pio,            only: file_desc_t, io_desc_t, var_desc_t, &
                             pio_double, pio_unlimited, pio_offset, &
                             pio_setdebuglevel, pio_setframe, &
                             pio_def_var, pio_def_dim, &
                             pio_inq_varid, &
                             pio_put_var, pio_get_var, &
                             pio_write_darray, pio_read_darray, &
                             pio_initdecomp, pio_freedecomp

  implicit none
  
  public read_restart_dynamics, write_restart_dynamics
  private
  
  real(r8), parameter ::  D0_0                    =  0.0_r8
  real(r8), parameter ::  D1E5                    =  1.0e5_r8
  
  public init_restart_dynamics
  integer, parameter :: namlen=16
  
  type restart_var_t
     real(r8), pointer :: v2d(:,:)
     real(r8), pointer :: v3d(:, :, :)
     real(r8), pointer :: v4d(:, :, :, :)
     real(r8), pointer :: v5d(:, :, :, :, :)
     real(r4), pointer :: v2dr4(:,:)
     real(r4), pointer :: v3dr4(:, :, :)
     real(r4), pointer :: v4dr4(:, :, :, :)
     real(r4), pointer :: v5dr4(:, :, :, :, :)
     
     type(var_desc_t), pointer  :: vdesc
     integer           :: ndims
     integer           :: timelevels
     character(len=namlen) :: name
  end type restart_var_t
  integer, parameter :: restartvarcnt = 6+pcnst

  type(var_desc_t) :: timedesc, tmass0desc

  type(restart_var_t) :: restartvars(restartvarcnt)
  logical ::  restart_varlistw_initialized=.false.
  logical ::  restart_varlistr_initialized=.false.

CONTAINS

  subroutine set_r_var(name, timelevels, index, v2, v3, v4, v5, v2r4, v3r4, v4r4, v5r4)
    use abortutils,      only: endrun


    character(len=*), intent(in) :: name
    integer, intent(in) :: timelevels, index
    real(r8), target, optional :: v2(:,:), v3(:,:,:), v4(:,:,:,:), v5(:,:,:,:,:)
    real(r4), target, optional :: v2r4(:,:), v3r4(:,:,:), v4r4(:,:,:,:), v5r4(:,:,:,:,:)

    restartvars(index)%name=name
    restartvars(index)%timelevels = timelevels

    if(present(v2)) then
       restartvars(index)%ndims = 2
       restartvars(index)%v2d => v2
    else if(present(v3)) then
       restartvars(index)%ndims = 3
       restartvars(index)%v3d => v3
    else if(present(v4)) then
       restartvars(index)%ndims = 4
       restartvars(index)%v4d => v4
    else if(present(v5)) then
       restartvars(index)%ndims = 5
       restartvars(index)%v5d => v5
    else if(present(v2r4)) then
       restartvars(index)%ndims = 2
       restartvars(index)%v2dr4 => v2r4
    else if(present(v3r4)) then
       restartvars(index)%ndims = 3
       restartvars(index)%v3dr4 => v3r4
    else if(present(v4r4)) then
       restartvars(index)%ndims = 4
       restartvars(index)%v4dr4 => v4r4
    else if(present(v5r4)) then
       restartvars(index)%ndims = 5
       restartvars(index)%v5dr4 => v5r4
    else
       call endrun('bad ndims in call to set_r_var')
    end if



    allocate(restartvars(index)%vdesc)

  end subroutine set_r_var
  
  subroutine init_restart_varlistw(dyn_out)
    use abortutils,      only: endrun
    use dyn_comp, only  : dyn_export_t
    use constituents, only : cnst_name

    type(dyn_export_t) :: dyn_out
    integer :: vcnt=1
    integer :: i, m


! Should only be called once
    if(restart_varlistw_initialized) return
    restart_varlistw_initialized=.true.
    call set_r_var('PHIS', 1, vcnt, v2=dyn_out%phis)

    vcnt=vcnt+1
    call set_r_var('U', 1, vcnt, v3=dyn_out%u3s)

    vcnt=vcnt+1
    call set_r_var('V', 1, vcnt, v3=dyn_out%v3s)

    vcnt=vcnt+1
    call set_r_var('DELP', 1, vcnt, v3=dyn_out%delp)

    vcnt=vcnt+1
    call set_r_var('PT', 1, vcnt, v3=dyn_out%pt)
    
    do m=1,pcnst
       vcnt=vcnt+1
       call set_r_var(cnst_name(m), 1, vcnt, v3=dyn_out%tracer(:,:,:,m) )
    end do

    vcnt=vcnt+1
    call set_r_var('PS', 1, vcnt, v2=dyn_out%ps )
    

    if(vcnt.ne.restartvarcnt) then
       write(iulog,*) 'vcnt= ',vcnt, ' restartvarcnt=',restartvarcnt
       call endrun('bad restartvarcnt')
    end if

  end subroutine init_restart_varlistw

  subroutine init_restart_varlistr(dyn_in)
    use abortutils,      only: endrun
    use dyn_comp, only  : dyn_import_t
    use constituents, only : cnst_name

    type(dyn_import_t) :: dyn_in
    integer :: vcnt=1
    integer :: i, m


! Should only be called once
    if(restart_varlistr_initialized) return
    restart_varlistr_initialized=.true.
    call set_r_var('PHIS', 1, vcnt, v2=dyn_in%phis)

    vcnt=vcnt+1
    call set_r_var('U', 1, vcnt, v3=dyn_in%u3s)

    vcnt=vcnt+1
    call set_r_var('V', 1, vcnt, v3=dyn_in%v3s)

    vcnt=vcnt+1
    call set_r_var('DELP', 1, vcnt, v3=dyn_in%delp)

    vcnt=vcnt+1
    call set_r_var('PT', 1, vcnt, v3=dyn_in%pt)
    
    do m=1,pcnst
       vcnt=vcnt+1
       call set_r_var(cnst_name(m), 1, vcnt, v3=dyn_in%tracer(:,:,:,m) )
    end do

    vcnt=vcnt+1
    call set_r_var('PS', 1, vcnt, v2=dyn_in%ps )
    

    if(vcnt.ne.restartvarcnt) then
       write(iulog,*) 'vcnt= ',vcnt, ' restartvarcnt=',restartvarcnt
       call endrun('bad restartvarcnt')
    end if

  end subroutine init_restart_varlistr



  subroutine init_restart_dynamics(File, hdimids, dyn_out)

    use dyn_comp, only: dyn_export_t
    use dyn_grid, only: get_horiz_grid_dim_d
    use hycoef,   only: init_restart_hycoef

    ! Input arguments
    type(File_desc_t),  intent(inout) :: File
    integer,            pointer       :: hdimids(:)
    type(Dyn_export_t), intent(in)    :: dyn_out

    integer :: vdimids(2)
    integer :: hdim1, hdim2
    character(len=namlen) :: name

    integer :: alldims(3), alldims2d(2)
    integer :: i, ierr, timelevels
    type(var_desc_t), pointer :: vdesc
    integer :: ndims
    
    call init_restart_hycoef(File, vdimids)

    call get_horiz_grid_dim_d(hdim1, hdim2)
    allocate(hdimids(2))
    ierr = PIO_Def_Dim(File, 'lon',hdim1, hdimids(1))
    ierr = PIO_Def_Dim(File, 'lat',hdim2, hdimids(2))

    ierr = PIO_Def_Var(File, 'time', pio_double, timedesc)

    ierr = PIO_Def_var(File, 'tmass0', pio_double, tmass0desc)

    alldims(1:2) = hdimids(1:2)
    alldims(3) = vdimids(1)

    alldims2d(1:2) = hdimids(1:2)

    call init_restart_varlistw(dyn_out)

    do i=1,restartvarcnt
    
       call get_restart_var(i, name, timelevels, ndims, vdesc)
       if(timelevels>1) then
          call endrun('not expecting timelevels>1 in fv dycore')
       else
          if(ndims==1) then
! broken i think
             ierr = PIO_Def_Var(File, name, pio_double, hdimids(2:2), vdesc)
          else if(ndims==2) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims2d(1:2), vdesc)
          else if(ndims==3) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims(1:3), vdesc)
          end if
       end if
    end do
#if ( defined OFFLINE_DYN )
      call write_met_restart( File )
#endif


  end subroutine init_restart_dynamics


  subroutine write_restart_dynamics (File, dyn_out)
    use cam_pio_utils, only : pio_subsystem
    use dyn_comp, only  : dyn_export_t
    use time_manager, only: get_curr_time
    use pmgrid,          only: plon, plat
    use dyn_grid, only : get_horiz_grid_dim_d
    use dyn_internal_state, only : get_dyn_state	
    use hycoef, only: write_restart_hycoef

    !
    ! Input arguments
    !
    type(File_desc_t), intent(inout) :: File     ! Unit number
    type(Dyn_export_t), intent(in) :: dyn_out ! Not used in eul dycore

    !
    ! Local workspace
    !
    type(io_desc_t) :: iodesc2d, iodesc3d
    integer :: ierr   ! error status
    integer :: ndcur, nscur
    real(r8) :: time, mold(1), null(0)
    integer :: i, s3d(1), s2d(1), ct
    integer(kind=pio_offset) :: t
    integer :: ndims, isize(1), timelevels
    type(var_desc_t), pointer :: vdesc
    integer :: hdim1, hdim2
    integer, pointer :: ldof(:)
    character(len=namlen) :: name
    logical :: use_transfer
    type (T_FVDYCORE_STATE), pointer :: dyn_state
    
    call write_restart_hycoef(File)

    !
    ! transfer is the fastest method to flatten the multidimensional arrays into the 1d needed by pio
    ! but it doesn't work correctly if the array is zero length...

    use_transfer = .true.
#if ( defined SPMD )
    dyn_state => get_dyn_state()
    if (dyn_state%grid%iam >= dyn_state%grid%npes_xy) then
       use_transfer = .false.
    end if       
#endif

    call get_curr_time(ndcur, nscur)
    call get_horiz_grid_dim_d(hdim1, hdim2)


    ldof => get_restart_decomp(hdim1, hdim2, plev)
    call pio_initdecomp(pio_subsystem, pio_double, (/hdim1, hdim2, plev/), ldof, iodesc3d)
    deallocate(ldof)
    ldof => get_restart_decomp(hdim1, hdim2, 1)
    call pio_initdecomp(pio_subsystem, pio_double, (/hdim1, hdim2/), ldof, iodesc2d)
    deallocate(ldof)


    ierr = pio_put_var(File, tmass0desc, (/tmass0/))

    time = ndcur+(real(nscur,kind=r8))/86400_r8
    ierr = pio_put_var(File,timedesc%varid, time)

    do i=1,restartvarcnt
       call get_restart_var(i, name, timelevels, ndims, vdesc)
       if(ndims==2) then
          if(use_transfer) then
             call pio_write_darray(File, vdesc, iodesc2d, transfer(restartvars(i)%v2d(:,:), mold), ierr)
          else
             call pio_write_darray(File, vdesc, iodesc2d, null, ierr)
          end if
       else if(ndims==3) then
          if(use_transfer) then
             call pio_write_darray(File, vdesc, iodesc3d, transfer(restartvars(i)%v3d(:,:,:), mold), ierr)
          else
             call pio_write_darray(File, vdesc, iodesc3d, null, ierr)
          end if
       end if
    end do
    call pio_freedecomp(File, iodesc2d)
    call pio_freedecomp(File, iodesc3d)

  end subroutine write_restart_dynamics


!
! Get the integer mapping of a variable in the dynamics decomp in memory.  
! The canonical ordering is as on the file. A 0 value indicates that the
! variable is not on the file (eg halo or boundary values)
!
  function get_restart_decomp(hdim1, hdim2, nlev) result(ldof)
    use dyn_grid, only : get_dyn_grid_parm
    integer, intent(in) :: hdim1, hdim2, nlev
    integer, pointer :: ldof(:)
    integer :: i, k, j
    integer :: lcnt

    integer :: bfirst, blast, ncols

    integer :: beglatxy, beglonxy, endlatxy, endlonxy


    beglonxy = get_dyn_grid_parm('beglonxy')
    endlonxy = get_dyn_grid_parm('endlonxy')

    beglatxy = get_dyn_grid_parm('beglatxy')
    endlatxy = get_dyn_grid_parm('endlatxy')

    lcnt=(endlatxy-beglatxy+1)*nlev*(endlonxy-beglonxy+1)
    allocate(ldof(lcnt))
    lcnt=0
    ldof(:)=0	

    do k=1,nlev
       do j=beglatxy,endlatxy
          do i=beglonxy, endlonxy
             lcnt=lcnt+1
             ldof(lcnt)=i+(j-(plat-hdim2+1))*hdim1+(k-1)*hdim1*hdim2
          end do
       end do
    end do

  end function get_restart_decomp


  subroutine get_restart_var(i,name, timelevels, ndims, vdesc)
    integer, intent(in) :: i
    character(len=namlen), intent(out) :: name
    integer, intent(out) :: ndims, timelevels
    type(var_desc_t), pointer :: vdesc

    name = restartvars(i)%name
    timelevels = restartvars(i)%timelevels
    ndims = restartvars(i)%ndims
    if(.not.associated(restartvars(i)%vdesc)) then
       allocate(restartvars(i)%vdesc)
    end if
    vdesc => restartvars(i)%vdesc
    call pio_setframe(vdesc,int(-1,kind=pio_offset))
  end subroutine get_restart_var

  subroutine read_restart_dynamics(File, dyn_in, dyn_out, NLFileName)
    use cam_pio_utils, only : pio_subsystem

    use pmgrid
    use dyn_comp, only : dyn_init
    use dyn_internal_state, only : get_dyn_state


    type(file_desc_t) :: File
    type(dyn_export_t) :: dyn_out
    type(dyn_import_t) :: dyn_in

    character(len=*), intent(in) :: NLFileName

    real(r8), allocatable :: tmp(:)
    integer ::  dims3d(3)
    integer :: ndims
    integer :: ierr
    character(len=namlen) :: name
    integer, pointer :: ldof(:)

    type (T_FVDYCORE_STATE), pointer :: dyn_state
    integer :: timelevels ! not used in fv
    type(var_desc_t), pointer :: vdesc
    type(io_desc_t) :: iodesc2d, iodesc3d
    integer :: s2d, s3d, i

#if ( defined OFFLINE_DYN )
    call read_met_restart( File )
#endif

    !
    ! Initialize the dynamics
    !
    dyn_state => get_dyn_state()

    call dyn_init(file, dyn_state, dyn_in, dyn_out, NLFileName )

    dims3d(1)=(endlonxy-beglonxy+1)
    dims3d(2)=(endlatxy-beglatxy+1)
    dims3d(3)=plev
    s2d = dims3d(1)*dims3d(2)
    s3d = s2d*dims3d(3)
    allocate(tmp(s3d))
    call init_restart_varlistr(dyn_in)

    ierr = PIO_Inq_varid(File, 'tmass0', tmass0desc)
    ierr = pio_get_var(File, tmass0desc, tmass0)


    ldof => get_restart_decomp(plon, plat, plev)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon, plat, plev/), ldof, iodesc3d)
    deallocate(ldof)
    ldof => get_restart_decomp(plon, plat, 1)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon, plat/), ldof, iodesc2d)
    deallocate(ldof)


    do i=1,restartvarcnt
       call get_restart_var(i, name, timelevels, ndims, vdesc)


       ierr = PIO_Inq_varid(File, name, vdesc)

       if(ndims==2) then
          call pio_read_darray(File, vdesc, iodesc2d, tmp(1:s2d), ierr)
          restartvars(i)%v2d(:,:) = reshape(tmp(1:s2d), dims3d(1:2))
       else if(ndims==3) then
          call pio_read_darray(File, restartvars(i)%vdesc, iodesc3d, tmp(1:s3d), ierr)
          restartvars(i)%v3d(:,:,:) = reshape(tmp(1:s3d), dims3d)
       end if
    end do
    deallocate(tmp)
    call pio_freedecomp(File, iodesc2d)
    call pio_freedecomp(File, iodesc3d)


  end subroutine read_restart_dynamics



end module restart_dynamics
