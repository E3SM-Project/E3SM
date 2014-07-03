module restart_dynamics

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use constituents,    only: pcnst
   use prognostics,     only: u3, v3, q3, qm1, div, dps, dpsl, dpsm, phis, phisl, &
                              phism, ps, urhs, vrhs, trhs, prhs, ql, qm, &
                              tl, tm, omga, n3, n3m1, parrsld, t3, etadot, &
                              initialize_prognostics
   use pmgrid,          only: plon, plat, plnlv, plevp, plev, beglat, endlat
   use spmd_utils,      only : masterproc
   use scanslt,         only: slt_alloc, qfcst, hw1, hw2, hw3, alpha
#ifdef SPMD
   use pspect,          only: psp
#endif
   use comspe,          only: lnpstar
   use abortutils,      only: endrun
   use sld_control_mod, only: tmass0
   use cam_logfile,     only: iulog
   use dyn_comp,        only: dyn_import_t, dyn_export_t
   use pio,             only: var_desc_t, file_desc_t, pio_double, pio_unlimited, &
         pio_def_var, pio_def_dim, io_desc_t, pio_offset, pio_put_var, &
         pio_write_darray, pio_setdebuglevel, pio_setframe, pio_freedecomp, &
         pio_initdecomp, pio_read_darray, pio_inq_varid, pio_get_var

   implicit none

   private
   save

!
! Public interfaces
!
  public :: read_restart_dynamics, init_restart_dynamics, write_restart_dynamics

  integer, parameter :: namlen=16

  type restart_var_t
     real(r8), pointer :: v1d(:) => null()
     real(r8), pointer :: v2d(:,:) => null()
     real(r8), pointer :: v3d(:, :, :) => null()
     real(r8), pointer :: v4d(:, :, :, :) => null()
     real(r8), pointer :: v5d(:, :, :, :, :) => null()

     type(var_desc_t), pointer  :: vdesc => null()
     integer           :: ndims
     integer           :: timelevels
     character(len=namlen) :: name
  end type restart_var_t
  integer, parameter :: restartvarcnt = 25
  type(var_desc_t) :: timedesc, tmass0desc, hw1desc, hw2desc, hw3desc, &
       alphadesc, lnpstardesc

  type(restart_var_t) :: restartvars(restartvarcnt)
  logical ::  restart_varlist_initialized=.false.


CONTAINS

  subroutine set_r_var(name, timelevels, index, v1, v2, v3, v4, v5)
    use abortutils,      only: endrun

    character(len=*), intent(in) :: name
    integer, intent(in) :: timelevels, index
    real(r8), target, optional :: v1(:), v2(:,:), v3(:,:,:), v4(:,:,:,:), v5(:,:,:,:,:)

    restartvars(index)%name=name
    restartvars(index)%timelevels = timelevels
    if(present(v1)) then
       restartvars(index)%ndims = 1
       restartvars(index)%v1d => v1
    else if(present(v2)) then
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
    else
       call endrun('bad ndims in call to set_r_var')
    end if
    allocate(restartvars(index)%vdesc)

  end subroutine set_r_var

  subroutine init_restart_varlist()
    use abortutils,      only: endrun
    use prognostics,     only: u3, v3, q3, qm1, div, dps, dpsl, dpsm, phis, phisl, &
         phism, ps, urhs, vrhs, trhs, prhs, ql, qm, &
         tl, tm, omga, parrsld, t3, etadot, ptimelevels
    use scanslt, only : qfcst
    integer :: vcnt=1

    ! Should only be called once
    if(restart_varlist_initialized) return
    restart_varlist_initialized=.true.

    call set_r_var('DIV', ptimelevels, vcnt, v4=div)

    vcnt=vcnt+1
    call set_r_var('URHS', 1, vcnt, v3=urhs)
    vcnt=vcnt+1
    call set_r_var('VRHS', 1, vcnt, v3=vrhs)
    vcnt=vcnt+1
    call set_r_var('TRHS', 1, vcnt, v3=trhs)
    vcnt=vcnt+1
    call set_r_var('PRHS', 1, vcnt, v3=prhs)


    vcnt=vcnt+1
    call set_r_var('DPSL', 1, vcnt, v2=dpsl)
    vcnt=vcnt+1
    call set_r_var('DPSM', 1, vcnt, v2=dpsm)

    vcnt=vcnt+1
    call set_r_var('QL', 1, vcnt, v3=ql)
    vcnt=vcnt+1
    call set_r_var('QM', 1, vcnt, v3=qm)
    vcnt=vcnt+1
    call set_r_var('DPS', 1, vcnt, v2=dps)
    vcnt=vcnt+1
    call set_r_var('PHIS', 1, vcnt, v2=phis)
    vcnt=vcnt+1
    call set_r_var('PHISL', 1, vcnt, v2=phisl)
    vcnt=vcnt+1
    call set_r_var('PHISM', 1, vcnt, v2=phism)


    vcnt=vcnt+1
    call set_r_var('TL', 1, vcnt, v3=tl)
    vcnt=vcnt+1
    call set_r_var('TM', 1, vcnt, v3=tm)
    vcnt=vcnt+1
    call set_r_var('OMEGA', 1, vcnt, v3=omga)

    vcnt=vcnt+1
    call set_r_var('U', ptimelevels, vcnt, v4=u3)

    vcnt=vcnt+1
    call set_r_var('V', ptimelevels, vcnt, v4=v3)

    vcnt=vcnt+1
    call set_r_var('T', ptimelevels, vcnt, v4=t3)

    vcnt=vcnt+1
    call set_r_var('PS', ptimelevels, vcnt, v3=ps)

    vcnt=vcnt+1
    call set_r_var('Q', ptimelevels, vcnt, v5=Q3 )

    vcnt=vcnt+1
    call set_r_var('QM1', 1, vcnt, v4=QM1 )
    vcnt=vcnt+1
    call set_r_var('QFCST', 1, vcnt, v4=Qfcst )


    vcnt=vcnt+1
    call set_r_var('PARRSLD', 1, vcnt, v3=parrsld )

    vcnt=vcnt+1
    call set_r_var('ETADOT', 1, vcnt, v3=etadot )


    if(vcnt.ne.restartvarcnt) then
       write(iulog,*) 'vcnt= ',vcnt, ' restartvarcnt=',restartvarcnt
       call endrun('bad restartvarcnt')
    end if

  end subroutine init_restart_varlist



  subroutine init_restart_dynamics(File, hdimids, dyn_out)

    use dyn_comp,        only: dyn_export_t
    use dyn_grid,        only: get_horiz_grid_dim_d
    use constituents,    only: pcnst
    use prognostics,     only: u3, v3, q3, qm1, div, dps, dpsl, dpsm, phis, phisl, &
         phism, ps, urhs, vrhs, trhs, prhs, ql, qm, &
         tl, tm, omga, parrsld, t3, etadot
    use pspect,          only: psp
    use hycoef,          only: init_restart_hycoef

    ! Input arguments
    type(File_desc_t),  intent(inout) :: File
    integer,            pointer       :: hdimids(:)
    type(Dyn_export_t), intent(in)    :: dyn_out

    character(len=namlen) :: name

    integer :: vdimids(2)
    integer :: alldims(4), alldims2d(3), qdims(5), pspdim
    integer :: timelevels_dimid, i, hdim1, hdim2, ierr
    type(var_desc_t), pointer :: vdesc
    integer :: ndims, timelevels

    call init_restart_hycoef(File, vdimids)

    call get_horiz_grid_dim_d(hdim1, hdim2)
    allocate(hdimids(2))
    ierr = PIO_Def_Dim(File, 'lon',hdim1, hdimids(1))

    ierr = PIO_Def_Dim(File, 'lat',hdim2, hdimids(2))

    ierr = PIO_Def_Dim(File,'timelevels',PIO_UNLIMITED,timelevels_dimid)

    ierr = PIO_Def_Dim(File,'pcnst',pcnst, qdims(4))

    ierr = PIO_Def_Dim(File,'psp',psp, pspdim)

    ierr = PIO_Def_Var(File, 'time', pio_double, (/timelevels_dimid/), timedesc)

    ierr = PIO_Def_var(File, 'tmass0', pio_double, tmass0desc)
    ierr = PIO_Def_var(File, 'hw1', pio_double, qdims(4:4), hw1desc)
    ierr = PIO_Def_var(File, 'hw2', pio_double, qdims(4:4), hw2desc)
    ierr = PIO_Def_var(File, 'hw3', pio_double, qdims(4:4), hw3desc)
    ierr = PIO_Def_var(File, 'alpha', pio_double, qdims(4:4), alphadesc)
    ierr = PIO_Def_var(File, 'lnpstar', pio_double, (/pspdim/), lnpstardesc)




    alldims(1:2) = hdimids(1:2)
    alldims(3) = vdimids(1)
    alldims(4) = timelevels_dimid

    alldims2d(1:2) = hdimids(1:2)
    alldims2d(3) = timelevels_dimid

    qdims(1:2) = hdimids(1:2)
    qdims(3) = vdimids(1)
    qdims(5) = timelevels_dimid

    call init_restart_varlist()

    do i=1,restartvarcnt

       call get_restart_var(i, name, timelevels, ndims, vdesc)
       if(timelevels>1) then
          if(ndims==3) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims2d, vdesc)
          else if(ndims==4) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims, vdesc)
          else if(ndims==5) then
             ierr = PIO_Def_Var(File, name, pio_double, qdims, vdesc)
          end if
       else
          if(ndims==1) then
             ! broken i think
             ierr = PIO_Def_Var(File, name, pio_double, hdimids(2:2), vdesc)
          else if(ndims==2) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims2d(1:2), vdesc)
          else if(ndims==3) then
             if(name.eq.'ETADOT') then
                ! this is the only plevp variable written
                alldims(3)=vdimids(2)
                ierr = PIO_Def_Var(File, name, pio_double, alldims(1:3), vdesc)
                alldims(3)=vdimids(1)
             else
                ierr = PIO_Def_Var(File, name, pio_double, alldims(1:3), vdesc)
             end if

          else if(ndims==4) then
             ierr = PIO_Def_Var(File, name, pio_double, qdims(1:4), vdesc)
          end if
       end if
    end do


  end subroutine init_restart_dynamics

  subroutine write_restart_dynamics (File, dyn_out)
    use dyn_comp,        only: dyn_export_t
    use constituents,    only: pcnst
    use prognostics,     only: n3, n3m1, ptimelevels
    use sld_control_mod, only: tmass0
    use scanslt,         only:  hw1, hw2, hw3, alpha   
    use comspe,          only: lnpstar
    use pmgrid,          only: plon, plat, plevp, plev, beglat, endlat
    use time_manager,    only: get_curr_time, get_step_size
    use cam_pio_utils,   only: pio_subsystem
    use hycoef,          only: write_restart_hycoef

    !
    ! Input arguments
    !
    type(File_desc_t), intent(inout) :: File     ! Unit number
    type(Dyn_export_t), intent(in) :: dyn_out ! Not used in sld dycore

    !
    ! Local workspace
    !
    integer :: ierr   ! error status
    integer :: ndcur, nscur
    real(r8) :: time, dtime, mold(1)
    integer :: i, s3d(1), s2d(1), ct
    integer(kind=pio_offset) :: t
    type(io_desc_t) :: iodesc3d, iodesc2d, iodesc4d, iodesc3dp1
    integer :: ndims, timelevels
    type(var_desc_t), pointer :: vdesc
    character(len=namlen) :: name
    integer, pointer :: ldof(:)
    !

    call write_restart_hycoef(File)

    call get_curr_time(ndcur, nscur)
    dtime = get_step_size()

    ldof => get_restart_decomp(plon, plat, plev)

    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat,plev/), ldof, iodesc3d)
    deallocate(ldof)

    ldof => get_restart_decomp(plon, plat, plevp)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat,plevp/), ldof, iodesc3dp1)
    deallocate(ldof)

    ldof => get_restart_decomp(plon, plat, plev*pcnst)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat,plev,pcnst/), ldof, iodesc4d)
    deallocate(ldof)

    ldof => get_restart_decomp(plon, plat, 1)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat/), ldof, iodesc2d)
    deallocate(ldof)

    ierr = pio_put_var(File, tmass0desc, (/tmass0/))

    ierr = pio_put_var(File, hw1desc, hw1)
    ierr = pio_put_var(File, hw2desc, hw2)

    ierr = pio_put_var(File, hw3desc, hw3)
    ierr = pio_put_var(File, alphadesc, alpha)
    ierr = pio_put_var(File, lnpstardesc, lnpstar)

    do t=1,ptimelevels
       time = ndcur+(real(nscur,kind=r8)+ (t-2)*dtime)/86400_r8
       ierr = pio_put_var(File,timedesc%varid, (/int(t)/), time)
    end do
    do i=1,restartvarcnt
       call get_restart_var(i, name, timelevels, ndims, vdesc)
       if(timelevels==1) then
          if(ndims==2) then
             call pio_write_darray(File, vdesc, iodesc2d, transfer(restartvars(i)%v2d(:,:), mold), ierr)
          else if(ndims==3) then
             if(name.eq.'ETADOT') then
                call pio_write_darray(File, vdesc, iodesc3dp1, transfer(restartvars(i)%v3d(:,:,:), mold), ierr)
             else
                call pio_write_darray(File, vdesc, iodesc3d, transfer(restartvars(i)%v3d(:,:,:), mold), ierr)
             end if
          else if(ndims==4) then
             call pio_write_darray(File, vdesc, iodesc4d, transfer(restartvars(i)%v4d(:,:,:,:), mold), ierr)
          end if
       else
          do t=1,timelevels
             if(t==1) ct=n3m1
             if(t==2) ct=n3

             call pio_setframe(vdesc, t)
             if(ndims==3) then
                call pio_write_darray(File, vdesc, iodesc2d, transfer(restartvars(i)%v3d(:,:,ct), mold), ierr)
             else if(ndims==4) then
                call pio_write_darray(File, vdesc, iodesc3d, transfer(restartvars(i)%v4d(:,:,:,ct), mold), ierr)
             else if(ndims==5) then
                call pio_write_darray(File, vdesc, iodesc4d, transfer(restartvars(i)%v5d(:,:,:,:,ct), mold), ierr)
             end if

          end do

       end if
    end do

    call pio_freedecomp(File,iodesc2d)
    call pio_freedecomp(File,iodesc3d)
    call pio_freedecomp(File,iodesc3dp1)
    call pio_freedecomp(File,iodesc4d)


    return
  end subroutine write_restart_dynamics

  function get_restart_decomp(hdim1, hdim2, nlev) result(ldof)
    use dyn_grid, only : get_dyn_grid_parm

    integer, intent(in) :: hdim1, hdim2, nlev
    integer, pointer :: ldof(:)
    integer :: i, k, j
    integer :: lcnt
    integer, allocatable :: gcols(:)

    integer :: beglatxy, beglonxy, endlatxy, endlonxy, plat


    beglonxy = get_dyn_grid_parm('beglonxy')
    endlonxy = get_dyn_grid_parm('endlonxy')
    beglatxy = get_dyn_grid_parm('beglatxy')
    endlatxy = get_dyn_grid_parm('endlatxy')

    plat = get_dyn_grid_parm('plat')
    
    
    lcnt=(endlatxy-beglatxy+1)*nlev*(endlonxy-beglonxy+1)

    allocate(ldof(lcnt))
    lcnt=0
    ldof(:)=0	
    do j=beglatxy,endlatxy
       do k=1,nlev
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
    call pio_setframe(vdesc, int(-1,pio_offset))

  end subroutine get_restart_var

  !#######################################################################

  subroutine read_restart_dynamics (File, dyn_in, dyn_out, NLFileName)
    use dyn_comp,        only: dyn_init, dyn_import_t, dyn_export_t
    use cam_pio_utils,   only: pio_subsystem
    use dyn_comp,        only: dyn_init
    use prognostics,     only: initialize_prognostics, n3, n3m1
    use pmgrid,          only: plon, plat, plevp, plev, beglat, endlat
    use constituents,    only: pcnst
    use sld_control_mod, only: tmass0
    use scanslt,         only: slt_alloc, hw1, hw2, hw3, alpha   
    use comspe,          only: lnpstar

    !
    ! Input arguments
    !
    type(file_desc_t), intent(inout) :: File     ! PIO file handle
    type(dyn_import_t) :: dyn_in    ! not used by this dycore, included for compatibility
    type(dyn_export_t) :: dyn_out ! not used by this dycore, included for compatibility    
    character(len=*), intent(in) :: NLFileName

    !
    ! Local workspace
    !

    integer :: dims4d(4), dims3d(3), dims2d(2), s2d, s3d, s4d
    real(r8), allocatable :: tmp(:)
    type(io_desc_t) :: iodesc3d, iodesc2d, iodesc4d, iodesc3dp1
    integer :: ndims
    type(var_desc_t), pointer :: vdesc
    character(len=namlen) :: name
    integer :: ierr
    integer :: timelevels
    integer :: i, ct
    integer(kind=pio_offset) :: t
    integer, pointer :: ldof(:)

    call dyn_init(file, NLFileName)

    call initialize_prognostics
    call slt_alloc()

    dims4d(1) = plon
    dims4d(2) = plev
    dims4d(3) = pcnst
    dims4d(4) = endlat-beglat+1
    s4d=dims4d(1)*dims4d(2)*dims4d(3)*dims4d(4)
    dims3d(1) = plon
    dims3d(2) = plev
    dims3d(3) = endlat-beglat+1
    s3d=dims3d(1)*dims3d(2)*dims3d(3)
    dims2d(1) = plon
    dims2d(2) = dims3d(3)
    s2d=dims2d(1)*dims2d(2)

    allocate(tmp(max(s4d,(s3d*plevp)/plev)))

    ldof => get_restart_decomp(plon, plat, plev)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat,plev/), ldof, iodesc3d)
    deallocate(ldof)

    ldof => get_restart_decomp(plon, plat, plevp)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat,plevp/), ldof, iodesc3dp1)
    deallocate(ldof)

    ldof => get_restart_decomp(plon, plat, plev*pcnst)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat,plev,pcnst/), ldof, iodesc4d)
    deallocate(ldof)

    ldof => get_restart_decomp(plon, plat, 1)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat/), ldof, iodesc2d)
    deallocate(ldof)

    ierr = PIO_Inq_varid(File, 'tmass0', tmass0desc)
    ierr = pio_get_var(File, tmass0desc, tmass0)

    ierr = PIO_Inq_varid(File, 'hw1', hw1desc)
    ierr = pio_get_var(File, hw1desc, hw1)
    ierr = PIO_Inq_varid(File, 'hw2', hw2desc)
    ierr = pio_get_var(File, hw2desc, hw2)
    ierr = PIO_Inq_varid(File, 'hw3', hw3desc)
    ierr = pio_get_var(File, hw3desc, hw3)
    ierr = PIO_Inq_varid(File,'alpha', alphadesc)
    ierr = pio_get_var(File, alphadesc, alpha)

    ierr = PIO_Inq_varid(File,'lnpstar', lnpstardesc)
    ierr = pio_get_var(File, lnpstardesc, lnpstar)

    call init_restart_varlist()

    do i=1,restartvarcnt
       call get_restart_var(i, name, timelevels, ndims, vdesc)

       ierr = PIO_Inq_varid(File, name, vdesc)
       if(timelevels == 1) then
          if(ndims==2) then
             call pio_read_darray(File, vdesc, iodesc2d, tmp(1:s2d), ierr)
             restartvars(i)%v2d(:,:) = reshape(tmp(1:s2d), dims2d)
          else if(ndims==3) then
             if(name.eq.'ETADOT') then
                s3d=s3d*plevp/plev
                dims3d(2)=plevp
                call pio_read_darray(File, restartvars(i)%vdesc, iodesc3dp1, tmp(1:s3d), ierr)
                restartvars(i)%v3d(:,:,:) = reshape(tmp(1:s3d), dims3d)
                s3d=s3d*plev/plevp
                dims3d(2)=plev
             else
                call pio_read_darray(File, restartvars(i)%vdesc, iodesc3d, tmp(1:s3d), ierr)
                restartvars(i)%v3d(:,:,:) = reshape(tmp(1:s3d), dims3d)
             end if
          else if(ndims==4) then
             call pio_read_darray(File, restartvars(i)%vdesc, iodesc4d, tmp(1:s4d), ierr)
             restartvars(i)%v4d(:,:,:,:) = reshape(tmp(1:s4d), dims4d)
          end if

       else
          do t=1,timelevels
             if(t==1) ct=n3m1
             if(t==2) ct=n3
             call pio_setframe(vdesc, t)
             if(ndims==3) then
                call pio_read_darray(File, vdesc, iodesc2d, tmp(1:s2d), ierr)
                restartvars(i)%v3d(:,:,ct) = reshape(tmp(1:s2d), dims2d)
             else if(ndims==4) then
                call pio_read_darray(File, vdesc, iodesc3d, tmp(1:s3d), ierr)
                restartvars(i)%v4d(:,:,:,ct) = reshape(tmp(1:s3d), dims3d)
             else if(ndims==5) then
                call pio_read_darray(File, vdesc, iodesc4d, tmp(1:s4d), ierr)
                restartvars(i)%v5d(:,:,:,:,ct) = reshape(tmp(1:s4d), dims4d)
             end if

          end do
       end if
    end do
    deallocate(tmp)
    call pio_freedecomp(File, iodesc2d)
    call pio_freedecomp(File, iodesc3d)
    call pio_freedecomp(File, iodesc3dp1)
    call pio_freedecomp(File, iodesc4d)

    return
  end subroutine read_restart_dynamics


end module restart_dynamics
