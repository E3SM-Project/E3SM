module nctopo_util_mod
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Driver for SE's hyper-viscsoity smoothing procedure
  !          used to create smoothed PHIS, SGH, SGH30 fields
  !
  !          This utility is not used during normal CAM simulations.
  !          It will be run if the user sets smooth_phis_numcycle>0 in the
  !          atm_in namelist, and adds PHIS_SM, SGH_SM and SGH30_SM to
  !          one of the history files. 
  !
  ! 
  ! Author:  M. Taylor (3/2011)
  ! 
  !-----------------------------------------------------------------------
  use cam_logfile, only : iulog
  use element_mod, only : element_t
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: iam
  use dimensions_mod,     only: nelemd, nlev, np, npsq
  implicit none
  private
  public nctopo_util_inidat, nctopo_util_driver


  real(r8),allocatable :: SGHdyn(:,:,:),SGH30dyn(:,:,:),PHISdyn(:,:,:)
  public sghdyn,sgh30dyn,phisdyn

contains



  subroutine nctopo_util_inidat(ncid_topo, elem)
    use control_mod,    only: smooth_phis_numcycle
    use parallel_mod,   only: par
    use bndry_mod,      only: bndry_exchangev
    use dof_mod,        only: putUniquePoints
    use edge_mod,       only: edgevpack, edgevunpack, InitEdgeBuffer, FreeEdgeBuffer
    use edgetype_mod,   only : EdgeBuffer_t
    use dyn_grid,       only: dyn_decomp
    use ncdio_atm,      only: infld
    use cam_abortutils, only: endrun
    use pio,            only: file_desc_t

    implicit none
    type(file_desc_t),intent(inout) :: ncid_topo
    type(element_t), pointer :: elem(:)

    real(r8), allocatable :: tmp(:,:)
    integer :: ie, i, j, indx
    character(len=40) :: fieldname
    logical :: found
    integer :: kptr
    type(EdgeBuffer_t) :: edge

    if (smooth_phis_numcycle==0) return

    if(iam > par%nprocs) then
       ! The special case of npes_se < npes_cam is not worth dealing with here
       call endrun('PHIS topo generation code code requires npes_se==npes_cam')
    end if

    allocate(tmp(npsq,nelemd))
    tmp = 0.0_r8

    allocate(PHISdyn(np,np,nelemd))
    PHISdyn = 0.0_r8
    allocate(SGHdyn(np,np,nelemd))
    SGHdyn = 0.0_r8
    allocate(SGH30dyn(np,np,nelemd))
    SGH30dyn = 0.0_r8


    fieldname = 'PHIS'
    if(par%masterproc  ) write(iulog,*) 'nctopo utility: reading PHIS:'
    call infld(fieldname, ncid_topo, 'ncol',                                  &
         1, npsq, 1, nelemd, tmp(:,:), found, gridname='GLL')
    if(.not. found) then
       call endrun('Could not find PHIS field on input datafile')
    end if
    do ie=1,nelemd
       indx = 1
       do j = 1, np
          do i = 1, np
             PHISdyn(i,j,ie) = tmp(indx,ie)
             indx = indx + 1
          end do
       end do
    end do

    fieldname = 'SGH'
    if(par%masterproc  ) write(iulog,*) 'nctopo utility: reading SGH:'
    call infld(fieldname, ncid_topo, 'ncol',                                  &
         1, npsq, 1, nelemd, tmp(:,:), found, gridname='GLL')
    if(.not. found) then
       call endrun('Could not find SGH field on input datafile')
    end if
    do ie=1,nelemd
       indx = 1
       do j = 1, np
          do i = 1, np
             SGHdyn(i,j,ie) = tmp(indx,ie)
             indx = indx + 1
          end do
       end do
    end do
    
    fieldname = 'SGH30'
    if(par%masterproc  ) write(iulog,*) 'nctopo utility: reading SGH30:'
    call infld(fieldname, ncid_topo, 'ncol',                                  &
         1, npsq, 1, nelemd, tmp(:,:), found, gridname='GLL')
    if(.not. found) then
       call endrun('Could not find SGH30 field on input datafile')
    end if
    do ie=1,nelemd
       indx = 1
       do j = 1, np
          do i = 1, np
             SGH30dyn(i,j,ie) = tmp(indx,ie)
             indx = indx + 1
          end do
       end do
    end do
    
    ! update non-unique points:
    call initEdgeBuffer(par, edge, elem, 3, numthreads_in=1)
    do ie=1,nelemd
       kptr=0
       call edgeVpack(edge, SGH30dyn(:,:,ie),1,kptr,ie)
       kptr=kptr+1
       call edgeVpack(edge, SGHdyn(:,:,ie),1,kptr,ie)
       kptr=kptr+1
       call edgeVpack(edge, PHISdyn(:,:,ie),1,kptr,ie)
    end do
    call bndry_exchangeV(par,edge)
    do ie=1,nelemd
       kptr=0
       call edgeVunpack(edge, SGH30dyn(:,:,ie),1,kptr,ie)
       kptr=kptr+1
       call edgeVunpack(edge, SGHdyn(:,:,ie),1,kptr,ie)
       kptr=kptr+1
       call edgeVunpack(edge, PHISdyn(:,:,ie),1,kptr,ie)
    end do
    call FreeEdgeBuffer(edge)
     
    
    deallocate(tmp)
    deallocate(PHISdyn)
    deallocate(SGHdyn)
    deallocate(SGH30dyn)

  end subroutine nctopo_util_inidat



  subroutine nctopo_util_driver(elem,hybrid,nets,nete)
    use prim_driver_mod,  only: smooth_topo_datasets
    use hybrid_mod,       only: hybrid_t
    use control_mod,      only: smooth_phis_numcycle

    type(element_t) :: elem(:)
    type(hybrid_t) :: hybrid
    integer :: nets,nete,i,j,ie
    real(r8) :: ftmp(npsq,1,1)

    if (smooth_phis_numcycle==0) return
    call smooth_topo_datasets(phisdyn,sghdyn,sgh30dyn,elem,hybrid,nets,nete)


  end subroutine 



end module nctopo_util_mod 
