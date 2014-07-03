module nctopo_util_mod
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Driver for SE's hyper-viscsoity smoothing procedure
  !          used to create smoothed PHIS, SGH, SGH30 fields
  !
  !          This utility is not used during normal CAM simulations.
  !          It will be run if the user sets smooth_phis_numcycle>0 in the
  !          atm_in namelist, and adds PHIS_SM, SGH_SM and SGH30_SM to
  !          on of the history files. 
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



  subroutine nctopo_util_inidat( ncid_topo, iodesc, elem)
    use control_mod,      only: smooth_phis_numcycle
    use parallel_mod,     only: par
    use bndry_mod,     only: bndry_exchangev
    use dof_mod, only           : putUniquePoints
    use edge_mod, only : edgevpack, edgevunpack, InitEdgeBuffer, FreeEdgeBuffer, EdgeBuffer_t
    use ncdio_atm, only : infld
    use abortutils,     only: endrun
    use pio, only : file_desc_t, io_desc_t, pio_double, pio_get_local_array_size, pio_freedecomp

    implicit none
    type(file_desc_t),intent(inout) :: ncid_topo
    type(io_desc_t),intent(inout) :: iodesc
    type(element_t), pointer :: elem(:)

    real(r8), allocatable :: tmp(:,:)
    integer :: tlncols, ig, ie, start, j, t, k
    character(len=40) :: fieldname
    logical :: found
    integer :: kptr
    type(EdgeBuffer_t) :: edge
    integer :: lsize, nets,nete

    if (smooth_phis_numcycle==0) return

    if(iam > par%nprocs) then
       ! The special case of npes_se < npes_cam is not worth dealing with here
       call endrun('PHIS topo generation code code requires npes_se==npes_cam')
    end if

    tlncols = pio_get_local_array_size(iodesc)	
    allocate(tmp(tlncols,1))

    allocate(PHISdyn(np,np,nelemd))
    allocate(SGHdyn(np,np,nelemd))
    allocate(SGH30dyn(np,np,nelemd))


    fieldname = 'PHIS'
    if(par%masterproc  ) write(iulog,*) 'nctopo utility: reading PHIS:'
    call infld(fieldname, ncid_topo, iodesc, tmp(:,1), found)
    if(.not. found) then
       call endrun('Could not find PHIS field on input datafile')
    end if
    start=1
    do ie=1,nelemd
       call putUniquePoints(elem(ie)%idxP, tmp(start:,1),PHISdyn(:,:,ie))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    fieldname = 'SGH'
    if(par%masterproc  ) write(iulog,*) 'nctopo utility: reading SGH:'
    call infld(fieldname, ncid_topo, iodesc, tmp(:,1), found)
    if(.not. found) then
       call endrun('Could not find SGH field on input datafile')
    end if
    start=1
    do ie=1,nelemd
       call putUniquePoints(elem(ie)%idxP, tmp(start:,1),SGHdyn(:,:,ie))
       start=start+elem(ie)%idxP%numUniquePts
    end do
    
    fieldname = 'SGH30'
    if(par%masterproc  ) write(iulog,*) 'nctopo utility: reading SGH30:'
    call infld(fieldname, ncid_topo, iodesc, tmp(:,1), found)
    if(.not. found) then
       call endrun('Could not find SGH30 field on input datafile')
    end if
    start=1
    do ie=1,nelemd
       call putUniquePoints(elem(ie)%idxP, tmp(start:,1),SGH30dyn(:,:,ie))
       start=start+elem(ie)%idxP%numUniquePts
    end do
    
    ! update non-unique points:
    call initEdgeBuffer(edge, 3)
    do ie=1,nelemd
       kptr=0
       call edgeVpack(edge, SGH30dyn(:,:,ie),1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVpack(edge, SGHdyn(:,:,ie),1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVpack(edge, PHISdyn(:,:,ie),1,kptr,elem(ie)%desc)
    end do
    call bndry_exchangeV(par,edge)
    do ie=1,nelemd
       kptr=0
       call edgeVunpack(edge, SGH30dyn(:,:,ie),1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVunpack(edge, SGHdyn(:,:,ie),1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVunpack(edge, PHISdyn(:,:,ie),1,kptr,elem(ie)%desc)
    end do
    call FreeEdgeBuffer(edge)
     
    
    deallocate(tmp)


  end subroutine 



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
