module od_common
!==========================================================================
! This module contains code common to different orographic drag 
! parameterizations.
! It includes 4 parts:
! orographic gravity wave drag (Xie et al.,2020), 
! flow-blocking drag (Xie et al.,2020),
! small-scale orographic gravity wave drag (Tsiringakis et al. 2017), 
! turbulent orographic form drag (Beljaars et al.,2004).
!==========================================================================
use shr_kind_mod,  only: i8 => shr_kind_i8, r8 => shr_kind_r8
use shr_sys_mod,   only: shr_sys_flush
use ppgrid,        only: pcols, pver, begchunk, endchunk
use cam_logfile,   only: iulog
use cam_abortutils,only: endrun
use spmd_utils,    only: masterproc
use pio,           only: file_desc_t
use phys_control,  only: use_od_ls, use_od_bl, use_od_ss
use physics_buffer,only: dtype_r8, physics_buffer_desc, pbuf_get_chunk
use physics_buffer,only: pbuf_get_index, pbuf_get_field, pbuf_add_field, pbuf_set_field

implicit none
private
save

! Public interface.
public :: oro_drag_readnl
public :: oro_drag_register
public :: oro_drag_init
public :: oro_drag_interface
public :: od_gsd,pblh_get_level_idx,grid_size

type(file_desc_t), pointer :: topo_file_ncid

! dimensions for topo shape data
integer, parameter :: ndir_asymmetry = 2+1 ! add 1 to avoid bug reading file - not sure why this happens
integer, parameter :: ndir_efflength = 180 ! 1-degree resolution with opposite directions mirrored

! pbuf indices for data read in from topo data file
integer :: oro_drag_convexity_idx = -1 ! Convexity
integer :: oro_drag_asymmetry_idx = -1 ! Asymmetry
integer :: oro_drag_efflength_idx = -1 ! Effective length
integer :: oro_drag_ribulk_idx    = -1 ! bulk richardson number (calculated in CLUBB)

!tunable parameter to the od schemes
real(r8),public, protected :: od_ls_ncleff = 3._r8 !tunable parameter for oGWD
real(r8),public, protected :: od_bl_ncd    = 3._r8 !tunable parameter for FBD
real(r8),public, protected :: od_ss_sncleff= 1._r8 !tunable parameter for sGWD

contains

!==========================================================================

subroutine oro_drag_readnl(nlfile)   
  
  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! Local variables
  integer :: unitn, ierr 
  character(len=*), parameter :: subname = 'oro_drag_readnl'

  ! More specific name for dc to prevent a name clash or confusion in the
  ! namelist.

  namelist /oro_drag_nl/ od_ls_ncleff, od_bl_ncd, od_ss_sncleff
  !---------------------------------------------------------------------
  !read oro_drag_nl only when use the od schemes
  if (use_od_ls.or.use_od_bl.or.use_od_ss) then
    if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'oro_drag_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, oro_drag_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun(subname // ':: ERROR reading namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)
    end if

    if (masterproc) write(iulog,*) "oro_drag_readnl od_ls_ncleff, od_bl_ncd, od_ss_sncleff ",od_ls_ncleff,od_bl_ncd,od_ss_sncleff

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(od_ls_ncleff,    1, mpir8,  0, mpicom)
  call mpibcast(od_bl_ncd,       1, mpir8,  0, mpicom)
  call mpibcast(od_ss_sncleff,   1, mpir8,  0, mpicom)
#endif
  !
  endif

end subroutine oro_drag_readnl

!==========================================================================

subroutine oro_drag_open_topo_file()
  use filenames,    only: bnd_topo
  use ioFileMod,    only: getfil
  use cam_pio_utils,only: cam_pio_openfile
  use pio,          only: pio_nowrite
  include 'netcdf.inc'
  !-----------------------------------------------------------------------
  character(len=256) :: bnd_topo_loc   ! filepath of topo file on local disk
  allocate(topo_file_ncid)
  call getfil(bnd_topo, bnd_topo_loc)
  call cam_pio_openfile(topo_file_ncid, bnd_topo_loc, PIO_NOWRITE)
end subroutine oro_drag_open_topo_file

!==========================================================================

subroutine oro_drag_close_topo_file
  use pio,          only: pio_closefile
  call pio_closefile(topo_file_ncid)
  deallocate(topo_file_ncid)
  nullify(topo_file_ncid)
end subroutine oro_drag_close_topo_file 

!==========================================================================

subroutine oro_drag_register()
  !-----------------------------------------------------------------------
  ! Register pbuf variables for orographic drag parameterizations
  !-----------------------------------------------------------------------
  ! create pbuf variables to hold oro drag data
  if (use_od_ls.or.use_od_bl) then
    call pbuf_add_field('oro_drag_convexity','physpkg',dtype_r8,(/pcols/),               oro_drag_convexity_idx)
    call pbuf_add_field('oro_drag_asymmetry','physpkg',dtype_r8,(/pcols,ndir_asymmetry/),oro_drag_asymmetry_idx)
    call pbuf_add_field('oro_drag_efflength','physpkg',dtype_r8,(/pcols,ndir_efflength/),oro_drag_efflength_idx)
  end if
  if (use_od_ss) then
    call pbuf_add_field('oro_drag_ribulk',   'physpkg',dtype_r8,(/pcols/),               oro_drag_ribulk_idx)
  end if

end subroutine oro_drag_register

!==========================================================================

subroutine oro_drag_init(pbuf2d)
  !-----------------------------------------------------------------------
  ! Initialization for orographic drag parameterizations
  !-----------------------------------------------------------------------
  use pio,                 only: file_desc_t
  use ncdio_atm,           only: infld
  use cam_grid_support,    only: cam_grid_check, cam_grid_get_decomp, cam_grid_id,cam_grid_get_dim_names
  use infnan,              only: nan, assignment(=)
  !-----------------------------------------------------------------------
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  !-----------------------------------------------------------------------
  logical :: found
  character(len=8) :: dim1name, dim2name
  character*11 :: subname='oro_drag_init'
  integer :: grid_id
  integer :: c

  real(r8), allocatable:: oro_drag_convexity_tmp(:,:)
  real(r8), allocatable:: oro_drag_asymmetry_tmp(:,:,:)
  real(r8), allocatable:: oro_drag_efflength_tmp(:,:,:)

  type(physics_buffer_desc), pointer :: pbuf_chunk(:) ! temporary pbuf pointer for single chunk
  !-----------------------------------------------------------------------
  if (.not.(use_od_ls.or.use_od_bl)) return

  grid_id = cam_grid_id('physgrid')
  if (.not. cam_grid_check(grid_id)) then
    call endrun(trim(subname)//': Internal error, no "physgrid" grid')
  end if

  ! Alocate variables for reading oro drag data
  allocate( oro_drag_convexity_tmp(pcols,begchunk:endchunk) )
  allocate( oro_drag_asymmetry_tmp(pcols,ndir_asymmetry,begchunk:endchunk) )
  allocate( oro_drag_efflength_tmp(pcols,ndir_efflength,begchunk:endchunk) )
  oro_drag_convexity_tmp(:,:)   = nan
  oro_drag_asymmetry_tmp(:,:,:) = nan
  oro_drag_efflength_tmp(:,:,:) = nan

  ! Read special orographic shape fields from topo file
  call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
  call oro_drag_open_topo_file()

  found=.false.
  call infld( 'OC', topo_file_ncid, dim1name, dim2name, 1, pcols, &
              begchunk, endchunk, oro_drag_convexity_tmp(:,:), found, gridname='physgrid')
  if(.not. found) call endrun('ERROR - oro_drag_init: topo file read error - OC')

  found=.false.
  call infld( 'OA', topo_file_ncid, dim1name, 'ndir_asymmetry', dim2name, 1, pcols, 1, ndir_asymmetry, &
              begchunk, endchunk, oro_drag_asymmetry_tmp(:,:,:), found, gridname='physgrid')
  if(.not. found) call endrun('ERROR - oro_drag_init: topo file read error - OA')

  found=.false.
  call infld( 'OL', topo_file_ncid, dim1name, 'ndir_efflength', dim2name, 1, pcols, 1, ndir_efflength, &
              begchunk, endchunk, oro_drag_efflength_tmp(:,:,:), found, gridname='physgrid')
  if(.not. found) call endrun('ERROR - oro_drag_init: topo file read error - OL')

  call oro_drag_close_topo_file()

  ! copy the oro drag data in pbuf
  do c=begchunk,endchunk
    pbuf_chunk => pbuf_get_chunk(pbuf2d, c)
    call pbuf_set_field(pbuf_chunk, oro_drag_convexity_idx, oro_drag_convexity_tmp(:,c) )
    call pbuf_set_field(pbuf_chunk, oro_drag_asymmetry_idx, oro_drag_asymmetry_tmp(:,:,c) )
    call pbuf_set_field(pbuf_chunk, oro_drag_efflength_idx, oro_drag_efflength_tmp(:,:,c) )
  end do

  deallocate(oro_drag_convexity_tmp)
  deallocate(oro_drag_asymmetry_tmp)
  deallocate(oro_drag_efflength_tmp)

end subroutine oro_drag_init
!==========================================================================

subroutine oro_drag_interface(state,    cam_in,   sgh,      pbuf,   dtime,  nm,     &
                              gwd_ls,   gwd_bl,   gwd_ss,   gwd_fd,                 &
                              od_ls_ncleff,       od_bl_ncd,od_ss_sncleff,          &
                              utgw,     vtgw,     ttgw,                             &
                              dtaux3_ls,dtauy3_ls,dtaux3_bl,dtauy3_bl,              &
                              dtaux3_ss,dtauy3_ss,dtaux3_fd,dtauy3_fd,              &
                              dusfc_ls, dvsfc_ls ,dusfc_bl, dvsfc_bl,               &
                              dusfc_ss, dvsfc_ss ,dusfc_fd, dvsfc_fd)
  use physics_types,  only: physics_state
  use camsrfexch,     only: cam_in_t
  use ppgrid,         only: pcols,pver,pverp
  use physconst,      only: gravit,rair,cpair,rh2o,zvir,pi
  use hycoef,         only: etamid
  !-----------------------------------------------------------------------
  type(physics_state), intent(in) :: state      ! physics state structure
  type(cam_in_t),      intent(in) :: cam_in
  real(r8), intent(in) :: sgh(pcols)
  type(physics_buffer_desc), pointer :: pbuf(:) ! Physics buffer
  real(r8), intent(in) :: dtime
  real(r8), intent(in) :: nm(state%ncol,pver)   ! midpoint Brunt-Vaisalla frequency
  !options for the 4 schemes
  logical , intent(in) :: gwd_ls
  logical , intent(in) :: gwd_bl
  logical , intent(in) :: gwd_ss
  logical , intent(in) :: gwd_fd
  !tunable parameter from namelist
  real(r8), intent(in) :: od_ls_ncleff
  real(r8), intent(in) :: od_bl_ncd
  real(r8), intent(in) :: od_ss_sncleff
  !vertical profile of the momentum tendencies
  real(r8), intent(out), optional :: utgw(state%ncol,pver)
  real(r8), intent(out), optional :: vtgw(state%ncol,pver)
  real(r8), intent(out), optional :: ttgw(state%ncol,pver)
  !output drag terms in 3D and surface
  real(r8), intent(out), optional :: dtaux3_ls(pcols,pver)
  real(r8), intent(out), optional :: dtauy3_ls(pcols,pver)
  real(r8), intent(out), optional :: dtaux3_bl(pcols,pver)
  real(r8), intent(out), optional :: dtauy3_bl(pcols,pver)
  real(r8), intent(out), optional :: dtaux3_ss(pcols,pver)
  real(r8), intent(out), optional :: dtauy3_ss(pcols,pver)
  real(r8), intent(out), optional :: dtaux3_fd(pcols,pver)
  real(r8), intent(out), optional :: dtauy3_fd(pcols,pver)
  real(r8), intent(out), optional :: dusfc_ls(pcols)
  real(r8), intent(out), optional :: dvsfc_ls(pcols)
  real(r8), intent(out), optional :: dusfc_bl(pcols)
  real(r8), intent(out), optional :: dvsfc_bl(pcols)
  real(r8), intent(out), optional :: dusfc_ss(pcols)
  real(r8), intent(out), optional :: dvsfc_ss(pcols)
  real(r8), intent(out), optional :: dusfc_fd(pcols)
  real(r8), intent(out), optional :: dvsfc_fd(pcols)
 
  real(r8) :: ztop(pcols,pver)             ! top interface height asl (m)
  real(r8) :: zbot(pcols,pver)             ! bottom interface height asl (m)
  real(r8) :: zmid(pcols,pver)             ! middle interface height asl (m)
  real(r8) :: dz(pcols,pver)               ! model layer height
  
  integer  :: pblh_idx     = 0
  integer  :: kpbl2d_in(pcols)
  integer  :: kpbl2d_reverse_in(pcols)
  real(r8), pointer :: pblh(:)
  real(r8) :: dx(pcols),dy(pcols)

  real(r8), pointer :: oro_drag_convexity(:)
  real(r8), pointer :: oro_drag_asymmetry(:,:)
  real(r8), pointer :: oro_drag_efflength(:,:)
  real(r8), pointer :: oro_drag_ribulk(:)                     ! pbuf pointer for bulk richardson number
  
  integer  :: ncol
  integer  :: i
  integer  :: k
  !-----------------------------------------------------------------------

  ncol=state%ncol
  !convert heights above surface to heights above sea level
  !obtain z,dz,dx,dy,and k for pblh
  kpbl2d_in=0_r8
  kpbl2d_reverse_in=0_r8
  ztop=0._r8
  zbot=0._r8
  zmid=0._r8
  dusfc_ls=0._r8
  dvsfc_ls=0._r8
  dusfc_bl=0._r8
  dvsfc_bl=0._r8
  dusfc_ss=0._r8
  dvsfc_ss=0._r8
  dusfc_fd=0._r8
  dvsfc_fd=0._r8
  dtaux3_ls=0._r8
  dtaux3_bl=0._r8
  dtauy3_ls=0._r8
  dtauy3_bl=0._r8
  dtaux3_ss=0._r8
  dtaux3_fd=0._r8
  dtauy3_ss=0._r8
  dtauy3_fd=0._r8

  do k=1,pver
    do i=1,ncol
      ! assign values for level top/bottom
      ztop(i,k)=state%zi(i,k)
      zbot(i,k)=state%zi(i,k+1)
    enddo
  end do

  !transform adding the pressure
  !transfer from surface to sea level
  do k=1,pver
    do i=1,ncol
      ztop(i,k)=ztop(i,k)+state%phis(i)/gravit
      zbot(i,k)=zbot(i,k)+state%phis(i)/gravit
      zmid(i,k)=state%zm(i,k)+state%phis(i)/gravit
      !dz is from bottom to top already for gw_drag
      dz(i,k)=ztop(i,k)-zbot(i,k)
    end do
  end do
  !get the layer index of pblh in layer for input in drag scheme
  pblh_idx = pbuf_get_index('pblh')
  call pbuf_get_field(pbuf, pblh_idx, pblh)
  do i=1,pcols
      kpbl2d_in(i)=pblh_get_level_idx(zbot(i,:)-(state%phis(i)/gravit),pblh(i))
      kpbl2d_reverse_in(i)=pverp-kpbl2d_in(i)!pverp-k
  end do

  call pbuf_get_field(pbuf, oro_drag_convexity_idx, oro_drag_convexity )
  call pbuf_get_field(pbuf, oro_drag_asymmetry_idx, oro_drag_asymmetry )
  call pbuf_get_field(pbuf, oro_drag_efflength_idx, oro_drag_efflength )
  call pbuf_get_field(pbuf, oro_drag_ribulk_idx,    oro_drag_ribulk)
  
  !get grid size for dx,dy
  call grid_size(state,dx,dy)

  !interface for orographic drag
  call od_gsd(u3d=state%u(:ncol,pver:1:-1),v3d=state%v(:ncol,pver:1:-1),t3d=state%t(:ncol,pver:1:-1),&
              qv3d=state%q(:ncol,pver:1:-1,1),p3d=state%pmid(:ncol,pver:1:-1),p3di=state%pint(:ncol,pver+1:1:-1),&
              pi3d=state%exner(:ncol,pver:1:-1),z=zbot(:ncol,pver:1:-1),&
              od_ls_ncleff=od_ls_ncleff,od_bl_ncd=od_bl_ncd,od_ss_sncleff=od_ss_sncleff,&
              rublten=utgw(:ncol,pver:1:-1),rvblten=vtgw(:ncol,pver:1:-1),rthblten=ttgw(:ncol,pver:1:-1),&
              dtaux3d_ls=dtaux3_ls(:ncol,pver:1:-1),dtauy3d_ls=dtauy3_ls(:ncol,pver:1:-1),&
              dtaux3d_bl=dtaux3_bl(:ncol,pver:1:-1),dtauy3d_bl=dtauy3_bl(:ncol,pver:1:-1),&
              dtaux3d_ss=dtaux3_ss(:ncol,pver:1:-1),dtauy3d_ss=dtauy3_ss(:ncol,pver:1:-1),&
              dtaux3d_fd=dtaux3_fd(:ncol,pver:1:-1),dtauy3d_fd=dtauy3_fd(:ncol,pver:1:-1),&
              dusfcg_ls=dusfc_ls(:ncol),dvsfcg_ls=dvsfc_ls(:ncol),&
              dusfcg_bl=dusfc_bl(:ncol),dvsfcg_bl=dvsfc_bl(:ncol),&
              dusfcg_ss=dusfc_ss(:ncol),dvsfcg_ss=dvsfc_ss(:ncol),&
              dusfcg_fd=dusfc_fd(:ncol),dvsfcg_fd=dvsfc_fd(:ncol),&
              xland=cam_in%landfrac,br=oro_drag_ribulk(:ncol),&
              var2d=sgh(:ncol),&
              oc12d=oro_drag_convexity(:ncol),&
              oa2d=oro_drag_asymmetry(:ncol,:),&
              ol2d=oro_drag_efflength(:ncol,:),&
              znu=etamid(pver:1:-1),dz=dz(:ncol,pver:1:-1),pblh=pblh(:ncol),&
              cp=cpair,g=gravit,rd=rair,rv=rh2o,ep1=zvir,pi=pi,bnvbg=nm(:ncol,pver:1:-1),&
              dt=dtime,dx=dx,dy=dy,&
              kpbl2d=kpbl2d_reverse_in,gwd_opt=0,&
              ids=1,ide=ncol,jds=0,jde=0,kds=1,kde=pver, &
              ims=1,ime=ncol,jms=0,jme=0,kms=1,kme=pver, &
              its=1,ite=ncol,jts=0,jte=0,kts=1,kte=pver, &
              gwd_ls=gwd_ls,gwd_bl=gwd_bl,gwd_ss=gwd_ss,gwd_fd=gwd_fd )

end subroutine oro_drag_interface

!==========================================================================

function pblh_get_level_idx(height_array,pblheight)
  implicit none
  real(r8),intent(in),dimension(pver) :: height_array
  real(r8),intent(in) :: pblheight
  integer :: pblh_get_level_idx
  !local 
  integer :: k
  logical :: found

  pblh_get_level_idx = -1
  found=.false.
  !get the pblh level index and return
  do k = 1, pver
    if((pblheight >= height_array(k+1).and.pblheight <height_array(k)))then
      pblh_get_level_idx =  k+1
      found=.true.
      return
    endif
  enddo

end function

!==========================================================================

subroutine dxygrid(dx,dy,theta_in,dxy)

  IMPLICIT NONE
  real(r8),intent(in) :: dx,dy,theta_in
  real(r8),intent(out):: dxy
  !local variables
  real(r8) :: rad,theta,theta1

  rad=4.0_r8*atan(1.0_r8)/180.0_r8
  theta1=MOD(theta_in,360._r8)
  !set negative axis into 0~360
  if (theta1.ge.-360._r8.and.theta1.lt.0._r8) then
    theta1=theta1+360._r8
  endif
  !in case the angle is not into the judgement
  theta=theta1
  !transform of angle into first quadrant
  if      (theta1.ge.  0._r8.and.theta1.lt. 90._r8) then
    theta=theta1
  else if (theta1.gt. 90._r8.and.theta1.lt.180._r8) then
    theta=(180._r8-theta1)
  else if (theta1.gt.180._r8.and.theta1.lt.270._r8) then
    theta=(theta1-180._r8)
  else if (theta1.gt.270._r8.and.theta1.lt.360._r8) then
    theta=(360._r8-theta1)
  else if (theta1.eq.90._r8.or.theta1.eq.270._r8) then
    theta=90._r8
  else if (theta1.eq.0._r8.or.theta1.eq.180._r8) then
    theta=0._r8
  endif

  !get dxy
  if   (theta.ge. 0._r8.and.theta.lt.atan2(dy,dx)/rad) then
    dxy=dx/cos(theta*rad)
  else if (theta.ge.atan2(dy,dx)/rad.and.theta.le.90._r8)then
    dxy=dy/sin(theta*rad)
  endif

end subroutine dxygrid

!==========================================================================

subroutine grid_size(state, grid_dx, grid_dy)
  ! Determine the size of the grid for each of the columns in state

  use phys_grid,       only: get_area_p
  use shr_const_mod,   only: shr_const_pi
  use physics_types,   only: physics_state
  use ppgrid,          only: pver, pverp, pcols

  type(physics_state), intent(in) :: state
  real(r8), intent(out)           :: grid_dx(pcols), grid_dy(pcols)   ! E3SM grid [m]

  real(r8), parameter :: earth_ellipsoid1 = 111132.92_r8 ! World Geodetic System 1984 (WGS84) 
                                                         ! first coefficient, meters per degree longitude at equator
  real(r8), parameter :: earth_ellipsoid2 = 559.82_r8 ! second expansion coefficient for WGS84 ellipsoid
  real(r8), parameter :: earth_ellipsoid3 = 1.175_r8 ! third expansion coefficient for WGS84 ellipsoid
  real(r8) :: mpdeglat, column_area, degree, lat_in_rad
  integer  :: i
  
  do i=1,state%ncol
    ! determine the column area in radians
    column_area = get_area_p(state%lchnk,i)
    ! convert to degrees
    degree = sqrt(column_area)*(180._r8/shr_const_pi)

    ! convert latitude to radians
    lat_in_rad = state%lat(i)*(shr_const_pi/180._r8)

    ! Now find meters per degree latitude
    ! Below equation finds distance between two points on an ellipsoid, derived from expansion
    !  taking into account ellipsoid using World Geodetic System (WGS84) reference 
    mpdeglat = earth_ellipsoid1 - earth_ellipsoid2 * cos(2._r8*lat_in_rad) + earth_ellipsoid3 * cos(4._r8*lat_in_rad)
    grid_dx(i) = mpdeglat * degree
    grid_dy(i) = grid_dx(i) ! Assume these are the same
  enddo
end subroutine grid_size

!==========================================================================

subroutine od_gsd(u3d,v3d,t3d,qv3d,p3d,p3di,pi3d,z,                            &
                  od_ls_ncleff,od_bl_ncd,od_ss_sncleff,                        &
                  rublten,rvblten,rthblten,                                    &
                  dtaux3d_ls,dtauy3d_ls,dtaux3d_bl,dtauy3d_bl,                 &
                  dtaux3d_ss,dtauy3d_ss,dtaux3d_fd,dtauy3d_fd,                 &
                  dusfcg_ls,dvsfcg_ls,dusfcg_bl,dvsfcg_bl,dusfcg_ss,dvsfcg_ss, &
                  dusfcg_fd,dvsfcg_fd,xland,br,                                &
                  var2d,oc12d,oa2d,ol2d,znu,znw,p_top,dz,pblh,                 &
                  cp,g,rd,rv,ep1,pi,bnvbg,                                     &
                  dt,dx,dy,kpbl2d,gwd_opt,                                     &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                  gwd_ls,gwd_bl,gwd_ss,gwd_fd)
  !-------------------------------------------------------------------------------
   implicit none
  !-------------------------------------------------------------------------------
  !                                                                       
  !-- u3d         3d u-velocity interpolated to theta points (m/s)
  !-- v3d         3d v-velocity interpolated to theta points (m/s)
  !-- t3d         temperature (k)
  !-- qv3d        3d water vapor mixing ratio (kg/kg)
  !-- p3d         3d pressure (pa)
  !-- p3di        3d pressure (pa) at interface level
  !-- pi3d        3d exner function (dimensionless)
  !-- rublten     u tendency due to pbl parameterization (m/s/s) 
  !-- rvblten     v tendency due to pbl parameterization (m/s/s)
  !-- rthblten    theta tendency due to pbl parameterization (K/s)
  !-- znu         eta values (sigma values)
  !-- cp          heat capacity at constant pressure for dry air (j/kg/k)
  !-- g           acceleration due to gravity (m/s^2)
  !-- rd          gas constant for dry air (j/kg/k)
  !-- z           height above sea level (m)
  !-- rv          gas constant for water vapor (j/kg/k)
  !-- dt          time step (s)
  !-- dx          model grid interval (m)
  !-- dz          height of model layers (m)
  !-- xland       land mask (1 for land, 2 for water)
  !-- br          bulk richardson number in surface layer
  !-- pblh        planetary boundary layer height (m)
  !-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
  !-- ids         start index for i in domain
  !-- ide         end index for i in domain
  !-- jds         start index for j in domain
  !-- jde         end index for j in domain
  !-- kds         start index for k in domain
  !-- kde         end index for k in domain
  !-- ims         start index for i in memory
  !-- ime         end index for i in memory
  !-- jms         start index for j in memory
  !-- jme         end index for j in memory
  !-- kms         start index for k in memory
  !-- kme         end index for k in memory
  !-- its         start index for i in tile
  !-- ite         end index for i in tile
  !-- jts         start index for j in tile
  !-- jte         end index for j in tile
  !-- kts         start index for k in tile
  !-- kte         end index for k in tile
  !-------------------------------------------------------------------------------
  integer,  intent(in)   ::      ids, ide, jds, jde, kds, kde
  integer,  intent(in)   ::      ims, ime, jms, jme, kms, kme
  integer,  intent(in)   ::      its, ite, jts, jte, kts, kte
  integer,  intent(in)   ::      gwd_opt
  real(r8), intent(in)   ::      cp,g,rd,rv,ep1,pi
  !input model grid length for the grid dx,dy
  real(r8), dimension(:),                   intent(in)  :: dx
  real(r8), dimension(:),                   intent(in)  :: dy
  !input atmospheric variables
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  :: qv3d
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  ::  p3d
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  :: pi3d
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  ::  t3d
  real(r8), dimension( ims:ime, kms:kme+1), intent(in)  :: p3di
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  ::  u3d
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  ::  v3d
  !input tunable parameters
  real(r8),                                 intent(in)  :: od_ls_ncleff
  real(r8),                                 intent(in)  :: od_bl_ncd
  real(r8),                                 intent(in)  :: od_ss_sncleff
  !logical variables for selection of 4 schemes
  logical,                                  intent(in)  :: gwd_ls
  logical,                                  intent(in)  :: gwd_bl
  logical,                                  intent(in)  :: gwd_ss
  logical,                                  intent(in)  :: gwd_fd
  !input variables
  integer,  dimension( ims:ime ),           intent(in)  :: kpbl2d
  real(r8), dimension( ims:ime ),           intent(in)  :: pblh, br, xland
  real(r8), dimension( ims:ime,kms:kme ),   intent(in)  ::          z
  real(r8), dimension( ims:ime,kms:kme ),   intent(in)  ::         dz
  !input topographic parameters
  real(r8), dimension( ims:ime ),           intent(in),   optional  :: var2d
  real(r8), dimension( ims:ime ),           intent(in),   optional  :: oc12d
  real(r8), dimension( ims:ime,ndir_efflength ),intent(in),   optional  :: ol2d
  real(r8), dimension( ims:ime,ndir_asymmetry ),intent(in),   optional  :: oa2d
  !input model parameters
  real(r8),                                 intent(in),   optional  :: dt
  real(r8),                                 intent(in),   optional  :: p_top
  real(r8), dimension(:),                   intent(in),   optional  :: znu
  real(r8), dimension(:),                   intent(in),   optional  :: znw
  !input bnv from gw_drag
  real(r8), dimension( ims:ime,kms:kme ),   intent(in),   optional  ::  bnvbg
  !output of vertical profile for momentum terms
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   rublten
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   rvblten
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   rthblten
  !output for 4 drag terms
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtaux3d_ls
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtauy3d_ls
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtaux3d_bl
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtauy3d_bl
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtaux3d_ss
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtauy3d_ss
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtaux3d_fd
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtauy3d_fd
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dusfcg_ls
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dvsfcg_ls
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dusfcg_bl
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dvsfcg_bl    
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dusfcg_ss
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dvsfcg_ss
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dusfcg_fd
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dvsfcg_fd
  !local drag terms
  real(r8), dimension( ims:ime, kms:kme )  ::  dtaux2d_ls
  real(r8), dimension( ims:ime, kms:kme )  ::  dtauy2d_ls
  real(r8), dimension( ims:ime, kms:kme )  ::  dtaux2d_bl
  real(r8), dimension( ims:ime, kms:kme )  ::  dtauy2d_bl
  real(r8), dimension( ims:ime, kms:kme )  ::  dtaux2d_ss
  real(r8), dimension( ims:ime, kms:kme )  ::  dtauy2d_ss
  real(r8), dimension( ims:ime, kms:kme )  ::  dtaux2d_fd
  real(r8), dimension( ims:ime, kms:kme )  ::  dtauy2d_fd
  real(r8), dimension( ims:ime ) ::  dusfc_ls
  real(r8), dimension( ims:ime ) ::  dvsfc_ls
  real(r8), dimension( ims:ime ) ::  dusfc_bl
  real(r8), dimension( ims:ime ) ::  dvsfc_bl
  real(r8), dimension( ims:ime ) ::  dusfc_ss
  real(r8), dimension( ims:ime ) ::  dvsfc_ss
  real(r8), dimension( ims:ime ) ::  dusfc_fd
  real(r8), dimension( ims:ime ) ::  dvsfc_fd
  !local variables
  real(r8),   dimension( its:ite, kts:kte )     ::  delprsi
  real(r8),   dimension( its:ite, kts:kte )     ::  pdh
  real(r8),   dimension( its:ite, kts:kte+1 )   ::  pdhi
  real(r8),   dimension( its:ite, ndir_asymmetry )  ::  oa4
  real(r8),   dimension( its:ite, ndir_efflength )  ::  ol4
  integer ::  i,j,k,kpblmax
  !determine the lowest level for planet boundary layer
  do k = kts,kte
    if( znu(k).gt.0.6_r8 ) kpblmax = k + 1
  enddo
  !get the interface pressure (hPa) for each level
  !for the last level over kte use mid-layer pressure instead
  do k = kts,kte+1
    do i = its,ite
      if(k.le.kte) pdh(i,k) = p3d(i,k)
        pdhi(i,k) = p3di(i,k)
    enddo
  enddo
  !get the layer pressure for each level
  do k = kts,kte
    do i = its,ite
      delprsi(i,k) = pdhi(i,k)-pdhi(i,k+1)
    enddo
  enddo
  !no need when there is no large drag
  if ( gwd_ls .or. gwd_bl ) then
    do i = its,ite
      oa4(i,:) = oa2d(i,:)
      ol4(i,:) = ol2d(i,:)
    enddo
  endif
  !call the od2d for calculatino of each grid
  call od2d(dudt=rublten(ims,kms),dvdt=rvblten(ims,kms)                   &
             ,dthdt=rthblten(ims,kms)                                       &
             ,ncleff=od_ls_ncleff,ncd=od_bl_ncd,sncleff=od_ss_sncleff                &
             ,dtaux2d_ls=dtaux2d_ls,dtauy2d_ls=dtauy2d_ls                   &
             ,dtaux2d_bl=dtaux2d_bl,dtauy2d_bl=dtauy2d_bl                   &
             ,dtaux2d_ss=dtaux2d_ss,dtauy2d_ss=dtauy2d_ss                   &
             ,dtaux2d_fd=dtaux2d_fd,dtauy2d_fd=dtauy2d_fd                   &
             ,u1=u3d(ims,kms),v1=v3d(ims,kms)                               &
             ,t1=t3d(ims,kms)                                               &
             ,q1=qv3d(ims,kms)                                              &
             ,del=delprsi(its,kts)                                          &
             ,prsi=pdhi(its,kts)                                            &
             ,prsl=pdh(its,kts),prslk=pi3d(ims,kms)                         &
             ,zl=z(ims,kms),rcl=1.0_r8                                      &
             ,xland1=xland(ims),br1=br(ims),hpbl=pblh(ims)                  &
             ,bnv_in=bnvbg(ims,kms)                                         &
             ,dz2=dz(ims,kms)                                               &
             ,kpblmax=kpblmax                                               &
             ,dusfc_ls=dusfc_ls,dvsfc_ls=dvsfc_ls                           &
             ,dusfc_bl=dusfc_bl,dvsfc_bl=dvsfc_bl                           &
             ,dusfc_ss=dusfc_ss,dvsfc_ss=dvsfc_ss                           &
             ,dusfc_fd=dusfc_fd,dvsfc_fd=dvsfc_fd                           &
             ,var=var2d(ims),oc1=oc12d(ims)                                 &
             ,oa4=oa4,ol4=ol4                                               &
             ,g=g,cp=cp,rd=rd,rv=rv,fv=ep1,pi=pi                            &
             ,dxmeter=dx,dymeter=dy,deltim=dt                               &
             ,kpbl=kpbl2d(ims)                                              &
             ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde             &
             ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme             &
             ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte             &
             ,gsd_gwd_ls=gwd_ls,gsd_gwd_bl=gwd_bl,gsd_gwd_ss=gwd_ss,gsd_gwd_fd=gwd_fd)
    !set the total stress output to each terms for the 4 drag schemes
    do i = its,ite
      dusfcg_ls(i) = dusfc_ls(i)
      dvsfcg_ls(i) = dvsfc_ls(i)
      dusfcg_bl(i) = dusfc_bl(i)
      dvsfcg_bl(i) = dvsfc_bl(i)
      dusfcg_ss(i) = dusfc_ss(i)
      dvsfcg_ss(i) = dvsfc_ss(i)
      dusfcg_fd(i) = dusfc_fd(i)
      dvsfcg_fd(i) = dvsfc_fd(i)
    enddo
    !set the 3D output tendencies to each terms for the 4 drag schemes
    dtaux3d_ls = dtaux2d_ls
    dtaux3d_bl = dtaux2d_bl
    dtauy3d_ls = dtauy2d_ls
    dtauy3d_bl = dtauy2d_bl
    dtaux3d_ss = dtaux2d_ss
    dtaux3d_fd = dtaux2d_fd
    dtauy3d_ss = dtauy2d_ss
    dtauy3d_fd = dtauy2d_fd

end subroutine od_gsd
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine od2d(dudt,dvdt,dthdt,ncleff,ncd,sncleff,                        &
                dtaux2d_ls,dtauy2d_ls,                                     &
                dtaux2d_bl,dtauy2d_bl,                                     &
                dtaux2d_ss,dtauy2d_ss,                                     &
                dtaux2d_fd,dtauy2d_fd,                                     &
                u1,v1,t1,q1,                                               &
                del,                                                       &
                prsi,prsl,prslk,zl,rcl,                                    &
                xland1,br1,hpbl,bnv_in,dz2,                                &
                kpblmax,dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,               &
                dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd,var,oc1,oa4,ol4,       &
                g,cp,rd,rv,fv,pi,dxmeter,dymeter,deltim,kpbl,              &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte,                                 &
                gsd_gwd_ls,gsd_gwd_bl,gsd_gwd_ss,gsd_gwd_fd)
  !  This code handles the time tendencies of u v due to the effect of mountain 
  !  induced gravity wave drag from sub-grid scale orography. It includes 4 parts:
  !  orographic gravity wave drag and flow-blocking drag (Xie et al.,2020),small-scale  
  !  orographic gravity wave drag (Tsiringakis et al. 2017), and turbulent orographic
  !  form drag (Beljaars et al.,2004).
  !
  !           Activation of each component is done by specifying the integer-parameters
  !           (defined below) to .true. (active) or .false. (inactive)
  !                    gsd_gwd_ls : large-scale
  !                    gsd_gwd_bl : blocking drag 
  !                    gsd_gwd_ss : small-scale gravity wave drag
  !                    gsd_gwd_fd : topographic form drag
  !
  !
  !        References:
  !        Xie et al. (2020), JAMES
  !        Tsiringakis et al. (2017), QJRMS
  !        Beljaars et al., (2004), QJRMS
  !-------------------------------------------------------------------------------
  !
  !  input                                                                
  !        dudt (ims:ime,kms:kme)  non-lin tendency for u wind component
  !        dvdt (ims:ime,kms:kme)  non-lin tendency for v wind component
  !        u1(ims:ime,kms:kme) zonal wind / sqrt(rcl)  m/sec  at t0-dt
  !        
  !        v1(ims:ime,kms:kme) meridional wind / sqrt(rcl) m/sec at t0-dt
  !        t1(ims:ime,kms:kme) temperature deg k at t0-dt
  !        q1(ims:ime,kms:kme) specific humidity at t0-dt
  !
  !        rcl     a scaling factor = reciprocal of square of cos(lat)
  !                for gmp.  rcl=1 if u1 and v1 are wind components.
  !        deltim  time step    secs                                       
  !        del(kts:kte)  positive increment of pressure across layer (pa)
  !                                                                       
  !  output
  !        dudt, dvdt    wind tendency due to gwdo
  !
  !-------------------------------------------------------------------------------
  implicit none
  !-------------------------------------------------------------------------------
  !intput parameters
  integer,  intent(in)         ::  ids,ide, jds,jde, kds,kde
  integer,  intent(in)         ::  ims,ime, jms,jme, kms,kme
  integer,  intent(in)         ::  its,ite, jts,jte, kts,kte
  !tunable parameter in oro_drag_nl, od_ls_ncleff,od_bl_ncd,od_ss_sncleff
  real(r8), intent(in) :: ncleff
  real(r8), intent(in) :: ncd
  real(r8), intent(in) :: sncleff
  !model timestep and other parameters
  real(r8), intent(in) ::  g,rd,rv,fv,cp,pi,deltim,rcl
  !input model grid length
  real(r8), dimension(:), intent(in)   ::  dxmeter
  real(r8), dimension(:), intent(in)   ::  dymeter
  !input topo variables 
  real(r8), dimension( ims:ime,ndir_asymmetry ), intent(in) ::  oa4
  real(r8), dimension( ims:ime,ndir_efflength ), intent(in) ::  ol4
  real(r8), dimension( ims:ime )           , intent(in) ::  var
  real(r8), dimension( ims:ime )           , intent(in) ::  oc1
  !input atmospheric variables
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  u1
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  v1
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  t1
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  q1
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  zl
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  prsl
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  prslk
  real(r8), dimension( ims:ime,kms:kme+1),   intent(in) ::  prsi
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  del
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  bnv_in
  !input pbl level index
  integer,dimension( ims:ime ),              intent(in) ::  kpbl
  integer,                                   intent(in) ::  kpblmax
  !input for small-scale ogwd
  real(r8), dimension( ims:ime ),            intent(in) :: br1
  real(r8), dimension( ims:ime ),            intent(in) :: hpbl
  real(r8), dimension( ims:ime ),            intent(in) :: xland1
  real(r8), dimension( ims:ime, kms:kme ),   intent(in) :: dz2
  !variables for open/close process
  logical,                                   intent(in) ::  gsd_gwd_ls
  logical,                                   intent(in) ::  gsd_gwd_bl
  logical,                                   intent(in) ::  gsd_gwd_ss
  logical,                                   intent(in) ::  gsd_gwd_fd
  !input/output tendencies
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dudt
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dvdt
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dthdt
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtaux2d_ls
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtauy2d_ls
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtaux2d_bl
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtauy2d_bl
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtaux2d_ss
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtauy2d_ss
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtaux2d_fd
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtauy2d_fd
  !input/output total column stress
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dusfc_ls
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dvsfc_ls              
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dusfc_bl
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dvsfc_bl              
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dusfc_ss
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dvsfc_ss              
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dusfc_fd
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dvsfc_fd
  !
  ! added for small-scale orographic wave drag
  !
  integer                                  :: kpbl2,kvar
  real(r8), dimension(its:ite,kts:kte)     :: utendwave,vtendwave,thx,thvx,za
  real(r8), dimension(its:ite)             :: govrth
  real(r8), dimension(its:ite,kts:kte+1)   :: zq
  real(r8)                                 :: tauwavex0,tauwavey0,bnrf,density,tvcon,hpbl2
  real(r8), parameter                      :: varmax = 200._r8
  ! Variables for scale-awareness:
  ! Small-scale GWD + turbulent form drag
  real(r8), parameter   :: dxmin_ss = 1000._r8, dxmax_ss = 12000._r8  ! min,max range of tapering (m)
  ! Large-scale GWD
  real(r8), parameter   :: dxmin_ls = 3000._r8, dxmax_ls = 13000._r8  ! min,max range of tapering (m)
  !Add y axis for taper consider
  real(r8), parameter   :: dymin_ls = 3000._r8, dymax_ls = 13000._r8  ! min,maxrange of tapering (m)
  real(r8), parameter   :: dymin_ss = 3000._r8, dymax_ss = 13000._r8  ! min,maxrange of tapering (m)
  real(r8)              :: ss_taper, ls_taper  ! small- and large-scale tapering factors (-)
  !
  ! added Beljaars orographic form drag
  real(r8), dimension( its:ite,kts:kte )   :: utendform
  real(r8), dimension( its:ite,kts:kte )   :: vtendform
  real(r8)                                 :: a1,a2,wsp
  ! critical richardson number for wave breaking : ! larger drag with larger value
  real(r8),parameter       ::  ric     = 1._r8!original 0.25, but 1. seems better at drag profile
  real(r8),parameter       ::  ric_rig  = 0.25_r8
  real(r8),parameter       ::  dw2min  = 1._r8
  real(r8),parameter       ::  rimin   = -100._r8
  real(r8),parameter       ::  bnv2min = 1.0e-5_r8
  real(r8),parameter       ::  efmin   = 0.0_r8
  real(r8),parameter       ::  efmax   = 10.0_r8
  real(r8),parameter       ::  xl      = 4.0e4_r8
  real(r8),parameter       ::  critac  = 1.0e-5_r8
  real(r8),parameter       ::  gmax    = 1._r8
  real(r8),parameter       ::  veleps  = 1.0_r8
  real(r8),parameter       ::  factop  = 0.5_r8
  real(r8),parameter       ::  frc     = 1.0_r8
  real(r8),parameter       ::  ce      = 0.8_r8
  real(r8),parameter       ::  cg      = 0.5_r8
  real(r8),parameter       ::  tndmax  = 400._r8 / 86400._r8 ! convert 400 m/s/day to m/s/s
  integer,parameter        ::  kpblmin = 2
  !number of direction for ogwd
  integer,parameter        ::  mdir=2*ndir_efflength
  !  variables for flow-blocking drag
  real(r8),parameter       ::  frmax  = 10._r8
  real(r8),parameter       ::  olmin  = 1.0e-5_r8
  real(r8),parameter       ::  odmin  = 0.1_r8
  real(r8),parameter       ::  odmax  = 10._r8
  real(r8),parameter       ::  erad   = 6371.315e+3_r8

  !local variables
  integer                  ::  i,j,k,lcap,lcapp1,nwd,idir
  integer                  ::  klcap,kp1,ikount,kk,nwd1!added nwd1
  real(r8)                 ::  rcs,rclcs,csg,fdir,cleff,cs,rcsks
  real(r8)                 ::  wdir,ti,rdz,temp,tem2,dw2,shr2,bvf2,rdelks
  real(r8)                 ::  wtkbj,tem,gfobnv,hd,fro,rim,temc,tem1,efact
  real(r8)                 ::  temv,dtaux,dtauy,eng0,eng1,theta,rad,wdir1
  !logical variables
  logical,dimension( its:ite )            ::  ldrag
  logical,dimension( its:ite )            ::  icrilv
  logical,dimension( its:ite )            ::  flag
  logical,dimension( its:ite )            ::  kloop1
  !parameters for taking over the input
  real(r8),dimension( its:ite )           :: taub
  real(r8),dimension( its:ite )           :: xn
  real(r8),dimension( its:ite )           :: yn
  real(r8),dimension( its:ite )           :: ubar
  real(r8),dimension( its:ite )           :: vbar
  real(r8),dimension( its:ite )           :: fr
  real(r8),dimension( its:ite )           :: ulow
  real(r8),dimension( its:ite )           :: rulow
  real(r8),dimension( its:ite )           :: bnv
  real(r8),dimension( its:ite )           :: oa1
  real(r8),dimension( its:ite )           :: ol
  real(r8),dimension( its:ite )           :: roll
  real(r8),dimension( its:ite )           :: dtfac
  real(r8),dimension( its:ite )           :: brvf
  real(r8),dimension( its:ite )           :: xlinv
  real(r8),dimension( its:ite )           :: delks
  real(r8),dimension( its:ite )           :: delks1
  real(r8),dimension( its:ite,kts:kte )   :: bnv2
  real(r8),dimension( its:ite,kts:kte )   :: usqj
  real(r8),dimension( its:ite,kts:kte )   :: taud_ls
  real(r8),dimension( its:ite,kts:kte )   :: taud_bl
  real(r8),dimension( its:ite,kts:kte )   :: ro
  real(r8),dimension( its:ite,kts:kte )   :: vtk
  real(r8),dimension( its:ite,kts:kte )   :: vtj
  real(r8),dimension( its:ite )           :: zlowtop
  real(r8),dimension( its:ite )           :: coefm
  real(r8),dimension( its:ite,kts:kte+1 ) :: taup
  real(r8),dimension( its:ite,kts:kte-1 ) :: velco
  !pbl related variables
  integer ,dimension( its:ite )           :: kbl
  integer ,dimension( its:ite )           :: klowtop
  integer ,dimension( its:ite )           :: komax
  !flow-blocking related variables
  integer                                 :: kblk
  real(r8)                                :: cd
  real(r8)                                :: zblk,tautem
  real(r8)                                :: pe,ke
  !grid related variables
  real(r8),dimension( its:ite )           :: delx
  real(r8),dimension( its:ite )           :: dely
  real(r8),dimension( its:ite )           :: dxy
  real(r8),dimension( its:ite )           :: dxyp
  real(r8),dimension( its:ite,ndir_efflength ):: dxy4
  real(r8),dimension( its:ite,ndir_efflength ):: dxy4p
  !topo parameters
  real(r8),dimension( its:ite )           :: olp
  real(r8),dimension( its:ite )           :: od
  real(r8),dimension( its:ite,kts:kte+1 ) :: taufb
  !readdata for low-level determination of ogwd
  real(r8), dimension( its:ite )          :: zl_hint
  real(r8)                                :: l1,l2,S
  logical                                 :: iint
  !open/close low-level momentum adjustment according to scorer parameter
  logical  :: scorer_on=.false.
  !
  !---- constants                                                         
  !                                                                       
  rcs    = sqrt(rcl)
  cs     = 1._r8 / sqrt(rcl)
  csg    = cs * g              
  lcap   = kte
  lcapp1 = lcap + 1
  fdir   = mdir / (2.0_r8*pi)
  !
  !--- calculate scale-aware tapering factors, currently set as no taper
  !
  ls_taper=1._r8
  ss_taper=1._r8
  !
  !--- calculate length of grid for flow-blocking drag
  !
  delx=dxmeter
  dely=dymeter  
  !
  !
  !-----initialize arrays
  dtaux = 0.0_r8
  dtauy = 0.0_r8
  do i = its,ite
    klowtop(i)    = 0
    kbl(i)        = 0
  enddo

  do i = its,ite
    xn(i)         = 0.0_r8
    yn(i)         = 0.0_r8
    ubar (i)      = 0.0_r8
    vbar (i)      = 0.0_r8
    roll (i)      = 0.0_r8
    taub (i)      = 0.0_r8
    oa1(i)        = 0.0_r8
    ol(i)         = 0.0_r8
    fr(i)         = 0.0_r8
    ulow (i)      = 0.0_r8
    dtfac(i)      = 1.0_r8
    ldrag(i)      = .false.
    icrilv(i)     = .false.
    flag(i)       = .true.
    zl_hint(i)    = 0.0_r8
  enddo
   
  do k = kts,kte
    do i = its,ite
      usqj(i,k) = 0.0_r8
      bnv2(i,k) = 0.0_r8
      vtj(i,k)  = 0.0_r8
      vtk(i,k)  = 0.0_r8
      taup(i,k) = 0.0_r8
      taud_ls(i,k) = 0.0_r8
      taud_bl(i,k) = 0.0_r8
      dtaux2d_ls(i,k)= 0.0_r8
      dtauy2d_ls(i,k)= 0.0_r8
      dtaux2d_bl(i,k)= 0.0_r8
      dtauy2d_bl(i,k)= 0.0_r8
      dtaux2d_ss(i,k)= 0.0_r8
      dtauy2d_ss(i,k)= 0.0_r8
      dtaux2d_fd(i,k)= 0.0_r8
      dtauy2d_fd(i,k)= 0.0_r8
    enddo
  enddo

  do i = its,ite
    dusfc_ls(i) = 0.0_r8
    dvsfc_ls(i) = 0.0_r8
    dusfc_bl(i) = 0.0_r8
    dvsfc_bl(i) = 0.0_r8
    dusfc_ss(i) = 0.0_r8
    dvsfc_ss(i) = 0.0_r8
    dusfc_fd(i) = 0.0_r8
    dvsfc_fd(i) = 0.0_r8
  enddo

  do i = its,ite
    taup(i,kte+1) = 0.0_r8
    xlinv(i)     = 1.0_r8/xl
  enddo
  !
  !  initialize array for flow-blocking drag
  !
  taufb(its:ite,kts:kte+1) = 0.0_r8
  komax(its:ite) = 0

  do k = kts,kte
    do i = its,ite
      vtj(i,k)  = t1(i,k)  * (1._r8+fv*q1(i,k))
      vtk(i,k)  = vtj(i,k) / prslk(i,k)
      ro(i,k)   = 1._r8/rd * prsl(i,k) / vtj(i,k) ! density kg/m**3
    enddo
  enddo
  !
  !  determine reference level: maximum of 2*var and pbl heights
  !
  do i = its,ite
    zlowtop(i) = 2._r8 * var(i)
  enddo

  do i = its,ite
    kloop1(i) = .true.
  enddo

  do k = kts+1,kte
    do i = its,ite
      if(kloop1(i).and.zl(i,k)-zl(i,1).ge.zlowtop(i)) then
        klowtop(i) = k+1
        kloop1(i)  = .false.
      endif
    enddo
  enddo

  do i = its,ite
    kbl(i)   = max(kpbl(i), klowtop(i))
    kbl(i)   = max(min(kbl(i),kpblmax),kpblmin)
  enddo
  !
  !  determine the level of maximum orographic height
  !
  komax(:) = kbl(:)
  ! 
  do i = its,ite
    delks(i)  = 1.0_r8 / (prsi(i,1) - prsi(i,kbl(i)))
    delks1(i) = 1.0_r8 / (prsl(i,1) - prsl(i,kbl(i)))
  enddo
  !
  !  compute low level averages within pbl
  !
  do k = kts,kpblmax
    do i = its,ite
      if (k.lt.kbl(i)) then
        rcsks   = rcs     * del(i,k) * delks(i)
        rdelks  = del(i,k)  * delks(i)
        ubar(i) = ubar(i) + rcsks  * u1(i,k)      ! pbl u  mean
        vbar(i) = vbar(i) + rcsks  * v1(i,k)      ! pbl v  mean
        roll(i) = roll(i) + rdelks * ro(i,k)      ! ro mean
      endif
    enddo
  enddo
  !
  ! For ls and bl only
  IF (gsd_gwd_ls.or.gsd_gwd_bl) then
    !     figure out low-level horizontal wind direction 
    !     order into a counterclockwise index instead
    !
    do i = its,ite
      wdir  = atan2(vbar(i),ubar(i)) + pi!changed into y/x
      wdir1 = wdir-pi
      if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
        nwd  = MOD(nint(fdir*wdir1),mdir) + 1
      else!(-pi,0)
        nwd  = MOD(nint(fdir*(wdir1+2._r8*pi)),mdir) + 1
      endif
      !turn backwords because start is pi
      !need turning
      rad    = 4.0_r8*atan(1.0_r8)/180.0_r8
      theta  = (real(nwd,kind=r8)-1._r8)*(360._r8/real(mdir,kind=r8))
      !select OA
      oa1(i) = oa4(i,1)*cos(theta*rad)+oa4(i,2)*sin(theta*rad)
      !select OL
      ol(i)  = ol4(i,MOD(nwd-1,int(mdir/2))+1)
      !calculate dxygrid, not so slow
      call dxygrid(dxmeter(i),dymeter(i),theta,dxy(i))
      !
      !----- compute orographic width along (ol) and perpendicular (olp)
      !----- the direction of wind
      !put wdir inside the (0,2*pi) section
      !changing pi/2 either way is perpendicular
      !wdir1=wdir-pi
      if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
        nwd1  = MOD(nint(fdir*(wdir1+pi/2._r8)),mdir) + 1
        olp(i)=ol4(i,MOD(nwd1-1,int(mdir/2))+1)
      else!(-pi,0)
        nwd1  = MOD(nint(fdir*(wdir1-pi/2._r8+2._r8*pi)),mdir) + 1
        olp(i)=ol4(i,MOD(nwd1-1,int(mdir/2))+1)
      endif
      theta=(real(nwd1,kind=r8)-1._r8)*(360._r8/real(mdir,kind=r8))
      call dxygrid(dxmeter(i),dymeter(i),theta,dxyp(i))
      !
      !
      !----- compute orographic direction (horizontal orographic aspect ratio)
      !
      od(i) = olp(i)/max(ol(i),olmin)
      od(i) = min(od(i),odmax)
      od(i) = max(od(i),odmin)
      !
      !----- compute length of grid in the along(dxy) and cross(dxyp) wind directions
      !
    enddo
  ENDIF
  !============================================
  ! END INITIALIZATION; BEGIN GWD CALCULATIONS:
  !============================================
  IF (gsd_gwd_ls.or.gsd_gwd_bl.and.(ls_taper .GT. 1.E-02) ) THEN 
    !
    !---  saving richardson number in usqj for migwdi                       
    !
    do k = kts,kte-1
      do i = its,ite
        ti        = 2.0_r8 / (t1(i,k)+t1(i,k+1))
        rdz       = 1._r8/(zl(i,k+1) - zl(i,k))
        tem1      = u1(i,k) - u1(i,k+1)
        tem2      = v1(i,k) - v1(i,k+1)
        dw2       = rcl*(tem1*tem1 + tem2*tem2)
        shr2      = max(dw2,dw2min) * rdz * rdz
        bvf2      = g*(g/cp+rdz*(vtj(i,k+1)-vtj(i,k))) * ti
        usqj(i,k) = max(bvf2/shr2,rimin)
        bnv2(i,k) = max(bnv_in(i,k)**2,bnv2min )
      enddo
    enddo
    !
    !----compute the "low level" or 1/3 wind magnitude (m/s)                
    !                                                                       
    do i = its,ite
      ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0_r8)
      rulow(i) = 1._r8/ulow(i)
    enddo

    do k = kts,kte-1
      do i = its,ite
        velco(i,k)  = (0.5_r8*rcs) * ((u1(i,k)+u1(i,k+1)) * ubar(i)                &
                                    + (v1(i,k)+v1(i,k+1)) * vbar(i))
        velco(i,k)  = velco(i,k) * rulow(i)
        if ((velco(i,k).lt.veleps) .and. (velco(i,k).gt.0._r8)) then
          velco(i,k) = veleps
        endif
      enddo
    enddo
    !                                                                       
    !  no drag when critical level in the base layer                        
    !                                                                       
    do i = its,ite
      ldrag(i) = velco(i,1).le.0._r8
    enddo
    !
    !  no drag when velco.lt.0                                               
    ! 
    do k = kpblmin,kpblmax
      do i = its,ite
        if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. velco(i,k).le.0._r8
      enddo
    enddo
    !                                                                       
    !  no drag when bnv2.lt.0                                               
    !                                                                       
    do k = kts,kpblmax
      do i = its,ite
        if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. bnv2(i,k).lt.0._r8
      enddo
    enddo
    !                                                                       
    !-----the low level weighted average ri is stored in usqj(1,1; im)      
    !-----the low level weighted average n**2 is stored in bnv2(1,1; im)    
    !---- this is called bnvl2 in phys_gwd_alpert_sub not bnv2                      
    !---- rdelks (del(k)/delks) vert ave factor so we can * instead of /    
    !                                                                       
    do i = its,ite
      wtkbj     = (prsl(i,1)-prsl(i,2)) * delks1(i)
      bnv2(i,1) = wtkbj * bnv2(i,1)
      usqj(i,1) = wtkbj * usqj(i,1)
    enddo
    
    do k = kpblmin,kpblmax
      do i = its,ite
        if (k .lt. kbl(i)) then
          rdelks    = (prsl(i,k)-prsl(i,k+1)) * delks1(i)
          bnv2(i,1) = bnv2(i,1) + bnv2(i,k) * rdelks
          usqj(i,1) = usqj(i,1) + usqj(i,k) * rdelks
        endif
      enddo
    enddo
                                                                         
    do i = its,ite
      ldrag(i) = ldrag(i) .or. bnv2(i,1).le.0.0_r8
      ldrag(i) = ldrag(i) .or. ulow(i)  .eq.1.0_r8
      ldrag(i) = ldrag(i) .or. var(i)   .le.0.0_r8
    enddo
    !                                                                       
    !  set all ri low level values to the low level value          
    !                                                                       
    do k = kpblmin,kpblmax
      do i = its,ite
        if (k .lt. kbl(i)) usqj(i,k) = usqj(i,1)
      enddo
    enddo

    do i = its,ite
      if (.not.ldrag(i))   then
        bnv(i) = sqrt( bnv2(i,1) )
        fr(i) = bnv(i)  * rulow(i) * 2._r8 * var(i) * od(i)
        fr(i) = min(fr(i),frmax)
        xn(i)  = ubar(i) * rulow(i)
        yn(i)  = vbar(i) * rulow(i)
      endif
    enddo
    !
    !  compute the base level stress and store it in taub
    !  calculate enhancement factor, number of mountains & aspect        
    !  ratio const. use simplified relationship between standard            
    !  deviation & critical hgt                                          
    !
    do i = its,ite
      if (.not. ldrag(i))   then
        !maintain (oa+2) greater than or equal to 0
        efact    = max(oa1(i)+2._r8,0._r8) ** (ce*fr(i)/frc)
        efact    = min(max(efact,efmin),efmax)
        ! cleff (effective grid length) is highly tunable parameter
        ! the bigger (smaller) value produce weaker (stronger) wave drag
        cleff    = sqrt(dxy(i)**2._r8 + dxyp(i)**2._r8)
        !tune the times of drag
        cleff    = (3._r8/ncleff) * max(dxmax_ls,cleff)
        coefm(i) = (1._r8 + ol(i)) ** (oa1(i)+1._r8)
        xlinv(i) = coefm(i) / cleff
        tem      = fr(i) * fr(i) * 1.!oc1(i)
        gfobnv   = gmax * tem / ((tem + cg)*bnv(i))
        !
        if (gsd_gwd_ls) then
           taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)                       &
                    * ulow(i) * gfobnv * efact
        else     ! We've gotten what we need for the blocking scheme
           taub(i) = 0.0_r8
        end if
      else
        taub(i) = 0.0_r8
        xn(i)   = 0.0_r8
        yn(i)   = 0.0_r8
      endif
    enddo

  ENDIF   ! (gsd_gwd_ls .eq. .true.).or.(gsd_gwd_bl .eq..true.)
  !=========================================================
  ! add small-scale wavedrag for stable boundary layer
  !=========================================================
  bnrf=0._r8
  tauwavex0=0._r8
  tauwavey0=0._r8
  utendwave=0._r8
  vtendwave=0._r8
  zq=0._r8

  IF (gsd_gwd_ss.and.(ss_taper.GT.1.E-02)) THEN
    !
    ! declaring potential temperature
    !
    do k = kts,kte
      do i = its,ite
        thx(i,k) = t1(i,k)/prslk(i,k)
      enddo
    enddo

    do k = kts,kte
      do i = its,ite
        tvcon = (1._r8+fv*q1(i,k))
        thvx(i,k) = thx(i,k)*tvcon
      enddo
    enddo
    !
    ! Defining layer height
    !
    do k = kts,kte
      do i = its,ite
        zq(i,k+1) = dz2(i,k)+zq(i,k)
      enddo
    enddo

    do k = kts,kte
      do i = its,ite
        za(i,k) = 0.5_r8*(zq(i,k)+zq(i,k+1))
      enddo
    enddo

    do i=its,ite
      hpbl2 = hpbl(i)+10._r8
      kpbl2 = kpbl(i)
      kvar = 1
      do k=kts+1,MAX(kpbl(i),kts+1)
        IF (za(i,k)>300._r8) then
          kpbl2 = k
          IF (k == kpbl(i)) then
            hpbl2 = hpbl(i)+10._r8
          ELSE
            hpbl2 = za(i,k)+10._r8
          ENDIF
          exit
        ENDIF
      enddo

      if(xland1(i).gt.0._r8 .and. 2._r8*var(i).le.hpbl(i))then
        if(br1(i).gt.0._r8 .and. thvx(i,kpbl2)-thvx(i,kts) > 0._r8)then
          cleff    = sqrt(dxy(i)**2_r8 + dxyp(i)**2_r8)
          cleff    = (2.0_r8/sncleff) * max(dxmax_ss,cleff)
          coefm(i) = (1._r8 + ol(i)) ** (oa1(i)+1._r8)
          xlinv(i) = coefm(i) / cleff
          govrth(i)=g/(0.5_r8*(thvx(i,kpbl2)+thvx(i,kts)))
          bnrf=sqrt(govrth(i)*(thvx(i,kpbl2)-thvx(i,kts))/hpbl2)

          if(abs(bnrf/u1(i,kpbl2)).gt.xlinv(i))then
            tauwavex0=0.5_r8*bnrf*xlinv(i)*(2._r8*MIN(var(i),varmax))**2_r8*ro(i,kvar)*u1(i,kvar)
            tauwavex0=tauwavex0*ss_taper   ! "Scale-awareness"
          else
            tauwavex0=0._r8
          endif

          if(abs(bnrf/v1(i,kpbl2)).gt.xlinv(i))then
            tauwavey0=0.5_r8*bnrf*xlinv(i)*(2._r8*MIN(var(i),varmax))**2._r8*ro(i,kvar)*v1(i,kvar)
            tauwavey0=tauwavey0*ss_taper   ! "Scale-awareness"
          else
            tauwavey0=0._r8
          endif

          do k=kts,kpbl(i) !MIN(kpbl2+1,kte-1)
            utendwave(i,k)=-1._r8*tauwavex0*2._r8*max((1._r8-za(i,k)/hpbl2),0._r8)/hpbl2
            vtendwave(i,k)=-1._r8*tauwavey0*2._r8*max((1._r8-za(i,k)/hpbl2),0._r8)/hpbl2
          enddo
        endif
      endif
    enddo ! end i loop

    do k = kts,kte
      do i = its,ite
        dudt(i,k)  = dudt(i,k) + utendwave(i,k)
        dvdt(i,k)  = dvdt(i,k) + vtendwave(i,k)
        dtaux2d_ss(i,k) = utendwave(i,k)
        dtauy2d_ss(i,k) = vtendwave(i,k)
        dusfc_ss(i) = dusfc_ss(i) + utendwave(i,k) * del(i,k)
        dvsfc_ss(i) = dvsfc_ss(i) + vtendwave(i,k) * del(i,k)
      enddo
    enddo

  ENDIF  ! end if gsd_gwd_ss == .true.
  !================================================================
  !add Beljaars et al. (2004, QJRMS, equ. 16) form drag:
  !================================================================
  IF (gsd_gwd_fd.and.(ss_taper.GT.1.E-02) ) THEN

    utendform=0._r8
    vtendform=0._r8
    zq=0._r8

    if (.not.gsd_gwd_ss.and.(ss_taper.GT.1.E-02) ) THEN
      ! Defining layer height. This is already done above is small-scale GWD is used
      do k = kts,kte
        do i = its,ite
          zq(i,k+1) = dz2(i,k)+zq(i,k)
        enddo
      enddo

      do k = kts,kte
        do i = its,ite
          za(i,k) = 0.5_r8*(zq(i,k)+zq(i,k+1))
        enddo
      enddo
    endif

    do i=its,ite
      if (xland1(i) .gt. 0..and.2._r8*var(i).gt.0) then
        ! refer to Beljaars (2004) eq.16.
        a1=0.00026615161_r8*var(i)**2_r8
        a2=a1*0.005363_r8
        do k=kts,kte
          wsp=SQRT(u1(i,k)**2_r8 + v1(i,k)**2_r8)
          ! refer to Beljaars (2004) eq.16.
          ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759 
          utendform(i,k)=-0.0759_r8*wsp*u1(i,k)* &
                         EXP(-(za(i,k)/1500._r8)**1.5_r8)*a2*za(i,k)**(-1.2_r8)*ss_taper
          vtendform(i,k)=-0.0759_r8*wsp*v1(i,k)* &
                         EXP(-(za(i,k)/1500._r8)**1.5_r8)*a2*za(i,k)**(-1.2_r8)*ss_taper
        enddo
      endif
    enddo

    do k = kts,kte
      do i = its,ite
        dudt(i,k)  = dudt(i,k) + utendform(i,k)
        dvdt(i,k)  = dvdt(i,k) + vtendform(i,k)
        !limit drag tendency 
        !some tendency is likely to even overturn the wind,
        !making wind reverse in 1 timestep and reverse again in next,
        !this limitation may help to make model stable,
        !and no more wind reversal due to drag,
        !which is suppose to decelerate, not accelerate
        utendform(i,k)  = sign(min(abs(utendform(i,k)),abs(u1(i,k))/deltim),utendform(i,k))
        vtendform(i,k)  = sign(min(abs(vtendform(i,k)),abs(v1(i,k))/deltim),vtendform(i,k))
        dtaux2d_fd(i,k) = utendform(i,k)
        dtauy2d_fd(i,k) = vtendform(i,k)
        dusfc_fd(i) = dusfc_fd(i) + utendform(i,k) * del(i,k)
        dvsfc_fd(i) = dvsfc_fd(i) + vtendform(i,k) * del(i,k)
      enddo
    enddo
  ENDIF  ! end if gsd_gwd_fd == .true.
  !=======================================================
  ! More for the large-scale gwd component
  !=======================================================
  IF (gsd_gwd_ls.and.(ls_taper.GT.1.E-02) ) THEN
    !                                                                       
    !   now compute vertical structure of the stress.
    !
    do k = kts,kpblmax
      do i = its,ite
        if (k .le. kbl(i)) taup(i,k) = taub(i)
      enddo
    enddo

    if (scorer_on) then
      !
      !determination of the interface height for scorer adjustment
      !
      do i=its,ite
        iint=.false.
        do k=kpblmin,kte-1
          if (k.gt.kbl(i).and.usqj(i,k)-usqj(i,k-1).lt.0.and.(.not.iint)) then
            iint=.true.
            zl_hint(i)=zl(i,k+1)
          endif
        enddo
      enddo
    endif

    do k = kpblmin, kte-1                   ! vertical level k loop!
      kp1 = k + 1
      do i = its,ite
        !
        !   unstablelayer if ri < ric
        !   unstable layer if upper air vel comp along surf vel <=0 (crit lay)
        !   at (u-c)=0. crit layer exists and bit vector should be set (.le.)
        !
        if (k .ge. kbl(i)) then
          !we modify the criteria for unstable layer
          !that the lv is critical under 0.25
          !while we keep wave breaking ric for
          !other larger lv
          icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric_rig)&
                                .or. (velco(i,k) .le. 0.0_r8)
          brvf(i)  = max(bnv2(i,k),bnv2min) ! brunt-vaisala frequency squared
          brvf(i)  = sqrt(brvf(i))          ! brunt-vaisala frequency
        endif
      enddo

      do i = its,ite
        if (k .ge. kbl(i) .and. (.not. ldrag(i)))   then
          if (.not.icrilv(i) .and. taup(i,k) .gt. 0.0_r8 ) then
            temv = 1.0_r8 / velco(i,k)
            tem1 = coefm(i)/(dxy(i)/ncleff)*(ro(i,kp1)+ro(i,k))*brvf(i)*velco(i,k)*0.5_r8
            hd   = sqrt(taup(i,k) / tem1)
            fro  = brvf(i) * hd * temv
            !
            !  rim is the minimum-richardson number by shutts (1985)
            !
            tem2   = sqrt(usqj(i,k))
            tem    = 1._r8 + tem2 * fro
            rim    = usqj(i,k) * (1._r8-fro) / (tem * tem)

            !
            !  check stability to employ the 'saturation hypothesis'
            !  of lindzen (1981) except at tropospheric downstream regions
            !
            if (rim .le. ric) then  ! saturation hypothesis!
              if ((oa1(i) .le. 0._r8).or.(kp1 .ge. kpblmin )) then
                temc = 2.0_r8 + 1.0_r8 / tem2
                hd   = velco(i,k) * (2.0_r8*sqrt(temc)-temc) / brvf(i)
                taup(i,kp1) = tem1 * hd * hd
                !
                !  taup is restricted to monotoncally decrease
                !  to avoid unexpected high taup in calculation
                !
                taup(i,kp1)=min(tem1*hd*hd,taup(i,k))
                !
                !  add vertical decrease at low level below hint (Kim and Doyle 2005)
                !  where Ri first decreases
                !
                if (scorer_on.and.k.gt.klowtop(i).and.zl(i,k).le.zl_hint(i).and.k.lt.kte-1) then
                        l1=(9.81_r8*bnv2(i,kp1)/velco(i,kp1)**2)
                        l2=(9.81_r8*bnv2(i,k)/velco(i,k)**2)
                        taup(i,kp1)=min(taup(i,k),taup(i,k)*(l1/l2),tem1*hd*hd)
                endif
              endif
            else                    ! no wavebreaking!
              taup(i,kp1) = taup(i,k)
            endif
          endif
        endif
      enddo
    enddo

    if(lcap.lt.kte) then
      do klcap = lcapp1,kte
        do i = its,ite
          taup(i,klcap) = prsi(i,klcap) / prsi(i,lcap) * taup(i,lcap)
        enddo
      enddo
    endif

  ENDIF !END LARGE-SCALE TAU CALCULATION
  !===============================================================
  !COMPUTE BLOCKING COMPONENT 
  !===============================================================
  IF (gsd_gwd_bl.and.(ls_taper .GT. 1.E-02)) THEN

    do i = its,ite
      if(.not.ldrag(i)) then
        !
        !------- determine the height of flow-blocking layer
        !
        kblk = 0
        pe = 0.0_r8

        do k = kte, kpblmin, -1
          if(kblk.eq.0 .and. k.le.komax(i)) then
            !flow block appears within the reference level
            !compare potential energy and kinetic energy
            !divided by g*ro is to turn del(pa) into height
            pe = pe + bnv2(i,k)*(zl(i,komax(i))-zl(i,k))*del(i,k)/g/ro(i,k)
            ke = 0.5_r8*((rcs*u1(i,k))**2._r8+(rcs*v1(i,k))**2._r8)
            !
            !---------- apply flow-blocking drag when pe >= ke 
            !
            if(pe.ge.ke) then
              kblk = k
              kblk = min(kblk,kbl(i))
              zblk = zl(i,kblk)-zl(i,kts)
            endif
          endif
        enddo

        if(kblk.ne.0) then
          !
          !--------- compute flow-blocking stress
          !

          !dxmax_ls is different than the usual one
          !because the taper is very different
          !dxy is a length scale mostly in the direction of the flow to the ridge
          !so it is good and not needed for an uneven grid area
          !ref Lott and Miller (1997) original scheme
          cd = max(2.0_r8-1.0_r8/od(i),0.0_r8)
          !
          !tuning of the drag magnitude
          cd=ncd*cd
          !
          taufb(i,kts) = 0.5_r8 * roll(i) * coefm(i) / max(dxmax_ls,dxy(i))**2 * cd * dxyp(i)   &
                         * olp(i) * zblk * ulow(i)**2
          !changed grid box area into dy*dy
          tautem = taufb(i,kts)/float(kblk-kts)
          do k = kts+1, kblk
            taufb(i,k) = taufb(i,k-1) - tautem
          enddo

          !
          !----------sum orographic GW stress and flow-blocking stress
          !
          !taup(i,:) = taup(i,:) + taufb(i,:)   ! Keep taup and taufb separate for now
        endif
      endif
    enddo

  ENDIF   ! end blocking drag
!===========================================================
  IF (gsd_gwd_ls.OR.gsd_gwd_bl.and.(ls_taper .GT. 1.E-02)) THEN
    !                                                                       
    !  calculate - (g)*d(tau)/d(pressure) and deceleration terms dtaux, dtauy
    !
    do k = kts,kte
      do i = its,ite
        taud_ls(i,k) = 1._r8 * (taup(i,k+1) - taup(i,k)) * csg / del(i,k)
        taud_bl(i,k) = 1._r8 * (taufb(i,k+1) - taufb(i,k)) * csg / del(i,k)
      enddo
    enddo
    !                                                                       
    !  limit de-acceleration (momentum deposition ) at top to 1/2 value 
    !  the idea is some stuff must go out the 'top'                     
    !                                                                       
    do klcap = lcap,kte
      do i = its,ite
        taud_ls(i,klcap) = taud_ls(i,klcap) * factop
        taud_bl(i,klcap) = taud_bl(i,klcap) * factop
      enddo
    enddo
    !                                                                       
    !  if the gravity wave drag would force a critical line             
    !  in the lower ksmm1 layers during the next deltim timestep,     
    !  then only apply drag until that critical line is reached.        
    !                                                                       
    do k = kts,kpblmax-1
      do i = its,ite
        if (k .le. kbl(i)) then
          if((taud_ls(i,k)+taud_bl(i,k)).ne.0._r8)                      &
            dtfac(i) = min(dtfac(i),abs(velco(i,k)                     &
                 /(deltim*rcs*(taud_ls(i,k)+taud_bl(i,k)))))
        endif
      enddo
    enddo

    do k = kts,kte
      do i = its,ite
        taud_ls(i,k)  = taud_ls(i,k) * dtfac(i) * ls_taper
        !apply limiter for ogwd
        !1.dudt < |c-u|/dt, so u-c cannot change sign(u^n+1 = u^n + du/dt * dt)
        !2.dudt<tndmax, eliminate ridiculous large tendency
        if (k.ne.kte) then!velco does not have top level value
          taud_ls(i,k)  = sign(min(abs(taud_ls(i,k)), 0.5*abs(velco(i,k))/deltim),taud_ls(i,k))
        endif
        taud_ls(i,k)  = sign(min(abs(taud_ls(i,k)),tndmax),taud_ls(i,k))
        taud_bl(i,k)  = taud_bl(i,k) * dtfac(i) * ls_taper
        dtaux2d_ls(i,k) = taud_ls(i,k) * xn(i)
        dtauy2d_ls(i,k) = taud_ls(i,k) * yn(i)
        dtaux2d_bl(i,k) = taud_bl(i,k) * xn(i)
        dtauy2d_bl(i,k) = taud_bl(i,k) * yn(i)
        dudt(i,k)  = dtaux2d_ls(i,k) + dtaux2d_bl(i,k) + dudt(i,k)
        dvdt(i,k)  = dtauy2d_ls(i,k) + dtauy2d_bl(i,k) + dvdt(i,k)
        dusfc_ls(i)  = dusfc_ls(i) + dtaux2d_ls(i,k) * del(i,k)
        dvsfc_ls(i)  = dvsfc_ls(i) + dtauy2d_ls(i,k) * del(i,k)
        dusfc_bl(i)  = dusfc_bl(i) + dtaux2d_bl(i,k) * del(i,k)
        dvsfc_bl(i)  = dvsfc_bl(i) + dtauy2d_bl(i,k) * del(i,k)
      enddo
    enddo

  ENDIF

  !  Finalize dusfc and dvsfc diagnoses
  do i = its,ite
    dusfc_ls(i) = (-1._r8/g*rcs) * dusfc_ls(i)
    dvsfc_ls(i) = (-1._r8/g*rcs) * dvsfc_ls(i)
    dusfc_bl(i) = (-1._r8/g*rcs) * dusfc_bl(i)
    dvsfc_bl(i) = (-1._r8/g*rcs) * dvsfc_bl(i)
    dusfc_ss(i) = (-1._r8/g*rcs) * dusfc_ss(i)
    dvsfc_ss(i) = (-1._r8/g*rcs) * dvsfc_ss(i)
    dusfc_fd(i) = (-1._r8/g*rcs) * dusfc_fd(i)
    dvsfc_fd(i) = (-1._r8/g*rcs) * dvsfc_fd(i)
  enddo

  return

end subroutine od2d

!==========================================================================


end module od_common
