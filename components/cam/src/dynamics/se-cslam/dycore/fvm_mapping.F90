!#define PCoM !replace PPM with PCoM for mass variables for fvm2phys and phys2fvm
!#define skip_high_order_fq_map !do mass and correlation preserving phys2fvm mapping but no high-order pre-mapping of fq
module fvm_mapping
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use dimensions_mod,         only: irecons_tracer
  use element_mod,            only: element_t
  use fvm_control_volume_mod, only: fvm_struct
  use perf_mod,       only: t_startf, t_stopf

  implicit none
  private

  public :: phys2dyn_forcings_fvm, dyn2phys, dyn2phys_vector, dyn2phys_all_vars,dyn2fvm_mass_vars
  public :: phys2dyn,fvm2dyn,fill_halo_phys
contains
  !
  ! map all mass variables from gll to fvm
  !
  subroutine phys2dyn_forcings_fvm(elem, fvm, hybrid,nets,nete,no_cslam, tl_f, tl_qdp)
    use dimensions_mod,         only: np, nc,nlev
    use dimensions_mod,         only: fv_nphys, nhc_phys,ntrac,nhc
    use dimensions_mod,         only: qsize_condensate_loading, qsize_condensate_loading_idx
    use fvm_control_volume_mod, only: n0_fvm
    use hybrid_mod,             only: hybrid_t
    use cam_abortutils,         only: endrun

    type (element_t), intent(inout):: elem(:)
    type(fvm_struct), intent(inout):: fvm(:)

    type (hybrid_t), intent(in)    :: hybrid  ! distributed parallel structure (shared)
    logical, intent(in)            :: no_cslam
    integer, intent(in)            :: nets, nete, tl_f, tl_qdp

    integer                                             :: ie, k, m_cnst,nq
    real (kind=r8), dimension(:,:,:,:,:)  , allocatable :: fld_phys, fld_gll, fld_fvm
    real (kind=r8), dimension(np,np,nlev,qsize_condensate_loading,nets:nete) :: qgll
    !
    ! for tensor product Lagrange interpolation
    !
    integer              :: nflds
    logical, allocatable :: llimiter(:)

    do ie=nets,nete
      do nq=1,qsize_condensate_loading
        qgll(:,:,:,nq,ie) = elem(ie)%state%Qdp(:,:,:,nq,tl_qdp)/elem(ie)%state%dp3d(:,:,:,tl_f)
      end do
    end do

    if (no_cslam) then
      call endrun("phys2dyn_forcings_fvm: no cslam case: NOT SUPPORTED")
    else if (nc.ne.fv_nphys) then
      !
      !***********************************************************
      !
      ! using cslam and different resolution physics grid
      !
      !***********************************************************
      !
      nflds = 4+ntrac
      allocate(fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,nflds,nets:nete))
      allocate(fld_gll(np,np,nlev,3,nets:nete))
      allocate(llimiter(nflds))
      fld_phys = -9.99E99_r8

      llimiter          = .false.
      do ie=nets,nete
        !
        ! pack fields that need to be interpolated
        !
        fld_phys(1:fv_nphys,1:fv_nphys,:,1,ie)       = fvm(ie)%ft(1:fv_nphys,1:fv_nphys,:)
        fld_phys(1:fv_nphys,1:fv_nphys,:,2,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,1,:)
        fld_phys(1:fv_nphys,1:fv_nphys,:,3,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,2,:)
        fld_phys(1:fv_nphys,1:fv_nphys,:,4,ie)       = fvm(ie)%dp_phys(1:fv_nphys,1:fv_nphys,:)
        do m_cnst=1,ntrac
          fld_phys(1:fv_nphys,1:fv_nphys,:,4+m_cnst,ie) = &
               fvm(ie)%fc_phys(1:fv_nphys,1:fv_nphys,:,m_cnst)
        end do
      end do
      call fill_halo_phys(elem,fld_phys,hybrid,nets,nete,nlev,nflds)
      !
      ! do mapping of fu,fv,ft
      !
      call phys2dyn(hybrid,elem,fld_phys(:,:,:,1:3,:),fld_gll(:,:,:,1:3,:),nets,nete,nlev,3,fvm,llimiter(1:3),2,.true.)
      do ie=nets,nete
        elem(ie)%derived%fT(:,:,:)   = fld_gll(:,:,:,1,ie)
        elem(ie)%derived%fM(:,:,1,:) = fld_gll(:,:,:,2,ie)
        elem(ie)%derived%fM(:,:,2,:) = fld_gll(:,:,:,3,ie)
      end do

      deallocate(fld_gll)
      !
      ! map fq from phys to fvm
      !
      do ie=nets,nete
         do k=1,nlev
!#ifdef just_map_q
!           do m_cnst=1,ntrac
!             call phys2fvm_scalar(ie,fvm(ie),fld_phys(:,:,k,4+m_cnst,ie),fvm(ie)%fc(:,:,k,m_cnst))
!             fvm(ie)%fc(:,:,k,m_cnst) = fvm(ie)%fc(:,:,k,m_cnst)*fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)
!           end do
!#else
           call phys2fvm(ie,k,fvm(ie),fld_phys(:,:,k,4,ie),&
                fld_phys(:,:,k,5:4+ntrac,ie),fvm(ie)%fc(:,:,k,1:ntrac),ntrac)
!#endif
         end do
       end do       
       !
       ! overwrite SE Q with cslam Q
       !
       nflds = qsize_condensate_loading
       allocate(fld_gll(np,np,nlev,nflds,nets:nete))
       allocate(fld_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,nflds,nets:nete))
       do ie=nets,nete
         !
         ! compute cslam updated Q value         
         do m_cnst=1,qsize_condensate_loading
           fld_fvm(1:nc,1:nc,:,m_cnst,ie) = fvm(ie)%c(1:nc,1:nc,:,qsize_condensate_loading_idx(m_cnst),n0_fvm)+&
                fvm(ie)%fc(1:nc,1:nc,:,qsize_condensate_loading_idx(m_cnst))/fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm)
         enddo
       end do
       llimiter(1:nflds) = .false.
       call fvm2dyn(elem,fld_fvm,fld_gll(:,:,:,1:nflds,:),hybrid,nets,nete,nlev,nflds,fvm,llimiter(1:nflds))
       !
       ! fld_gll now holds q cslam value on gll grid
       !
       ! convert fld_gll to increment (q_new-q_old)
       !
       do ie=nets,nete
         do m_cnst=1,qsize_condensate_loading
           elem(ie)%derived%fq(:,:,:,m_cnst)   =&
                fld_gll(:,:,:,m_cnst,ie)-qgll(:,:,:,m_cnst,ie)
         end do
       end do
       deallocate(fld_fvm)
     else
       !
       !
       !*****************************************************************************************
       !
       ! using cslam with same physics grid resolution as cslam resolution
       !
       !*****************************************************************************************
       !
       nflds = 3+qsize_condensate_loading
       allocate(fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,nflds,nets:nete))
       allocate(fld_gll(np,np,nlev,nflds,nets:nete))
       allocate(llimiter(nflds))
       llimiter(1:nflds) = .false.
       do ie=nets,nete
         !
         ! pack fields that need to be interpolated
         !
         fld_phys(1:fv_nphys,1:fv_nphys,:,1,ie)       = fvm(ie)%ft(1:fv_nphys,1:fv_nphys,:)
         fld_phys(1:fv_nphys,1:fv_nphys,:,2,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,1,:)
         fld_phys(1:fv_nphys,1:fv_nphys,:,3,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,2,:)
         !
         ! compute cslam mixing ratio with physics update
         !
         do m_cnst=1,qsize_condensate_loading
           do k=1,nlev
             fld_phys(1:fv_nphys,1:fv_nphys,k,m_cnst+3,ie) = &
                  fvm(ie)%c(1:fv_nphys,1:fv_nphys,k,qsize_condensate_loading_idx(m_cnst),n0_fvm)+&
                  fvm(ie)%fc_phys(1:fv_nphys,1:fv_nphys,k,qsize_condensate_loading_idx(m_cnst))
           end do
         end do
       end do
       !
       ! do mapping
       !
       call phys2dyn(hybrid,elem,fld_phys,fld_gll,nets,nete,nlev,nflds,fvm,llimiter,2)
       do ie=nets,nete
         elem(ie)%derived%fT(:,:,:)   = fld_gll(:,:,:,1,ie)
         elem(ie)%derived%fM(:,:,1,:) = fld_gll(:,:,:,2,ie)
         elem(ie)%derived%fM(:,:,2,:) = fld_gll(:,:,:,3,ie)
       end do
       do ie=nets,nete
         do m_cnst=1,qsize_condensate_loading
           !
           ! convert fq so that it will effectively overwrite SE q with CSLAM q
           !
           elem(ie)%derived%fq(:,:,:,m_cnst) = fld_gll(:,:,:,m_cnst+3,ie)-&
                qgll(:,:,:,m_cnst,ie)
         end do
         do m_cnst = 1,ntrac
           fvm(ie)%fc(1:nc,1:nc,:,m_cnst) = fvm(ie)%fc_phys(1:nc,1:nc,:,m_cnst)*fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm)
         end do
       end do
     end if
     deallocate(fld_phys,llimiter,fld_gll)
  end subroutine phys2dyn_forcings_fvm

  subroutine fvm2dyn(elem,fld_fvm,fld_gll,hybrid,nets,nete,numlev,num_flds,fvm,llimiter)
    use dimensions_mod, only: np, nhc, nc
    use hybrid_mod    , only: hybrid_t

    use bndry_mod     , only: ghost_exchange
    use edge_mod      , only: initghostbuffer, freeghostbuffer,ghostpack,ghostunpack
    use edgetype_mod  , only: edgebuffer_t


    integer              , intent(in)    :: nets,nete,num_flds,numlev
    real (kind=r8), intent(inout) :: fld_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,numlev,num_flds,nets:nete)
    real (kind=r8), intent(out)   :: fld_gll(np,np,numlev,num_flds,nets:nete)
    type (hybrid_t)      , intent(in)    :: hybrid
    type (element_t)     , intent(in)    :: elem(nets:nete)
    type(fvm_struct)     , intent(in)    :: fvm(nets:nete)
    logical              , intent(in)    :: llimiter(num_flds)

    integer               :: ie, iwidth
    type (edgeBuffer_t)   :: cellghostbuf
    !
    !*********************************************
    !
    ! halo exchange
    !
    !*********************************************
    !
    call t_startf('fvm2dyn:initbuffer')
    call initghostbuffer(hybrid%par,cellghostbuf,elem,numlev*num_flds,nhc,nc)
    call t_stopf('fvm2dyn:initbuffer')
    do ie=nets,nete
       call ghostpack(cellghostbuf, fld_fvm(:,:,:,:,ie),numlev*num_flds,0,ie)
    end do
    call ghost_exchange(hybrid,cellghostbuf)
    do ie=nets,nete
       call ghostunpack(cellghostbuf, fld_fvm(:,:,:,:,ie),numlev*num_flds,0,ie)
    end do
    call freeghostbuffer(cellghostbuf)
    !
    ! mapping
    !
    iwidth=2
!    iwidth=1
    do ie=nets,nete
      call tensor_lagrange_interp(fvm(ie)%cubeboundary,np,nc,nhc,numlev,num_flds,fld_fvm(:,:,:,:,ie),&
           fld_gll(:,:,:,:,ie),llimiter,iwidth,fvm(ie)%norm_elem_coord)
    end do
  end subroutine fvm2dyn


  subroutine fill_halo_phys(elem,fld_phys,hybrid,nets,nete,num_lev,num_flds)
    use dimensions_mod, only: nhc_phys, fv_nphys
    use hybrid_mod    , only: hybrid_t
    use bndry_mod     , only: ghost_exchange
    use edge_mod      , only: initghostbuffer, freeghostbuffer, ghostpack, ghostunpack
    use edgetype_mod  , only: edgebuffer_t


    integer              , intent(in)    :: nets,nete,num_lev,num_flds
    real (kind=r8), intent(inout) :: fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,num_lev,num_flds, &
         nets:nete)
    type (hybrid_t)      , intent(in)    :: hybrid  ! distributed parallel structure (shared)
    type (element_t)     , intent(inout) :: elem(:)

    integer                 :: ie
    type (edgeBuffer_t)   :: cellghostbuf
    !
    !*********************************************
    !
    ! halo exchange
    !
    !*********************************************
    !
    call t_startf('fvm:fill_halo_phys')
    call t_startf('fvm:fill_halo_phys:initbuffer')
    call initghostbuffer(hybrid%par,cellghostbuf,elem,num_lev*num_flds,nhc_phys,fv_nphys)
    call t_stopf('fvm:fill_halo_phys:initbuffer')
    do ie=nets,nete
       call ghostpack(cellghostbuf, fld_phys(:,:,:,:,ie),num_lev*num_flds,0,ie)
    end do
    call ghost_exchange(hybrid,cellghostbuf)
    do ie=nets,nete
       call ghostunpack(cellghostbuf, fld_phys(:,:,:,:,ie),num_lev*num_flds,0,ie)
    end do
    call freeghostbuffer(cellghostbuf)
    !
    call t_stopf('fvm:fill_halo_phys')
  end subroutine fill_halo_phys
  !
  ! must call fill_halo_phys before calling this subroutine
  !
  subroutine phys2dyn(hybrid,elem,fld_phys,fld_gll,nets,nete,num_lev,num_flds,fvm,llimiter,istart_vector,halo_filled)
    use dimensions_mod, only: np, nhc_phys, fv_nphys
    use hybrid_mod, only : hybrid_t
    type (hybrid_t), intent(in)   :: hybrid  ! distributed parallel structure (shared)
    integer       , intent(in)    :: nets,nete,num_flds,num_lev
    real (kind=r8), intent(inout) :: fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,num_lev,num_flds, &
         nets:nete)
    real (kind=r8), intent(out)   :: fld_gll(np,np,num_lev,num_flds,nets:nete)
    type (element_t)     , intent(inout) :: elem(:)
    type(fvm_struct)     , intent(in)    :: fvm(:)
    integer, optional    , intent(in)    :: istart_vector
    logical              , intent(in)    :: llimiter(num_flds)
    logical, optional    , intent(in)    :: halo_filled

    integer                 :: i, j, ie, k, iwidth
    real (kind=r8)   :: v1,v2

    if (present(halo_filled)) then
      if (.not.halo_filled) call fill_halo_phys(elem,fld_phys,hybrid,nets,nete,num_lev,num_flds)
    else
      call fill_halo_phys(elem,fld_phys,hybrid,nets,nete,num_lev,num_flds)
    end if
    if (present(istart_vector)) then
      do ie=nets,nete
        do k=1,num_lev
          do j=1-nhc_phys,fv_nphys+nhc_phys
            do i=1-nhc_phys,fv_nphys+nhc_phys
              !
              ! convert lat-lon vectors to contra-variant gnomonic
              !
              v1 = fld_phys(i,j,k,istart_vector  ,ie)
              v2 = fld_phys(i,j,k,istart_vector+1,ie)
              fld_phys(i,j,k,istart_vector  ,ie)=fvm(ie)%Dinv_physgrid(i,j,1,1)*v1 + fvm(ie)%Dinv_physgrid(i,j,1,2)*v2
              fld_phys(i,j,k,istart_vector+1,ie)=fvm(ie)%Dinv_physgrid(i,j,2,1)*v1 + fvm(ie)%Dinv_physgrid(i,j,2,2)*v2
            end do
          end do
        end do
      end do
    end if
    !
    ! mapping
    !
    iwidth=2
!    iwidth=1
    if (fv_nphys==1) iwidth=1
    do ie=nets,nete
      call tensor_lagrange_interp(fvm(ie)%cubeboundary,np,fv_nphys,nhc_phys,num_lev,num_flds,fld_phys(:,:,:,:,ie),&
           fld_gll(:,:,:,:,ie),llimiter,iwidth,fvm(ie)%norm_elem_coord_physgrid)
    end do

    if (present(istart_vector)) then
      !
      ! convert contra-variant to lat-lon
      !
      do ie=nets,nete
        do k=1,num_lev
          do j=1,np
            do i=1,np
              v1 = fld_gll(i,j,k,istart_vector  ,ie)
              v2 = fld_gll(i,j,k,istart_vector+1,ie)
              fld_gll(i,j,k,istart_vector  ,ie) = elem(ie)%D(i,j,1,1)*v1 + elem(ie)%D(i,j,1,2)*v2
              fld_gll(i,j,k,istart_vector+1,ie) = elem(ie)%D(i,j,2,1)*v1 + elem(ie)%D(i,j,2,2)*v2
            end do
          end do
        end do
      end do
    end if
  end subroutine phys2dyn
  !
  ! map all mass variables from gll to fvm
  !
  subroutine dyn2fvm_mass_vars(dp_gll,ps_gll,q_gll,&
       dp_fvm,ps_fvm,q_fvm,num_trac,metdet,inv_area)
    use dimensions_mod, only: np, nc,nlev
    integer, intent(in) :: num_trac
    real (kind=r8), dimension(np,np,nlev)         , intent(in) :: dp_gll
    real (kind=r8), dimension(np,np,nlev,num_trac), intent(in) :: q_gll
    real (kind=r8), dimension(np,np)              , intent(in) :: ps_gll


    real (kind=r8), dimension(nc,nc,nlev)         , intent(inout) :: dp_fvm
    real (kind=r8), dimension(nc,nc,nlev,num_trac), intent(inout) :: q_fvm
    real (kind=r8), dimension(nc,nc)     , intent(inout)        :: ps_fvm
    real (kind=r8), dimension(nc,nc)     , intent(out)          :: inv_area

    real (kind=r8), intent(in)           :: metdet(np,np)

    real (kind=r8) :: se_area_sphere(nc,nc), tmp(np,np)
    real (kind=r8) :: inv_darea_dp_fvm(nc,nc)
    integer        :: k,m_cnst

    tmp = 1.0_r8
    se_area_sphere = dyn2fvm(tmp,metdet)
    inv_area = 1.0_r8/se_area_sphere

    ps_fvm(:,:) = dyn2fvm(ps_gll,metdet,inv_area)
    do k=1,nlev
      dp_fvm(:,:,k) = dyn2fvm(dp_gll(:,:,k),metdet,inv_area)
      inv_darea_dp_fvm = inv_area/dp_fvm(:,:,k)
      do m_cnst=1,num_trac
        q_fvm(:,:,k,m_cnst) = &
             dyn2fvm(q_gll(:,:,k,m_cnst)*dp_gll(:,:,k),metdet,&
             inv_darea_dp_fvm,q_gll(:,:,k,m_cnst))
      end do
    end do
  end subroutine dyn2fvm_mass_vars
  !
  ! this subroutine assumes that the fvm halo has already been filled
  ! (if nc/=fv_nphys)
  !
  subroutine dyn2phys_all_vars(ie,dp_gll,T_gll,omega_gll,&
       dp_fvm,q_fvm,num_trac,metdet,fvm,ptop,&
       dp3d_phys,ps_phys,q_phys,T_phys,omega_phys,phis_phys)
    use dimensions_mod, only: np, nc,nlev,fv_nphys,nhc
    use dp_mapping,     only: nphys_pts

    integer, intent(in) :: ie,num_trac
    real (kind=r8), dimension(np,np,nlev)             , intent(in) :: dp_gll,T_gll,omega_gll

    real (kind=r8), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev)         , intent(inout) :: dp_fvm
    real (kind=r8), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,num_trac), intent(inout) :: q_fvm
    type(fvm_struct)                                                         , intent(in)    :: fvm

    real (kind=r8), intent(in)           :: metdet(np,np)
    real (kind=r8), intent(in)           :: ptop

    real (kind=r8), dimension(nphys_pts)               , intent(out) :: ps_phys,phis_phys
    real (kind=r8), dimension(nphys_pts,nlev)          , intent(out) :: dp3d_phys,T_phys,omega_phys
    real (kind=r8), dimension(nphys_pts,nlev,num_trac) , intent(out) :: q_phys


    real (kind=r8) :: tmp(np,np)
    real (kind=r8), dimension(fv_nphys,fv_nphys)          :: inv_area,inv_darea_dp_phys,se_area_sphere,dp3d_tmp
    real (kind=r8), dimension(fv_nphys,fv_nphys)          :: dp_phys_tmp
    real (kind=r8), dimension(fv_nphys,fv_nphys,num_trac) :: q_phys_tmp
    real (kind=r8), dimension(nc,nc)                      :: inv_darea_dp_fvm
    integer :: k,m_cnst

    tmp = 1.0_r8
    se_area_sphere = dyn2phys(tmp,metdet)
    inv_area = 1.0_r8/se_area_sphere
    phis_phys(:) = RESHAPE(fvm%phis_physgrid,SHAPE(phis_phys(:)))!not used !!!!!!

    ps_phys = ptop

    if (nc.ne.fv_nphys) then
      tmp = 1.0_r8
      do k=1,nlev
        inv_darea_dp_fvm = dyn2fvm(dp_gll(:,:,k),metdet)
        inv_darea_dp_fvm = 1.0_r8/inv_darea_dp_fvm

        T_phys(:,k) = RESHAPE(dyn2phys(T_gll(:,:,k),metdet,inv_area),SHAPE(T_phys(:,k)))
        Omega_phys(:,k) = RESHAPE(dyn2phys(Omega_gll(:,:,k),metdet,inv_area),SHAPE(Omega_phys(:,k)))

        call fvm2phys(ie,fvm,dp_fvm(:,:,k),dp_phys_tmp,q_fvm(:,:,k,:),q_phys_tmp,num_trac)
        
        dp3d_phys(:,k) = RESHAPE(dp_phys_tmp,SHAPE(dp3d_phys(:,k)))
        ps_phys(:) = ps_phys(:)+RESHAPE(dp_phys_tmp,SHAPE(ps_phys(:)))
        do m_cnst=1,num_trac
          q_phys(:,k,m_cnst) = RESHAPE(q_phys_tmp(:,:,m_cnst),SHAPE(q_phys(:,k,m_cnst)))
        end do
      end do
    else
      do k=1,nlev
        dp3d_tmp       = dyn2phys(dp_gll(:,:,k),metdet,inv_area)
        inv_darea_dp_phys = inv_area/dp3d_tmp
        T_phys(:,k) = RESHAPE(dyn2phys(T_gll(:,:,k)*dp_gll(:,:,k),metdet,&
             inv_darea_dp_phys),SHAPE(T_phys(:,k)))
        Omega_phys(:,k) = RESHAPE(dyn2phys(Omega_gll(:,:,k),metdet,inv_area),SHAPE(Omega_phys(:,k)))
        !
        ! no mapping needed - just copy fields into physics structure
        !
        dp3d_phys(:,k) = RESHAPE(dp_fvm(1:nc,1:nc,k),SHAPE(dp3d_phys(:,k)))
        ps_phys(:) = ps_phys(:)+RESHAPE(dp_fvm(1:nc,1:nc,k),SHAPE(ps_phys(:)))
        do m_cnst=1,num_trac
          q_phys(:,k,m_cnst) = RESHAPE(q_fvm(1:nc,1:nc,k,m_cnst),SHAPE(q_phys(:,k,m_cnst)))
        end do
      end do
    end if
  end subroutine dyn2phys_all_vars


  function dyn2phys(qdp_gll,metdet,inv_dp_darea_phys) result(qdp_phys)
    use dimensions_mod, only: np, nc, fv_nphys
    use derivative_mod, only: subcell_integration
    real (kind=r8), intent(in)           :: qdp_gll(np,np)
    real (kind=r8)                       :: qdp_phys(fv_nphys,fv_nphys)
    real (kind=r8), intent(in)           :: metdet(np,np)
    real (kind=r8), intent(in), optional :: inv_dp_darea_phys(fv_nphys,fv_nphys)
    integer :: i,j
    real (kind=r8) :: min_val, max_val

    call subcell_integration(qdp_gll(:,:), np, fv_nphys, metdet,qdp_phys,nc.ne.fv_nphys)
    if (present(inv_dp_darea_phys)) &
      qdp_phys = qdp_phys*inv_dp_darea_phys ! convert qdp to q
  end function dyn2phys


  function dyn2fvm(qdp_gll,metdet,inv_dp_darea_phys,q_gll) result(qdp_phys)
    use dimensions_mod, only: np, nc
    use derivative_mod, only: subcell_integration
    real (kind=r8), intent(in)           :: qdp_gll(np,np)
    real (kind=r8)                       :: qdp_phys(nc,nc)
    real (kind=r8), intent(in)           :: metdet(np,np)
    real (kind=r8), intent(in), optional :: inv_dp_darea_phys(nc,nc)
    real (kind=r8), intent(in), optional :: q_gll(np,np)
    integer :: i,j
    real (kind=r8) :: min_val, max_val

    call subcell_integration(qdp_gll(:,:), np, nc, metdet,qdp_phys)
    if (present(inv_dp_darea_phys)) then
      !
      ! convert qdp to q
      !
      qdp_phys = qdp_phys*inv_dp_darea_phys
      !
      ! simple limiter
      !
      if (present(q_gll)) then
        do j = 1, nc
          do i = 1, nc
            !
            ! simple limiter: only coded for nc=3 and np4
            !
            min_val = minval(q_gll(i:i+1,j:j+1))
            max_val = maxval(q_gll(i:i+1,j:j+1))
            qdp_phys(i,j) = max(min_val,min(max_val,qdp_phys(i,j)))
          end do
        end do
      end if
    end if
  end function dyn2fvm

  function dyn2phys_vector(v_gll,elem) result(v_phys)
    use dimensions_mod, only: np, nlev,  fv_nphys
    use interpolate_mod,only: interpdata_t,interpolate_2d,interpolate_t
    use cube_mod       ,only: dmap
    use control_mod    ,only: cubed_sphere_map

    type (interpdata_t):: interpdata
    type (element_t), intent(in)   :: elem
    type (interpolate_t) , target :: interp_p
    real (kind=r8), intent(in)           :: v_gll(np,np,2,nlev)
    real (kind=r8)                       :: v_phys(fv_nphys*fv_nphys,2,nlev)

    integer :: i,j,k

    ! Local variables
    real (kind=r8)    ::  fld_contra(np,np,2,nlev) ! vector field

    real (kind=r8)    ::  v1,v2
    real (kind=r8)    ::  D(2,2,fv_nphys*fv_nphys)   ! derivative of gnomonic mapping
    !
    ! this could be done at initialization and does not need to be repeated
    !
    call setup_interpdata_for_gll_to_phys_vec_mapping(interpdata, interp_p)
    ! convert to contra
    do k=1,nlev
      do j=1,np
        do i=1,np
          ! latlon->contra
          fld_contra(i,j,1,k) = elem%Dinv(i,j,1,1)*v_gll(i,j,1,k) + elem%Dinv(i,j,1,2)*v_gll(i,j,2,k)
          fld_contra(i,j,2,k) = elem%Dinv(i,j,2,1)*v_gll(i,j,1,k) + elem%Dinv(i,j,2,2)*v_gll(i,j,2,k)
        enddo
      enddo
    end do

    do k=1,nlev
      do i=1,interpdata%n_interp
        v_phys(i,1,k)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp_p,np)
        v_phys(i,2,k)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp_p,np)
      end do
    end do
    do i=1,interpdata%n_interp
      ! convert fld from contra->latlon
      call dmap(D(:,:,i),interpdata%interp_xy(i)%x,interpdata%interp_xy(i)%y,&
           elem%corners3D,cubed_sphere_map,elem%corners,elem%u2qmap,elem%facenum)
    end do
    do k=1,nlev
      do i=1,interpdata%n_interp
        ! convert fld from contra->latlon
        v1 = v_phys(i,1,k)
        v2 = v_phys(i,2,k)

        v_phys(i,1,k)=D(1,1,i)*v1 + D(1,2,i)*v2
        v_phys(i,2,k)=D(2,1,i)*v1 + D(2,2,i)*v2
      end do
    end do
  end function dyn2phys_vector

  subroutine setup_interpdata_for_gll_to_phys_vec_mapping(interpdata,interp_p)
    !
    ! initialize interpolation data structures to interpolate to phys grid
    ! using interpolate_mod subroutines
    !
    use interpolate_mod, only: interpolate_t, interpdata_t, interpolate_create
    use dimensions_mod, only : np
    use quadrature_mod, only : quadrature_t, gausslobatto
    use dimensions_mod, only : fv_nphys
    type (interpdata_t)  , intent(out)         :: interpdata
    type (interpolate_t) , intent(out), target :: interp_p

    ! local
    type (quadrature_t)   :: gp_quadrature
    integer i,j,ioff,ngrid
    real (kind=r8) ::  dx

    ngrid = fv_nphys*fv_nphys
    interpdata%n_interp=ngrid
    !
    ! initialize interpolation stuff related to basis functions
    !
    gp_quadrature = gausslobatto(np)
    call interpolate_create(gp_quadrature,interp_p)
    allocate(interpdata%interp_xy(ngrid))
    allocate(interpdata%ilat(ngrid) )
    allocate(interpdata%ilon(ngrid) )
    !
    !WARNING: THIS CODE INTERFERES WITH LAT-LON OUTPUT
    !         OF REGULAR SE IF nc>0
    !
    ioff=1
    dx = 2.0_r8/dble(fv_nphys)
    do j=1,fv_nphys
      do i=1,fv_nphys
        interpdata%interp_xy(ioff)%x = -1_r8+(i-0.5_r8)*dx
        interpdata%interp_xy(ioff)%y = -1_r8+(j-0.5_r8)*dx
        interpdata%ilon(ioff) = i
        interpdata%ilat(ioff) = j
        ioff=ioff+1
      enddo
    enddo
  end subroutine setup_interpdata_for_gll_to_phys_vec_mapping


  function lagrange_1d(src_grid,src_val,ngrid,dst_point,iwidth) result(val)
    integer              , intent(in)  :: ngrid,iwidth
    real (kind=r8), intent(in)  :: src_grid(ngrid), src_val(ngrid)
    real (kind=r8)              :: val

    real (kind=r8), intent(in)  :: dst_point

    integer :: iref, j,k
    real (kind=r8)              :: w(ngrid)

    if (dst_point.LE.src_grid(1)) then
      iref=1
    else
      iref=1
      do while (dst_point>src_grid(iref))
        iref = iref + 1
        if (iref>ngrid) then
          exit
        end if
      end do
      iref=iref-1
    end if

    iref=MIN(MAX(iref,iwidth),ngrid-iwidth)

    w = 1.0_r8
    do j=iref-(iwidth-1),iref+iwidth
      do k=iref-(iwidth-1),iref+iwidth
        if (k.ne.j) then
          w(j)=w(j)*(dst_point-src_grid(k))/(src_grid(j)-src_grid(k))
        end if
      end do
    end do

    val=0.0_r8
    do j=iref-(iwidth-1),iref+iwidth
      val=val+w(j)*src_val(j)
    end do
  end function lagrange_1d

  subroutine tensor_lagrange_interp(cubeboundary,np,nc,nhc,num_lev,nflds,psi,interp_value,llimiter,iwidth,norm_elem_coord)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    implicit none

    integer              , intent(in)    :: cubeboundary,nc, np, iwidth,nhc,num_lev,nflds
    logical              , intent(in)    :: llimiter(nflds)                   !apply limiter
    real (kind=r8), intent(inout) :: psi(1-nhc:nc+nhc,1-nhc:nc+nhc,num_lev,nflds) !fvm grid values with filled halo
    real (kind=r8), intent(out)   :: interp_value(np,np,num_lev,nflds)            !interpolated field
    real (kind=r8), intent(in)    :: norm_elem_coord(2,1-nhc:nc+nhc,1-nhc:nc+nhc)
    integer :: which_nc_cell(np)

    real (kind=r8):: dx,gll_points(np)
    real (kind=r8):: nc_points(1-nc:nc+nc)

    real (kind=r8):: value(1-iwidth:nc+iwidth)
    real (kind=r8):: val_tmp(1-nhc:nc+nhc,1-nhc:nc+nhc)

    real (kind=r8):: min_value(np,np,num_lev,nflds), max_value(np,np,num_lev,nflds)

    integer :: imin(1-nhc:nc+nhc), imax(1-nhc:nc+nhc)
    integer :: k,i,j,isearch,igll,jgll,jrow,h,irow,itr

    gll_points(1) = -1.0_r8
    gll_points(2) = -sqrt(1.0_r8/5.0_r8)
    gll_points(3) =  sqrt(1.0_r8/5.0_r8)
    gll_points(4) =  1.0_r8

    dx = 2_r8/dble(nc)
    do k=1-nc,2*nc
      nc_points(k) = -1.0_r8+dx*0.5_r8+dble(k-1)*dx
    end do
    !
    ! find fvm point surrounding gll points for simple limiter
    !
    do k=1,np
      do isearch=0,nc+1
        if (nc_points(isearch)<gll_points(k).and.nc_points(isearch+1).ge.gll_points(k)) exit
      end do
      which_nc_cell(k)=isearch
    end do
    do itr=1,nflds
      if (llimiter(itr)) then
        !
        ! fill non-existent halo cells for limiter
        !
        if (cubeboundary>4) then
          h=1
          select case(cubeboundary)
          case (nwest)
            psi(0,nc+h  ,:,itr) = psi(1-h,nc  ,:,itr)
            psi(1-h,nc+1,:,itr) = psi(1  ,nc+h,:,itr)
          case (swest)
            psi(1-h,0,:,itr) = psi(1,1-h,:,itr)
            psi(0,1-h,:,itr) = psi(1-h,1,:,itr)
          case (seast)
            psi(nc+h,0,:,itr) = psi(nc,1-h,:,itr)
            psi(nc+1,1-h,:,itr) = psi(nc+h,1,:,itr)
          case (neast)
            psi(nc+h,nc+1,:,itr) = psi(nc,nc+h,:,itr)
            psi(nc+1,nc+h,:,itr) = psi(nc+h,nc,:,itr)
          end select
        end if
        do k=1,num_lev
          do j=1,np
            do i=1,np
              max_value(i,j,k,itr) = max(&
                   psi(which_nc_cell(i)  ,which_nc_cell(j)  ,k,itr),&
                   psi(which_nc_cell(i)+1,which_nc_cell(j)  ,k,itr),&
                   psi(which_nc_cell(i)  ,which_nc_cell(j)+1,k,itr),&
                   psi(which_nc_cell(i)+1,which_nc_cell(j)+1,k,itr) &
                   )
              min_value(i,j,k,itr) = min(&
                   psi(which_nc_cell(i)  ,which_nc_cell(j)  ,k,itr),&
                   psi(which_nc_cell(i)+1,which_nc_cell(j)  ,k,itr),&
                   psi(which_nc_cell(i)  ,which_nc_cell(j)+1,k,itr),&
                   psi(which_nc_cell(i)+1,which_nc_cell(j)+1,k,itr) &
                   )
            end do
          end do
        end do
      end if
    end do

    imin=1-nhc
    imax=nc+nhc
    !
    ! special corner treatment
    !
    if (cubeboundary==swest) then
      do itr=1,nflds
        do k=1,num_lev
          do jrow=1,nc+iwidth
            !
            ! cubic along constant x (i=irow) in west halo to fvm points in halo
            !
            do irow=1-iwidth,0
              val_tmp(irow,jrow) = lagrange_1d(norm_elem_coord(2,irow,1:nc+nhc),psi(irow,1:nc+nhc,k,itr),nc+nhc,&
                   norm_elem_coord(2,1,jrow),iwidth)
            end do
          end do
          psi(1-iwidth:0,1:nc+iwidth,k,itr) = val_tmp(1-iwidth:0,1:nc+iwidth)
        enddo
      end do
      imin(1-nhc:0) = 1
    end if
    if (cubeboundary==nwest) then
      do itr=1,nflds
        do k=1,num_lev
          do jrow=1-iwidth,nc
            !
            ! cubic along constant x (i=irow) in west halo to fvm points in halo
            !
            do irow=1-iwidth,0
              val_tmp(irow,jrow) = lagrange_1d(norm_elem_coord(2,irow,1-nhc:nc),psi(irow,1-nhc:nc,k,itr),nc+nhc,&
                   norm_elem_coord(2,1,jrow),iwidth)
            end do
          end do
          psi(1-iwidth:0,1-iwidth:nc,k,itr) = val_tmp(1-iwidth:0,1-iwidth:nc)
        end do
      end do
      imin(nc+1:nc+nhc) = 1
    end if

    if (cubeboundary==seast) then
      do itr=1,nflds
        do k=1,num_lev
          do jrow=1,nc+iwidth
            value=0.0_r8
            !
            ! cubic along constant y in ease halo to fvm points in halo
            !
            do irow=nc+1,nc+iwidth
              val_tmp(irow,jrow) = lagrange_1d(norm_elem_coord(2,irow,1:nc+nhc),psi(irow,1:nc+nhc,k,itr),nc+nhc,&
                   norm_elem_coord(2,1,jrow),iwidth)
            end do
          end do
          psi(nc+1:nc+iwidth,1:nc+iwidth,k,itr) = val_tmp(nc+1:nc+iwidth,1:nc+iwidth)
        end do
      end do
      imax(1-nhc:0) = nc
    end if

    if (cubeboundary==neast) then
      do itr=1,nflds
        do k=1,num_lev
          do jrow=1-iwidth,nc
            !
            ! cubic along constant y in ease halo to fvm points in halo
            !
            do irow=nc+1,nc+iwidth
              val_tmp(irow,jrow) = lagrange_1d(norm_elem_coord(2,irow,1-nhc:nc),psi(irow,1-nhc:nc,k,itr),nc+nhc,&
                   norm_elem_coord(2,1,jrow),iwidth)
            end do
          end do
          psi(nc+1:nc+iwidth,1-iwidth:nc,k,itr) = val_tmp(nc+1:nc+iwidth,1-iwidth:nc)
        end do
      end do
      imax(nc+1:nc+nhc) = nc
    end if
    !
    ! mapping
    !
    !
    if (cubeboundary==0.or.cubeboundary==north.or.cubeboundary==south.or.&
         cubeboundary==swest.or.cubeboundary==nwest.or.&
         cubeboundary==seast.or.cubeboundary==neast) then
      do itr=1,nflds
        do k=1,num_lev
          do igll=1,np
            !
            ! cubic along constant y (j=jrow)
            !
            do jrow=1-iwidth,nc+iwidth
              value(jrow) = lagrange_1d(norm_elem_coord(1,imin(jrow):imax(jrow),jrow),psi(imin(jrow):imax(jrow),jrow,k,itr),&
                   imax(jrow)-imin(jrow)+1,gll_points(igll),iwidth)
            end do
            do jgll=1,np
              interp_value(igll,jgll,k,itr) = lagrange_1d(norm_elem_coord(2,1,1-iwidth:nc+iwidth),value,nc+2*iwidth,&
                   gll_points(jgll),iwidth)
            end do
          end do
        end do
      end do
    else if (cubeboundary==east.or.cubeboundary==west) then
      do itr=1,nflds
        do k=1,num_lev
          do jgll=1,np
            !
            ! cubic along constant x (i=irow)
            !
            do irow=1-iwidth,nc+iwidth
              value(irow) = lagrange_1d(norm_elem_coord(2,irow,1-nhc:nc+nhc),psi(irow,1-nhc:nc+nhc,k,itr),nc+2*nhc,&
                   gll_points(jgll),iwidth)
            end do
            do igll=1,np
              interp_value(igll,jgll,k,itr) = lagrange_1d(norm_elem_coord(1,1-iwidth:nc+iwidth,1),value,nc+2*iwidth,&
                   gll_points(igll),iwidth)
            end do
          end do
        end do
      end do
    end if
    do itr=1,nflds
      if (llimiter(itr)) then
        do k=1,num_lev
          do j=1,np
            do i=1,np
              interp_value(i,j,k,itr)=max(min_value(i,j,k,itr),min(max_value(i,j,k,itr),interp_value(i,j,k,itr)))
            end do
          enddo
        end do
      end if
    end do
  end subroutine tensor_lagrange_interp


  subroutine fvm2phys(ie,fvm,dp_fvm,dp_phys,q_fvm,q_phys,num_trac)
    use dimensions_mod, only: nc,nhc,fv_nphys
    !
    ! weights must be initialized in fvm2phys_init before using these functions
    !
    use dp_mapping, only: weights_all_fvm2phys, weights_eul_index_all_fvm2phys
    use dp_mapping, only: weights_lgr_index_all_fvm2phys, jall_fvm2phys

    type(fvm_struct)     , intent(in)           :: fvm
    integer              , intent(in)           :: ie
    integer              , intent(in)           :: num_trac
    real (kind=r8), intent(inout)        :: dp_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,1)
    real (kind=r8), intent(out)          :: dp_phys(fv_nphys,fv_nphys)

    real (kind=r8), intent(inout)        :: q_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,num_trac)
    real (kind=r8), intent(out)          :: q_phys(fv_nphys,fv_nphys,num_trac)

    real (kind=r8)                       :: recons    (irecons_tracer,1:nc,1:nc,1)
    real (kind=r8)                       :: recons_q  (irecons_tracer,1:nc,1:nc, &
                                            num_trac)

    real (kind=r8)                       :: recons_tmp(irecons_tracer)

    logical                              :: llimiter(1),llimiter_q(num_trac)
    integer                              :: h,jx,jy,jdx,jdy,m_cnst
    real (kind=r8)                       :: dp_phys_inv(fv_nphys,fv_nphys),dp_tmp

    llimiter=.false.
    call get_fvm_recons(fvm,dp_fvm,recons,1,llimiter)

    dp_phys = 0.0_r8
    do h=1,jall_fvm2phys(ie)
       jx  = weights_lgr_index_all_fvm2phys(h,1,ie)
       jy  = weights_lgr_index_all_fvm2phys(h,2,ie)
       jdx = weights_eul_index_all_fvm2phys(h,1,ie)
       jdy = weights_eul_index_all_fvm2phys(h,2,ie)
!#ifdef PCoM
!       dp_phys(jx,jy) = dp_phys(jx,jy) + weights_all_fvm2phys(h,1,ie)*dp_fvm(jdx,jdy,1)
!#else
       dp_phys(jx,jy) = dp_phys(jx,jy) + SUM(weights_all_fvm2phys(h,:,ie)*recons(:,jdx,jdy,1))
!#endif
    end do

    llimiter_q=.true.
    call get_fvm_recons(fvm,q_fvm,recons_q,num_trac,llimiter_q)
    !
    ! q-dp coupling as described in equation (55) in Appendinx B of
    ! Nair and Lauritzen, 2010: A Class of Deformational Flow Test Cases for Linear Transport Problems on the Sphere.
    ! J. Comput. Phys.: Vol. 229, Issue 23, pp. 8868-8887, DOI:10.1016/j.jcp.2010.08.014.
    !
    q_phys = 0.0_r8
    do h=1,jall_fvm2phys(ie)
       jx  = weights_lgr_index_all_fvm2phys(h,1,ie)
       jy  = weights_lgr_index_all_fvm2phys(h,2,ie)
       jdx = weights_eul_index_all_fvm2phys(h,1,ie)
       jdy = weights_eul_index_all_fvm2phys(h,2,ie)
       recons_tmp    = recons(:,jdx,jdy,1)
       recons_tmp(1) = recons(1,jdx,jdy,1)-dp_fvm(jdx,jdy,1)
       dp_tmp = SUM(weights_all_fvm2phys(h,:,ie)*recons_tmp(:))!does not need to be recomputed for each tracer
       do m_cnst=1,num_trac
!#ifdef PCoM
!         q_phys(jx,jy,m_cnst) = q_phys(jx,jy,m_cnst) + &
!              weights_all_fvm2phys(h,1,ie)*q_fvm(jdx,jdy,m_cnst)*dp_fvm(jdx,jdy,1)
!#else
         q_phys(jx,jy,m_cnst) = q_phys(jx,jy,m_cnst) + &
              dp_fvm(jdx,jdy,1)*SUM(weights_all_fvm2phys(h,:,ie)*recons_q(:,jdx,jdy,m_cnst))+&
              q_fvm(jdx,jdy,m_cnst) *dp_tmp
!#endif
       end do
    end do
    !
    ! convert to mixing ratio
    !
    dp_phys_inv=1.0_r8/dp_phys
    do m_cnst=1,num_trac
      q_phys(:,:,m_cnst) = q_phys(:,:,m_cnst)*dp_phys_inv(:,:)
    end do
    dp_phys = dp_phys/fvm%area_sphere_physgrid
  end subroutine fvm2phys

  subroutine phys2fvm(ie,k,fvm,dp_phys,fq_phys,fqdp_fvm,num_trac)
    use dimensions_mod, only: nhc_phys,fv_nphys,nc,nhc
    use fvm_control_volume_mod, only: n0_fvm

    integer              , intent(in)           :: ie,k
    type(fvm_struct)     , intent(inout)        :: fvm
    integer              , intent(in)           :: num_trac
    real (kind=r8), intent(inout)        :: dp_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys)
    real (kind=r8), intent(inout)        :: fq_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,num_trac)
    real (kind=r8), intent(out)          :: fqdp_fvm (nc,nc,num_trac)


    integer                              :: h,jx,jy,jdx,jdy,m_cnst

    real(kind=r8), dimension(fv_nphys,fv_nphys) :: phys_cdp_max, phys_cdp_min

    integer, parameter :: max_overlap = 4 !max number of mass overlap areas between phys and fvm grids
    integer      , dimension(fv_nphys,fv_nphys)               :: num_overlap
    integer      , dimension(2,max_overlap,fv_nphys,fv_nphys) :: overlap_idx
    real(kind=r8), dimension(max_overlap,fv_nphys,fv_nphys)   :: air_mass_overlap,dq_min_overlap,dq_max_overlap
    real(kind=r8), dimension(max_overlap,fv_nphys,fv_nphys)   :: phys_air_mass_overlap,fq_phys_overlap
    real(kind=r8), dimension(max_overlap,fv_nphys,fv_nphys)   :: q_overlap,overlap_area
    real(kind=r8), dimension(fv_nphys,fv_nphys)               :: q_phys
    integer       :: num
    real(kind=r8) :: tmp,dq_star,sum_dq_min,sum_dq_max

    real(kind=r8) :: mass_phys, mass_star
    real(kind=r8) :: min_patch,max_patch
    real (kind=r8):: q_prev,mass_forcing,mass_forcing_phys,inv_sum_dq 

    call get_dp_overlap(ie,k,fvm,max_overlap,air_mass_overlap,num_overlap,overlap_idx,overlap_area)
    do m_cnst=1,num_trac
      fqdp_fvm(:,:,m_cnst) = 0.0_r8
      call get_q_overlap(ie,k,fvm,max_overlap,air_mass_overlap,&
           fvm%c(1-nhc:nc+nhc,1-nhc:nc+nhc,k,m_cnst,n0_fvm),q_overlap,num_overlap,1,&
           dp_phys(1:fv_nphys,1:fv_nphys),q_phys)
!      call get_fq_overlap(ie,fvm,fvm%dp_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,k),&
      call get_fq_overlap(ie,fvm,dp_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys),&
           fq_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,m_cnst),max_overlap,&
           phys_air_mass_overlap,fq_phys_overlap,1)

      min_patch = MINVAL(fvm%c(0:nc+1,0:nc+1,k,m_cnst,n0_fvm))
      max_patch = MAXVAL(fvm%c(0:nc+1,0:nc+1,k,m_cnst,n0_fvm))
      do jy=1,fv_nphys
        do jx=1,fv_nphys
          num = num_overlap(jx,jy)
          tmp = q_phys(jx,jy)+fq_phys(jx,jy,m_cnst) !updated physics grid mixing ratio
          phys_cdp_max(jx,jy)= MAX(MAX(MAXVAL(q_overlap(1:num,jx,jy)),max_patch),tmp)
          phys_cdp_min(jx,jy)= MIN(MIN(MINVAL(q_overlap(1:num,jx,jy)),min_patch),tmp)
          !
          ! add high-order fq change when it does not violate monotonicity
          !
          mass_forcing_phys = 0.0_r8
          do h=1,num
            jdx = overlap_idx(1,h,jx,jy); jdy = overlap_idx(2,h,jx,jy)
            q_prev = q_overlap(h,jx,jy)   
!#ifndef skip_high_order_fq_map
            q_overlap(h,jx,jy) = q_overlap(h,jx,jy)+fq_phys_overlap(h,jx,jy)
            if (fq_phys_overlap(h,jx,jy)<0.0_r8) then
              q_overlap(h,jx,jy) = MAX(q_overlap(h,jx,jy),phys_cdp_min(jx,jy))            
            elseif (fq_phys_overlap(h,jx,jy)>0.0_r8) then
              q_overlap(h,jx,jy) = MIN(q_overlap(h,jx,jy),phys_cdp_max(jx,jy))
            end if
            mass_forcing = (q_overlap(h,jx,jy)-q_prev)*air_mass_overlap(h,jx,jy)
            mass_forcing_phys = mass_forcing_phys + mass_forcing
            fqdp_fvm(jdx,jdy,m_cnst) = fqdp_fvm(jdx,jdy,m_cnst)+mass_forcing 
!#endif
            !
            ! prepare for mass fixing algorithm
            !       
            dq_min_overlap(h,jx,jy)   = q_overlap(h,jx,jy)-phys_cdp_min(jx,jy)
            dq_max_overlap  (h,jx,jy) = q_overlap(h,jx,jy)-phys_cdp_max(jx,jy)
          end do
          fq_phys(jx,jy,m_cnst) = fq_phys(jx,jy,m_cnst)-mass_forcing_phys/&
               (dp_phys(jx,jy)*fvm%area_sphere_physgrid(jx,jy))
        end do
      end do
      !
      ! let physics mass tendency remove excess mass (as defined above) first proportional to how much is availabe
      ! 
      do jy=1,fv_nphys
        do jx=1,fv_nphys
          !
          ! total mass change from physics on physics grid
          !
          mass_phys = fq_phys(jx,jy,m_cnst)*fvm%area_sphere_physgrid(jx,jy)*dp_phys(jx,jy)  
          num = num_overlap(jx,jy)

          if (fq_phys(jx,jy,m_cnst)<0.0_r8) then
            sum_dq_min = SUM(dq_min_overlap(1:num,jx,jy))
            if (sum_dq_min>1.0E-14_r8) then
              inv_sum_dq = 1.0_r8/sum_dq_min
              mass_star = SUM(dq_min_overlap(1:num,jx,jy)*air_mass_overlap(1:num,jx,jy))
              tmp = mass_star*inv_sum_dq
              dq_star = mass_phys/tmp
              dq_star = MAX(-1.0_r8,dq_star*inv_sum_dq)
              do h=1,num
                jdx = overlap_idx(1,h,jx,jy); jdy = overlap_idx(2,h,jx,jy)
                fqdp_fvm(jdx,jdy,m_cnst) = fqdp_fvm(jdx,jdy,m_cnst)&
                     +dq_star*dq_min_overlap(h,jx,jy)*air_mass_overlap(h,jx,jy)
              end do
            end if
          end if
          if (fq_phys(jx,jy,m_cnst)>0.0_r8) then
            sum_dq_max = SUM(dq_max_overlap(1:num,jx,jy))
            if (sum_dq_max<-1.0E-14_r8) then
              inv_sum_dq = 1.0_r8/sum_dq_max
              mass_star = SUM(dq_max_overlap(1:num,jx,jy)*air_mass_overlap(1:num,jx,jy))
              tmp = mass_star*inv_sum_dq
              dq_star = mass_phys/tmp
              dq_star = MIN(1.0_r8,dq_star*inv_sum_dq)
              do h=1,num
                jdx = overlap_idx(1,h,jx,jy); jdy = overlap_idx(2,h,jx,jy)
                fqdp_fvm(jdx,jdy,m_cnst) = fqdp_fvm(jdx,jdy,m_cnst)&
                    +dq_star*dq_max_overlap(h,jx,jy)*air_mass_overlap(h,jx,jy)
              end do
            end if
          end if
        end do
      end do
      !
      ! convert to mass per unit area
      !
      fqdp_fvm(:,:,m_cnst) = fqdp_fvm(:,:,m_cnst)*fvm%inv_area_sphere(:,:)
    end do
  end subroutine phys2fvm



  subroutine get_dp_overlap(ie,k,fvm,max_overlap,air_mass_overlap,num_overlap,overlap_idx,overlap_area)
    use dimensions_mod, only: nc,nhr,nhc,fv_nphys
    use fvm_control_volume_mod, only: n0_fvm
    !
    ! weights must be initialized in fvm2phys_init before using these functions
    !
    use dp_mapping, only: weights_all_fvm2phys, weights_eul_index_all_fvm2phys
    use dp_mapping, only: weights_lgr_index_all_fvm2phys, jall_fvm2phys
    !
    ! setting nhe=0 because we do not need reconstruction outside of element
    !
    integer, parameter :: nh = nhr!+(nhe-1) ! = 2 (nhr=2; nhe_local=1),! = 3 (nhr=2; nhe_local=2)

    type(fvm_struct)                                         , intent(inout):: fvm
    integer                                                  , intent(in)   :: ie, k, max_overlap
    integer, dimension(fv_nphys,fv_nphys)                    , intent(out)  :: num_overlap
    integer      , dimension(2,max_overlap,fv_nphys,fv_nphys), intent(out)  :: overlap_idx
    real(kind=r8), dimension(max_overlap,fv_nphys,fv_nphys)  , intent(out)  :: air_mass_overlap,overlap_area


    real (kind=r8)                       :: recons    (irecons_tracer,nc,nc)
    logical                              :: llimiter(1)
    integer                              :: h,jx,jy,jdx,jdy,idx
    llimiter=.false.
    call get_fvm_recons(fvm,fvm%dp_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,k,n0_fvm),recons,1,llimiter)

    num_overlap(:,:) = 0
    do h=1,jall_fvm2phys(ie)
       jx  = weights_lgr_index_all_fvm2phys(h,1,ie); jy  = weights_lgr_index_all_fvm2phys(h,2,ie)       
       jdx = weights_eul_index_all_fvm2phys(h,1,ie); jdy = weights_eul_index_all_fvm2phys(h,2,ie)
       num_overlap(jx,jy) = num_overlap(jx,jy)+1
       idx = num_overlap(jx,jy)
       overlap_idx(1,idx,jx,jy) = jdx; overlap_idx(2,idx,jx,jy) = jdy;
       overlap_area(idx,jx,jy)  = weights_all_fvm2phys(h,1,ie)
       air_mass_overlap(idx,jx,jy) = SUM(weights_all_fvm2phys(h,:,ie)*recons(:,jdx,jdy))
!#ifdef PCoM
!       air_mass_overlap(idx,jx,jy) = fvm%dp_fvm(jdx,jdy,k,n0_fvm)*weights_all_fvm2phys(h,1,ie)!PCoM
!#endif
    end do
  end subroutine get_dp_overlap

  subroutine get_q_overlap(ie,k,fvm,max_overlap,air_mass_overlap,q_fvm,q_overlap,num_overlap,num_trac,dp_phys,q_phys)
    use fvm_control_volume_mod, only: n0_fvm
    use dimensions_mod, only: nc,nhr,nhc,fv_nphys
    !
    ! weights must be initialized in fvm2phys_init before using these functions
    !
    use dp_mapping, only: weights_all_fvm2phys, weights_eul_index_all_fvm2phys
    use dp_mapping, only: weights_lgr_index_all_fvm2phys, jall_fvm2phys
    !
    ! setting nhe=0 because we do not need reconstruction outside of element
    !
    integer, parameter :: nhe_local=0
    integer, parameter :: nh = nhr!+(nhe-1) ! = 2 (nhr=2; nhe_local=1),! = 3 (nhr=2; nhe_local=2)

    type(fvm_struct)     , intent(inout)        :: fvm
    integer              , intent(in)           :: ie, k, max_overlap
    integer              , intent(in)           :: num_trac

    real(kind=r8), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,1:num_trac)      :: q_fvm
    real(kind=r8), dimension(max_overlap,fv_nphys,fv_nphys,num_trac)    :: q_overlap
    real(kind=r8), dimension(max_overlap,fv_nphys,fv_nphys), intent(in) :: air_mass_overlap
    integer      , dimension(fv_nphys,fv_nphys), intent(out)            :: num_overlap
    real(kind=r8), dimension(fv_nphys,fv_nphys),             intent(in) :: dp_phys
    real(kind=r8), dimension(fv_nphys,fv_nphys, num_trac),   intent(out):: q_phys

    real (kind=r8)                       :: recons_q  (irecons_tracer,1:nc,1:nc,num_trac)
    logical                              :: llimiter_q(num_trac)
    integer                              :: h,jx,jy,jdx,jdy,m_cnst,idx
    real (kind=r8)                       :: dp_tmp, dp_fvm_tmp, tmp
    real (kind=r8)                       :: dp_phys_inv(fv_nphys,fv_nphys)

    llimiter_q=.true.
    call get_fvm_recons(fvm,q_fvm,recons_q,num_trac,llimiter_q)
    num_overlap(:,:) = 0
    q_phys = 0.0_r8
    do h=1,jall_fvm2phys(ie)
       jx  = weights_lgr_index_all_fvm2phys(h,1,ie); jy  = weights_lgr_index_all_fvm2phys(h,2,ie)       
       jdx = weights_eul_index_all_fvm2phys(h,1,ie); jdy = weights_eul_index_all_fvm2phys(h,2,ie)

       num_overlap(jx,jy) = num_overlap(jx,jy)+1
       idx = num_overlap(jx,jy)

       dp_fvm_tmp = fvm%dp_fvm(jdx,jdy,k,n0_fvm)
       dp_tmp = air_mass_overlap(idx,jx,jy)-dp_fvm_tmp*weights_all_fvm2phys(h,1,ie)
       do m_cnst=1,num_trac
         tmp = dp_fvm_tmp*SUM(weights_all_fvm2phys(h,:,ie)*recons_q(:,jdx,jdy,m_cnst))+q_fvm(jdx,jdy,m_cnst)*dp_tmp
!#ifdef PCoM
!         tmp = dp_fvm_tmp*weights_all_fvm2phys(h,1,ie)*q_fvm(jdx,jdy,m_cnst)
!#endif
         q_overlap(idx,jx,jy,m_cnst) = tmp/air_mass_overlap(idx,jx,jy)
         q_phys(jx,jy,m_cnst) = q_phys(jx,jy,m_cnst)+tmp
       end do
     end do
     !
     ! q_phys holds mass - convert to mixing ratio
     !
     dp_phys_inv = 1.0_r8/(dp_phys*fvm%area_sphere_physgrid)
     do m_cnst=1,num_trac
       q_phys(:,:,m_cnst) = q_phys(:,:,m_cnst)*dp_phys_inv
     end do
   end subroutine get_q_overlap

  subroutine get_q_fvm(ie,fvm,dp_phys,dp_fvm,q_phys,dp_q_fvm,num_trac,return_mixing_ratio)
    use dimensions_mod, only: fv_nphys,nhc_phys,fv_nphys,nhe_phys,nc
    !
    ! weights must be initialized in phys2fvm_init before using this function
    !
    use dp_mapping, only: weights_all_phys2fvm, weights_eul_index_all_phys2fvm
    use dp_mapping, only: weights_lgr_index_all_phys2fvm, jall_phys2fvm

    integer              , intent(in)           :: ie
    type(fvm_struct)     , intent(in)           :: fvm
    integer              , intent(in)           :: num_trac
    logical              , intent(in)           :: return_mixing_ratio
    real (kind=r8), intent(inout)        :: dp_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,1)
    real (kind=r8), intent(out)          :: dp_fvm(nc,nc)

    real (kind=r8), intent(inout)        :: q_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,num_trac)
    real (kind=r8), intent(out)          :: dp_q_fvm (nc,nc,num_trac)

    real (kind=r8)                       :: recons    (irecons_tracer,1-nhe_phys:fv_nphys+nhe_phys,&
         1-nhe_phys:fv_nphys+nhe_phys,1)
    real (kind=r8)                       :: recons_q  (irecons_tracer,1-nhe_phys:fv_nphys+nhe_phys,&
         1-nhe_phys:fv_nphys+nhe_phys,num_trac)
    real (kind=r8)                       :: recons_tmp(irecons_tracer)

    logical                              :: llimiter(1),llimiter_q(num_trac)
    integer                              :: h,jx,jy,jdx,jdy,m_cnst
    real (kind=r8)                       :: dp_fvm_inv(nc,nc), dp_tmp

    llimiter=.false.
    call get_physgrid_recons(fvm,dp_phys,recons,1,llimiter)

    dp_fvm = 0.0_r8
    do h=1,jall_phys2fvm(ie)
       jx  = weights_lgr_index_all_phys2fvm(h,1,ie)
       jy  = weights_lgr_index_all_phys2fvm(h,2,ie)
       jdx = weights_eul_index_all_phys2fvm(h,1,ie)
       jdy = weights_eul_index_all_phys2fvm(h,2,ie)
!#ifdef PCoM
!       dp_fvm(jx,jy) = dp_fvm(jx,jy) + weights_all_phys2fvm(h,1,ie)*dp_phys(jdx,jdy,1)
!#else
       dp_fvm(jx,jy) = dp_fvm(jx,jy) + SUM(weights_all_phys2fvm(h,:,ie)*recons(:,jdx,jdy,1))
!#endif
    end do

    llimiter_q=.true.
    call get_physgrid_recons(fvm,q_phys,recons_q,num_trac,llimiter_q)
    !
    ! q-dp coupling as described in equation (55) in Appendinx B of
    ! Nair and Lauritzen, 2010: A Class of Deformational Flow Test Cases for Linear Transport Problems on the Sphere.
    ! J. Comput. Phys.: Vol. 229, Issue 23, pp. 8868-8887, DOI:10.1016/j.jcp.2010.08.014.
    !
    dp_q_fvm = 0.0_r8
    do h=1,jall_phys2fvm(ie)
       jx  = weights_lgr_index_all_phys2fvm(h,1,ie)
       jy  = weights_lgr_index_all_phys2fvm(h,2,ie)
       jdx = weights_eul_index_all_phys2fvm(h,1,ie)
       jdy = weights_eul_index_all_phys2fvm(h,2,ie)
       recons_tmp    = recons(:,jdx,jdy,1)
       recons_tmp(1) = recons(1,jdx,jdy,1)-dp_phys(jdx,jdy,1)
       dp_tmp = SUM(weights_all_phys2fvm(h,:,ie)*recons_tmp(:))
       do m_cnst=1,num_trac
!#ifdef PCoM
!         dp_q_fvm(jx,jy,m_cnst) = dp_q_fvm(jx,jy,m_cnst) + &
!              dp_phys(jdx,jdy,1)*weights_all_phys2fvm(h,1,ie)*q_phys(jdx,jdy,m_cnst)
!#else
         dp_q_fvm(jx,jy,m_cnst) = dp_q_fvm(jx,jy,m_cnst) + &
              dp_phys(jdx,jdy,1)*SUM(weights_all_phys2fvm(h,:,ie)*recons_q(:,jdx,jdy,m_cnst))+&
              q_phys(jdx,jdy,m_cnst) *dp_tmp
!#endif
       end do
    end do
    !
    ! convert to mixing ratio
    !
    if (return_mixing_ratio) then
      dp_fvm_inv=1.0_r8/dp_fvm
      do m_cnst=1,num_trac
        dp_q_fvm(:,:,m_cnst) = dp_q_fvm(:,:,m_cnst)*dp_fvm_inv(:,:)
      end do
    end if
  end subroutine get_q_fvm

  subroutine get_fq_overlap(ie,fvm,dp_phys,fq_phys,max_overlap,phys_air_mass_overlap,fq_phys_overlap,num_trac)
    use dimensions_mod, only: fv_nphys,nhc_phys
    use dp_mapping, only: weights_all_fvm2phys,weights_lgr_index_all_fvm2phys, jall_fvm2phys

    integer              , intent(in)           :: ie
    type(fvm_struct)     , intent(in)           :: fvm
    integer              , intent(in)           :: num_trac, max_overlap
    real (kind=r8), intent(inout)        :: dp_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,1)
    real(kind=r8), dimension(max_overlap,fv_nphys,fv_nphys,num_trac),intent(out) :: fq_phys_overlap
    real(kind=r8), dimension(max_overlap,fv_nphys,fv_nphys)         ,intent(out) :: phys_air_mass_overlap

    real (kind=r8), intent(inout)        :: fq_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,num_trac)

    real (kind=r8)                       :: recons    (irecons_tracer,1:fv_nphys,fv_nphys,1)
    real (kind=r8)                       :: recons_q  (irecons_tracer,fv_nphys,fv_nphys,num_trac)

    logical                              :: llimiter(1),llimiter_q(num_trac)
    integer                              :: h,jx,jy,m_cnst
    real (kind=r8)                       :: dp_tmp, dp_phys_tmp
    integer:: idx
    integer, dimension(fv_nphys,fv_nphys):: num_overlap

    llimiter=.false.
    call get_physgrid_recons(fvm,dp_phys,recons,1,llimiter)
    num_overlap(:,:) = 0
    do h=1,jall_fvm2phys(ie)
       jx  = weights_lgr_index_all_fvm2phys(h,1,ie)
       jy  = weights_lgr_index_all_fvm2phys(h,2,ie)

       num_overlap(jx,jy) = num_overlap(jx,jy)+1
       idx = num_overlap(jx,jy)
!#ifdef PCoM
!       phys_air_mass_overlap(idx,jx,jy) = weights_all_fvm2phys(h,1,ie)*dp_phys(jx,jy,1)
!#else
       phys_air_mass_overlap(idx,jx,jy) = SUM(weights_all_fvm2phys(h,:,ie)*recons(:,jx,jy,1))
!#endif
    end do

    llimiter_q=.true.
    call get_physgrid_recons(fvm,fq_phys,recons_q,num_trac,llimiter_q)
    !
    ! q-dp coupling as described in equation (55) in Appendinx B of
    ! Nair and Lauritzen, 2010: A Class of Deformational Flow Test Cases for Linear Transport Problems on the Sphere.
    ! J. Comput. Phys.: Vol. 229, Issue 23, pp. 8868-8887, DOI:10.1016/j.jcp.2010.08.014.
    !
    num_overlap(:,:) = 0
    do h=1,jall_fvm2phys(ie)
       jx  = weights_lgr_index_all_fvm2phys(h,1,ie)
       jy  = weights_lgr_index_all_fvm2phys(h,2,ie)
       num_overlap(jx,jy) = num_overlap(jx,jy)+1
       idx = num_overlap(jx,jy)
       dp_phys_tmp = dp_phys(jx,jy,1)
       dp_tmp = phys_air_mass_overlap(idx,jx,jy)-dp_phys_tmp*weights_all_fvm2phys(h,1,ie)
       do m_cnst=1,num_trac
!#ifdef PCoM
!         fq_phys_overlap(idx,jx,jy,m_cnst) = dp_phys_tmp*weights_all_fvm2phys(h,1,ie)*fq_phys(jx,jy,m_cnst)
!#else
         fq_phys_overlap(idx,jx,jy,m_cnst) = &
              (dp_phys_tmp*SUM(weights_all_fvm2phys(h,:,ie)*recons_q(:,jx,jy,m_cnst))+&
              fq_phys(jx,jy,m_cnst)*dp_tmp)/phys_air_mass_overlap(idx,jx,jy)
!#endif
       end do
     end do
  end subroutine get_fq_overlap

  subroutine phys2fvm_scalar(ie,fvm,q_phys,dp_q_fvm)
    use dimensions_mod, only: fv_nphys,nhc_phys,nc
    !
    ! weights must be initialized in phys2fvm_init before using this function
    !
    use dp_mapping, only: weights_all_phys2fvm, weights_eul_index_all_phys2fvm
    use dp_mapping, only: weights_lgr_index_all_phys2fvm, jall_phys2fvm

    integer              , intent(in)    :: ie
    type(fvm_struct)     , intent(in)    :: fvm

    real (kind=r8), intent(inout)        :: q_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys)
    real (kind=r8), intent(out)          :: dp_q_fvm (nc,nc)

    real (kind=r8)                       :: recons_q(irecons_tracer,fv_nphys,fv_nphys)

    logical                              :: llimiter(1)
    integer                              :: h,jx,jy,jdx,jdy
    real (kind=r8)                       :: dp_fvm_inv(nc,nc)

    llimiter=.true.
    call get_physgrid_recons(fvm,q_phys,recons_q,1,llimiter)

    dp_q_fvm = 0.0_r8
    do h=1,jall_phys2fvm(ie)
       jx  = weights_lgr_index_all_phys2fvm(h,1,ie)
       jy  = weights_lgr_index_all_phys2fvm(h,2,ie)
       jdx = weights_eul_index_all_phys2fvm(h,1,ie)
       jdy = weights_eul_index_all_phys2fvm(h,2,ie)
       dp_q_fvm(jx,jy) = dp_q_fvm(jx,jy) + &
            SUM(weights_all_phys2fvm(h,:,ie)*recons_q(:,jdx,jdy))
    end do

    dp_fvm_inv=1.0_r8/fvm%area_sphere(:,:)
    dp_q_fvm(:,:) = dp_q_fvm(:,:)*dp_fvm_inv(:,:)
  end subroutine phys2fvm_scalar

  subroutine get_physgrid_recons(fvm,field_phys,recons_phys,num_trac,llimiter)
    use dimensions_mod, only: fv_nphys,nhr_phys,nhc_phys,ns_phys
    use fvm_reconstruction_mod, only: reconstruction
    type(fvm_struct), intent(in)           :: fvm
    integer,          intent(in)           :: num_trac
    real (kind=r8),   intent(inout)        :: field_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,num_trac)
    real (kind=r8),   intent(out)          :: recons_phys(irecons_tracer,1:fv_nphys,1:fv_nphys,num_trac)
    logical,          intent(in)           :: llimiter(num_trac)

    integer, dimension(3)                  :: jx_min_local, jx_max_local, jy_min_local, jy_max_local

    jx_min_local(1) = 1            ; jx_max_local(1) = fv_nphys+1
    jy_min_local(1) = 1            ; jy_max_local(1) = fv_nphys+1
    jx_min_local(2) = 0            ; jx_max_local(2) = -1
    jy_min_local(2) = 0            ; jy_max_local(2) = -1
    jx_min_local(3) = 0            ; jx_max_local(3) = -1
    jy_min_local(3) = 0            ; jy_max_local(3) = -1

    call reconstruction(field_phys,recons_phys,irecons_tracer,llimiter,num_trac,&
       fv_nphys,0,nhr_phys,nhc_phys,nhr_phys,ns_phys,nhr_phys,&
       jx_min_local,jx_max_local,jy_min_local,jy_max_local,&
       fvm%cubeboundary,fvm%halo_interp_weight_physgrid(1:ns_phys,1-nhr_phys:fv_nphys+nhr_phys,1:nhr_phys,:),&
       fvm%ibase_physgrid(1-nhr_phys:fv_nphys+nhr_phys,1:nhr_phys,:),&
       fvm%spherecentroid_physgrid(:,1:fv_nphys,1:fv_nphys),&
       fvm%recons_metrics_physgrid(:,1:fv_nphys,1:fv_nphys),&
       fvm%recons_metrics_integral_physgrid(:,1:fv_nphys,1:fv_nphys)    ,&
       fvm%rot_matrix_physgrid,&
       fvm%centroid_stretch_physgrid(1:7,1:fv_nphys,1:fv_nphys),&
       fvm%vertex_recons_weights_physgrid(1:irecons_tracer-1,:,1:fv_nphys,1:fv_nphys),&
       fvm%vtx_cart_physgrid(:,:,1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys))
  end subroutine get_physgrid_recons

  subroutine get_fvm_recons(fvm,field_fvm,recons_fvm,num_trac,llimiter)
    use dimensions_mod, only: nc,nhr,nhc,ns
    use fvm_reconstruction_mod, only: reconstruction

    type(fvm_struct), intent(in)   :: fvm
    integer,          intent(in)   :: num_trac
    logical,          intent(in)   :: llimiter(num_trac)
    real (kind=r8),   intent(inout):: field_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,num_trac)
    real (kind=r8),   intent(out)  :: recons_fvm(irecons_tracer,nc,nc,num_trac)
    integer                        :: jx_min_local(3), jx_max_local(3), jy_min_local(3), jy_max_local(3)

    jx_min_local(1) = 1            ; jx_max_local(1) = nc+1
    jy_min_local(1) = 1            ; jy_max_local(1) = nc+1
    jx_min_local(2) = 0            ; jx_max_local(2) = -1
    jy_min_local(2) = 0            ; jy_max_local(2) = -1
    jx_min_local(3) = 0            ; jx_max_local(3) = -1
    jy_min_local(3) = 0            ; jy_max_local(3) = -1

    call reconstruction(field_fvm,recons_fvm,irecons_tracer,&
         llimiter,num_trac,nc,0,nhr,nhc,nhr,ns,nhr,&
         jx_min_local,jx_max_local,jy_min_local,jy_max_local,&
         fvm%cubeboundary,fvm%halo_interp_weight(1:ns,1-nhr:nc+nhr,1:nhr,:),fvm%ibase(1-nhr:nc+nhr,1:nhr,:),&
         fvm%spherecentroid(:,1:nc,1:nc),&
         fvm%recons_metrics(:,1:nc,1:nc),&
         fvm%recons_metrics_integral(:,1:nc,1:nc)    ,&
         fvm%rot_matrix,fvm%centroid_stretch(1:7,1:nc,1:nc),&
         fvm%vertex_recons_weights(1:irecons_tracer-1,:,1:nc,1:nc),&
         fvm%vtx_cart(:,:,1-nhc:nc+nhc,1-nhc:nc+nhc))
  end subroutine get_fvm_recons
end module fvm_mapping
