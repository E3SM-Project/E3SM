!-----------------------------------------------------------------------------!
!MODULE FVM_MOD-----------------------------------------------------CE-for FVM!
! FVM_MOD File for the fvm project in HOMME                                   !
! Author: Christoph Erath                                                     !
! Date: 25.January 2011                                                       !
! MAIN module to run fvm on HOMME                                             !
! 14.November 2011: reorganisation done                                       !
! 7.Februar 2012: cslam_run and cslam_runair                                  !
!-----------------------------------------------------------------------------!

module fvm_mod    
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use edge_mod,               only: initghostbuffer, freeghostbuffer, ghostpack, ghostunpack
  use edgetype_mod,           only: edgebuffer_t
  use bndry_mod,              only: ghost_exchange

  use element_mod,            only: element_t
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod,             only: hybrid_t
  
  implicit none
  private
  save
  
  type (EdgeBuffer_t)                         :: edgeveloc
  type (EdgeBuffer_t), public  :: ghostBufQnhc, ghostBufQ1, ghostBufFlux

  interface fill_halo_fvm
     module procedure fill_halo_fvm_noprealloc
     module procedure fill_halo_fvm_prealloc
  end interface


  public :: edgeveloc, fvm_init1,fvm_init2, fill_halo_fvm, fvm_pg_init,fvm_init3,fill_halo_and_extend_panel

contains

  subroutine fill_halo_fvm_noprealloc(elem,fvm,hybrid,nets,nete,tnp0,ndepth,kmin,kmax)
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    use dimensions_mod, only: nc, ntrac
    implicit none
    type (element_t),intent(inout)            :: elem(:)
    type (fvm_struct),intent(inout)           :: fvm(:)
    type (hybrid_t),intent(in)                :: hybrid

    type (edgeBuffer_t)                      :: cellghostbuf

    integer,intent(in)                        :: nets,nete,tnp0,ndepth,kmin,kmax
    integer                                   :: ie,i1,i2,num_levels
    !
    !

    call t_startf('FVM:initbuf')
    i1=1-ndepth
    i2=nc+ndepth
    num_levels = kmax-kmin+1
    call initghostbuffer(hybrid%par,cellghostbuf,elem,num_levels*(ntrac+1),ndepth,nc)
    call t_stopf('FVM:initbuf')
    call t_startf('FVM:pack')
    do ie=nets,nete
       call ghostpack(cellghostbuf, fvm(ie)%dp_fvm(i1:i2,i1:i2,kmin:kmax,tnp0),num_levels,      0,ie)
       call ghostpack(cellghostbuf, fvm(ie)%c(i1:i2,i1:i2,kmin:kmax,:,tnp0)   ,num_levels*ntrac,num_levels,ie)
    end do
    call t_stopf('FVM:pack')
    call t_startf('FVM:Communication')
    call ghost_exchange(hybrid,cellghostbuf)
    call t_stopf('FVM:Communication')
    !-----------------------------------------------------------------------------------!                        
    call t_startf('FVM:Unpack')
    do ie=nets,nete
       call ghostunpack(cellghostbuf, fvm(ie)%dp_fvm(i1:i2,i1:i2,kmin:kmax,tnp0),num_levels      ,0,ie)
       call ghostunpack(cellghostbuf, fvm(ie)%c(i1:i2,i1:i2,kmin:kmax,:,tnp0),   num_levels*ntrac,num_levels,ie)
    enddo
    call t_stopf('FVM:Unpack')
    call t_startf('FVM:freebuf')
    call freeghostbuffer(cellghostbuf)
    call t_stopf('FVM:freebuf')
  end subroutine fill_halo_fvm_noprealloc

subroutine fill_halo_fvm_prealloc(cellghostbuf,elem,fvm,hybrid,nets,nete,tnp0,ndepth,kmin,kmax)
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    use dimensions_mod, only: nc, ntrac
    implicit none
    type (EdgeBuffer_t), intent(inout)       :: cellghostbuf
    type (element_t),intent(inout)            :: elem(:)
    type (fvm_struct),intent(inout)           :: fvm(:)
    type (hybrid_t),intent(in)                :: hybrid


    integer,intent(in)                        :: nets,nete,tnp0,ndepth,kmin,kmax
    integer                                   :: ie,i1,i2,num_levels
    !
    !

!    call t_startf('FVM:initbuf')
    i1=1-ndepth
    i2=nc+ndepth
    num_levels = kmax-kmin+1
    call t_startf('FVM:pack')
    do ie=nets,nete
       call ghostpack(cellghostbuf, fvm(ie)%dp_fvm(i1:i2,i1:i2,kmin:kmax,tnp0),num_levels,      0,ie)
       call ghostpack(cellghostbuf, fvm(ie)%c(i1:i2,i1:i2,kmin:kmax,:,tnp0) ,num_levels*ntrac,num_levels,ie)
    end do
    call t_stopf('FVM:pack')
    call t_startf('FVM:Communication')
    call ghost_exchange(hybrid,cellghostbuf)
    call t_stopf('FVM:Communication')
    !-----------------------------------------------------------------------------------!                        
    call t_startf('FVM:Unpack')
    do ie=nets,nete
       call ghostunpack(cellghostbuf, fvm(ie)%dp_fvm(i1:i2,i1:i2,kmin:kmax,tnp0),num_levels      ,0,ie)
       call ghostunpack(cellghostbuf, fvm(ie)%c(i1:i2,i1:i2,kmin:kmax,:,tnp0), num_levels*ntrac,num_levels,ie)
    enddo
    call t_stopf('FVM:Unpack')

  end subroutine fill_halo_fvm_prealloc

 subroutine PrintArray(i1,i2,array)
   ! debug routine potentially called from any MPI rank
   integer :: i1,i2
   real(kind=r8) :: array(i1:i2,i1:i2)
   integer :: sz,i,ub

   sz = size(array,dim=1)

   if (sz == 9) then 
     do i=i2,i1,-1
        write(6,9) array(-2,i),array(-1,i), array(0,i), &
                   array( 1,i), array(2,i), array(3,i), &
                   array( 4,i), array(5,i), array(6,i)
     enddo 
   endif

 9      format('|',9(f10.1,'|'))


 end subroutine


  subroutine fill_halo_and_extend_panel(elem,fvm,fld,hybrid,nets,nete,nphys,nhcc, ndepth,numlev,num_flds,lfill_halo,lextend_panel)
    use hybrid_mod,             only: hybrid_t
    use edge_mod,               only: initghostbuffer, freeghostbuffer, ghostpack, ghostunpack

    use fvm_reconstruction_mod, only: extend_panel_interpolate
    use cam_abortutils,         only: endrun
    use dimensions_mod,         only: fv_nphys,nhr,nhr_phys,nhc,nhc_phys,ns,ns_phys,nhe_phys,nc
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL

    integer              , intent(in)    :: nets,nete,nphys,ndepth,numlev,num_flds,nhcc
    real (kind=r8)       , intent(inout) :: fld(1-nhcc:nphys+nhcc,1-nhcc:nphys+nhcc,numlev,num_flds,nets:nete)
    type (hybrid_t)      , intent(in)    :: hybrid  ! distributed parallel structure (shared)
    type (element_t)     , intent(inout) :: elem(:)
    type(fvm_struct)     , intent(in)    :: fvm(:)
    logical              , intent(in)    :: lfill_halo,lextend_panel
!    real (kind=r8)       , allocatable   :: ftmp(:,:)
!    real (kind=r8)               :: ftmp(1-nhcc:nphys+nhcc,1-nhcc:nphys+nhcc,numlev,num_flds,nets:nete)
    real (kind=r8), allocatable  :: fld_tmp(:,:)

    integer                 :: ie,k,itr,nht_phys,nh_phys
    type (edgeBuffer_t)   :: cellghostbuf
    
    if (lfill_halo) then
      !
      !*********************************************
      !
      ! halo exchange
      !
      !*********************************************
      !
      call t_startf('fill_halo_and_extend_panel initbuffer')
      call initghostbuffer(hybrid%par,cellghostbuf,elem,numlev*num_flds,nhcc,nphys)
      call t_stopf('fill_halo_and_extend_panel initbuffer')
      do ie=nets,nete
        call ghostpack(cellghostbuf, fld(:,:,:,:,ie),numlev*num_flds,0,ie)
      end do
      call ghost_exchange(hybrid,cellghostbuf)
      do ie=nets,nete
        call ghostunpack(cellghostbuf, fld(:,:,:,:,ie),numlev*num_flds,0,ie)
      end do
      call freeghostbuffer(cellghostbuf)
    end if
    if (lextend_panel) then
      !
      !*********************************************
      !
      ! extend panel
      !
      !*********************************************
      !
      if (nphys==fv_nphys) then
        if (ndepth>nhr_phys) &
             call endrun("fill_halo_and_extend_panel: ndepth>nhr_phys")
        nht_phys = nhe_phys+nhr_phys
        nh_phys  = nhr_phys
        allocate(fld_tmp(1-nht_phys:nphys+nht_phys,1-nht_phys:nphys+nht_phys))
        do ie=nets,nete
          do itr=1,num_flds
            do k=1,numlev
              call extend_panel_interpolate(fv_nphys,nhc_phys,nhr_phys,nht_phys,ns_phys,nh_phys,&
                   fld(:,:,k,itr,ie),fvm(ie)%cubeboundary,&
                   fvm(ie)%halo_interp_weight_physgrid(1:ns_phys,1-nh_phys:fv_nphys+nh_phys,1:nhr_phys,:),&
                   fvm(ie)%ibase_physgrid(1-nh_phys:fv_nphys+nh_phys,1:nhr_phys,:),&
                   fld_tmp)
              fld(1-ndepth:nphys+ndepth,1-ndepth:nphys+ndepth,k,itr,ie) = fld_tmp(1-ndepth:nphys+ndepth,1-ndepth:nphys+ndepth)
            end do
          end do
        end do
        deallocate(fld_tmp)
      else if (nphys==nc) then
        if (ndepth>nhr) &
             call endrun("fill_halo_and_extend_panel: ndepth>nhr")
        nhe_phys= 0
        nht_phys= nhe_phys+nhr
        nh_phys = nhr
        allocate(fld_tmp(1-nht_phys:nphys+nht_phys,1-nht_phys:nphys+nht_phys))
        do ie=nets,nete
          do itr=1,num_flds
            do k=1,numlev
              call extend_panel_interpolate(nc,nhc,nhr,nht_phys,ns,nh_phys,&
                   fld(:,:,k,itr,ie),fvm(ie)%cubeboundary,&
                   fvm(ie)%halo_interp_weight(1:ns,1-nh_phys:nc+nh_phys,1:nhr,:),&
                   fvm(ie)%ibase(1-nh_phys:nc+nh_phys,1:nhr,:),&
                   fld_tmp)
              fld(1-ndepth:nphys+ndepth,1-ndepth:nphys+ndepth,k,itr,ie) = fld_tmp(1-ndepth:nphys+ndepth,1-ndepth:nphys+ndepth)
            end do
          end do
        end do
        deallocate(fld_tmp)
      else
        call endrun("fill_halo_and_extend_panel: resolution not supported")
      end if
    end if
  end subroutine fill_halo_and_extend_panel

  
  ! initialize global buffers shared by all threads
  subroutine fvm_init1(par,elem)
    use parallel_mod,           only: parallel_t
    use cam_abortutils,         only: endrun
    use cam_logfile,            only: iulog
    use control_mod,            only: tracer_transport_type, rsplit
    use control_mod,            only: TRACERTRANSPORT_CONSISTENT_SE_FVM
    use fvm_control_volume_mod, only: n0_fvm, np1_fvm, fvm_supercycling
    use dimensions_mod,         only: qsize, qsize_d
    use dimensions_mod,         only: nc,nhe, nhc, nlev,ntrac, ntrac_d,ns, nhr

    type (parallel_t) :: par
    type (element_t),intent(inout)            :: elem(:)
    !
    ! initialize fvm time-levels
    !
    n0_fvm  = 1
    np1_fvm = 2
    !
    if (ntrac>0) then
      if (par%masterproc) then 
        write(iulog,*) "                                           "
        write(iulog,*) "|-----------------------------------------|"
        write(iulog,*) "| FVM tracer transport scheme information |"
        write(iulog,*) "|-----------------------------------------|"
        write(iulog,*) "                                           "
      end if
      if (tracer_transport_type == TRACERTRANSPORT_CONSISTENT_SE_FVM) then
        if (par%masterproc) then 
          write(iulog,*) "Running consistent SE-CSLAM, Lauritzen et al. (2017, MWR)."
          write(iulog,*) "CSLAM = Conservative Semi-LAgrangian Multi-tracer scheme"
          write(iulog,*) "Lauritzen et al., (2010), J. Comput. Phys."
          write(iulog,*) "  "
        end if
      end if

      if (ntrac>ntrac_d) &
           call endrun("PARAMETER ERROR for fvm: ntrac > ntrac_d")

            if (qsize>0.and.mod(rsplit,fvm_supercycling).ne.0) then
        if (par%masterproc) then
          write(iulog,*)'cannot supercycle fvm tracers with respect to se tracers'
          write(iulog,*)'with this choice of rsplit =',rsplit
          write(iulog,*)'rsplit must be a multiple of fvm_supercycling=',fvm_supercycling
          call endrun("PARAMETER ERROR for fvm: mod(rsplit,fvm_supercycling<>0")
        end if
      endif
      
      
      if (par%masterproc) then 
        write(iulog,*) "                                            "
        write(iulog,*) "Done Tracer transport scheme information    "
        write(iulog,*) "                                            "
      end if
    end if

      
    if (par%masterproc) write(iulog,*) "fvm resolution is nc*nc in each element: nc = ",nc
    if (par%masterproc) write(iulog,*)'ntrac,ntrac_d=',ntrac,ntrac_d      
    if (par%masterproc) write(iulog,*)'qsize,qsize_d=',qsize,qsize_d

    
    if (nc<3) then
      if (par%masterproc) then 
        write(iulog,*) "NUMBER OF CELLS ERROR for fvm: Number of cells parameter"
        write(iulog,*) "parameter nc at least 3 (nc>=3), nc*nc cells per element. This is"
        write(iulog,*) "needed for the cubic reconstruction, which is only implemented yet! STOP"
      endif
      call endrun("stopping")
    end if
    
    if (par%masterproc) then
      write(iulog,*) "  "
      if (ns==1) then
        write(iulog,*) "ns==1: using no interpolation for mapping cell averages values across edges"
        write(iulog,*) "Note: this is not a recommended setting - large errors at panel edges!"
      else if (ns==2) then
        write(iulog,*) "ns==2: using linear interpolation for mapping cell averages values across edges"
        write(iulog,*) "Note that ns=4 is default CSLAM setting used in Lauritzen et al. (2010)"
        write(iulog,*) "so this option is slightly less accurate (but the stencil is smaller near panel edges!)"
        
      else if (ns==3) then
        write(iulog,*) "ns==3: using quadratic interpolation for mapping cell averages values across edges"
        write(iulog,*) "Note that ns=4 is default CSLAM setting used in Lauritzen et al. (2010)"
        write(iulog,*) "so this option is slightly less accurate (but the stencil is smaller near panel edges!)"
      else if (ns==4) then
        write(iulog,*) "ns==4: using cubic interpolation for mapping cell averages values across edges"
        write(iulog,*) "This is default CSLAM setting used in Lauritzen et al. (2010)"
      else 
        write(iulog,*) "Not a tested value for ns but it should work! You choose ns = ",ns
      end if
      
      !       if (ns.NE.3) then
      !         write(*,*) "In fvm_reconstruction_mod function matmul_w has been hard-coded for ns=3 for performance"
      !         write(*,*) "Revert to general code - outcommented above"
      !         call endrun("stopping")
      !       end if
    end if
    
    if (MOD(ns,2)==0.and.nhr+(nhe-1)+ns/2>nc+nc) then
      write(iulog,*) "to run this combination of ns and nhr you need to increase nc to ",nhr+ns/2+nhe-1
      write(iulog,*) "You choose (ns,nhr,nc,nhe)=",ns,nhr,nc,nhe
      call endrun("stopping")
    end if
    if (MOD(ns,2)==1.and.nhr+(ns-1)/2+(nhe-1)>nc+nc) then
      write(iulog,*) "to run this combination of ns and nhr you need to increase nc to ",nhr+(ns-1)/2+nhe-1
      write(iulog,*) "You choose (ns,nhr,nc,nhe)=",ns,nhr,nc,nhe
      call endrun("stopping")
    end if
    
    if (nc==3.and.ns.ne.3) then
      if (par%masterproc) then
        write(iulog,*) "Recommended setting for nc=3 is ns=3 (linear interpolation in halo)"
        write(iulog,*) "You choose ns=",ns
        write(iulog,*) "Goto dimensions_mod to change value of ns"
        write(iulog,*) "or outcomment call haltmop below (i.e. you know what you are doing!)"
      endif
      call endrun("stopping")
    end if
    
    if (nc==4.and.ns.ne.4) then
      if (par%masterproc) then
        write(iulog,*) "Recommended setting for nc=4 is ns=4 (cubic interpolation in halo)"
        write(iulog,*) "You choose ns=",ns
        write(iulog,*) "Goto dimensions_mod to change value of ns"
        write(iulog,*) "or outcomment call haltmop below (i.e. you know what you are doing!)"
      endif
      call endrun("stopping")
    end if
    
    if (nhe .ne. 1) then
      if (par%masterproc) then
        write(iulog,*) "PARAMETER ERROR for fvm: Number of halo zone for the extended"
        write(iulog,*) "element nhe has to be 1, only this is available now! STOP!"
      endif
      call endrun("stopping")
    end if
    
  end subroutine fvm_init1
  
  
  
  
  
  ! initialization that can be done in threaded regions
  subroutine fvm_init2(elem,fvm,hybrid,nets,nete)
    use fvm_control_volume_mod, only: fvm_mesh,fvm_set_cubeboundary,n0_fvm,np1_fvm
    use bndry_mod,              only: compute_ghost_corner_orientation
    use dimensions_mod,         only: nlev, nc, nhc, nhe, ntrac, ntrac_d, np
    use hycoef,                 only: hyai, hybi, ps0
    use derivative_mod,         only: subcell_integration
    
    type (fvm_struct) :: fvm(:)
    type (element_t)  :: elem(:)
    type (hybrid_t)   :: hybrid
    integer           :: ie,nets,nete,k
    real(kind=r8)     :: one(np,np)

    one = 1.0_r8
    do ie=nets,nete
      do k = 1, nlev
        fvm(ie)%dp_ref(k)         = ( hyai(k+1) - hyai(k) )*ps0 + ( hybi(k+1) - hybi(k) )*ps0
        fvm(ie)%dp_ref_inverse(k) = 1.0_r8/fvm(ie)%dp_ref(k)
      end do
    end do

    n0_fvm  = 1 !in case no cslam but physgrid
    np1_fvm = 2
    
    call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
    ! run some tests:
    !    call test_ghost(hybrid,elem,nets,nete)
    
    do ie=nets,nete
      call fvm_set_cubeboundary(elem(ie),fvm(ie))
      call fvm_mesh(elem(ie),fvm(ie))
      fvm(ie)%inv_area_sphere    = 1.0_r8/fvm(ie)%area_sphere
      !
      ! compute CSLAM areas consistent with SE area (at 1 degree they can be up to 
      ! 1E-6 different than the correct spherical areas used in CSLAM)
      !
      call subcell_integration(one, np, nc, elem(ie)%metdet,fvm(ie)%inv_se_area_sphere)
      fvm(ie)%inv_se_area_sphere = 1.0_r8/fvm(ie)%inv_se_area_sphere

      fvm(ie)%fc(:,:,:,:) = 0.0_r8
      fvm(ie)%fm(:,:,:,:) = 0.0_r8
      fvm(ie)%ft(:,:,:  ) = 0.0_r8
    enddo
    ! Need to allocate ghostBufQnhc after compute_ghost_corner_orientation because it 
    ! changes the values for reverse
    call initghostbuffer(hybrid%par,ghostBufQnhc,elem,nlev*(ntrac+1),nhc,nc)
    call initghostbuffer(hybrid%par,ghostBufQ1,elem,nlev*(ntrac+1),1,nc)
    call initghostbuffer(hybrid%par,ghostBufFlux,elem,4*nlev,nhe,nc)

  end subroutine fvm_init2

  
  subroutine fvm_init3(elem,fvm,hybrid,nets,nete,irecons)
    use control_mod     ,       only: neast, nwest, seast, swest
    use fvm_analytic_mod,       only: compute_reconstruct_matrix
    use dimensions_mod  ,       only: fv_nphys
    use dimensions_mod,         only: nlev, nc, nhe, nlev, ntrac, ntrac_d,nhc
    use coordinate_systems_mod, only: cartesian2D_t,cartesian3D_t
    use fvm_control_volume_mod, only: n0_fvm
    use coordinate_systems_mod, only: cubedsphere2cart, cart2cubedsphere
    implicit none
    type (element_t) ,intent(inout)  :: elem(:)
    type (fvm_struct),intent(inout)  :: fvm(:) 
    type (hybrid_t)  ,intent(in)     :: hybrid                      
    integer          ,intent(in)     :: nets,nete,irecons
    !
    type (edgeBuffer_t)     :: cellghostbuf
    integer                 :: ie, ixy, ivertex, i, j,istart,itot,ishft,imin,imax
    integer, dimension(2,4) :: unit_vec
    integer                 :: rot90_matrix(2,2), iside

    type (cartesian2D_t)                :: tmpgnom
    type (cartesian2D_t)                :: gnom
    type(cartesian3D_t)                 :: tmpcart3d

    if (ntrac>0.and.nc.ne.fv_nphys) then
      !
      ! fill the fvm halo for mapping in d_p_coupling if
      ! physics grid resolution is different than fvm resolution
      !
      call fill_halo_fvm(elem,fvm,hybrid,nets,nete,n0_fvm,nhc,1,nlev)
    end if


    imin=1-nhc
    imax=nc+nhc
    !
    ! fill halo start
    !
    itot=9+irecons-1+2
    call initghostbuffer(hybrid%par,cellghostbuf,elem,itot,nhc,nc)
    do ie=nets,nete
      istart = 0
      call ghostpack(cellghostbuf, fvm(ie)%norm_elem_coord(1,:,:),1,istart,ie)
      istart = istart+1
      call ghostpack(cellghostbuf, fvm(ie)%norm_elem_coord(2,:,:),1,istart,ie)
      istart = istart+1        
      do ixy=1,2
        do ivertex=1,4
          call ghostpack(cellghostbuf, fvm(ie)%vtx_cart(ivertex,ixy,:,:) ,1,istart,ie)
          istart = istart+1
        end do
      end do
      call ghostpack(cellghostbuf, fvm(ie)%flux_orient(1,:,:) ,1,istart,ie)
      do ixy=1,irecons-1
        istart=istart+1
        call ghostpack(cellghostbuf, fvm(ie)%spherecentroid(ixy,:,:) ,1,istart,ie)
      end do
    end do
    call ghost_exchange(hybrid,cellghostbuf)
    do ie=nets,nete
      istart = 0
      call ghostunpack(cellghostbuf, fvm(ie)%norm_elem_coord(1,:,:),1,istart,ie)
      istart = istart+1
      call ghostunpack(cellghostbuf, fvm(ie)%norm_elem_coord(2,:,:),1,istart,ie)
      istart = istart+1        
      do ixy=1,2
        do ivertex=1,4
          call ghostunpack(cellghostbuf, fvm(ie)%vtx_cart(ivertex,ixy,:,:) ,1,istart,ie)
          istart = istart+1
        end do
      end do
      call ghostunpack(cellghostbuf, fvm(ie)%flux_orient(1,:,:) ,1,istart,ie)
      do ixy=1,irecons-1
        istart=istart+1
        call ghostunpack(cellghostbuf, fvm(ie)%spherecentroid(ixy,:,:) ,1,istart,ie)
      end do
    enddo
    call freeghostbuffer(cellghostbuf)    
    !
    ! indicator for non-existing cells 
    ! set vtx_cart to corner value in non-existent cells
    !
    do ie=nets,nete
       if (fvm(ie)%cubeboundary==nwest) then
         fvm(ie)%flux_orient     (:  ,1-nhc      :0     ,nc      +1 :nc      +nhc      ) = -1
         fvm(ie)%spherecentroid  (:,    1-nhc      :0     ,nc      +1 :nc      +nhc    ) = -1e5_r8
         fvm(ie)%vtx_cart(:,1,1-nhc:0     ,nc+1 :nc+nhc) = fvm(ie)%vtx_cart(4,1,1,nc)
         fvm(ie)%vtx_cart(:,2,1-nhc:0     ,nc+1 :nc+nhc) = fvm(ie)%vtx_cart(4,2,1,nc)
       else if (fvm(ie)%cubeboundary==swest) then
         fvm(ie)%flux_orient     (:,1-nhc      :0     ,1-nhc      :0   ) = -1
         fvm(ie)%spherecentroid  (:,1-nhc      :0     ,1-nhc      :0   ) = -1e5_r8
         fvm(ie)%vtx_cart(:,1,1-nhc:0     ,1-nhc:0     ) = fvm(ie)%vtx_cart(1,1,1,1)
         fvm(ie)%vtx_cart(:,2,1-nhc:0     ,1-nhc:0     ) = fvm(ie)%vtx_cart(1,2,1,1)
       else if (fvm(ie)%cubeboundary==neast) then
         fvm(ie)%flux_orient     (:,nc      +1 :nc      +nhc      ,nc      +1 :nc      +nhc    ) = -1
         fvm(ie)%spherecentroid  (:,nc      +1 :nc      +nhc      ,nc      +1 :nc      +nhc    ) = -1e5_r8
         fvm(ie)%vtx_cart(:,1,nc+1 :nc+nhc,nc+1 :nc+nhc) = fvm(ie)%vtx_cart(3,1,nc,nc)
         fvm(ie)%vtx_cart(:,2,nc+1 :nc+nhc,nc+1 :nc+nhc) = fvm(ie)%vtx_cart(3,2,nc,nc)
       else if (fvm(ie)%cubeboundary==seast) then
         fvm(ie)%flux_orient     (:,nc      +1 :nc      +nhc      ,1-nhc      :0   ) = -1
         fvm(ie)%spherecentroid  (:,nc      +1 :nc      +nhc      ,1-nhc      :0   ) = -1e5_r8
         fvm(ie)%vtx_cart(:,1,nc+1 :nc+nhc,1-nhc:0     ) = fvm(ie)%vtx_cart(2,1,nc,1)
         fvm(ie)%vtx_cart(:,2,nc+1 :nc+nhc,1-nhc:0     ) = fvm(ie)%vtx_cart(2,2,nc,1)
       end if
     end do
     
     !
     ! set vectors for perpendicular flux vector
     !
     rot90_matrix(1,1) = 0; rot90_matrix(2,1) =  1 !counter-clockwise rotation matrix
     rot90_matrix(1,2) =-1; rot90_matrix(2,2) =  0 !counter-clockwise rotation matrix 
     
     iside = 1
     unit_vec(1,iside) = 0 !x-component of displacement vector for side 1
     unit_vec(2,iside) = 1 !y-component of displacement vector for side 1
     
     do iside=2,4
       unit_vec(:,iside) = MATMUL(rot90_matrix(:,:),unit_vec(:,iside-1))
     end do
     
     !
     ! fill halo done
     !
     !-------------------------------
     
     do ie=nets,nete
       fvm(ie)%displ_max = 0.0_r8
       do j=imin,imax
         do i=imin,imax
           !
           ! rotate gnomonic coordinate vector
           !
           !           fvm(ie)%norm_elem_coord(:,i,j) = MATMUL(fvm(ie)%rot_matrix(:,:,i,j),fvm(ie)%norm_elem_coord(:,i,j))
           !    
           ishft = NINT(fvm(ie)%flux_orient(2,i,j))
           do ixy=1,2
             !
             ! rotate coordinates if needed through permutation
             !
             fvm(ie)%vtx_cart(1:4,ixy,i,j) = cshift(fvm(ie)%vtx_cart(1:4,ixy,i,j),shift=ishft)
             fvm(ie)%flux_vec        (ixy,i,j,1:4) = cshift(unit_vec                (ixy,1:4    ),shift=ishft)
             !
             ! set flux vector to zero in non-existent cells (corner halo) 
             !
             fvm(ie)%flux_vec        (ixy,i,j,1:4) = fvm(ie)%ifct(i,j)*fvm(ie)%flux_vec(ixy,i,j,1:4)
             
             iside=1
             fvm(ie)%displ_max(i,j,iside) = fvm(ie)%displ_max(i,j,iside)+&
                  ABS(fvm(ie)%vtx_cart(4,ixy,i,j)-fvm(ie)%vtx_cart(1,ixy,i,j))
             iside=2
             fvm(ie)%displ_max(i,j,iside) = fvm(ie)%displ_max(i,j,iside)+&
                  ABS(fvm(ie)%vtx_cart(1,ixy,i,j)-fvm(ie)%vtx_cart(2,ixy,i,j))
             iside=3
             fvm(ie)%displ_max(i,j,iside) = fvm(ie)%displ_max(i,j,iside)+&
                  ABS(fvm(ie)%vtx_cart(2,ixy,i,j)-fvm(ie)%vtx_cart(3,ixy,i,j))
             iside=4
             fvm(ie)%displ_max(i,j,iside) = fvm(ie)%displ_max(i,j,iside)+&
                  ABS(fvm(ie)%vtx_cart(2,ixy,i,j)-fvm(ie)%vtx_cart(1,ixy,i,j))
           end do
         end do
       end do
     end do
     !
     ! pre-compute derived metric terms used for integration, polynomial
     ! evaluation at fvm cell vertices, etc.
     !    
     do ie=nets,nete
       call compute_reconstruct_matrix(nc,nhe,nhc,irecons,fvm(ie)%dalpha,fvm(ie)%dbeta,&
           fvm(ie)%spherecentroid,fvm(ie)%vtx_cart,fvm(ie)%centroid_stretch,&
           fvm(ie)%vertex_recons_weights,fvm(ie)%recons_metrics,fvm(ie)%recons_metrics_integral)
     end do
      !
      ! create a normalized element coordinate system with a halo
      !    
     do ie=nets,nete
       do j=1-nhc,nc+nhc
         do i=1-nhc,nc+nhc
           !
           ! only compute for physically existent cells
           !
           if (fvm(ie)%ifct(i,j)>0) then
             gnom%x = fvm(ie)%norm_elem_coord(1,i,j)
             gnom%y = fvm(ie)%norm_elem_coord(2,i,j)
             !
             ! coordinate transform only necessary for points on another panel
             !
             if (NINT(fvm(ie)%flux_orient(1,1,1)).NE.NINT(fvm(ie)%flux_orient(1,i,j))) then
               tmpcart3d=cubedsphere2cart(gnom,NINT(fvm(ie)%flux_orient(1,i,j)))
               tmpgnom=cart2cubedsphere(tmpcart3d,NINT(fvm(ie)%flux_orient(1,1,1)))
             else
               tmpgnom%x = fvm(ie)%norm_elem_coord(1,i,j)
               tmpgnom%y = fvm(ie)%norm_elem_coord(2,i,j)
             end if
             !
             ! convert to element normalized coordinates
             !
             fvm(ie)%norm_elem_coord(1,i,j) =(tmpgnom%x-elem(ie)%corners(1)%x)/&
                  (0.5_r8*dble(nc)*fvm(ie)%dalpha)-1.0_r8
             fvm(ie)%norm_elem_coord(2,i,j) =(tmpgnom%y-elem(ie)%corners(1)%y)/&
                  (0.5_r8*dble(nc)*fvm(ie)%dalpha)-1.0_r8
           else
             fvm(ie)%norm_elem_coord(1,i,j) = 1D9
             fvm(ie)%norm_elem_coord(2,i,j) = 1D9
           end if
         end do
       end do
     end do

   end subroutine fvm_init3
  

  subroutine fvm_pg_init(elem, fvm, hybrid, nets, nete,irecons)
    use coordinate_systems_mod, only : cartesian2D_t,cartesian3D_t
    use control_mod, only : neast, nwest, seast, swest
    use coordinate_systems_mod, only : cubedsphere2cart, cart2cubedsphere
    use dimensions_mod, only: fv_nphys, nhe_phys,nhc_phys
    use dimensions_mod, only: ntrac_d
    use cube_mod       ,only: dmap
    use control_mod    ,only: cubed_sphere_map
    use fvm_analytic_mod, only: compute_reconstruct_matrix

    type (element_t) , intent(in)    :: elem(:)
    type (fvm_struct), intent(inout) :: fvm(:)
    type (hybrid_t)  , intent(in)    :: hybrid

    type (cartesian2D_t)                :: gnom
    type(cartesian3D_t)                 :: tmpcart3d
    type (cartesian2D_t)                :: tmpgnom


    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete,irecons  ! ending thread element number   (private)

    ! ==================================
    ! Local variables
    ! ==================================

    integer                 :: ie, ixy, ivertex, i, j,istart,itot,ishft,imin,imax
    integer, dimension(2,4) :: unit_vec
    integer                 :: rot90_matrix(2,2), iside

    type (edgeBuffer_t)                     :: cellghostbuf

    ! D is derivative of gnomonic mapping
    real (kind=r8)              :: D(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,2,2)
    real (kind=r8)              :: detD,x1,x2

    if (fv_nphys>0) then
      !
      ! do the same as fvm_init3 for the metric terms of physgrid
      !
      imin=1-nhc_phys
      imax=fv_nphys+nhc_phys
      !
      ! fill halo start
      !
      itot=9+irecons-1+2
      call initghostbuffer(hybrid%par,cellghostbuf,elem,itot,nhc_phys,fv_nphys)
      do ie=nets,nete
        istart = 0
        call ghostpack(cellghostbuf, fvm(ie)%norm_elem_coord_physgrid(1,:,:),1,istart,ie)
        istart = istart+1
        call ghostpack(cellghostbuf, fvm(ie)%norm_elem_coord_physgrid(2,:,:),1,istart,ie)
        istart = istart+1
        do ixy=1,2
          do ivertex=1,4
            call ghostpack(cellghostbuf, fvm(ie)%vtx_cart_physgrid(ivertex,ixy,:,:) ,1,istart,ie)
            istart = istart+1
          end do
        end do
        call ghostpack(cellghostbuf, fvm(ie)%flux_orient_physgrid(1,:,:) ,1,istart,ie)
        do ixy=1,irecons-1
          istart=istart+1
          call ghostpack(cellghostbuf, fvm(ie)%spherecentroid_physgrid(ixy,:,:) ,1,istart,ie)
        end do
      end do
      call ghost_exchange(hybrid,cellghostbuf)
      do ie=nets,nete
        istart = 0
        call ghostunpack(cellghostbuf, fvm(ie)%norm_elem_coord_physgrid(1,:,:),1,istart,ie)
        istart = istart+1
        call ghostunpack(cellghostbuf, fvm(ie)%norm_elem_coord_physgrid(2,:,:),1,istart,ie)
        istart = istart+1
        do ixy=1,2
          do ivertex=1,4
            call ghostunpack(cellghostbuf, fvm(ie)%vtx_cart_physgrid(ivertex,ixy,:,:) ,1,istart,ie)
            istart = istart+1
          end do
        end do
        call ghostunpack(cellghostbuf, fvm(ie)%flux_orient_physgrid(1,:,:) ,1,istart,ie)
        do ixy=1,irecons-1
          istart=istart+1
          call ghostunpack(cellghostbuf, fvm(ie)%spherecentroid_physgrid(ixy,:,:) ,1,istart,ie)
        end do
      enddo
      call freeghostbuffer(cellghostbuf)    
      !
      ! indicator for non-existing cells 
      ! set vtx_cart to corner value in non-existent cells
      !
      do ie=nets,nete
        if (fvm(ie)%cubeboundary==nwest) then
          fvm(ie)%flux_orient_physgrid   (:  ,1-nhc_phys      :0     ,fv_nphys      +1 :fv_nphys      +nhc_phys    ) = -1
          fvm(ie)%spherecentroid_physgrid(:,  1-nhc_phys      :0     ,fv_nphys      +1 :fv_nphys      +nhc_phys    ) = -1e5_r8
          fvm(ie)%vtx_cart_physgrid(:,1,1-nhc_phys:0     ,fv_nphys+1 :fv_nphys+nhc_phys) = &
               fvm(ie)%vtx_cart_physgrid(4,1,1,fv_nphys)
          fvm(ie)%vtx_cart_physgrid(:,2,1-nhc_phys:0     ,fv_nphys+1 :fv_nphys+nhc_phys) = &
               fvm(ie)%vtx_cart_physgrid(4,2,1,fv_nphys)
        else if (fvm(ie)%cubeboundary==swest) then
          fvm(ie)%flux_orient_physgrid   (:,1-nhc_phys      :0     ,1-nhc_phys      :0   ) = -1
          fvm(ie)%spherecentroid_physgrid(:,1-nhc_phys      :0     ,1-nhc_phys      :0   ) = -1e5_r8
          fvm(ie)%vtx_cart_physgrid(:,1,1-nhc_phys:0     ,1-nhc_phys:0     ) = fvm(ie)%vtx_cart_physgrid(1,1,1,1)
          fvm(ie)%vtx_cart_physgrid(:,2,1-nhc_phys:0     ,1-nhc_phys:0     ) = fvm(ie)%vtx_cart_physgrid(1,2,1,1)
        else if (fvm(ie)%cubeboundary==neast) then
          fvm(ie)%flux_orient_physgrid   (:,fv_nphys      +1 :fv_nphys      +nhc_phys      , &
               fv_nphys      +1 :fv_nphys      +nhc_phys      ) = -1
          fvm(ie)%spherecentroid_physgrid(:,fv_nphys      +1 :fv_nphys      +nhc_phys      , &
               fv_nphys      +1 :fv_nphys      +nhc_phys      ) = -1e5_r8
          fvm(ie)%vtx_cart_physgrid(:,1,fv_nphys+1 :fv_nphys+nhc_phys,fv_nphys+1 :fv_nphys+nhc_phys) = &
               fvm(ie)%vtx_cart_physgrid(3,1,fv_nphys,fv_nphys)
          fvm(ie)%vtx_cart_physgrid(:,2,fv_nphys+1 :fv_nphys+nhc_phys,fv_nphys+1 :fv_nphys+nhc_phys) = &
               fvm(ie)%vtx_cart_physgrid(3,2,fv_nphys,fv_nphys)
        else if (fvm(ie)%cubeboundary==seast) then
          fvm(ie)%flux_orient_physgrid   (:,fv_nphys      +1 :fv_nphys      +nhc_phys      ,1-nhc_phys      :0   ) = -1
          fvm(ie)%spherecentroid_physgrid(:,fv_nphys      +1 :fv_nphys      +nhc_phys      ,1-nhc_phys      :0   ) = -1e5_r8
          fvm(ie)%vtx_cart_physgrid(:,1,fv_nphys+1 :fv_nphys+nhc_phys,1-nhc_phys:0     ) = &
               fvm(ie)%vtx_cart_physgrid(2,1,fv_nphys,1)
          fvm(ie)%vtx_cart_physgrid(:,2,fv_nphys+1 :fv_nphys+nhc_phys,1-nhc_phys:0     ) = &
               fvm(ie)%vtx_cart_physgrid(2,2,fv_nphys,1)
        end if
      end do
      
      !
      ! set vectors for perpendicular flux vector
      !
      rot90_matrix(1,1) = 0; rot90_matrix(2,1) =  1 !counter-clockwise rotation matrix
      rot90_matrix(1,2) =-1; rot90_matrix(2,2) =  0 !counter-clockwise rotation matrix 
      
      iside = 1
      unit_vec(1,iside) = 0 !x-component of displacement vector for side 1
      unit_vec(2,iside) = 1 !y-component of displacement vector for side 1
      
      do iside=2,4
        unit_vec(:,iside) = MATMUL(rot90_matrix(:,:),unit_vec(:,iside-1))
      end do
      
      !
      ! fill halo done
      !
      !-------------------------------
      
      do ie=nets,nete
        do j=imin,imax
          do i=imin,imax
            !
            ! rotate gnomonic coordinate vector
            !    
            ishft = NINT(fvm(ie)%flux_orient_physgrid(2,i,j))
            do ixy=1,2
              !
              ! rotate coordinates if needed through permutation
              !
              fvm(ie)%vtx_cart_physgrid(1:4,ixy,i,j) = cshift(fvm(ie)%vtx_cart_physgrid(1:4,ixy,i,j),shift=ishft)
            end do
          end do
        end do
      end do
      !
      ! pre-compute derived metric terms used for integration, polynomial
      ! evaluation at fvm cell vertices, etc.
      !    
      do ie=nets,nete
        call compute_reconstruct_matrix(fv_nphys,nhe_phys,nhc_phys,irecons,fvm(ie)%dalpha_physgrid,fvm(ie)%dbeta_physgrid,&
             fvm(ie)%spherecentroid_physgrid,fvm(ie)%vtx_cart_physgrid,fvm(ie)%centroid_stretch_physgrid,&
             fvm(ie)%vertex_recons_weights_physgrid,fvm(ie)%recons_metrics_physgrid,fvm(ie)%recons_metrics_integral_physgrid)
      end do      
      !
      ! code specific for physgrid
      !
      !
      ! create a normalized element coordinate system with a halo
      !    
      do ie=nets,nete
        do j=1-nhc_phys,fv_nphys+nhc_phys
          do i=1-nhc_phys,fv_nphys+nhc_phys
            !
            ! only compute for physically existent cells
            !
            if (fvm(ie)%ifct_physgrid(i,j)>0) then
              gnom%x = fvm(ie)%norm_elem_coord_physgrid(1,i,j)
              gnom%y = fvm(ie)%norm_elem_coord_physgrid(2,i,j)
              !
              ! coordinate transform only necessary for points on another panel
              !
              if (NINT(fvm(ie)%flux_orient_physgrid(1,1,1)).NE.NINT(fvm(ie)%flux_orient_physgrid(1,i,j))) then
                tmpcart3d=cubedsphere2cart(gnom,NINT(fvm(ie)%flux_orient_physgrid(1,i,j)))
                tmpgnom=cart2cubedsphere(tmpcart3d,NINT(fvm(ie)%flux_orient_physgrid(1,1,1)))
              else
                tmpgnom%x = fvm(ie)%norm_elem_coord_physgrid(1,i,j)
                tmpgnom%y = fvm(ie)%norm_elem_coord_physgrid(2,i,j)
              end if
              !
              ! convert to element normalized coordinates
              !
              fvm(ie)%norm_elem_coord_physgrid(1,i,j) =(tmpgnom%x-elem(ie)%corners(1)%x)/&
                   (0.5_r8*dble(fv_nphys)*fvm(ie)%dalpha_physgrid)-1.0_r8
              fvm(ie)%norm_elem_coord_physgrid(2,i,j) =(tmpgnom%y-elem(ie)%corners(1)%y)/&
                   (0.5_r8*dble(fv_nphys)*fvm(ie)%dalpha_physgrid)-1.0_r8
            else
              fvm(ie)%norm_elem_coord_physgrid(1,i,j) = 1D9
              fvm(ie)%norm_elem_coord_physgrid(2,i,j) = 1D9
            end if
          end do
        end do
      end do
      !
      ! compute Dinv
      !
      do ie=nets,nete
        do j=1-nhc_phys,fv_nphys+nhc_phys
          do i=1-nhc_phys,fv_nphys+nhc_phys
            x1 = fvm(ie)%norm_elem_coord_physgrid(1,i,j)
            x2 = fvm(ie)%norm_elem_coord_physgrid(2,i,j)
            call Dmap(D(i,j,:,:),x1,x2,elem(ie)%corners3D,cubed_sphere_map,elem(ie)%corners,elem(ie)%u2qmap,elem(ie)%facenum)
            detD = D(i,j,1,1)*D(i,j,2,2) - D(i,j,1,2)*D(i,j,2,1)      
            
            fvm(ie)%Dinv_physgrid(i,j,1,1) =  D(i,j,2,2)/detD
            fvm(ie)%Dinv_physgrid(i,j,1,2) = -D(i,j,1,2)/detD
            fvm(ie)%Dinv_physgrid(i,j,2,1) = -D(i,j,2,1)/detD
            fvm(ie)%Dinv_physgrid(i,j,2,2) =  D(i,j,1,1)/detD
          end do
        end do
      end do            
    end if

  end subroutine fvm_pg_init


end module fvm_mod
