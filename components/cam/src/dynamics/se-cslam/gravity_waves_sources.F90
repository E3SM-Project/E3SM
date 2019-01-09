module gravity_waves_sources
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: np,nlev
  use edgetype_mod,   only: EdgeBuffer_t
  use element_mod,    only: element_t
  use hybrid_mod,     only: hybrid_t
  use shr_kind_mod,   only: r8 => shr_kind_r8

  implicit none
  private
  save

  !! gravity_waves_sources created by S Santos, 10 Aug 2011
  !!
  !! gws_src_fnct starts parallel environment and computes frontogenesis
  !!   for use by WACCM (via dp_coupling)

  public  :: gws_src_fnct
  public  :: gws_init
  private :: compute_frontogenesis

  type (EdgeBuffer_t) :: edge3
  type (derivative_t)   :: deriv
  real(r8) :: psurf_ref

!----------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------

  subroutine gws_init(elem)
    use parallel_mod, only   : par
    use edge_mod, only       : initEdgeBuffer
    use hycoef, only         : hypi
    use pmgrid, only         : plev
    use thread_mod, only     : horz_num_threads
    implicit none

    ! Elem will be needed for future updates to edge code
    type(element_t), pointer :: elem(:)

    ! Set up variables similar to dyn_comp and prim_driver_mod initializations
    call initEdgeBuffer(par, edge3, elem, 3*nlev,nthreads=horz_num_threads)

    psurf_ref = hypi(plev+1)

  end subroutine gws_init

  subroutine gws_src_fnct(elem, tl, tlq, frontgf, frontga,nphys)
    use derivative_mod, only  : derivinit
    use dimensions_mod, only  : npsq, nelemd
    use dof_mod, only         : UniquePoints
    use hybrid_mod, only      : config_thread_region, get_loop_ranges
    use parallel_mod, only    : par
    use ppgrid, only          : pver
    use thread_mod, only      : horz_num_threads
    use dimensions_mod, only  : fv_nphys
    implicit none
    type (element_t), intent(inout), dimension(:) :: elem
    integer, intent(in)          :: tl, nphys, tlq
    real (kind=r8), intent(out) :: frontgf(nphys*nphys,pver,nelemd)
    real (kind=r8), intent(out) :: frontga(nphys*nphys,pver,nelemd)

    ! Local variables
    type (hybrid_t) :: hybrid
    integer :: nets, nete, ithr, ncols, ie
    real(kind=r8), allocatable  ::  frontgf_thr(:,:,:,:)
    real(kind=r8), allocatable  ::  frontga_thr(:,:,:,:)

    ! This does not need to be a thread private data-structure
    call derivinit(deriv)
    !$OMP PARALLEL NUM_THREADS(horz_num_threads),  DEFAULT(SHARED), PRIVATE(nets,nete,hybrid,ie,ncols,frontgf_thr,frontga_thr)
    hybrid = config_thread_region(par,'horizontal')
!JMD    hybrid = config_thread_region(par,'serial')
    call get_loop_ranges(hybrid,ibeg=nets,iend=nete)

    allocate(frontgf_thr(nphys,nphys,nlev,nets:nete))
    allocate(frontga_thr(nphys,nphys,nlev,nets:nete))    
    call compute_frontogenesis(frontgf_thr,frontga_thr,tl,tlq,elem,deriv,hybrid,nets,nete,nphys)
    if (fv_nphys>0) then
      do ie=nets,nete
        frontgf(:,:,ie) = RESHAPE(frontgf_thr(:,:,:,ie),(/nphys*nphys,nlev/))
        frontga(:,:,ie) = RESHAPE(frontga_thr(:,:,:,ie),(/nphys*nphys,nlev/))
      end do
    else
      do ie=nets,nete
        ncols = elem(ie)%idxP%NumUniquePts
        call UniquePoints(elem(ie)%idxP, nlev, frontgf_thr(:,:,:,ie), frontgf(1:ncols,:,ie))
        call UniquePoints(elem(ie)%idxP, nlev, frontga_thr(:,:,:,ie), frontga(1:ncols,:,ie))
      end do
    end if
    deallocate(frontga_thr)
    deallocate(frontgf_thr)
    !$OMP END PARALLEL

  end subroutine gws_src_fnct

  subroutine compute_frontogenesis(frontgf,frontga,tl,tlq,elem,ederiv,hybrid,nets,nete,nphys)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute frontogenesis function F
  !   F =  -gradth dot C
  ! with:
  !   theta  = potential temperature
  !   gradth = grad(theta)
  !   C = ( gradth dot grad ) U
  !
  ! Original by Mark Taylor, July 2011
  ! Change by Santos, 10 Aug 2011:
  ! Integrated into gravity_waves_sources module, several arguments made global
  !  to prevent repeated allocation/initialization
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use physconst,      only: cappa
    use derivative_mod, only: gradient_sphere, ugradv_sphere
    use edge_mod,       only: edgevpack, edgevunpack
    use bndry_mod,      only: bndry_exchange
    use dyn_grid,       only: hvcoord
    use dimensions_mod, only: fv_nphys,ntrac
    use dimensions_mod, only: qsize_condensate_loading_idx_gll,qsize_condensate_loading
    use fvm_mapping,    only: dyn2phys_vector,dyn2phys
    
    type(hybrid_t),     intent(in)            :: hybrid
    type(element_t),    intent(inout), target :: elem(:)
    type(derivative_t), intent(in)            :: ederiv
    integer,            intent(in)            :: nets,nete,nphys
    integer,            intent(in)            :: tl,tlq
    real(r8),           intent(out)           :: frontgf(nphys,nphys,nlev,nets:nete)
    real(r8),           intent(out)           :: frontga(nphys,nphys,nlev,nets:nete)

    ! local
    real(r8) :: area_inv(fv_nphys,fv_nphys), tmp(np,np)
    real(r8) :: uv_tmp(fv_nphys*fv_nphys,2,nlev)
    real(r8) :: frontgf_gll(np,np,nlev,nets:nete)
    real(r8) :: frontga_gll(np,np,nlev,nets:nete)
    integer  :: k,kptr,i,j,ie,component,h,nq,m_cnst
    real(r8) :: gradth(np,np,2,nlev,nets:nete) ! grad(theta)
    real(r8) :: p(np,np)                       ! pressure at mid points
    real(r8) :: pint(np,np)                    ! pressure at interface points
    real(r8) :: theta(np,np)                   ! potential temperature at mid points
    real(r8) :: C(np,np,2), sum_water(np,np)

    do ie=nets,nete
      ! pressure at model top
      pint(:,:) = hvcoord%hyai(1) 
      do k=1,nlev
        ! moist pressure at mid points
        sum_water(:,:) = 1.0_r8
        do nq=1,qsize_condensate_loading
          m_cnst = qsize_condensate_loading_idx_gll(nq)
          !
          ! make sure Q is updated
          !
          sum_water(:,:) = sum_water(:,:) + elem(ie)%state%Qdp(:,:,k,m_cnst,tlq)/elem(ie)%state%dp3d(:,:,k,tl)
        end do
        p(:,:) = pint(:,:) + 0.5_r8*sum_water(:,:)*elem(ie)%state%dp3d(:,:,k,tl)
        ! moist pressure at interface for next iteration
        pint(:,:) = pint(:,:)+elem(ie)%state%dp3d(:,:,k,tl)
        !
        theta(:,:) = elem(ie)%state%T(:,:,k,tl)*(psurf_ref / p(:,:))**cappa
        ! gradth(:,:,:,k,ie) = gradient_sphere(theta,ederiv,elem(ie)%Dinv)        
        call gradient_sphere(theta,ederiv,elem(ie)%Dinv,gradth(:,:,:,k,ie))        
        ! compute C = (grad(theta) dot grad ) u
        C(:,:,:) = ugradv_sphere(gradth(:,:,:,k,ie), elem(ie)%state%v(:,:,:,k,tl),ederiv,elem(ie))        
        ! gradth dot C
        frontgf_gll(:,:,k,ie) = -( C(:,:,1)*gradth(:,:,1,k,ie) +  C(:,:,2)*gradth(:,:,2,k,ie)  )        
        ! apply mass matrix
        gradth(:,:,1,k,ie)=gradth(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
        gradth(:,:,2,k,ie)=gradth(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
        frontgf_gll(:,:,k,ie)=frontgf_gll(:,:,k,ie)*elem(ie)%spheremp(:,:)        
      enddo
      ! pack
      call edgeVpack(edge3, frontgf_gll(:,:,:,ie),nlev,0,ie)
      call edgeVpack(edge3, gradth(:,:,:,:,ie),2*nlev,nlev,ie)
    enddo
    call bndry_exchange(hybrid,edge3,location='compute_frontogenesis')
    do ie=nets,nete
      call edgeVunpack(edge3, frontgf_gll(:,:,:,ie),nlev,0,ie)
      call edgeVunpack(edge3, gradth(:,:,:,:,ie),2*nlev,nlev,ie)
      ! apply inverse mass matrix,
      do k=1,nlev
        gradth(:,:,1,k,ie)=gradth(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
        gradth(:,:,2,k,ie)=gradth(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        frontgf_gll(:,:,k,ie)=frontgf_gll(:,:,k,ie)*elem(ie)%rspheremp(:,:)        
      end do
      if (fv_nphys>0) then
        uv_tmp(:,:,:) = dyn2phys_vector(gradth(:,:,:,:,ie),elem(ie))
        do k=1,nlev
          h=0
          do j=1,fv_nphys
            do i=1,fv_nphys
              h=h+1
              frontga(i,j,k,ie) = atan2 ( uv_tmp(h,2,k) , uv_tmp(h,1,k) + 1.e-10_r8 )
            end do
          end do
        end do
        !
        ! compute inverse physgrid area for mapping of scaler
        !
        tmp = 1.0_r8
        area_inv = dyn2phys(tmp,elem(ie)%metdet)
        area_inv = 1.0_r8/area_inv
        do k=1,nlev
          frontgf(:,:,k,ie) = dyn2phys(frontgf_gll(:,:,k,ie),elem(ie)%metdet,area_inv)
        end do        
      else
        do k=1,nlev
          frontgf(:,:,k,ie)=frontgf_gll(:,:,k,ie)
          ! Frontogenesis angle
          frontga(:,:,k,ie) = atan2 ( gradth(:,:,2,k,ie) , gradth(:,:,1,k,ie) + 1.e-10_r8 )
        end do
      end if
    enddo
  end subroutine compute_frontogenesis


end module gravity_waves_sources
