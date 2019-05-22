module gravity_waves_sources
  use derivative_mod, only : derivative_t
  use dimensions_mod, only : np,nlev
  use edgetype_mod, only       : EdgeBuffer_t
  use element_mod, only    : element_t
  use hybrid_mod, only     : hybrid_t
  use kinds, only          : real_kind
  use shr_kind_mod, only   : r8 => shr_kind_r8
  use thread_mod, only   : hthreads

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
  type (derivative_t), allocatable   :: deriv(:)
  real(r8) :: psurf_ref

CONTAINS
  !-------------------------------------------------------------------------------------------------
  subroutine gws_init(elem)
    use edge_mod, only       : initEdgeBuffer
    use hycoef, only         : hypi
    use pmgrid, only         : plev
    use parallel_mod, only    : par
    implicit none
    type (element_t), intent(inout), dimension(:) :: elem
    !---------------------------------------------------------------------------
    ! Set up variables similar to dyn_comp and prim_driver_mod initializations
    call initEdgeBuffer(par,edge3,elem,3*nlev)
    allocate(deriv(0:hthreads-1))
    
    psurf_ref = hypi(plev+1)

  end subroutine gws_init
  !-------------------------------------------------------------------------------------------------
  subroutine gws_src_fnct(elem, tl, nphys, frontgf, frontga)
    use derivative_mod, only  : derivinit
    use dimensions_mod, only  : npsq, nelemd
    use dof_mod, only         : UniquePoints
    use dyn_comp, only        : dom_mt
    use hybrid_mod, only      : hybrid_create
    use parallel_mod, only    : par
    use ppgrid, only          : pver
    use thread_mod, only      : omp_get_thread_num
    use dyn_grid, only        : fv_nphys
    implicit none
    type (element_t),      intent(inout), dimension(:) :: elem
    integer,               intent(in   ) :: tl
    integer,               intent(in   ) :: nphys
    real (kind=real_kind), intent(out  ) :: frontgf(nphys*nphys,pver,nelemd)
    real (kind=real_kind), intent(out  ) :: frontga(nphys*nphys,pver,nelemd)
    
    ! Local variables
    type (hybrid_t) :: hybrid
    integer :: nets, nete, ithr, ncols, ie, k
    real(kind=real_kind), allocatable  ::  frontgf_thr(:,:,:,:)
    real(kind=real_kind), allocatable  ::  frontga_thr(:,:,:,:)
    !---------------------------------------------------------------------------
    !$OMP PARALLEL NUM_THREADS(hthreads), DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid,ie,ncols,frontgf_thr,frontga_thr)
    ithr = omp_get_thread_num()
    nets = dom_mt(ithr)%start
    nete = dom_mt(ithr)%end
    hybrid = hybrid_create(par,ithr,hthreads)
    call derivinit(deriv(hybrid%ithr))
    allocate(frontgf_thr(nphys,nphys,nlev,nets:nete))
    allocate(frontga_thr(nphys,nphys,nlev,nets:nete))
    call compute_frontogenesis(frontgf_thr,frontga_thr,tl,elem,deriv(hybrid%ithr),hybrid,nets,nete,nphys)
    if (fv_nphys>0) then
      do ie = nets,nete
        do k = 1,nlev
          frontgf(:,k,ie) = RESHAPE( frontgf_thr(:,:,k,ie), (/nphys*nphys/))
          frontga(:,k,ie) = RESHAPE( frontga_thr(:,:,k,ie), (/nphys*nphys/))
        end do
      end do
    else
      do ie = nets,nete
        ncols = elem(ie)%idxP%NumUniquePts
        call UniquePoints(elem(ie)%idxP, nlev, frontgf_thr(:,:,:,ie), frontgf(1:ncols,:,ie))
        call UniquePoints(elem(ie)%idxP, nlev, frontga_thr(:,:,:,ie), frontga(1:ncols,:,ie))
      end do
    end if ! fv_nphys>0
    deallocate(frontga_thr)
    deallocate(frontgf_thr)
    !$OMP END PARALLEL

  end subroutine gws_src_fnct
  !-------------------------------------------------------------------------------------------------
  subroutine compute_frontogenesis(frontgf,frontga,tl,elem,ederiv,hybrid,nets,nete,nphys)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  ! Change by Aaron Donahue, April 2017:
  !   Fixed bug where boundary information was called for processors not associated
  !   with dynamics when dyn_npes<npes
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use physical_constants, only: kappa
    use derivative_mod,     only: gradient_sphere, ugradv_sphere
    use edge_mod,           only: edgevpack, edgevunpack
    use bndry_mod,          only: bndry_exchangev
    use dyn_comp,           only: hvcoord
    use spmd_utils,         only: iam 
    use parallel_mod,       only: par 
    use dyn_grid,           only: fv_nphys
    use derivative_mod,     only: subcell_integration
    use element_ops,        only : get_temperature
    implicit none
    type(hybrid_t),        intent(in   ) :: hybrid
    type(element_t),target,intent(inout) :: elem(:)
    type(derivative_t),    intent(in   ) :: ederiv
    integer,               intent(in   ) :: nets,nete,nphys
    integer,               intent(in   ) :: tl ! timelevel to use
    real(kind=real_kind),  intent(out  ) :: frontgf(nphys,nphys,nlev,nets:nete)
    real(kind=real_kind),  intent(out  ) :: frontga(nphys,nphys,nlev,nets:nete)
  
    ! local
    integer :: k,kptr,i,j,ie,component
    real(kind=real_kind) :: frontgf_gll(np,np,nlev,nets:nete)
    real(kind=real_kind) :: gradth_gll(np,np,2,nlev,nets:nete)  ! grad(theta)
    real(kind=real_kind) :: gradth_fv(fv_nphys,fv_nphys,2)      ! grad(theta)
    real(kind=real_kind) :: p(np,np)        ! pressure at mid points
    real(kind=real_kind) :: theta(np,np)    ! potential temperature at mid points
    real(kind=real_kind) :: temperature(np,np,nlev)  ! Temperature
    real(kind=real_kind) :: C(np,np,2)     
    real(r8), dimension(np,np)             :: tmp_area
    real(r8), dimension(fv_nphys,fv_nphys) :: inv_area
    !---------------------------------------------------------------------------

    do ie = nets,nete
      do k = 1,nlev
        
        p(:,:) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,tl)
        call get_temperature(elem(ie),temperature,hvcoord,tl)

        ! potential temperature: theta = T (p/p0)^kappa
        theta(:,:) = temperature(:,:,k)*(psurf_ref / p(:,:))**kappa

        ! pressure at mid points
        p(:,:) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,tl)
        
        gradth_gll(:,:,:,k,ie) = gradient_sphere(theta,ederiv,elem(ie)%Dinv)
        
        ! compute C = (grad(theta) dot grad ) u
        C(:,:,:) = ugradv_sphere(gradth_gll(:,:,:,k,ie), elem(ie)%state%v(:,:,:,k,tl),ederiv,elem(ie))
        
        ! gradth dot C
        frontgf_gll(:,:,k,ie) = -( C(:,:,1)*gradth_gll(:,:,1,k,ie) + C(:,:,2)*gradth_gll(:,:,2,k,ie)  )
        
        ! apply mass matrix
        gradth_gll(:,:,1,k,ie)=gradth_gll(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
        gradth_gll(:,:,2,k,ie)=gradth_gll(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
        frontgf_gll(:,:,k,ie)=frontgf_gll(:,:,k,ie)*elem(ie)%spheremp(:,:)
          
      end do ! k
      ! pack
      call edgeVpack(edge3, frontgf_gll(:,:,:,ie),nlev,0,ie)
      call edgeVpack(edge3, gradth_gll(:,:,:,:,ie),2*nlev,nlev,ie)
    end do ! ie

    ! Boundary exchange
    if (par%dynproc) call bndry_exchangeV(hybrid,edge3)

    do ie = nets,nete
      ! unpack
      call edgeVunpack(edge3, frontgf_gll(:,:,:,ie),nlev,0,ie)
      call edgeVunpack(edge3, gradth_gll(:,:,:,:,ie),2*nlev,nlev,ie)
      ! apply inverse mass matrix
      do k = 1,nlev
        gradth_gll(:,:,1,k,ie) = gradth_gll(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
        gradth_gll(:,:,2,k,ie) = gradth_gll(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        frontgf_gll(:,:,k,ie) = frontgf_gll(:,:,k,ie)*elem(ie)%rspheremp(:,:)
      end do
      if (fv_nphys>0) then
        tmp_area(:,:) = 1.0_r8
        inv_area(:,:) = 1.0_r8/subcell_integration(tmp_area,np,fv_nphys,elem(ie)%metdet(:,:))
        do k = 1,nlev
          
          frontgf(:,:,k,ie) = subcell_integration(frontgf_gll(:,:,k,ie),      &
                              np, fv_nphys, elem(ie)%metdet(:,:) ) * inv_area
          
          gradth_fv(:,:,1) = subcell_integration(gradth_gll(:,:,1,k,ie),      &
                             np, fv_nphys, elem(ie)%metdet(:,:) ) * inv_area

          gradth_fv(:,:,2) = subcell_integration(gradth_gll(:,:,2,k,ie),      &
                             np, fv_nphys, elem(ie)%metdet(:,:) ) * inv_area

          do j=1,fv_nphys
            do i=1,fv_nphys
              frontga(i,j,k,ie) = atan2 ( gradth_fv(i,j,2) , &
                                          gradth_fv(i,j,1) + 1.e-10_r8 )
            end do
          end do

        end do ! k
      else
        do k=1,nlev
          frontgf(:,:,k,ie) = frontgf_gll(:,:,k,ie)
          ! Frontogenesis angle
          frontga(:,:,k,ie) = atan2 ( gradth_gll(:,:,2,k,ie) , &
                                      gradth_gll(:,:,1,k,ie) + 1.e-10_real_kind )
        end do
      end if ! fv_nphys>0
    end do ! ie

  end subroutine compute_frontogenesis
  !-------------------------------------------------------------------------------------------------
end module gravity_waves_sources
