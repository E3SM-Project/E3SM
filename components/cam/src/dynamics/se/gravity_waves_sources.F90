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
  
  real(r8) :: psurf_ref

CONTAINS
  !-------------------------------------------------------------------------------------------------
  subroutine gws_init(elem)
    use hycoef, only         : hypi
    use pmgrid, only         : plev
    use parallel_mod, only    : par
    implicit none
    type (element_t), intent(inout), dimension(:) :: elem
    !---------------------------------------------------------------------------
    ! Set up variables similar to dyn_comp and prim_driver_mod initializations
    
    psurf_ref = hypi(plev+1)

  end subroutine gws_init
  !-------------------------------------------------------------------------------------------------
  subroutine gws_src_fnct(elem, tl, nphys, frontgf, frontga)
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
    allocate(frontgf_thr(nphys,nphys,nlev,nets:nete))
    allocate(frontga_thr(nphys,nphys,nlev,nets:nete))
    call compute_frontogenesis(frontgf_thr,frontga_thr,tl,elem,hybrid,nets,nete,nphys)
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
  subroutine compute_frontogenesis(frontgf,frontga,tl,elem,hybrid,nets,nete,nphys)
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
    use edge_mod,           only : edge_g, edgevpack_nlyr, edgevunpack_nlyr
    use bndry_mod,          only: bndry_exchangev
    use dyn_comp,           only: hvcoord
    use spmd_utils,         only: iam 
    use parallel_mod,       only: par 
    use element_ops,        only : get_temperature
    use dyn_grid,           only: fv_nphys
    use prim_driver_mod, only     : deriv1
    use element_ops,        only : get_temperature
    use gllfvremap_mod, only  : gfr_g2f_scalar, gfr_g2f_vector
    implicit none
    type(hybrid_t),        intent(in   ) :: hybrid
    type(element_t),target,intent(inout) :: elem(:)
    integer,               intent(in   ) :: nets,nete,nphys
    integer,               intent(in   ) :: tl ! timelevel to use
    real(kind=real_kind),  intent(out  ) :: frontgf(nphys,nphys,nlev,nets:nete)
    real(kind=real_kind),  intent(out  ) :: frontga(nphys,nphys,nlev,nets:nete)
  
    ! local
    integer :: k,kptr,i,j,ie,component
    real(kind=real_kind) :: frontgf_gll(np,np,nlev,nets:nete)
    real(kind=real_kind) :: gradth_gll(np,np,2,nlev,nets:nete)  ! grad(theta)
    real(kind=real_kind) :: p(np,np)        ! pressure at mid points
    real(kind=real_kind) :: theta(np,np)    ! potential temperature at mid points
    real(kind=real_kind) :: temperature(np,np,nlev)  ! Temperature
    real(kind=real_kind) :: C(np,np,2), wf1(nphys*nphys,nlev), wf2(nphys*nphys,nlev)
    !---------------------------------------------------------------------------

    do ie = nets,nete
      do k = 1,nlev
        ! pressure at mid points
        p(:,:) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,tl)
        ! potential temperature: theta = T (p/p0)^kappa
        call get_temperature(elem(ie),temperature,hvcoord,tl)
        theta(:,:) = temperature(:,:,k)*(psurf_ref / p(:,:))**kappa
        gradth_gll(:,:,:,k,ie) = gradient_sphere(theta,deriv1,elem(ie)%Dinv)
        ! compute C = (grad(theta) dot grad ) u
        C(:,:,:) = ugradv_sphere(gradth_gll(:,:,:,k,ie), elem(ie)%state%v(:,:,:,k,tl),deriv1,elem(ie))
        ! gradth_gll dot C
        frontgf_gll(:,:,k,ie) = -( C(:,:,1)*gradth_gll(:,:,1,k,ie) + C(:,:,2)*gradth_gll(:,:,2,k,ie)  )
        ! apply mass matrix
        gradth_gll(:,:,1,k,ie) = gradth_gll(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
        gradth_gll(:,:,2,k,ie) = gradth_gll(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
        frontgf_gll(:,:,k,ie)  = frontgf_gll(:,:,k,ie)*elem(ie)%spheremp(:,:)
      end do ! k
      ! pack
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, frontgf_gll(:,:,:,ie),nlev,0,3*nlev)
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, gradth_gll(:,:,:,:,ie),2*nlev,nlev,3*nlev)

    end do ! ie

    ! Boundary exchange
    if (par%dynproc) call bndry_exchangeV(hybrid,edge_g)

    do ie = nets,nete
      ! unpack
     call edgeVunpack_nlyr(edge_g, elem(ie)%desc,frontgf_gll(:,:,:,ie),nlev,0,3*nlev)
     call edgeVunpack_nlyr(edge_g, elem(ie)%desc,gradth_gll(:,:,:,:,ie),2*nlev,nlev,3*nlev)
      ! apply inverse mass matrix
      do k = 1,nlev
        gradth_gll(:,:,1,k,ie) = gradth_gll(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
        gradth_gll(:,:,2,k,ie) = gradth_gll(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        frontgf_gll(:,:,k,ie) = frontgf_gll(:,:,k,ie)*elem(ie)%rspheremp(:,:)
        if (fv_nphys == 0) then ! physics on GLL nodes
          frontgf(:,:,k,ie) = frontgf_gll(:,:,k,ie)
          ! Frontogenesis angle
          frontga(:,:,k,ie) = atan2(gradth_gll(:,:,2,k,ie) , &
                                    gradth_gll(:,:,1,k,ie) + 1.e-10_real_kind )
        end if
      end do ! k
      if (fv_nphys > 0) then
        call gfr_g2f_scalar(ie, elem(ie)%metdet, frontgf_gll(:,:,:,ie), wf1)
        frontgf(:,:,:,ie) = reshape(wf1, (/nphys,nphys,nlev/))
        call gfr_g2f_vector(ie, elem, &
             gradth_gll(:,:,1,:,ie), gradth_gll(:,:,2,:,ie), &
             wf1, wf2)
        frontga(:,:,:,ie) = reshape( &
             atan2(wf2, wf1 + 1.e-10_real_kind), &
             (/nphys,nphys,nlev/))
      end if
    end do ! ie

  end subroutine compute_frontogenesis
  !-------------------------------------------------------------------------------------------------
end module gravity_waves_sources
