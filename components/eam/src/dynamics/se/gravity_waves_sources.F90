module gravity_waves_sources
  use derivative_mod, only : derivative_t
  use dimensions_mod, only : np,nlev,nlevp
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
  subroutine gws_src_fnct(elem, tl, nphys, use_fgf_pgrad_correction, use_fgf_zgrad_correction, frontgf, frontga)
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
    logical,               intent(in   ) :: use_fgf_pgrad_correction
    logical,               intent(in   ) :: use_fgf_zgrad_correction
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
    call compute_frontogenesis( frontgf_thr, frontga_thr, tl, &
                                use_fgf_pgrad_correction, &
                                use_fgf_zgrad_correction, &
                                elem, hybrid, nets, nete, nphys )
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
  subroutine compute_frontogenesis( frontgf, frontga, tl, &
                                    use_fgf_pgrad_correction, &
                                    use_fgf_zgrad_correction, &
                                    elem, hybrid, nets, nete, nphys )
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute frontogenesis function F
  !   F =  -gradth dot C
  ! with:
  !   theta  = potential temperature
  !   gradth = grad(theta)
  !   C = ( gradth dot grad ) U
  ! 
  ! For more details on the frontal GW scheme and fronotogenesis source calculation:
  !   Charron, M., and E. Manzini, 2002: Gravity Waves from Fronts:
  !   Parameterization and Middle Atmosphere Response in a General
  !   Circulation Model. J. Atmos. Sci., 59, 923–941
  !
  ! Original by Mark Taylor, July 2011
  ! Change by Santos, 10 Aug 2011:
  !   Integrated into gravity_waves_sources module, several arguments made global
  !   to prevent repeated allocation/initialization
  ! Change by Aaron Donahue, April 2017:
  !   Fixed bug where boundary information was called for processors not associated
  !   with dynamics when dyn_npes<npes
  ! Change by Walter Hannah and Wandi Yu, Oct 2025:
  !   added pressure gradient correction options for consistency with original C&M paper
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use physical_constants, only: kappa, g
    use derivative_mod,     only: gradient_sphere, ugradv_sphere
    use edge_mod,           only: edge_g, edgevpack_nlyr, edgevunpack_nlyr
    use bndry_mod,          only: bndry_exchangev
    use dyn_comp,           only: hvcoord
    use spmd_utils,         only: iam 
    use parallel_mod,       only: par 
    use element_ops,        only: get_temperature, get_hydro_pressure_i
#ifdef MODEL_THETA_L
    use element_ops,        only: get_phi
#else
    use element_ops,        only: get_phi, get_phi_i
#endif
    use dyn_grid,           only: fv_nphys
    use prim_driver_mod,    only: deriv1
    use gllfvremap_mod,     only: gfr_g2f_scalar, gfr_g2f_vector
    implicit none
    type(hybrid_t),        intent(in   ) :: hybrid
    type(element_t),target,intent(inout) :: elem(:)
    integer,               intent(in   ) :: nets,nete,nphys
    integer,               intent(in   ) :: tl ! timelevel to use
    logical,               intent(in   ) :: use_fgf_pgrad_correction
    logical,               intent(in   ) :: use_fgf_zgrad_correction
    real(kind=real_kind),  intent(out  ) :: frontgf(nphys,nphys,nlev,nets:nete)
    real(kind=real_kind),  intent(out  ) :: frontga(nphys,nphys,nlev,nets:nete)
  
    ! local
    integer :: k,kptr,i,j,ie,component
    real(kind=real_kind) :: frontgf_gll(np,np,nlev,nets:nete)
    real(kind=real_kind) :: gradth_gll(np,np,2,nlev,nets:nete)  ! grad(theta)
    real(kind=real_kind) :: zint(np,np,nlevp)           ! interface geopotential
    real(kind=real_kind) :: zmid(np,np,nlev)            ! mid-point geopotential
    real(kind=real_kind) :: pint(np,np,nlevp)           ! interface hydrostatic pressure
    real(kind=real_kind) :: pmid(np,np,nlev)            ! mid-point hydrostatic pressure
    real(kind=real_kind) :: temperature(np,np,nlev)     ! Temperature
    real(kind=real_kind) :: C(np,np,2), wf1(nphys*nphys,nlev), wf2(nphys*nphys,nlev)

    real(kind=real_kind) :: theta(np,np,nlev)           ! potential temperature at mid points
    real(kind=real_kind) :: dum_grad(np,np,2)           ! temporary variable for spherical gradient calcualtions
    real(kind=real_kind) :: dum_cart(np,np,3,nlev)      ! temporary variable for cartesian gradient calcualtions

    ! variables needed for eta to pressure surface correction
    real(kind=real_kind) :: grad_vert_gll(np,np,2,nlev) ! grad(vertical coordinate) - either pressure or geopotential
    real(kind=real_kind) :: theta_dvert(np,np,nlev)     ! d(theta)/dp for eta to pressure surface correction
    real(kind=real_kind) :: dum_cart_dvert(np,np,3,nlev)! vertical pressure derivative of dum_cart

    !---------------------------------------------------------------------------
    !
    ! For a vector velocity "v", a tensor "grad(v)", and a vector "grad(theta)",
    ! this loop computes the vector "grad(theta)*grad(v)"
    !
    ! Representing the tensor "grad(v)" in spherical coordinates is difficult.
    ! This routine avoids this by computing a mathematically equivalent form 
    ! using a mixture of Cartesian and spherical coordinates
    !
    ! This routine is a modified version of derivative_mod.F90:ugradv_sphere()
    ! in that the grad(v) term is modified to compute grad_z(v) (or grad_p(v))
    ! - the gradient on z-surfaces expressed in terms of the gradient on model
    ! surfaces and a vertical geopotential gradient.
    !
    ! The old version only computed gradients on model surfaces, which creates
    ! issues around topograpy. This is address with use_fgf_pgrad_correction=.true.
    !
    ! First, v is represented in cartesian coordinates  v(c) for c=1,2,3
    ! For each v(c), we compute its gradient on z (or p) surfaces via:
    !    grad(v(c)) - d(v(c))/dz grad(z)
    ! Each of these gradients is represented in *spherical* coordinates (i=1,2)
    !
    ! We then dot each of these vectors with grad(theta).  This dot product is
    ! computed in spherical coordinates.  The end result is dum_cart(c),
    ! for c=1,2,3 These three scalars are the three Cartesian coefficients of
    ! the vector "grad(theta)*grad(v)"
    !
    ! This Cartesian vector is then transformed back to spherical coordinates
    !
    !---------------------------------------------------------------------------

    do ie = nets,nete

      if (use_fgf_pgrad_correction .or. use_fgf_zgrad_correction) then

        ! compute pressure at interfaces and mid-points
        call get_hydro_pressure_i(pint,elem(ie)%state%dp3d(:,:,:,tl),hvcoord)

        call get_temperature(elem(ie),temperature,hvcoord,tl)

        ! calculate pmid and potential temperature: theta = T (p/p0)^kappa
        do k = 1,nlev
          pmid(:,:,k) = pint(:,:,k) + elem(ie)%state%dp3d(:,:,k,tl)/2
          theta(:,:,k) = temperature(:,:,k)*(psurf_ref / pmid(:,:,k))**kappa
        end do

        if (use_fgf_pgrad_correction) then
          ! compute d(theta)/dp
          call compute_vertical_derivative(pint,theta,theta_dvert)
        end if

        if (use_fgf_zgrad_correction) then
          ! compute geopotential
#ifdef MODEL_THETA_L
          call get_phi(elem(ie), zmid, zint, hvcoord, tl)
#else
          call get_phi(elem(ie), zmid, hvcoord, tl, -1)
          call get_phi_i(elem(ie), zint, hvcoord, tl, -1)
#endif
          ! compute d(theta)/dz
          call compute_vertical_derivative(zint,theta,theta_dvert)
        end if

        do k = 1,nlev
          gradth_gll(:,:,:,k,ie) = gradient_sphere(theta(:,:,k),deriv1,elem(ie)%Dinv)
          if (use_fgf_pgrad_correction) grad_vert_gll(:,:,:,k) = gradient_sphere(pmid(:,:,k),deriv1,elem(ie)%Dinv)
          if (use_fgf_zgrad_correction) grad_vert_gll(:,:,:,k) = gradient_sphere(zmid(:,:,k),deriv1,elem(ie)%Dinv)
          do component=1,2
            gradth_gll(:,:,component,k,ie) = gradth_gll(:,:,component,k,ie) - theta_dvert(:,:,k) * grad_vert_gll(:,:,component,k)
          end do
        end do

        do k = 1,nlev
          ! latlon -> cartesian - Summing along the third dimension is a sum over components for each point
          do component=1,3
            dum_cart(:,:,component,k)=sum( elem(ie)%vec_sphere2cart(:,:,component,:) * elem(ie)%state%v(:,:,:,k,tl) ,3)
          end do
        end do

        do component=1,3
          if (use_fgf_pgrad_correction) call compute_vertical_derivative(pint,dum_cart(:,:,component,:),dum_cart_dvert(:,:,component,:))
          if (use_fgf_zgrad_correction) call compute_vertical_derivative(zint,dum_cart(:,:,component,:),dum_cart_dvert(:,:,component,:))
        end do

        do k = 1,nlev
          ! Do ugradv on the cartesian components - Dot u with the gradient of each component
          do component=1,3
            dum_grad(:,:,:) = gradient_sphere(dum_cart(:,:,component,k),deriv1,elem(ie)%Dinv)
            do i=1,2
              dum_grad(:,:,i) = dum_grad(:,:,i) - dum_cart_dvert(:,:,component,k) * grad_vert_gll(:,:,i,k)
            end do
            dum_cart(:,:,component,k) = sum( gradth_gll(:,:,:,k,ie) * dum_grad ,3)
          enddo
          ! cartesian -> latlon - vec_sphere2cart is its own pseudoinverse.
          do i=1,2
            C(:,:,i) = sum(dum_cart(:,:,:,k)*elem(ie)%vec_sphere2cart(:,:,:,i), 3)
          end do
          ! gradth_gll dot C
          frontgf_gll(:,:,k,ie) = -( C(:,:,1)*gradth_gll(:,:,1,k,ie) + C(:,:,2)*gradth_gll(:,:,2,k,ie)  )
          ! apply mass matrix
          gradth_gll(:,:,1,k,ie) = gradth_gll(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
          gradth_gll(:,:,2,k,ie) = gradth_gll(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
          frontgf_gll(:,:,k,ie)  = frontgf_gll(:,:,k,ie)*elem(ie)%spheremp(:,:)
        end do ! k

      else ! .not. use_fgf_pgrad_correction

        do k = 1,nlev
          ! pressure at mid points - this pressure preserves the old behavior of E3SMv3 and prior
          pmid(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,tl)
          ! potential temperature: theta = T (p/p0)^kappa
          call get_temperature(elem(ie),temperature,hvcoord,tl)
          theta(:,:,k) = temperature(:,:,k)*(psurf_ref / pmid(:,:,k))**kappa
          gradth_gll(:,:,:,k,ie) = gradient_sphere(theta(:,:,k),deriv1,elem(ie)%Dinv)
          ! compute C = (grad(theta) dot grad ) u
          C(:,:,:) = ugradv_sphere(gradth_gll(:,:,:,k,ie), elem(ie)%state%v(:,:,:,k,tl),deriv1,elem(ie))
          ! gradth_gll dot C
          frontgf_gll(:,:,k,ie) = -( C(:,:,1)*gradth_gll(:,:,1,k,ie) + C(:,:,2)*gradth_gll(:,:,2,k,ie)  )
          ! apply mass matrix
          gradth_gll(:,:,1,k,ie) = gradth_gll(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
          gradth_gll(:,:,2,k,ie) = gradth_gll(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
          frontgf_gll(:,:,k,ie)  = frontgf_gll(:,:,k,ie)*elem(ie)%spheremp(:,:)
        end do ! k

      end if ! use_fgf_pgrad_correction

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
  subroutine compute_vertical_derivative(vert_int, data_mid, ddata_dvert)
    !---------------------------------------------------------------------------
    real(kind=real_kind),   intent(in ) :: vert_int(np,np,nlev) ! vertical coord on interfaces (i.e. pint or zint)
    real(kind=real_kind),   intent(in ) :: data_mid(np,np,nlev) ! input data in mid-points
    real(kind=real_kind),   intent(out) :: ddata_dvert(np,np,nlev) ! vertical derivative of data
    !---------------------------------------------------------------------------
    integer :: k
    real(kind=real_kind) :: vert_above(np,np) ! pressure interpolated to interface above the current k mid-point
    real(kind=real_kind) :: vert_below(np,np) ! pressure interpolated to interface below the current k mid-point
    real(kind=real_kind) :: data_above(np,np) ! data interpolated to interface above the current k mid-point
    real(kind=real_kind) :: data_below(np,np) ! data interpolated to interface below the current k mid-point
    !---------------------------------------------------------------------------
    do k = 1,nlev
      if (k==1) then
        data_above = data_mid(:,:,k)
        data_below = ( data_mid(:,:,k+1) + data_mid(:,:,k) ) / 2.0 ! interpolate to interface k+1
        vert_above = ( vert_int(:,:,k+1) + vert_int(:,:,k) ) / 2.0 ! interpolate to mid-point k
        vert_below = vert_int(:,:,k+1)
      elseif (k==nlev) then
        data_above = ( data_mid(:,:,k-1) + data_mid(:,:,k) ) / 2.0 ! interpolate to interface
        data_below = data_mid(:,:,k)
        vert_above = vert_int(:,:,k)
        vert_below = ( vert_int(:,:,k+1) + vert_int(:,:,k) ) / 2.0 ! interpolate to mid-point k
      else
        data_above = ( data_mid(:,:,k-1) + data_mid(:,:,k) ) / 2.0 ! interpolate to interface k
        data_below = ( data_mid(:,:,k+1) + data_mid(:,:,k) ) / 2.0 ! interpolate to interface k+1
        vert_above = vert_int(:,:,k)
        vert_below = vert_int(:,:,k+1)
      end if
      ddata_dvert(:,:,k) = ( data_above - data_below ) / ( vert_above - vert_below )
    end do
  end subroutine compute_vertical_derivative
  !-------------------------------------------------------------------------------------------------
end module gravity_waves_sources
