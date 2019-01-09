module global_norms_mod

  use shr_kind_mod, only: r8=>shr_kind_r8
  use cam_logfile,  only: iulog
  use edgetype_mod, only: EdgeBuffer_t

  implicit none
  private
  save

  public :: l1_snorm
  public :: l2_snorm
  public :: linf_snorm

  public :: l1_vnorm
  public :: l2_vnorm
  public :: linf_vnorm

  public :: print_cfl
  public :: test_global_integral
  public :: global_integral
  public :: global_integrals_general
  public :: wrap_repro_sum

  private :: global_maximum
  type (EdgeBuffer_t), private :: edgebuf

contains


  subroutine global_integrals(elem, h,hybrid,npts,num_flds,nets,nete,I_sphere)
    use hybrid_mod,     only: hybrid_t
    use element_mod,    only: element_t
    use dimensions_mod, only: np, nelemd
    use physconst,      only: pi
    use parallel_mod,   only: global_shared_buf, global_shared_sum

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete,num_flds
    real (kind=r8), intent(in) :: h(npts,npts,num_flds,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=r8) :: I_sphere(num_flds)
    
    real (kind=r8) :: I_priv
    real (kind=r8) :: I_shared
    common /gblintcom/I_shared
    !
    ! Local variables
    !
    integer :: ie,j,i,q

    real (kind=r8) :: da
    real (kind=r8) :: J_tmp(nets:nete,num_flds)
    !
    ! This algorithm is independent of thread count and task count.
    ! This is a requirement of consistancy checking in cam.
    !
    J_tmp = 0.0_r8

!JMD    print *,'global_integral: before loop'
    do ie=nets,nete
      do q=1,num_flds
        do j=1,np
          do i=1,np
            da = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
            J_tmp(ie,q) = J_tmp(ie,q) + da*h(i,j,q,ie)
          end do
        end do
      end do
    end do
    do ie=nets,nete
      global_shared_buf(ie,1:num_flds) = J_tmp(ie,:)
    enddo
    !JMD    print *,'global_integral: before wrap_repro_sum'
    call wrap_repro_sum(nvars=num_flds, comm=hybrid%par%comm)
    !JMD    print *,'global_integral: after wrap_repro_sum'
    I_sphere(:) =global_shared_sum(1:num_flds) /(4.0_r8*PI)
  end subroutine global_integrals

  subroutine global_integrals_general(h,hybrid,npts,da,num_flds,nets,nete,I_sphere)
    use hybrid_mod,     only: hybrid_t
    use dimensions_mod, only: nc, nelemd
    use physconst,      only: pi
    use parallel_mod,   only: global_shared_buf, global_shared_sum

    integer,         intent(in) :: npts,nets,nete,num_flds
    real (kind=r8),  intent(in) :: h(npts,npts,num_flds,nets:nete)
    type (hybrid_t), intent(in) :: hybrid
    real (kind=r8),  intent(in) :: da(npts,npts,nets:nete)

    real (kind=r8) :: I_sphere(num_flds)
    
    real (kind=r8) :: I_priv
    real (kind=r8) :: I_shared
    common /gblintcom/I_shared
    !
    ! Local variables
    !
    integer :: ie,j,i,q

    real (kind=r8) :: J_tmp(nets:nete,num_flds)
    !
    ! This algorithm is independent of thread count and task count.
    ! This is a requirement of consistancy checking in cam.
    !
    J_tmp = 0.0_r8

!JMD    print *,'global_integral: before loop'
    do ie=nets,nete
      do q=1,num_flds
        do j=1,npts
          do i=1,npts
            J_tmp(ie,q) = J_tmp(ie,q) + da(i,j,ie)*h(i,j,q,ie)
          end do
        end do
      end do
    end do
    do ie=nets,nete
      global_shared_buf(ie,1:num_flds) = J_tmp(ie,:)
    enddo
    !JMD    print *,'global_integral: before wrap_repro_sum'
    call wrap_repro_sum(nvars=num_flds, comm=hybrid%par%comm)
    !JMD    print *,'global_integral: after wrap_repro_sum'
    I_sphere(:) =global_shared_sum(1:num_flds) /(4.0_r8*PI)
  end subroutine global_integrals_general


  ! ================================
  ! global_integral:
  !
  ! eq 81 in Williamson, et. al. p 218
  ! for spectral elements
  !
  ! ================================
  ! --------------------------
  function global_integral(elem, h,hybrid,npts,nets,nete) result(I_sphere)
    use hybrid_mod,     only: hybrid_t
    use element_mod,    only: element_t
    use dimensions_mod, only: np, nelemd
    use physconst,      only: pi
    use parallel_mod,   only: global_shared_buf, global_shared_sum

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: h(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=r8) :: I_sphere

    real (kind=r8) :: I_priv
    real (kind=r8) :: I_shared
    common /gblintcom/I_shared

    ! Local variables

    integer :: ie,j,i
    real(kind=r8) :: I_tmp(1)

    real (kind=r8) :: da
    real (kind=r8) :: J_tmp(nets:nete)
!
! This algorythm is independent of thread count and task count.
! This is a requirement of consistancy checking in cam.
!
    J_tmp = 0.0_r8

!JMD    print *,'global_integral: before loop'
       do ie=nets,nete
          do j=1,np
             do i=1,np
                da = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
                J_tmp(ie) = J_tmp(ie) + da*h(i,j,ie)
             end do
          end do
       end do
    do ie=nets,nete
      global_shared_buf(ie,1) = J_tmp(ie)
    enddo
!JMD    print *,'global_integral: before wrap_repro_sum'
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
!JMD    print *,'global_integral: after wrap_repro_sum'
    I_tmp = global_shared_sum(1)
!JMD    print *,'global_integral: after global_shared_sum'
    I_sphere = I_tmp(1)/(4.0_r8*PI)

  end function global_integral

  ! ================================
  ! test_global_integral:
  !
  ! test that the global integral of
  ! the area of the sphere is 1.
  !
  ! ================================

  subroutine test_global_integral(elem,hybrid,nets,nete,mindxout)
    use hybrid_mod,     only: hybrid_t
    use element_mod,    only: element_t
    use dimensions_mod, only: np,ne, nelem, nelemd
    use mesh_mod,       only: MeshUseMeshFile
    use reduction_mod,  only: ParallelMin,ParallelMax
    use physconst,      only: pi, ra, rearth
    use parallel_mod,   only: global_shared_buf, global_shared_sum

    type(element_t), intent(inout)         :: elem(:)
    integer,         intent(in)            :: nets,nete
    type(hybrid_t),  intent(in)            :: hybrid

    real(kind=r8),   intent(out), optional :: mindxout

    ! Local variables
    real(kind=r8) :: I_sphere
    real(kind=r8) :: h(np,np,nets:nete)
    ! Element statisics
    real(kind=r8) :: min_area,max_area,avg_area, max_ratio
    real(kind=r8) :: min_min_dx, max_min_dx, avg_min_dx
    real(kind=r8) :: min_normDinv, max_normDinv
    real(kind=r8) :: min_len
    integer       :: ie


    h(:,:,nets:nete)=1.0_r8

    ! Calculate surface area by integrating 1.0d0 over sphere and dividing by 4*PI
    ! (Should be 1)
    I_sphere = global_integral(elem, h(:,:,nets:nete),hybrid,np,nets,nete)

    min_area=1d99
    max_area=0
    avg_area=0_r8

    max_ratio = 0

    min_normDinv=1d99
    max_normDinv=0

    min_min_dx=1d99
    max_min_dx=0
    avg_min_dx=0_r8

    do ie=nets,nete

       elem(ie)%area = sum(elem(ie)%spheremp(:,:))
       min_area=min(min_area,elem(ie)%area)
       max_area=max(max_area,elem(ie)%area)

       min_normDinv = min(min_normDinv,elem(ie)%normDinv)
       max_normDinv = max(max_normDinv,elem(ie)%normDinv)

       max_ratio   = max(max_ratio,elem(ie)%dx_long/elem(ie)%dx_short)


       min_min_dx = min(min_min_dx,elem(ie)%dx_short)
       max_min_dx = max(max_min_dx,elem(ie)%dx_short)


       global_shared_buf(ie,1) = elem(ie)%area
       global_shared_buf(ie,2) = elem(ie)%dx_short

    enddo

    min_area=ParallelMin(min_area,hybrid)
    max_area=ParallelMax(max_area,hybrid)

    min_normDinv=ParallelMin(min_normDinv,hybrid)
    max_normDinv=ParallelMax(max_normDinv,hybrid)

    max_ratio=ParallelMax(max_ratio,hybrid)

    min_min_dx=ParallelMin(min_min_dx,hybrid)
    max_min_dx=ParallelMax(max_min_dx,hybrid)

    call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)

    avg_area = global_shared_sum(1)/dble(nelem)
    avg_min_dx = global_shared_sum(2)/dble(nelem)

    ! Physical units for area
    min_area = min_area*rearth*rearth/1000000._r8
    max_area = max_area*rearth*rearth/1000000._r8
    avg_area = avg_area*rearth*rearth/1000000._r8


    ! for an equation du/dt = i c u, leapfrog is stable for |c u dt| < 1
    ! Consider a gravity wave at the equator, c=340m/s
    ! u = exp(i kmax x/ a ) with x = longitude,  and kmax =  pi a / dx,
    ! u = exp(i pi x / dx ),   so du/dt = c du/dx becomes du/dt = i c pi/dx u
    ! stable for dt < dx/(c*pi)
    ! CAM 26 level AMIP simulation: max gravity wave speed 341.75 m/s
    if (hybrid%masterthread) then
       write(iulog,* )""
       write(iulog,* )"Running Global Integral Diagnostic..."
       write(iulog,*)"Area of unit sphere is",I_sphere
       write(iulog,*)"Should be 1.0 to round off..."
       write(iulog,'(a,f9.3)') 'Element area:  max/min',(max_area/min_area)
       if (.not.MeshUseMeshFile) then
           write(iulog,'(a,f6.3,f8.2)') "Average equatorial node spacing (deg, km) = ", &
                dble(90)/dble(ne*(np-1)), PI*rearth/(2000.0d0*dble(ne*(np-1)))
       end if
       write(iulog,'(a,2f9.3)') 'norm of Dinv (min, max): ', min_normDinv, max_normDinv
       write(iulog,'(a,1f8.2)') 'Max Dinv-based element distortion: ', max_ratio
       write(iulog,'(a,3f8.2)') 'dx based on Dinv svd:          ave,min,max = ', avg_min_dx, min_min_dx, max_min_dx
       write(iulog,'(a,3f8.2)') "dx based on sqrt element area: ave,min,max = ", &
                sqrt(avg_area)/(np-1),sqrt(min_area)/(np-1),sqrt(max_area)/(np-1)
    end if

    if(present(mindxout)) then
        ! min_len now based on norm(Dinv)
        min_len = 0.002d0*rearth/(dble(np-1)*max_normDinv)
        mindxout=1000_r8*min_len
    end if

  end subroutine test_global_integral


!------------------------------------------------------------------------------------

  ! ================================
  ! print_cfl:
  !
  ! Calculate / output CFL info
  ! (both advective and based on
  ! viscosity or hyperviscosity)
  !
  ! ================================

  subroutine print_cfl(elem,hybrid,nets,nete,dtnu)
!
!   estimate various CFL limits
!   also, for variable resolution viscosity coefficient, make sure
!   worse viscosity CFL (given by dtnu) is not violated by reducing
!   viscosity coefficient in regions where CFL is violated
!
    use hybrid_mod,     only: hybrid_t, PrintHybrid
    use element_mod,    only: element_t
    use dimensions_mod, only: np,ne,nelem,nelemd,nc,nhe,qsize,ntrac
    use quadrature_mod, only: gausslobatto, quadrature_t

    use reduction_mod,  only: ParallelMin,ParallelMax
    use physconst,      only: ra, rearth, pi
    use control_mod,    only: nu, nu_q, nu_div, nu_top, &
                            fine_ne, rk_stage_user, max_hypervis_courant
    use control_mod,    only: tstep_type, hypervis_power, hypervis_scaling
    use cam_abortutils, only: endrun
    use parallel_mod,   only: global_shared_buf, global_shared_sum
    use edge_mod,       only: initedgebuffer, FreeEdgeBuffer, edgeVpack, edgeVunpack
    use bndry_mod,      only: bndry_exchange
    use time_mod,       only: tstep

    type(element_t)      , intent(inout) :: elem(:)
    integer              , intent(in) :: nets,nete
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8), intent(in) :: dtnu

    ! Element statisics
    real (kind=r8) :: min_max_dx,max_unif_dx   ! used for normalizing scalar HV
    real (kind=r8) :: max_normDinv  ! used for CFL
    real (kind=r8) :: min_hypervis, max_hypervis, avg_hypervis, stable_hv
    real (kind=r8) :: normDinv_hypervis
    real (kind=r8) :: x, y, noreast, nw, se, sw
    real (kind=r8), dimension(np,np,nets:nete) :: zeta
    real (kind=r8) :: lambda_max, lambda_vis, min_gw, lambda
    integer :: ie,corner, i, j, rowind, colind
    type (quadrature_t)    :: gp


    ! Eigenvalues calculated by folks at UMich (Paul U & Jared W)
    select case (np)
        case (2)
            lambda_max = 0.5d0
            lambda_vis = 0.0d0  ! need to compute this
        case (3)
            lambda_max = 1.5d0
            lambda_vis = 12.0d0
        case (4)
            lambda_max = 2.74d0
            lambda_vis = 30.0d0
        case (5)
            lambda_max = 4.18d0
            lambda_vis = 91.6742d0
        case (6)
            lambda_max = 5.86d0
            lambda_vis = 190.1176d0
        case (7)
            lambda_max = 7.79d0
            lambda_vis = 374.7788d0
        case (8)
            lambda_max = 10.0d0
            lambda_vis = 652.3015d0
        case DEFAULT
            lambda_max = 0.0d0
            lambda_vis = 0.0d0
    end select

    if ((lambda_max.eq.0d0).and.(hybrid%masterthread)) then
        print*, "lambda_max not calculated for NP = ",np
        print*, "Estimate of gravity wave timestep will be incorrect"
    end if
    if ((lambda_vis.eq.0d0).and.(hybrid%masterthread)) then
        print*, "lambda_vis not calculated for NP = ",np
        print*, "Estimate of viscous CFLs will be incorrect"
    end if

    do ie=nets,nete
      elem(ie)%variable_hyperviscosity = 1.0_r8
    end do

    gp=gausslobatto(np)
    min_gw = minval(gp%weights)

    max_normDinv=0
    min_max_dx=1d99
    do ie=nets,nete
        max_normDinv = max(max_normDinv,elem(ie)%normDinv)
        min_max_dx = min(min_max_dx,elem(ie)%dx_long)
    enddo
    max_normDinv=ParallelMax(max_normDinv,hybrid)
    min_max_dx=ParallelMin(min_max_dx,hybrid)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SCALAR, RESOLUTION-AWARE HYPERVISCOSITY
!  this block of code initializes the variable_hyperviscsoity() array
!  based on largest length scale in each element and user specified scaling
!  it then limits the coefficient if the user specifed a max CFL
!  this limiting is based on the smallest length scale of each element
!  since that controls the CFL.
!  Mike Levy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (hypervis_power /= 0) then

        min_hypervis = 1d99
        max_hypervis = 0
        avg_hypervis = 0


        max_unif_dx = min_max_dx  ! use this for average resolution, unless:
!       viscosity in namelist specified for smallest element:
        if (fine_ne>0) then
           ! viscosity in namelist specified for regions with a resolution
           ! equivilant to a uniform grid with ne=fine_ne
           if (np /= 4 ) call endrun('ERROR: setting fine_ne only supported with NP=4')
           max_unif_dx = (111.28_r8*30)/dble(fine_ne)   ! in km
        endif

!
! note: if L = eigenvalue of metinv, then associated length scale (km) is
! dx = 1.0d0/( sqrt(L)*0.5d0*dble(np-1)*ra*1000.0d0)
!
!       for viscosity *tensor*, we take at each point:
!            nu1 = nu*(dx1/max_unif_dx)**3.2      dx1 associated with eigenvalue 1
!            nu2 = nu*(dx2/max_unif_dx)**3.2      dx2 associated with eigenvalue 2
!       with this approach:
!          - with this formula, no need to adjust for CFL violations
!          - if nu comes from a 3.2 scaling that is stable for coarse and fine resolutions,
!            this formulat will be stable.
!          - gives the correct answer in long skinny rectangles:
!            large viscosity in the long direction, small viscosity in the short direction
!
!

        normDinv_hypervis = 0

        do ie=nets,nete
           ! variable viscosity based on map from ulatlon -> ucontra

           ! dx_long
           elem(ie)%variable_hyperviscosity = sqrt((elem(ie)%dx_long/max_unif_dx) ** hypervis_power)
           elem(ie)%hv_courant = dtnu*(elem(ie)%variable_hyperviscosity(1,1)**2) * &
                (lambda_vis**2) * ((ra*elem(ie)%normDinv)**4)

            ! Check to see if this is stable
            if (elem(ie)%hv_courant.gt.max_hypervis_courant) then
                stable_hv = sqrt( max_hypervis_courant / &
                     (  dtnu * (lambda_vis)**2 * (ra*elem(ie)%normDinv)**4 ) )

#if 0
         ! Useful print statements for debugging the adjustments to hypervis
                 print*, "Adjusting hypervis on elem ", elem(ie)%GlobalId
                 print*, "From ", nu*elem(ie)%variable_hyperviscosity(1,1)**2, " to ", nu*stable_hv
                 print*, "Difference = ", nu*(/elem(ie)%variable_hyperviscosity(1,1)**2-stable_hv/)
                 print*, "Factor of ", elem(ie)%variable_hyperviscosity(1,1)**2/stable_hv
                 print*, " "
#endif

!                make sure that: elem(ie)%hv_courant <=  max_hypervis_courant
                elem(ie)%variable_hyperviscosity = stable_hv
                elem(ie)%hv_courant = dtnu*(stable_hv**2) * (lambda_vis)**2 * (ra*elem(ie)%normDinv)**4
            end if
            normDinv_hypervis = max(normDinv_hypervis, elem(ie)%hv_courant/dtnu)

            min_hypervis = min(min_hypervis, elem(ie)%variable_hyperviscosity(1,1))
            max_hypervis = max(max_hypervis, elem(ie)%variable_hyperviscosity(1,1))
            global_shared_buf(ie,1) = elem(ie)%variable_hyperviscosity(1,1)
        end do

        min_hypervis = ParallelMin(min_hypervis, hybrid)
        max_hypervis = ParallelMax(max_hypervis, hybrid)
        call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
        avg_hypervis = global_shared_sum(1)/dble(nelem)

        normDinv_hypervis = ParallelMax(normDinv_hypervis, hybrid)

        ! apply DSS (aka assembly procedure) to variable_hyperviscosity (makes continuous)
        call initEdgeBuffer(hybrid%par,edgebuf,elem,1)
        do ie=nets,nete
            zeta(:,:,ie) = elem(ie)%variable_hyperviscosity(:,:)*elem(ie)%spheremp(:,:)
            call edgeVpack(edgebuf,zeta(1,1,ie),1,0,ie)
        end do
        call bndry_exchange(hybrid,edgebuf)
        do ie=nets,nete
            call edgeVunpack(edgebuf,zeta(1,1,ie),1,0,ie)
            elem(ie)%variable_hyperviscosity(:,:) = zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
        end do
        call FreeEdgeBuffer(edgebuf)

        ! replace hypervis w/ bilinear based on continuous corner values
        do ie=nets,nete
            noreast = elem(ie)%variable_hyperviscosity(np,np)
            nw = elem(ie)%variable_hyperviscosity(1,np)
            se = elem(ie)%variable_hyperviscosity(np,1)
            sw = elem(ie)%variable_hyperviscosity(1,1)
            do i=1,np
                x = gp%points(i)
                do j=1,np
                    y = gp%points(j)
                    elem(ie)%variable_hyperviscosity(i,j) = 0.25d0*( &
                                            (1.0d0-x)*(1.0d0-y)*sw + &
                                            (1.0d0-x)*(y+1.0d0)*nw + &
                                            (x+1.0d0)*(1.0d0-y)*se + &
                                            (x+1.0d0)*(y+1.0d0)*noreast)
                end do
            end do
        end do
    else  if (hypervis_scaling/=0) then
       ! tensorHV.  New eigenvalues are the eigenvalues of the tensor V
       ! formulas here must match what is in cube_mod.F90
       ! for tensorHV, we scale out the rearth dependency
       lambda = max_normDinv**2
       normDinv_hypervis = (lambda_vis**2) * (max_normDinv**4) * &
            (lambda**(-hypervis_scaling/2) )
    else
       ! constant coefficient formula:
       normDinv_hypervis = (lambda_vis**2) * (ra*max_normDinv)**4
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  TENSOR, RESOLUTION-AWARE HYPERVISCOSITY
!  The tensorVisc() array is computed in cube_mod.F90
!  this block of code will DSS it so the tensor if C0
!  and also make it bilinear in each element.
!  Oksana Guba
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (hypervis_scaling /= 0) then

    call initEdgeBuffer(hybrid%par,edgebuf,elem,1)
    do rowind=1,2
      do colind=1,2
    do ie=nets,nete
      zeta(:,:,ie) = elem(ie)%tensorVisc(:,:,rowind,colind)*elem(ie)%spheremp(:,:)
      call edgeVpack(edgebuf,zeta(1,1,ie),1,0,ie)
    end do

    call bndry_exchange(hybrid,edgebuf)
    do ie=nets,nete
      call edgeVunpack(edgebuf,zeta(1,1,ie),1,0,ie)
          elem(ie)%tensorVisc(:,:,rowind,colind) = zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
    end do
      enddo !rowind
    enddo !colind
    call FreeEdgeBuffer(edgebuf)

!IF BILINEAR MAP OF V NEEDED

    do rowind=1,2
      do colind=1,2
    ! replace hypervis w/ bilinear based on continuous corner values
    do ie=nets,nete
      noreast = elem(ie)%tensorVisc(np,np,rowind,colind)
      nw = elem(ie)%tensorVisc(1,np,rowind,colind)
      se = elem(ie)%tensorVisc(np,1,rowind,colind)
      sw = elem(ie)%tensorVisc(1,1,rowind,colind)
      do i=1,np
        x = gp%points(i)
        do j=1,np
        y = gp%points(j)
        elem(ie)%tensorVisc(i,j,rowind,colind) = 0.25d0*( &
                    (1.0d0-x)*(1.0d0-y)*sw + &
                    (1.0d0-x)*(y+1.0d0)*nw + &
                    (x+1.0d0)*(1.0d0-y)*se + &
                    (x+1.0d0)*(y+1.0d0)*noreast)
        end do
      end do
    end do
      enddo !rowind
    enddo !colind

    endif
    deallocate(gp%points)
    deallocate(gp%weights)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (hybrid%masterthread) then
       write(iulog,'(a,f10.2)') 'CFL estimates in terms of S=time step stability region'
       write(iulog,'(a,f10.2)') '(i.e. advection w/leapfrog: S=1, viscosity w/forward Euler: S=2)'
       if (rk_stage_user>0) then
          write(iulog,'(a,f10.2,a)') 'SSP preservation (120m/s) RKSSP euler step dt  < S *', &
               min_gw/(120.0d0*max_normDinv*ra),'s'
       endif
       if (qsize>0) &
          write(iulog,'(a,f10.2,a)') 'Stability: advective (120m/s)   dt_tracer < S *',&
               1/(120.0d0*max_normDinv*lambda_max*ra),'s'
       if (ntrac>0) then
          !
          ! rough estimate of Courant number limted time-step:
          !
          ! U_max*dt_tracer/dx < nhe
          !
          ! where U_max=120 m/s and dx = 360 degrees/(4*ne*nc) = (2*pi*Rearth m)/(4*ne*nc)
          !
          write(iulog,'(a,f10.2,a)') "Stability (fvm Courant number): advective (120m/s)   dt_tracer < ",&
               dble(nhe)*(2.0_r8*pi*Rearth/dble(4.0_r8*ne*nc))/120.0e0_r8,'s'
          write(iulog,*) "(note that fvm stability is also limited by flow deformation - Lipschitz criterion!)"
       end if
       write(iulog,'(a,f10.2,a)') 'Stability: advective (120m/s)   dt_tracer < S *', &
                                   1/(120.0d0*max_normDinv*lambda_max*ra),'s'
       write(iulog,'(a,f10.2,a)') 'Stability: gravity wave(342m/s)   dt_dyn  < S *', &
                                   1/(342.0d0*max_normDinv*lambda_max*ra),'s'
       if (nu>0) then
!          if (hypervis_order==1) then
!              write(iulog,'(a,f10.2,a)') 'Stability: viscosity dt < S *',1/(((ra*max_normDinv)**2)*lambda_vis),'s'
!          endif
!          if (hypervis_order==2) then
             ! counrant number = dtnu*normDinv_hypervis  < S
             !  dt < S  1/nu*normDinv
             write(iulog,'(a,f10.2,a)') "Stability: nu_q   hyperviscosity dt < S *", 1/(nu_q*normDinv_hypervis),'s'
             write(iulog,'(a,f10.2,a)') "Stability: nu_vor hyperviscosity dt < S *", 1/(nu*normDinv_hypervis),'s'

             write(iulog,'(a,f10.2,a)') "Stability: nu_div hyperviscosity dt < S *", 1/(nu_div*normDinv_hypervis),'s'
!          endif
       endif
       if(nu_top>0) then
          write(iulog,'(a,f10.2,a)') 'TOP3 viscosity CFL: dt < S*', &
                                  1.0d0/(4*nu_top*((ra*max_normDinv)**2)*lambda_vis),'s'
       end if
      if (hypervis_power /= 0) then
        write(iulog,'(a,3e11.4)')'Hyperviscosity (dynamics): ave,min,max = ', &
                                  nu*(/avg_hypervis**2,min_hypervis**2,max_hypervis**2/)
!         print*, 'fine_ne = ', fine_ne
!         print*, 'Using max_unif_dx = ', max_unif_dx
      end if
      write(iulog,*) 'tstep_type = ',tstep_type
    end if

  end subroutine print_cfl

  ! ================================
  ! global_maximum:
  !
  ! Find global maximum on sphere
  !
  ! ================================

  function global_maximum(h,hybrid,npts,nets,nete) result(Max_sphere)

    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : red_max, pmax_mt

    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: h(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=r8) :: Max_sphere

    ! Local variables

    real (kind=r8) :: redp(1)

    Max_sphere = MAXVAL(h(:,:,nets:nete))

    redp(1) = Max_sphere
    call pmax_mt(red_max,redp,1,hybrid)
    Max_sphere = red_max%buf(1)

  end function global_maximum

  ! ==========================================================
  ! l1_snorm:
  !
  ! computes the l1 norm per Williamson et al, p. 218 eq(8)
  ! for a scalar quantity
  ! ===========================================================

  function l1_snorm(elem, h,ht,hybrid,npts,nets,nete) result(l1)

    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: l1

    ! Local variables

    real (kind=r8) :: dhabs(npts,npts,nets:nete)
    real (kind=r8) :: htabs(npts,npts,nets:nete)
    real (kind=r8) :: dhabs_int
    real (kind=r8) :: htabs_int
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dhabs(i,j,ie) = ABS(h(i,j,ie)-ht(i,j,ie))
             htabs(i,j,ie) = ABS(ht(i,j,ie))
          end do
       end do
    end do

    dhabs_int = global_integral(elem, dhabs(:,:,nets:nete),hybrid,npts,nets,nete)
    htabs_int = global_integral(elem, htabs(:,:,nets:nete),hybrid,npts,nets,nete)

    l1 = dhabs_int/htabs_int

  end function l1_snorm

  ! ===========================================================
  ! l1_vnorm:
  !
  ! computes the l1 norm per Williamson et al, p. 218 eq(97),
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function l1_vnorm(elem, v,vt,hybrid,npts,nets,nete) result(l1)
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: l1

    ! Local variables

    real (kind=r8), dimension(:,:,:,:), pointer :: met
    real (kind=r8) :: dvsq(npts,npts,nets:nete)
    real (kind=r8) :: vtsq(npts,npts,nets:nete)
    real (kind=r8) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: dv1,dv2
    real (kind=r8) :: vt1,vt2
    real (kind=r8) :: dvsq_int
    real (kind=r8) :: vtsq_int

    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met
       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = SQRT(dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2)
             vtsq(i,j,ie) = SQRT(vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2)

          end do
       end do
    end do

    dvsq_int = global_integral(elem, dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_int = global_integral(elem, vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    l1 = dvsq_int/vtsq_int

  end function l1_vnorm

  ! ==========================================================
  ! l2_snorm:
  !
  ! computes the l2 norm per Williamson et al, p. 218 eq(83)
  ! for a scalar quantity on the pressure grid.
  !
  ! ===========================================================

  function l2_snorm(elem, h,ht,hybrid,npts,nets,nete) result(l2)
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t), intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: l2

    ! Local variables

    real (kind=r8) :: dh2(npts,npts,nets:nete)
    real (kind=r8) :: ht2(npts,npts,nets:nete)
    real (kind=r8) :: dh2_int
    real (kind=r8) :: ht2_int
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dh2(i,j,ie)=(h(i,j,ie)-ht(i,j,ie))**2
             ht2(i,j,ie)=ht(i,j,ie)**2
          end do
       end do
    end do

    dh2_int = global_integral(elem,dh2(:,:,nets:nete),hybrid,npts,nets,nete)
    ht2_int = global_integral(elem,ht2(:,:,nets:nete),hybrid,npts,nets,nete)

    l2 = SQRT(dh2_int)/SQRT(ht2_int)

  end function l2_snorm

  ! ==========================================================
  ! l2_vnorm:
  !
  ! computes the l2 norm per Williamson et al, p. 219 eq(98)
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function l2_vnorm(elem, v,vt,hybrid,npts,nets,nete) result(l2)
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: l2

    ! Local variables

    real (kind=r8), dimension(:,:,:,:), pointer :: met
    real (kind=r8) :: dvsq(npts,npts,nets:nete)
    real (kind=r8) :: vtsq(npts,npts,nets:nete)
    real (kind=r8) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: dv1,dv2
    real (kind=r8) :: vt1,vt2
    real (kind=r8) :: dvsq_int
    real (kind=r8) :: vtsq_int
    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met
       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2
             vtsq(i,j,ie) = vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2

          end do
       end do
    end do

    dvsq_int = global_integral(elem, dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_int = global_integral(elem, vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    l2 = SQRT(dvsq_int)/SQRT(vtsq_int)

  end function l2_vnorm

  ! ==========================================================
  ! linf_snorm:
  !
  ! computes the l infinity norm per Williamson et al, p. 218 eq(84)
  ! for a scalar quantity on the pressure grid...
  !
  ! ===========================================================

  function linf_snorm(h,ht,hybrid,npts,nets,nete) result(linf)
    use hybrid_mod, only : hybrid_t
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: linf

    ! Local variables

    real (kind=r8) :: dhabs(npts,npts,nets:nete)
    real (kind=r8) :: htabs(npts,npts,nets:nete)
    real (kind=r8) :: dhabs_max
    real (kind=r8) :: htabs_max
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dhabs(i,j,ie)=ABS(h(i,j,ie)-ht(i,j,ie))
             htabs(i,j,ie)=ABS(ht(i,j,ie))
          end do
       end do
    end do

    dhabs_max = global_maximum(dhabs(:,:,nets:nete),hybrid,npts,nets,nete)
    htabs_max = global_maximum(htabs(:,:,nets:nete),hybrid,npts,nets,nete)

    linf = dhabs_max/htabs_max

  end function linf_snorm


  ! ==========================================================
  ! linf_vnorm:
  !
  ! computes the linf norm per Williamson et al, p. 218 eq(99),
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function linf_vnorm(elem,v,vt,hybrid,npts,nets,nete) result(linf)
    use hybrid_mod, only : hybrid_t
    use element_mod, only : element_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: linf

    ! Local variables

    real (kind=r8), dimension(:,:,:,:), pointer :: met
    real (kind=r8) :: dvsq(npts,npts,nets:nete)
    real (kind=r8) :: vtsq(npts,npts,nets:nete)
    real (kind=r8) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: dv1,dv2
    real (kind=r8) :: vt1,vt2
    real (kind=r8) :: dvsq_max
    real (kind=r8) :: vtsq_max
    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met

       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = SQRT(dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2)
             vtsq(i,j,ie) = SQRT(vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2)

          end do
       end do
    end do

    dvsq_max = global_maximum(dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_max = global_maximum(vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    linf = dvsq_max/vtsq_max

  end function linf_vnorm

  subroutine wrap_repro_sum (nvars, comm, nsize)
    use dimensions_mod,   only: nelemd
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc
    use cam_abortutils,   only: endrun
    use parallel_mod,     only: global_shared_buf, global_shared_sum, nrepro_vars

    integer :: nvars            !  number of variables to be summed (cannot exceed nrepro_vars)
    integer :: comm             !  mpi communicator
    integer, optional :: nsize  !  local buffer size (defaults to nelemd - number of elements in mpi task)

    integer nsize_use

    if (present(nsize)) then
       nsize_use = nsize
    else
       nsize_use = nelemd
    endif
    if (nvars .gt. nrepro_vars) call endrun('ERROR: repro_sum_buffer_size exceeded')

! Repro_sum contains its own OpenMP, so only one thread should call it (AAM)

!$OMP BARRIER
!$OMP MASTER

    call repro_sum(global_shared_buf, global_shared_sum, nsize_use, nelemd, nvars, commid=comm)


!$OMP END MASTER
!$OMP BARRIER

    end subroutine wrap_repro_sum

end module global_norms_mod
