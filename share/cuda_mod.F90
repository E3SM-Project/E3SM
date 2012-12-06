!This is where all of the PGI CUDA FORTRAN code will go, and these routines will be called from prim_advection_mod.
!This is compiled regardless, but PGI-specific calls are always wrapped in the _ACCEL ifdefs that are automagically
!activated when -Mcuda is specified during compilation with a PGI compiler. Thus, it will be ignored unless explicitly
!activated by the user
!
!As a general rule, all of the routines in here will be called within a threaded context (assuming ELEMENT_OPENMP is not
!deifned), and therefore, we enforce BARRIERS, MASTERS, and SINGLES from within these routines rather than outside them.
!This is to minimize the visible code impacts on the existing CPU code.

! Please pay attention to this all caps passive aggresive banner.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!                     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  STATUS INCOMPLETE  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  DO NOT USE YET     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  UNTIL THIS BANNER  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  IS REMOVED         !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!                     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cuda_mod
#ifdef _ACCEL
!Put everything CUDA-specific in here so it doesn't get compiled without -Mcuda enabled on a PGI compiler
  use cudafor
  use kinds          , only : real_kind
  use dimensions_mod, only: np,nlevp,nlev,qsize,qsize_d,max_corner_elem,max_neigh_edges,nelemd
  use element_mod, only: timelevels
  use edge_mod, only: EdgeBuffer_t
  implicit none
  private

  !First listed are all externally accescible routines
  public :: cuda_mod_init
  public :: euler_step_cuda 

  type(EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdvQ2, edgeAdvDSS
  real(kind=real_kind), allocatable :: qmin(:,:,:), qmax(:,:,:)
  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1



contains



  !The point of this is to initialize any data required in other routines of this module as well
  !as to run one initial CUDA kernel just to get those overheads out of the way so that subsequent
  !timing routines are accurage.
  subroutine cuda_mod_init()
    use edge_mod      , only: initEdgeBuffer
    implicit none
    integer     :: ierr, ie
    type (dim3) :: griddim,blockdim
!$OMP BARRIER
!$OMP MASTER
    write(*,*) "cuda_mod_init"

    blockdim = dim3(1,1,1)
    griddim  = dim3(1,1,1)
    call warmup <<< griddim , blockdim >>> ( ie )
    ierr = cudaThreadSynchronize()

!$OMP END MASTER
    write(*,*) __LINE__
    call initEdgeBuffer(edgeAdv   ,qsize*nlev  )
    write(*,*) __LINE__
    call initEdgeBuffer(edgeAdvDSS,      nlev  )
    write(*,*) __LINE__
    call initEdgeBuffer(edgeAdvQ2 ,qsize*nlev*2)
    write(*,*) __LINE__
    call initEdgeBuffer(edgeAdvQ3 ,qsize*nlev*3)
    write(*,*) __LINE__
!$OMP MASTER

    allocate(qmin(nlev,qsize,nelemd))
    allocate(qmax(nlev,qsize,nelemd))

    write(*,*)"done cuda_mod_init"
!$OMP END MASTER
!$OMP BARRIER
  end subroutine cuda_mod_init

  !Meaningless kernel just to get initial kernel overheads out of the way.
  attributes(global) subroutine warmup(a)
    integer,value :: a
    a = 2.0 * a
  end subroutine warmup




  subroutine euler_step_cuda( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.  
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds             , only: real_kind
  use dimensions_mod    , only: np, npdg, nlev, qsize
  use hybrid_mod        , only: hybrid_t
  use element_mod       , only: element_t
  use derivative_mod    , only: derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod          , only: edgevpack, edgevunpack
  use bndry_mod         , only: bndry_exchangev
  use hybvcoord_mod     , only: hvcoord_t
  use control_mod       , only:  nu_q, nu_p, limiter_option
  use perf_mod          , only: t_startf, t_stopf  ! _EXTERNAL
  use viscosity_mod     , only: biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  ! local
  real(kind=real_kind), dimension(np,np                       ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np,np,2                     ) :: gradQ
  real(kind=real_kind), dimension(np,np,2,nlev                ) :: Vstar
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: Qtens
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: dp,dp_star
  real(kind=real_kind), dimension(np,np  ,nlev,qsize,nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=real_kind) :: dp0
  integer :: ie,q,i,j,k
  integer :: rhs_viss = 0

#if (defined ELEMENT_OPENMP)
write(*,*) 'ERROR: Do not use ELEMENT_OPENMP and CUDA FORTRAN'
stop
#endif
! call t_barrierf('sync_euler_step', hybrid%par%comm)
  call t_startf('euler_step')

  if ( DSSopt /= DSSno_var ) then
    do ie = nets , nete
      if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
      if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
      if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
      do k = 1 , nlev
        DSSvar(:,:,k) = elem(ie)%spheremp(:,:)*DSSvar(:,:,k) 
      enddo
      call edgeVpack(edgeAdvDSS,DSSvar(:,:,1:nlev),nlev,0,elem(ie)%desc)
    enddo
    call bndry_exchangeV(hybrid,edgeAdvDSS)
    do ie = nets , nete
      if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
      if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
      if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
      call edgeVunpack(edgeAdvDSS,DSSvar(:,:,1:nlev),nlev,0,elem(ie)%desc)
      do k = 1 , nlev
        DSSvar(:,:,k)=DSSvar(:,:,k)*elem(ie)%rspheremp(:,:)
      enddo
    enddo
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0
  if ( limiter_option == 8 .or. nu_p > 0 ) then
    ! for limiter=0,4 or 8 we will apply dissipation in the RHS,
    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity 
    !     
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
    !        apply dissipation to Q (not Qdp) to preserve Q=1
    !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)                
    !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    do ie = nets , nete
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
      do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
        dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k) 
        do q = 1 , qsize
          Qtens_biharmonic(:,:,k,q,ie) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp)/dp(:,:,k)
        enddo
      enddo
    enddo

    ! compute element qmin/qmax
    if ( rhs_multiplier == 0 ) then
      do ie = nets , nete
        do k = 1 , nlev    
          do q = 1 , qsize
            qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
            qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
          enddo
        enddo
      enddo
      ! update qmin/qmax based on neighbor data for lim8
      if ( limiter_option == 8 ) call neighbor_minmax(elem,hybrid,edgeAdvQ2,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    endif

    ! lets just reuse the old neighbor min/max, but update based on local data
    if ( rhs_multiplier == 1 ) then
      do ie = nets , nete
        do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
          do q = 1 , qsize
            qmin(k,q,ie)=min(qmin(k,q,ie),minval(Qtens_biharmonic(:,:,k,q,ie)))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
            qmax(k,q,ie)=max(qmax(k,q,ie),maxval(Qtens_biharmonic(:,:,k,q,ie)))
          enddo
        enddo
      enddo
    endif

    ! get niew min/max values, and also compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! compute element qmin/qmax  
      do ie = nets , nete
        do k = 1  ,nlev    
          do q = 1 , qsize
            qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
            qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
          enddo
        enddo
      enddo
      ! two scalings depending on nu_p:
      ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
      ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
      if ( nu_p > 0 ) then
        do ie = nets , nete
          do k = 1 , nlev    
            dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
            dpdiss(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%derived%psdiss_ave(:,:)
            do q = 1 , qsize
              ! NOTE: divide by dp0 since we multiply by dp0 below
              Qtens_biharmonic(:,:,k,q,ie)=Qtens_biharmonic(:,:,k,q,ie)*dpdiss(:,:)/dp0
            enddo
          enddo
        enddo
      endif
      if ( limiter_option == 8 ) then
        ! biharmonic and update neighbor min/max
        call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , deriv , edgeAdvQ3 , hybrid , nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
      else
        ! regular biharmonic, no need to updat emin/max
        call biharmonic_wk_scalar( elem , qtens_biharmonic , deriv , edgeAdv , hybrid , nets , nete )
      endif
      do ie = nets , nete
        do k = 1 , nlev    !  Loop inversion (AAM)
          dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
          do q = 1 , qsize
            ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
            qtens_biharmonic(:,:,k,q,ie) = -rhs_viss*dt*nu_q*dp0*Qtens_biharmonic(:,:,k,q,ie) / elem(ie)%spheremp(:,:)
          enddo
        enddo
      enddo
    endif
  endif  ! compute biharmonic mixing term and qmin/qmax


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie = nets , nete
    ! Compute velocity used to advance Qdp 
    do k = 1 , nlev    !  Loop index added (AAM)
      ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
      ! but that's ok because rhs_multiplier=0 on the first stage:
      dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(:,:,k) 
      Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k) / dp(:,:,k)
      Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k) / dp(:,:,k)
    enddo

    ! advance Qdp
    do q = 1 , qsize
      do k = 1 , nlev  !  dp_star used as temporary instead of divdp (AAM)
        ! div( U dp Q), 
        gradQ(:,:,1) = Vstar(:,:,1,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        gradQ(:,:,2) = Vstar(:,:,2,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        dp_star(:,:,k) = divergence_sphere( gradQ , deriv , elem(ie) )
        Qtens(:,:,k) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) - dt * dp_star(:,:,k)
        ! optionally add in hyperviscosity computed above:
        if ( rhs_viss /= 0 ) Qtens(:,:,k) = Qtens(:,:,k) + Qtens_biharmonic(:,:,k,q,ie)
      enddo
         
      if ( limiter_option == 8 ) then
        do k = 1 , nlev  ! Loop index added (AAM)
          ! UN-DSS'ed dp at timelevel n0+1:  
          dp_star(:,:,k) = dp(:,:,k) - dt * elem(ie)%derived%divdp(:,:,k)  
          if ( nu_p > 0 .and. rhs_viss /= 0 ) then
            ! add contribution from UN-DSS'ed PS dissipation
            dpdiss(:,:) = ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * elem(ie)%derived%psdiss_biharmonic(:,:)
            dp_star(:,:,k) = dp_star(:,:,k) - rhs_viss * dt * nu_q * dpdiss(:,:) / elem(ie)%spheremp(:,:)
          endif
        enddo
        ! apply limiter to Q = Qtens / dp_star 
        call limiter_optim_iter_full( Qtens(:,:,:) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , qmax(:,q,ie) , dp_star(:,:,:) )
      endif

      ! apply mass matrix, overwrite np1 with solution:
      ! dont do this earlier, since we allow np1 to be the same as n0
      ! and we dont want to overwrite n0 until we are done using it
      do k = 1 , nlev
        elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%spheremp(:,:) * Qtens(:,:,k) 
      enddo

      if ( limiter_option == 4 ) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! sign-preserving limiter, applied after mass matrix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,np1_qdp) , hvcoord ) 
      endif
    enddo

    call edgeVpack(edgeAdv , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
  enddo

  call bndry_exchangeV( hybrid , edgeAdv )

  do ie = nets , nete
    call edgeVunpack( edgeAdv , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
    do q = 1 , qsize
      do k = 1 , nlev    !  Potential loop inversion (AAM)
        elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
      enddo
    enddo
  enddo
  call t_stopf('euler_step')
end subroutine euler_step_cuda




subroutine limiter2d_zero(Q,hvcoord)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  use hybvcoord_mod     , only: hvcoord_t
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np,nlev)
  type (hvcoord_t)     , intent(in   ) :: hvcoord

  ! local
  real (kind=real_kind) :: dp(np,np)
  real (kind=real_kind) :: mass,mass_new,ml
  integer i,j,k

  do k = nlev , 1 , -1
    mass = 0
    do j = 1 , np
      do i = 1 , np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        ml = Q(i,j,k)
        mass = mass + ml
      enddo
    enddo

    ! negative mass.  so reduce all postive values to zero 
    ! then increase negative values as much as possible
    if ( mass < 0 ) Q(:,:,k) = -Q(:,:,k) 
    mass_new = 0
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          ml = Q(i,j,k)
          mass_new = mass_new + ml
        endif
      enddo
    enddo

    ! now scale the all positive values to restore mass
    if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass) / mass_new
    if ( mass     < 0 ) Q(:,:,k) = -Q(:,:,k) 
  enddo
end subroutine limiter2d_zero




subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,dpmass)
  !THIS IS A NEW VERSION OF LIM8, POTENTIALLY FASTER BECAUSE INCORPORATES KNOWLEDGE FROM
  !PREVIOUS ITERATIONS
  
  !The idea here is the following: We need to find a grid field which is closest
  !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
  !So, first we find values which do not satisfy constraints and bring these values
  !to a closest constraint. This way we introduce some mass change (addmass),
  !so, we redistribute addmass in the way that l2 error is smallest. 
  !This redistribution might violate constraints thus, we do a few iterations. 
  use kinds         , only : real_kind
  use dimensions_mod, only : np, np, nlev
  real (kind=real_kind), dimension(np,np,nlev), intent(inout)            :: ptens
  real (kind=real_kind), dimension(np,np     ), intent(in   )            :: sphweights
  real (kind=real_kind), dimension(      nlev), intent(inout)            :: minp
  real (kind=real_kind), dimension(      nlev), intent(inout)            :: maxp
  real (kind=real_kind), dimension(np,np,nlev), intent(in   ), optional  :: dpmass

  real (kind=real_kind), dimension(np,np,nlev) :: weights
  real (kind=real_kind), dimension(np,np     ) :: ptens_mass
  integer  k1, k, i, j, iter, i1, i2
  integer :: whois_neg(np*np), whois_pos(np*np), neg_counter, pos_counter
  real (kind=real_kind) :: addmass, weightssum, mass
  real (kind=real_kind) :: x(np*np),c(np*np)
  real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch
  real (kind=real_kind) :: tol_limiter = 1e-15
  integer, parameter :: maxiter = 5

  do k = 1 , nlev
    weights(:,:,k) = sphweights(:,:) * dpmass(:,:,k)
    ptens(:,:,k) = ptens(:,:,k) / dpmass(:,:,k)
  enddo

  do k = 1 , nlev
    k1 = 1
    do i = 1 , np
      do j = 1 , np
        c(k1) = weights(i,j,k)
        x(k1) = ptens(i,j,k)
        k1 = k1 + 1
      enddo
    enddo

    mass = sum(c*x)

    ! relax constraints to ensure limiter has a solution:
    ! This is only needed if runnign with the SSP CFL>1 or 
    ! due to roundoff errors
    if( (mass / sum(c)) < minp(k) ) then
      minp(k) = mass / sum(c)
    endif
    if( (mass / sum(c)) > maxp(k) ) then
      maxp(k) = mass / sum(c)
    endif

    addmass = 0.0d0
    pos_counter = 0;
    neg_counter = 0;
    
    ! apply constraints, compute change in mass caused by constraints 
    do k1 = 1 , np*np
      if ( ( x(k1) >= maxp(k) ) ) then
        addmass = addmass + ( x(k1) - maxp(k) ) * c(k1)
        x(k1) = maxp(k)
        whois_pos(k1) = -1
      else
        pos_counter = pos_counter+1;
        whois_pos(pos_counter) = k1;
      endif
      if ( ( x(k1) <= minp(k) ) ) then
        addmass = addmass - ( minp(k) - x(k1) ) * c(k1)
        x(k1) = minp(k)
        whois_neg(k1) = -1
      else
        neg_counter = neg_counter+1;
        whois_neg(neg_counter) = k1;
      endif
    enddo
    
    ! iterate to find field that satifies constraints and is l2-norm closest to original 
    weightssum = 0.0d0
    if ( addmass > 0 ) then
      do i2 = 1 , maxIter
        weightssum = 0.0
        do k1 = 1 , pos_counter
          i1 = whois_pos(k1)
          weightssum = weightssum + c(i1)
          al_pos(i1) = maxp(k) - x(i1)
        enddo
        
        if( ( pos_counter > 0 ) .and. ( addmass > tol_limiter * abs(mass) ) ) then
          do k1 = 1 , pos_counter
            i1 = whois_pos(k1)
            howmuch = addmass / weightssum
            if ( howmuch > al_pos(i1) ) then
              howmuch = al_pos(i1)
              whois_pos(k1) = -1
            endif
            addmass = addmass - howmuch * c(i1)
            weightssum = weightssum - c(i1)
            x(i1) = x(i1) + howmuch
          enddo
          !now sort whois_pos and get a new number for pos_counter
          !here neg_counter and whois_neg serve as temp vars
          neg_counter = pos_counter
          whois_neg = whois_pos
          whois_pos = -1
          pos_counter = 0
          do k1 = 1 , neg_counter
            if ( whois_neg(k1) .ne. -1 ) then
              pos_counter = pos_counter+1
              whois_pos(pos_counter) = whois_neg(k1)
            endif
          enddo
        else
          exit
        endif
      enddo
    else
       do i2 = 1 , maxIter
         weightssum = 0.0
         do k1 = 1 , neg_counter
           i1 = whois_neg(k1)
           weightssum = weightssum + c(i1)
           al_neg(i1) = x(i1) - minp(k)
         enddo
         
         if ( ( neg_counter > 0 ) .and. ( (-addmass) > tol_limiter * abs(mass) ) ) then
           do k1 = 1 , neg_counter
             i1 = whois_neg(k1)
             howmuch = -addmass / weightssum
             if ( howmuch > al_neg(i1) ) then
               howmuch = al_neg(i1)
               whois_neg(k1) = -1
             endif
             addmass = addmass + howmuch * c(i1)
             weightssum = weightssum - c(i1)
             x(i1) = x(i1) - howmuch
           enddo
           !now sort whois_pos and get a new number for pos_counter
           !here pos_counter and whois_pos serve as temp vars
           pos_counter = neg_counter
           whois_pos = whois_neg
           whois_neg = -1
           neg_counter = 0
           do k1 = 1 , pos_counter
             if ( whois_pos(k1) .ne. -1 ) then
               neg_counter = neg_counter+1
               whois_neg(neg_counter) = whois_pos(k1)
             endif
           enddo
         else
           exit
         endif
       enddo
    endif
    
    k1 = 1
    do i = 1 , np
      do j = 1 , np
        ptens(i,j,k) = x(k1)
        k1 = k1+1
      enddo
    enddo
    
 enddo
 
 do k = 1 , nlev
   ptens(:,:,k) = ptens(:,:,k) * dpmass(:,:,k)
 enddo
end subroutine limiter_optim_iter_full





#endif
end module cuda_mod


