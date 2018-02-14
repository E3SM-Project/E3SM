!
! Vertical spectral-element operators for PESE and NHSE dynamics
!_______________________________________________________________________
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module vertical_se

  use control_mod,		only: vanalytic, vtop
  use derivative_mod, only: derivative_t,  divergence_sphere, gradient_sphere, vorticity_sphere
	use dimensions_mod, only: np, nlev, nlevp
  use element_mod,		only: element_t
  use element_state,  only: elem_state_t, derived_state_t
	use hybrid_mod,			only: hybrid_t
  use hybvcoord_mod,	only: hvcoord_t, set_layer_locations
  use kinds,					only: rl => real_kind, dd => longdouble_kind
	use physical_constants, only : pi => dd_pi
	use quadrature_mod, only: quadrature_t, gausslobatto, legendre

	implicit none

  ! dimension for vertical se coordinates
  integer, parameter, public :: npv = 5         ! number of points per vertical element
!integer, parameter, public :: npv = 30         ! number of points per vertical element

  integer, parameter, public :: nev = nlev/npv  ! number of vertical elements
  integer, parameter, public :: n_unique = nev*(npv-1)+1

	!_____________________________________________________________________
	type, public :: velem_t

		! data structure for a vertical element in general vertical coordinates s

    integer :: kt, kb   ! top and bottom indices
    real(rl):: t , b    ! top and bottom coordinates

  end type

	!_____________________________________________________________________
	type, public :: solver_args

		! data structure containing arguments commonly used by solvers

		integer	  					:: np1,nm1,n0																		! time level indices
		integer  						:: qn0																					! time level used for virtual temperature
    real(rl)            :: dt                                           ! dynamics timestep
    type (hybrid_t)			:: hybrid																				! hybrid omp/mpi structure
		type (hvcoord_t)		:: hvcoord																			! hybrid vertical coordinate structure
		integer							:: nets,nete																		! start and end element indices
		type (derivative_t)	:: deriv																				! horizontal derivative data struct
    logical             :: compute_diagnostics
    real (rl)           :: eta_ave_w
	end type

  real(rl), dimension(npv,npv) :: ddn,ddn2,ddn3,ddn4                    ! derivatives w.r.t eta
  real(rl), dimension(npv,npv) :: ddn_1,ddn_n                           ! eta integration with upper bc, lower bc
  real(rl), dimension(npv,npv) :: M,LU_1,LU_n                           ! mass, LU decomposition matrices

  type (quadrature_t)	:: gll																						! gll nodes for 1 vertical element
  real(rl) :: eta_t, eta_b                                              ! position of bottom and top of the column
  real(rl) :: ds_deta, deta_ds												 								  ! metric terms for linear map
  real(rl) :: eta(nlev)                                                 ! array of all eta values in the column
  real(rl) :: ddn_hyam(nlev),  ddn_hybm(nlev)                           ! vertical derivatives of hybrid coefficients
  real(rl) :: ddn_hyai(nlevp), ddn_hybi(nlevp)                          ! vertical derivatives of hybrid coefficients

  integer  :: ipiv_1(npv),ipiv_n(npv)                                   ! pivot indices for vertical integration

  type(velem_t), allocatable :: ev(:)                                   ! array of vertical element data structures
  real(rl) :: elem_height                                               ! height of vertical elements in s coords

  real(rl), allocatable:: M_interp(:,:), L_interp(:)                    ! vertical interpolation matrix and levels

	CONTAINS

   !_____________________________________________________________________
  function pack_solver_args(np1,nm1,n0,qn0,dt,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,eta_ave_w) result(a)

    ! convenience routine for packing data into solver_args structure

    integer,							intent(in)		:: np1,nm1,n0										! time indices
    integer,							intent(in)		:: qn0													! time level used for virtual temperature
    real*8,               intent(in)    :: dt                           ! dynamics timestep
    type (hvcoord_t),			intent(in)    :: hvcoord											! hybrid vertical coord data struct
    type (hybrid_t),			intent(in)		:: hybrid												! mpi/omp data struct
    integer,							intent(in)		:: nets,nete										! start and end element indices
    type (derivative_t),	intent(in)		:: deriv												! horizontal derivative data struct
    logical,              intent(in)    :: compute_diagnostics          ! TODO: enable energy diagnostics
    real (rl),            intent(in)    :: eta_ave_w                    ! TODO: enable qsplit
    type(solver_args) :: a                                              ! resultant arg data structure

    a%np1     = np1;
    a%nm1     = nm1;
    a%n0      = n0;
    a%qn0     = qn0;
    a%dt      = dt
    a%hvcoord	= hvcoord;
    a%hybrid	= hybrid;
    a%nets    = nets;
    a%nete		= nete;
    a%deriv		= deriv;
    a%compute_diagnostics = compute_diagnostics
    a%eta_ave_w = eta_ave_w

  end function

	!_____________________________________________________________________
	real(rl) elemental function s2eta(s_coord, eta_t, eta_b)

		! map reference-element coordinate to vertical coordinate s

		real(dd), intent(in) :: s_coord                                     ! coorindate in ref element space
    real(rl), intent(in) :: eta_t,eta_b                                 ! top and bottom coord of element
    real(rl) :: deta_ds                                                 ! metric term
    deta_ds = (eta_b - eta_t)/2.0_dd
    s2eta = (s_coord + 1.0_dd)*deta_ds + eta_t

	end function

	!_____________________________________________________________________
	real(dd) elemental function eta2s(eta_coord, eta_t,eta_b)

		! map eta coordinate to reference-element coordinate

		real(rl), intent(in) :: eta_coord
    real(rl), intent(in) :: eta_t,eta_b                                 ! top and bottom coord of element
    real(rl) :: ds_deta                                                 ! metric term
    ds_deta = 2.0_dd/(eta_b-eta_t)
		eta2s = (eta_coord - eta_t)*ds_deta - 1.0_dd

	end function

  !_____________________________________________________________________
  function mass_matrix(gll,n) result(M)

    ! get diagonal mass matrix

    type (quadrature_t),	intent(in)	:: gll														! gll nodes and weights
    integer,							intent(in)  :: n															! number of points in 1d element
    real(rl) :: M(n,n)                                                  ! resultant mass matrix
    integer :: i
    M = 0
    forall(i=1:n) M(i,i) = gll%weights(i)

  end function

	!_____________________________________________________________________
  function barycentric_weights(gll,n) result(weights)

    ! compute barycentric weights from node points (kopriva2009 p.76)

    type (quadrature_t),	intent(in)	:: gll														! gll nodes and weights
    integer,							intent(in)  :: n															! number of points in 1d element
    real(rl)                          :: weights(n)                     ! result
    integer :: j,k

    weights = 1.0d0
    do j=2,n
      do k=1,j-1
        weights(k) = weights(k)*(gll%points(k)-gll%points(j))
        weights(j) = weights(j)*(gll%points(j)-gll%points(k))
      enddo
    enddo
    weights = 1.0d0/weights

  end function

  !_____________________________________________________________________
  function first_derivative_matrix(gll,n) result(D)

		! build 1st derivative matrix

    type (quadrature_t),	intent(in)	:: gll														! gll nodes and weights
		integer,							intent(in)  :: n															! number of points in 1d element
    real(rl) :: w(n)                                                    ! barycentric weights
		real(rl) :: D(n,n)																									! resultant derivative matrix
    integer  :: i,j                                                     ! loop indices

    w = barycentric_weights(gll,n)
    D = 0
    do i = 1,n
      do j = 1,n
        if(j .ne. i) then
          D(i,j) = w(j)/w(i) * 1.0d0/(gll%points(i)-gll%points(j))
          D(i,i) = D(i,i) - D(i,j)
        endif
      enddo
    enddo
  end function

  !_____________________________________________________________________
  function column_sum(f) result(f_sum)

    ! get quadrature of f over the vertical column

    real (rl), intent(in) :: f(np,np,nlev)
    real (rl) :: f_sum(np,np)
    integer		:: kt,k,l

    f_sum = 0.0d0

    do l=1,nev
      kt = ev(l)%kt
      do k=1,npv; f_sum = f_sum + deta_ds * gll%weights(k) * f(:,:,(kt-1)+k); enddo
    enddo

  end function

  !_____________________________________________________________________
  function compress(f) result(fout)

    ! remove duplicate values at vertical edge nodes

    real(rl), intent(in) :: f(nlev)
    real(rl) :: fout(n_unique)
    integer :: l,i,j, ni
    ni = npv-1 ! number of interiod nodes per element

    i=1; j=1
    do l=1,nev-1
      fout(i:i+ni-1) = f(j:j+ni-1)
      i=i+ni; j=j+npv
    enddo
    fout(i:i+npv-1)=f(j:j+npv-1)

  end function


  !_____________________________________________________________________
  function decompress(f) result(fout)

    ! restore duplicate values at vertical edge nodes

    real(rl), intent(in) :: f(n_unique)
    real(rl) :: fout(nlev)
    integer :: l,i,j
    i=1; j=1
    do l=1,nev
      fout(i:i+npv-1)=f(j:j+npv-1)
      i=i+npv; j=j+npv-1
    enddo

  end function
  !_____________________________________________________________________
  subroutine precompute_LU()

    ! pre-compute LU decomposition to make integration faster

    integer  :: info																									  ! status flag: 0=success

		! get LU decomposition of vertical integration matrix ddn
		LU_1 = ddn_1
    LU_n = ddn_n

    call DGETRF(npv,npv,LU_1,npv,ipiv_1,info)
		if(info .ne. 0) then																								! halt if routine failed
      print *,"LU_1=",LU_1
			print *,"least squares DGETRF info = ",info; stop
		endif

    call DGETRF(npv,npv,LU_n,npv,ipiv_n,info)
		if(info .ne. 0) then																								! halt if routine failed
      print *,"LU_n=",LU_n
			print *,"least squares DGETRF info = ",info; stop
		endif

  end subroutine

#if 0
	!_____________________________________________________________________
	function least_squares_integral(B,D,LU,ipiv) result(x)

		! solve D^T D x = D^T B where D is diff matrix, x is anti-deriv of 1d curve B

		real(rl), intent(in) :: B(npv)																			! 1d field to integrate
    real(rl), intent(in) :: D(npv,npv)                                  ! derivative matrix
    real(rl), intent(in):: LU(npv,npv)                                  ! LU factorization
    integer,  intent(in):: ipiv(npv)                                    ! pivot indices

    real(rl) :: x   (npv)                                               ! resulting 1d integral
		integer  :: info																									  ! status flag: 0=success
		logical  :: LU_initialized = .false.
    real(rl) :: Dtr(npv,npv)                                            !transpose of A

		! apply LU factorization to solve for definite integral function
    Dtr = transpose(D)
		x = matmul(Dtr,B)

		call DGETRS('N',npv,1,LU,npv,ipiv,x,npv,info)                       ! solve for x
		if(info .ne. 0) then																								! halt if routine failed
			print *,"integrate DGETRS info = ",info; stop
		endif

	end function
#endif

	!_____________________________________________________________________
	function solve_LU(B,LU,ipiv,NRHS) result(x)

		! LU x = B were LU is decomposed version of matrix A

		real(rl), intent(in):: B(npv,NRHS)														  					! 1d field to integrate
    real(rl), intent(in):: LU(npv,npv)                                  ! LU factorization
    integer,  intent(in):: ipiv(npv)                                    ! pivot indices
    integer,  intent(in):: NRHS                                         ! pivot indices

    real(rl) :: x(npv,NRHS)                                                  ! resulting 1d integral
		integer  :: info																									  ! status flag: 0=success

		x = B
		call DGETRS('N',npv,NRHS,LU,npv,ipiv,x,npv,info)                       ! solve for x
		if(info .ne. 0) then																								! halt if routine failed
			print *,"integrate DGETRS info = ",info; stop
		endif

	end function


  !_____________________________________________________________________
	function least_squares_integral(B,A) result(x)

    ! to solve  A x =     B using least squares
    ! where A is diff matrix, x is anti-deriv of 1d curve B
		! solve [A^T A] x = [A^T B]

		real(rl), intent(in) :: B(npv)																			! 1d field to integrate
    real(rl), intent(in) :: A(npv,npv)                                  ! derivative matrix

    integer, parameter :: M=npv, N=npv, NRHS=1, LDA=npv, LDB=npv
    integer, parameter :: NB=32, LWORK = M + M*NB


    real(rl) :: x(M), WORK(LWORK)
		integer  :: INFO																									  ! status flag: 0=success

    call DGELS('N',M,N,NRHS,A,LDA,B,LDB,WORK,LWORK,INFO)

		if(info .ne. 0) then																								! halt if routine failed
			print *,"least_squares_integral: DGELS info = ",info; stop
		endif

  end function

  !_____________________________________________________________________
  function solve(A,B,M) result (x)

    ! solve A x = B for x, using DGESV

		real(rl), intent(in) :: B(M)                                        ! 1d field to integrate
    real(rl), intent(in) :: A(M,M)                                      ! derivative matrix
    integer,  intent(in) :: M                                           ! nrows = ncols

    integer  :: IPIV(M), INFO
    real(rl) :: x(M), B_(M), A_(M,M)
    B_=B; A_=A

    !call DGESV(N, NRHS, A_, LDA, IPIV, B_, LDB, INFO)
    call  DGESV(M, 1   , A_, M  , IPIV, B_, M  , INFO)

		if(info .ne. 0) then																								! halt if routine failed
			print *,"solve: DGESV info = ",info;
		endif

    x = B_

  end function

  !_____________________________________________________________________
	function eta_integral_from_1(f, bc_1) result(x)

		! Get anti-derivative wrt vertical eta coordinates

		real (rl), intent(in) :: f(np,np,nlev)														  ! field to differentiate
    real (rl), intent(in) :: bc_1(np,np)                                ! boundary condition at n

    integer, parameter :: NRHS=np*np
		real (rl) :: x(np,np,nlev),bc(np,np)
    real (rl) :: B(npv,NRHS),R(npv,NRHS)
    integer		:: kt,kb,l,i,j

    bc = bc_1

    do l=1,nev,1
      kt = ev(l)%kt; kb = ev(l)%kb
      B           = reshape( f(:,:,kt:kb), (/npv,NRHS/), order=(/2,1/) )
      B(1,:)      = reshape(bc,(/NRHS/))
      R           = solve_LU(B,LU_1,ipiv_1,NRHS)
      x(:,:,kt:kb)= reshape( R, (/np,np,npv/), order=(/3,1,2/))
      bc = x(:,:,kb)
    enddo

	end function

  !_____________________________________________________________________
	function eta_integral_from_n(f, bc_n) result(x)

		! Get anti-derivative wrt vertical eta coordinates

		real (rl), intent(in) :: f(np,np,nlev)														  ! field to differentiate
    real (rl), intent(in) :: bc_n(np,np)                                ! boundary condition at n

    integer, parameter :: NRHS=np*np
		real (rl) :: x(np,np,nlev),bc(np,np)
    real (rl) :: B(npv,NRHS),R(npv,NRHS)
    integer		:: kt,kb,l,i,j

    bc = bc_n

    do l=nev,1,-1
      kt = ev(l)%kt; kb = ev(l)%kb
      B           = reshape( f(:,:,kt:kb), (/npv,NRHS/), order=(/2,1/) )
      B(npv,:)    = reshape(bc,(/NRHS/))
      R           = solve_LU(B,LU_n,ipiv_n,NRHS)
      x(:,:,kt:kb)= reshape( R, (/np,np,npv/), order=(/3,1,2/))
      bc = x(:,:,kt)
    enddo

	end function

  !_____________________________________________________________________
	function eta_derivative_1d(f) result(deriv)

		! get vertical derivative spanning entire column

		real (rl), intent(in) :: f(nlev)                                    ! field to differentiate
		real (rl) :: deriv(nlev)
    integer		:: kt,kb,l,i,j

    do l=1,nev
      kt = ev(l)%kt; kb = ev(l)%kb;
      deriv(kt:kb)  = matmul(ddn,f(kt:kb))
    enddo
    call vertical_dss_1d(deriv)

	end function

  !_____________________________________________________________________
	function eta_derivative(f) result(deriv)

		! get vertical derivative spanning entire column

		real (rl), intent(in) :: f(np,np,nlev)														  ! field to differentiate

    integer, parameter :: NRHS=np*np
		real (rl) :: deriv(np,np,nlev), B(npv,NRHS), R(npv,NRHS)
		integer		:: i,j,kt,kb,l

    do l=1,nev
      kt = ev(l)%kt; kb = ev(l)%kb;
      B = reshape( f(:,:,kt:kb), (/npv,NRHS/), order=(/2,1/) )
      R = matmul(ddn, B )
      deriv(:,:,kt:kb) = reshape( R, (/np,np,npv/), order=(/3,1,2/) )
    enddo
    call vertical_dss(deriv)

	end function

  !_____________________________________________________________________
	function eta_derivative_pow(f,x) result(deriv)

		! apply eta deriv matrix x times to entire column

		real (rl), intent(in) :: f(np,np,nlev)														  ! field to differentiate
    integer,   intent(in) :: x																					! exponent of the deriv

		real (rl) :: deriv(np,np,nlev)
		integer		:: i,j,m,kt,kb,l

    deriv = f

    do m=1,x
        deriv = eta_derivative(deriv)
    enddo

	end function

  !_____________________________________________________________________
  function evenly_spaced_eta_coords(ni) result(s)
    integer :: ni ! number of points
    real(rl) :: s(ni)
    integer :: i

    forall(i=2:ni-1) s(i) = eta_t + (eta_b-eta_t)/(ni-1)*(i-1)
    s(1)  = eta_t
    s(ni) = eta_b

  end function

  !_____________________________________________________________________
  function advection(f,u,v,eta_dot,a, elem) result (f_adv)

    ! compute advection terms adv(f) = v*grad(f) + eta_dot df/deta

    real(rl), intent(in) :: f(np,np,nlev)                               ! scalar field to advect
    real(rl), intent(in) :: u(np,np,nlev)                               ! horizontal velocity 1
    real(rl), intent(in) :: v(np,np,nlev)                               ! horizontal velocity 2
    real(rl), intent(in) :: eta_dot(np,np,nlev)                         ! vertical velocity component
    type(solver_args), intent(in):: a                                   ! solver args
    type (element_t),  intent(in), target :: elem                       ! element to operate upon

    real(rl) :: f_adv   (np,np,nlev)                                    ! resultant advection terms
    real(rl) :: df_deta (np,np,nlev)                                    ! vertical deriv of f
    real(rl) :: grad_f  (np,np,2,nlev)                                  ! horizontal derivs of f
    integer :: k

    ! get horizontal derivatives of f
    do k=1,nlev; grad_f(:,:,:,k) = gradient_sphere(f(:,:,k),a%deriv, elem%Dinv); enddo

    ! get vertical derivative of f
    df_deta = eta_derivative(f)

    ! compute and return advective terms
    f_adv = u*grad_f(:,:,1,:) + v*grad_f(:,:,2,:) + eta_dot * df_deta
!call vertical_dss(f_adv)

  end function

  !_____________________________________________________________________
  function advection2(f,df_deta,u,v,eta_dot,a, elem) result (f_adv)

    ! compute advection terms adv(f) = v*grad(f) + eta_dot df/deta

    real(rl), intent(in) :: f(np,np,nlev)                               ! scalar field to advect
    real(rl), intent(in) :: df_deta(np,np,nlev)                         ! scalar field to advect
    real(rl), intent(in) :: u(np,np,nlev)                               ! horizontal velocity 1
    real(rl), intent(in) :: v(np,np,nlev)                               ! horizontal velocity 2
    real(rl), intent(in) :: eta_dot(np,np,nlev)                         ! vertical velocity component
    type(solver_args), intent(in):: a                                   ! solver args
    type (element_t),  intent(in), target :: elem                       ! element to operate upon

    real(rl) :: f_adv   (np,np,nlev)                                    ! resultant advection terms
    real(rl) :: grad_f  (np,np,2,nlev)                                  ! horizontal derivs of f
    integer :: k

    ! get horizontal derivatives of f
    do k=1,nlev; grad_f(:,:,:,k) = gradient_sphere(f(:,:,k),a%deriv, elem%Dinv); enddo

    ! compute and return advective terms
    f_adv = u*grad_f(:,:,1,:) + v*grad_f(:,:,2,:) + eta_dot * df_deta
!call vertical_dss(f_adv)

  end function

    !_____________________________________________________________________
  function self_advection(u,v,eta_dot,a, elem) result (v_adv)

    ! compute self advection of horiztonal velocity
    ! v*grad(v) = 1/2 grad(v^2) + (curl u) cross u  + eta_dot dv/deta

    real(rl), intent(in) :: u(np,np,nlev)                               ! horizontal velocity 1
    real(rl), intent(in) :: v(np,np,nlev)                               ! horizontal velocity 2
    real(rl), intent(in) :: eta_dot(np,np,nlev)                         ! vertical velocity component
    type(solver_args), intent(in):: a                                   ! solver args
    type (element_t), intent(in), target :: elem                        ! element to operate upon

    real(rl) :: v_adv(np,np,2,nlev)                                     ! advection terms
    real(rl) :: vsqr (np,np,nlev)                                       ! velocity^2
    real(rl) :: vort (np,np,nlev)                                       ! vorticity

    real(rl) :: grad_vsqr(np,np,2,nlev)                                 ! gradient of velocity^2
    real(rl) :: vort_cross_v(np,np,2,nlev)                              ! vorticity cross v
    real(rl) :: v2d    (np,np,2,nlev)                                   ! (u,v) velocity vector
    real(rl) :: du_deta(np,np,nlev)                                     ! vertical deriv of u
    real(rl) :: dv_deta(np,np,nlev)                                     ! vertical deriv of v
    integer :: k

    ! get vertical derivative of u and v
    du_deta = eta_derivative(u)
    dv_deta = eta_derivative(v)

    ! get grad v^2 term
    vsqr = u*u + v*v
    do k=1,nlev; grad_vsqr(:,:,:,k) = gradient_sphere(vsqr(:,:,k),a%deriv, elem%Dinv); enddo

    ! get curl u terms
    v2d(:,:,1,:) = u; v2d(:,:,2,:)=v
    do k=1,nlev; vort(:,:,k)= vorticity_sphere(v2d(:,:,:,k),a%deriv,elem); enddo
    vort_cross_v(:,:,1,:) = -v*vort
    vort_cross_v(:,:,2,:) = +u*vort

    ! compute and return advective terms
    v_adv(:,:,1,:) = 0.5d0*grad_vsqr(:,:,1,:) + vort_cross_v(:,:,1,:) + eta_dot*du_deta
    v_adv(:,:,2,:) = 0.5d0*grad_vsqr(:,:,2,:) + vort_cross_v(:,:,2,:) + eta_dot*dv_deta
!call vertical_dss(v_adv(:,:,1,:))
!call vertical_dss(v_adv(:,:,2,:))

  end function

  !_____________________________________________________________________
  function self_advection2(u,v,du_dn,dv_dn,eta_dot,a,elem) result (v_adv)

    ! compute self advection of horiztonal velocity
    ! v*grad(v) = 1/2 grad(v^2) + (curl u) cross u  + eta_dot dv/deta

    real(rl), intent(in) :: u(np,np,nlev)                               ! u vel
    real(rl), intent(in) :: v(np,np,nlev)                               ! v vel
    real(rl), intent(in) :: du_dn(np,np,nlev)                           ! u velocity gradient
    real(rl), intent(in) :: dv_dn(np,np,nlev)                           ! v velcoity gradient
    real(rl), intent(in) :: eta_dot(np,np,nlev)                         ! vertical velocity component
    type(solver_args), intent(in):: a                                   ! solver args
    type (element_t), intent(in), target :: elem                        ! element to operate upon

    real(rl) :: v_adv(np,np,2,nlev)                                     ! advection terms
    real(rl) :: vsqr (np,np,nlev)                                       ! velocity^2
    real(rl) :: vort (np,np,nlev)                                       ! vorticity

    real(rl) :: grad_vsqr(np,np,2,nlev)                                 ! gradient of velocity^2
    real(rl) :: vort_cross_v(np,np,2,nlev)                              ! vorticity cross v
    real(rl) :: v2d    (np,np,2,nlev)                                   ! (u,v) velocity vector
    integer :: k

    ! get grad v^2 term
    vsqr = u*u + v*v
    do k=1,nlev; grad_vsqr(:,:,:,k) = gradient_sphere(vsqr(:,:,k),a%deriv, elem%Dinv); enddo

    ! get curl u terms
    v2d(:,:,1,:) = u; v2d(:,:,2,:)=v
    do k=1,nlev; vort(:,:,k)= vorticity_sphere(v2d(:,:,:,k),a%deriv,elem); enddo
    vort_cross_v(:,:,1,:) = -v*vort
    vort_cross_v(:,:,2,:) = +u*vort

    ! compute and return advective terms
    v_adv(:,:,1,:) = 0.5d0*grad_vsqr(:,:,1,:) + vort_cross_v(:,:,1,:) + eta_dot*du_dn
    v_adv(:,:,2,:) = 0.5d0*grad_vsqr(:,:,2,:) + vort_cross_v(:,:,2,:) + eta_dot*dv_dn
!call vertical_dss(v_adv(:,:,1,:))
!call vertical_dss(v_adv(:,:,2,:))

  end function

  !_____________________________________________________________________
  subroutine vertical_dss(f)

    ! Apply direct-stiffness-summation (DSS) at vertical element edges

    real(rl), intent(inout) :: f(np,np,nlev)
    integer :: kt,kb, iv
    real(rl):: f_star(np,np)

    do iv=1,nev-1

      kb = ev(iv  )%kb; ! bottom of element iv
      kt = ev(iv+1)%kt; ! top    of element iv+1

      f_star    = 0.5d0*(f(:,:,kt) + f(:,:,kb))
      f(:,:,kt) = f_star;
      f(:,:,kb) = f_star

    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine vertical_dss_1d(f)

    ! Apply direct-stiffness-summation (DSS) at vertical element edges

    real(rl), intent(inout) :: f(nlev)
    integer :: kt,kb, iv
    real(rl):: f_star

    do iv=1,nev-1

      kb = ev(iv  )%kb; ! bottom of element iv
      kt = ev(iv+1)%kt; ! top    of element iv+1

      f_star  = 0.5d0*(f(kt) + f(kb))
      f(kt)   = f_star;
      f(kb)   = f_star

    enddo

  end subroutine

  !_____________________________________________________________________
  function lagrange_polynomial(s,j,gll) result(L)

    ! get jth Lagrange interpolating-polynomial at points s
    ! L_j(s) = PROD(i/=j) (s-x_i)/(x_j-x_i)°

    real(rl),             intent(in) :: s(:)                              ! s points in range [-1,1]
    integer,              intent(in) :: j                                 ! index of polynomial
    type (quadrature_t),  intent(in) :: gll                               ! gll nodes and weights

    real(rl)  :: L(size(s))
    real(rl)  :: D
    integer   :: i,npts
    npts = size(gll%points)

    L = 1.0d0
    D = 1.0d0

    do i=1,j-1
      L = L * (s             - gll%points(i))
      D = D * (gll%points(j) - gll%points(i))
    enddo

    do i=j+1,npts
      L = L * (s             - gll%points(i))
      D = D * (gll%points(j) - gll%points(i))
    enddo

    L = L/D

  end function

  !_____________________________________________________________________________
  subroutine init_vertical_interp(levels, ni)

    ! Get vertically interpolated levels

    real(rl), intent(in) :: levels(ni)  ! vector of vertical interp coordinates
    integer,  intent(in) :: ni          ! number of interpolated levels

    ! allocate and store vertical interpolation matrix

    if(.not. allocated(M_interp)) then

      allocate(L_interp(ni))
      allocate(M_interp(ni,nlev))

      L_interp = levels
      M_interp = vertical_interp_matrix(levels,ni)
    endif

  end subroutine

  !_____________________________________________________________________________
  function v_interpolate(f,ni) result(f_i)

    real(rl), intent(in) :: f(nlev)
    integer,  intent(in) :: ni
    real(rl) :: f_i(ni)

    f_i = matmul(M_interp,f)

  end function

  !_____________________________________________________________________
  function vertical_interp_matrix(s, ni) result(M)

    ! get matrix mapping nlev values to vertical interpolation coordinates s_i
    ! note:
    real(rl), intent(in) :: s(:)    ! vector of vertical interp coordinates
    integer :: ni                   ! number of points in array s_i

    real(rl):: M(ni,nlev )          ! result, interpolation matrix

    integer :: i,j,k                  ! loop index
    integer :: v_ind(ni)            ! vertical element indices
    real(rl):: ref(ni)              ! ref coord locations
    real(rl):: ds
    real(rl):: val(1)

    ! initialize matrix to 0
    M     = 0
    v_ind =-1
    ref   =-2

    ! assign interpolation coefficients for each point s(i)
    do i=1,ni

      ! locate element and get position in reference element
      do j=1,nev
        ! check to see if point is inside vertical element j
        if( s(i)<=ev(j)%b .and. s(i)>=ev(j)%t ) then

          ! store element index
          v_ind(i) = j

          ! find reference coordinate location
          ref(i)   = eta2s(s(i), ev(j)%t, ev(j)%b)

          ! assign matrix coefficients for row(i)
          do k=1,npv
            val =lagrange_polynomial(ref(i:i),k,gll)
            M(i, ev(j)%kt-1 + k) = val(1)
          enddo

          exit ! exit loop, to make sure we assign only one element

        endif
      enddo

      ! output for debugging
      !print *,"s(i)=",s(i)," v_ind(i)=",v_ind(i)," ref(i)=",ref(i)

    enddo

   ! do j=1,nev
   !   print *, "ev(j)%b=",ev(j)%b ," ev(j)%t=", ev(j)%t
   ! enddo

    if(any(v_ind==-1)) stop 'vertical_interp_matrix: point outside vertical mesh'

  end function

 !_____________________________________________________________________
  subroutine neumann_bcs_bottom_elem(f,dfdn)

    ! Neumann bc derived without explicit surface terms

    real(rl), intent(inout) :: f(np,np,npv)
    real(rl), intent(in)    :: dfdn(np,np)

    real(rl) :: r_k
    integer  :: i,j
    logical  :: mask(npv)

    mask(:)   = .true.
    mask(npv)= .false.

    ! compute unknown value at index k
    do j=1,np
      do i=1,np

        ! get rhs of Neumann boundary condition
        r_k = dfdn(i,j) - sum( ddn(npv,:)*f(i,j,:), MASK=mask)

        ! solve for unknown value
        f(i,j,npv) = r_k/ddn(npv,npv)

      enddo
    enddo
  end subroutine

  !_____________________________________________________________________
  subroutine neumann_bcs_top_elem(f,dfdn)

    ! Neumann bc derived without explicit surface terms

    real(rl), intent(inout) :: f(np,np,npv)
    real(rl), intent(in)    :: dfdn(np,np)

    real(rl) :: r_k
    integer  :: i,j
    logical  :: mask(npv)

    mask(:)= .true.
    mask(1)= .false.

    ! compute unknown value at index k
    do j=1,np
      do i=1,np

        ! get rhs of Neumann boundary condition
        r_k = dfdn(i,j) - sum( ddn(1,:)*f(i,j,:), MASK=mask)

        ! solve for unknown value
        f(i,j,1) = r_k/ddn(1,1)

      enddo
    enddo
  end subroutine

  !_____________________________________________________________________
  subroutine neumann_bcs_bottom(f,dfdn)

    ! Neumann bc derived without explicit surface terms

    real(rl), intent(inout) :: f(np,np,nlev)
    real(rl), intent(in)    :: dfdn(np,np)
    call neumann_bcs_bottom_elem(f(:,:,1+(nlev-npv):nlev),dfdn)
call vertical_dss(f)

  end subroutine

  !_____________________________________________________________________
  subroutine neumann_bcs_top(f,dfdn)

    ! Neumann bc derived without explicit surface terms

    real(rl), intent(inout) :: f(np,np,nlev)
    real(rl), intent(in)    :: dfdn(np,np)
    call neumann_bcs_top_elem(f(:,:,1:npv),dfdn)
call vertical_dss(f)

  end subroutine
  !_____________________________________________________________________
  subroutine neumann_bcs_top_bottom(f,dfdn_t,dfdn_b)

    ! assumes nelemv > 1

    real(rl), intent(inout) :: f(np,np,nlev)
    real(rl), intent(in)    :: dfdn_t(np,np)
    real(rl), intent(in)    :: dfdn_b(np,np)

    call neumann_bcs_top   (f,dfdn_t)
    call neumann_bcs_bottom(f,dfdn_b)

  end subroutine

  !_____________________________________________________________________
  subroutine set_hvcoeffs_from_etai(hvcoord,hybrid)

    ! Set hybrid coefficients consistent with analytical etai levels

    type(hybrid_t),     intent(in)            :: hybrid
    type(hvcoord_t),    intent(inout)         :: hvcoord
    real(rl)  :: hyai,hybi, eta, tmp
    integer   :: k
    real(rl), parameter :: c = 3.0 ! 3.5 is best

    if(hybrid%par%masterproc) print *,"set_hvcoeffs_from_etai"

    ! set hybrid coefficients as suggested in dcmip test-case document
    do k=1,nlevp
      eta   = hvcoord%etai(k)

      hybi = 0.0
      tmp  = ((eta-vtop)/(1.0-vtop))
      if(tmp>0.0) hybi = tmp**c

      hyai  = eta - hybi
      hvcoord%hyai(k) = hyai
      hvcoord%hybi(k) = hybi
      if(hybrid%par%masterproc) print *,k,': etai = ',eta,' Ai = ',hyai,' Bi = ',hybi
    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine append_scalar_to_file(filename,scalar,erase)

    character*(*),    intent(in) :: filename
    real(rl),         intent(in) :: scalar
    logical,          intent(in) :: erase

    if(erase) then; open(unit=12, file=adjustl(filename),action='write',status='replace')
    else          ; open(unit=12, file=adjustl(filename),action='write',status='old',position='append')
    endif

    write(12,*) scalar
    close(unit=12)

  end subroutine

  !_____________________________________________________________________
  subroutine write_scalar_field_1d(filename,field)
    character*(*), intent(in) :: filename
    real(rl), intent(in) :: field(:)

    open(unit=12, file=adjustl(filename),action='write',status='replace')
    write(12,*) field

  end subroutine

  !_____________________________________________________________________
  subroutine write_scalar_field_2d(filename,field)
    character*(*), intent(in) :: filename
    real(rl), intent(in) :: field(:,:)
    integer i

    open(unit=12, file=adjustl(filename),action='write',status='replace')
    do i=1,size(field,1)
      write(12,*) field(i,:)
    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine write_scalar_field_3d(filename,field)
    character*(*), intent(in) :: filename
    real(rl), intent(in) :: field(:,:,:)
    integer i,j

    open(unit=12, file=adjustl(filename),action='write',status='replace')
    do i=1,size(field,1)
      do j=1,size(field,2)
          write(12,*) field (i,j,:)
      enddo
    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine make_vertical_mesh(hybrid, hvcoord)

    type (hybrid_t), intent(in)     :: hybrid                           ! mpi/omp data struct
    type (hvcoord_t),intent(inout)	:: hvcoord                          ! hybrid vertical coord data struct

    integer   :: ie, k, iv, k1, k2
    real(rl)  :: eta1, eta2
    real(rl)  :: D(npv,npv)                                           ! resultant derivative matrix

    if (hybrid%masterthread) print *,"make vertical mesh:"

    ! ensure nlev is an integer multple of npv
    !nev = nlev/npv
    if( MOD(nlev,npv) /= 0) then
      print *,"error: nlev=",nlev,"must be a multple of npv=",npv
      stop
    endif
    if (hybrid%masterthread) print *,"nev=",nev," npv=",npv

    ! set extent of vertical domain
    eta_b = 1.0d0             ! set bottom coord
    eta_t = hvcoord%etam(1)   ! set top coord from test or file

    ! if analytic vcoords get vtop from namelist
    if(vanalytic==1) eta_t = vtop
    if(hybrid%masterthread)  print *, "eta_b=",eta_b, "eta_t=",eta_t

    ! store gll nodes and weights
    gll = gausslobatto(npv)

    ! init linear map from s to eta coordinates
    elem_height = (eta_b - eta_t)/nev
    ds_deta     = 2.0_dd/elem_height
    deta_ds     = 1.0_dd/ds_deta
    if (hybrid%masterthread)  print *, "elem_height = ",elem_height

    ! set matrix operators
    D     = ds_deta*first_derivative_matrix(gll,npv)
    ddn   = D                       ! get 1st deriv matrix
    ddn2  = matmul(ddn,ddn)         ! get 2nd deriv matrix
    ddn3  = matmul(ddn,ddn2)        ! get 3rd deriv matrix
    ddn4  = matmul(ddn,ddn3)        ! get 4th deriv matrix
    M     = mass_matrix(gll,npv)

    ddn_1 = ddn; ddn_1(1  ,:)=0.0_rl; ddn_1(1  ,1  )=1.0_rl
    ddn_n = ddn; ddn_n(npv,:)=0.0_rl; ddn_n(npv,npv)=1.0_rl

    ! store LU decomposition needed for least-squares integration
    call precompute_LU()

    ! allocate vertical element array
    allocate ( ev(nev) )

    ! construct vertical mesh
    do iv=1,nev

      eta1 = eta_t + elem_height*(iv-1)   ! get eta at element top
      eta2 = eta_b - elem_height*(nev-iv)     ! get eta at element bottom
      k1   = 1  + npv*(iv-1)             ! get index at element top
      k2   = k1 + npv-1                  ! get index at element bottom

      ! store coordinates and indices in a data structure
      ev(iv)%t  = eta1;  ev(iv)%b  = eta2;
      ev(iv)%kt = k1;    ev(iv)%kb = k2;

      ! set eta levels at gll node locations
      eta(k1:k2) = s2eta(gll%points, eta1, eta2)
      if (hybrid%masterthread)  print *,"eta1 = ",eta1," eta2 = ",eta2

    enddo

    ! store vertical coords in hvcoord struct
    hvcoord%etai(2:nlev+1)= eta
    hvcoord%etai(1)       = eta(1)/2.0d0
    call set_hvcoeffs_from_etai(hvcoord,hybrid)

    hvcoord%etam  = eta
    hvcoord%hyam  = hvcoord%hyai(2:nlev+1)
    hvcoord%hybm  = hvcoord%hybi(2:nlev+1)
    call set_layer_locations(hvcoord, .true., hybrid%masterthread)

    ! compute vertical derivatives of hybrid coefficients 
    ddn_hyam = eta_derivative_1d(hvcoord%hyam)
    ddn_hybm = eta_derivative_1d(hvcoord%hybm)
    ddn_hyai(2:nlevp) = ddn_hyam; ddn_hyai(1)=0
    ddn_hybi(2:nlevp) = ddn_hybm; ddn_hybi(1)=0

    call test_vertical_operators(hybrid)

  end subroutine

  !_____________________________________________________________________
	subroutine test_vertical_operators(hybrid)

      integer, parameter :: ni = 200 ! number of interpolation pts
      logical, parameter :: write_files = .false.

			type (hybrid_t), intent(in)	:: hybrid

			real(rl), dimension(np,np,nlev) :: f,f2,f3,f4, f_deriv,f_integral, f_integral2

			real(rl), dimension(np,np,nlev) :: err_deriv, err_integral, err_integral2
			real(rl), dimension(np,np,nlev) :: comp, comp2, err_comp, err_comp2
			real(rl), dimension(np,np,nlev) :: c,s,d1,d2,d3,d4
      real(rl), dimension(np,np)	    :: zeros

      real(rl), dimension(nlev)::f1     ! 1d function
      real(rl), dimension(ni  )::f1_i   ! 1d interpolation function

      real(rl), dimension(ni)         :: s_interp
      real(rl), dimension(ni,nlev)    :: M_interp

			real(rl) :: max_err_deriv, max_err_integral,max_err_integral2, max_err_comp, max_err_comp2
			real(rl) :: e1,e2,e3,e4
			real(rl) :: frq
			integer :: i,j,k,l

      if (hybrid%masterthread) print *,"test_vertical_operators"
      zeros = 0.0d0

      ! let f = e^eta * something
      do j=1,np; do i=1,np; do k=1,nlev;
         f(i,j,k) = exp(eta(k))*cos(1.0d0*i/np)*cos(1.0d0*j/np)
        !  f(i,j,k) = cos(5.0d0*eta(k))
      enddo; enddo; enddo;

			! compute derivative of f wrt eta
      f_deriv       = eta_derivative(f)
			err_deriv			= f_deriv - f
			max_err_deriv = maxval(abs(err_deriv))

			! compute integral of f wrt eta from top
      f_integral      = eta_integral_from_1(f, f(:,:,1))
      err_integral		= f_integral - f

			max_err_integral= maxval(abs(err_integral))

      ! get interpolated vertical values
      s_interp = evenly_spaced_eta_coords(ni)
      M_interp = vertical_interp_matrix(s_interp, ni)
      do k=1,nlev; f1(k)=cos(7.0*pi*eta(k)); enddo
      f1_i = matmul(M_interp,f1)


      ! compute integral of f wrt eta from bottom
      f_integral2      = eta_integral_from_n(f, f(:,:,nlev))
      err_integral2		 = f_integral2 - f
      max_err_integral2= maxval(abs(err_integral2))

			! compose deriv and inverse matrices
      comp      = eta_integral_from_1( f_deriv, f(:,:,1))                ! integral of deriv
      comp2			= eta_derivative( f_integral )                     ! deriv of integral
			err_comp  = comp  - f
			err_comp2 = comp2 - f

			max_err_comp  = maxval(abs(err_comp))
			max_err_comp2 = maxval(abs(err_comp2))

      !if (hybrid%masterthread) then
     ! 		print *, "f(1,1,:)           =", f(1,1,:)
     !     print *, "f_integral(1,1,:)  =", f_integral(1,1,:)
     ! 		print *, "f_deriv(1,1,:)     =", f_deriv(1,1,:)
     ! endif

     if(write_files) then
        call write_scalar_field_1d("eta.txt",         eta(:))
        call write_scalar_field_1d("f.txt",           f(1,1,:))
        call write_scalar_field_1d("f_deriv.txt",     f_deriv(1,1,:))
        call write_scalar_field_1d("f_integral.txt",  f_integral(2,2,:))
        call write_scalar_field_1d("f_integral2.txt", f_integral2(1,1,:))
        call write_scalar_field_1d("comp.txt",        comp(1,1,:))
        call write_scalar_field_1d("comp2.txt",       comp2(1,1,:))
        call write_scalar_field_1d("f1.txt",          f1)
        call write_scalar_field_1d("f1_i.txt",        f1_i)
        call write_scalar_field_1d("s_i.txt",         s_interp)
      endif

			if (hybrid%masterthread) then
				print *,"max error, vertical deriv     of exp(eta)=",max_err_deriv
				print *,"max error, vertical integral  of exp(eta)=",max_err_integral
        print *,"max error, vertical integral2 of exp(eta)=",max_err_integral2
				print *,"max error, integral of derivative       =",max_err_comp
				print *,"max error, derivative of integral       =",max_err_comp2
			endif

			! test higher vertical derivatives
			frq = 1.0_rl
			forall(k=1:nlev) s(:,:,k) = sin(frq*eta(k))
			forall(k=1:nlev) c(:,:,k) = cos(frq*eta(k))

			! compute 1st, 2nd, 3rd, and 4th deriv of sin(eta) wrt eta
			d1 = eta_derivative_pow(s,1)
			d2 = eta_derivative_pow(s,2)
			d3 = eta_derivative_pow(s,3)
			d4 = eta_derivative_pow(s,4)

			! compute max error in higher vertical derivs
			e1 = maxval(abs( d1 - frq**1 * c))
			e2 = maxval(abs( d2 + frq**2 * s))
			e3 = maxval(abs( d3 + frq**3 * c))
			e4 = maxval(abs( d4 - frq**4 * s))

			if (hybrid%masterthread) then
				print *,"max error, D   cos'(eta)=",e1
				print *,"max error, D^2 cos'(eta)=",e2
				print *,"max error, D^3 cos'(eta)=",e3
				print *,"max error, D^4 cos'(eta)=",e4
			endif

      ! test solution of Neumann boundary conditions
      f2 = f; f3 = f; f4 = f

      ! compute neumann bcs at top and bottom of the column
      !call neumann_bcs_top_bottom(f2,f_deriv(:,:,1),f_deriv(:,:,nlev))

			!if (hybrid%masterthread) then
      !  print *,"Check neumann_bcs_top_bottom"
      !  print *,"f2_t - f_t =",f2(:,:,1)-f(:,:,1)
      !  print *,"f2_b - f_b =",f2(:,:,nlev)-f(:,:,nlev)
			!endif

      ! apply neuman bc at top of column
      call neumann_bcs_top(f3,f_deriv(:,:,1))

      ! apply neumann bc at bottom of column
      call neumann_bcs_bottom(f4,f_deriv(:,:,nlev))

      if (hybrid%masterthread) then
        print *,"Check neumann_bcs at index 1 and n"
        print *,"max error f3 - f =",maxval(abs(f3(:,:,1)-f(:,:,1)))
        print *,"max error f4 - f =",maxval(abs(f4(:,:,nlev)-f(:,:,nlev) ))
			endif

	end subroutine

 end module

