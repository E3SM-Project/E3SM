#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_mod
  ! ------------------------------
  use kinds, only : real_kind, long_kind, int_kind
  ! ------------------------------
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2D_t, cartesian3D_t, distance
  !--------------------------
#ifdef _FVM
  use dimensions_mod, only : np, npsq, nlev, nlevp, qsize_d, max_neigh_edges, nc, nhc, ntrac_d
#else
  use dimensions_mod, only : np, npsq, nlev, nlevp, qsize_d, max_neigh_edges
#endif
  !--------------------------
  use edge_mod, only : edgedescriptor_t, rotation_t
  !--------------------------
  use gridgraph_mod, only : gridvertex_t
  !--------------------------
  implicit none
  private
  integer, public, parameter :: timelevels=3

#ifdef _PRIM

  type, public :: elem_state_t
!
!    note: variables must match that in prim_restart_mod.F90
!
! prognostic variables
!
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)   ! velocity                                            1
     real (kind=real_kind) :: T(np,np,nlev,timelevels)     ! temperature                                         2
     real (kind=real_kind) :: ps_v(np,np,timelevels)       ! surface pressure                         4
     real (kind=real_kind) :: lnps(np,np,timelevels)       ! log surface pressure                                3
     real (kind=real_kind) :: phis(np,np)         ! surface geopotential (prescribed)                            5

     ! vetically lagrangian code advects dp instead of ps_v
     real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)  ! delta p on levels

     !tracers always use 2 level time scheme: 
     real (kind=real_kind) :: Q(np,np,nlev,qsize_d)  ! Tracer concentration
     real (kind=real_kind) :: Qdp(np,np,nlev,qsize_d,2)  ! Tracer mass           


  end type elem_state_t
  ! number of variables in the state struct 
  ! used for standalone HOMME restart file
  integer(kind=int_kind),public,parameter :: StateComponents=8


  type, public :: derived_state_t


!    storage needed when subcycling tracers/dynamics 
     real (kind=real_kind) :: vn0(np,np,2,nlev) ! velocity used for SE tracer advection
                                                ! if compute_mean_flux=1:  time average of U*dp
                                                !    compute_mean_flux=0:  U at tracer time t
     real (kind=real_kind) :: vstar(np,np,2,nlev) ! velocity on Lagrangian surfaces
!     real (kind=real_kind) :: psdiss_biharmonic(np,np)  ! mean ps dissipation tendency, if nu_p>0
!     real (kind=real_kind) :: psdiss_ave(np,np)         ! mean ps used to compute psdiss_tens
     real (kind=real_kind) :: dpdiss_biharmonic(np,np,nlev)  ! mean dp dissipation tendency, if nu_p>0
     real (kind=real_kind) :: dpdiss_ave(np,np,nlev)         ! mean dp used to compute psdiss_tens

! diagnostic variables computed during explicit timestep
     real (kind=real_kind) :: phi(np,np,nlev)     ! geopotential                      
     real (kind=real_kind) :: omega_p(np,np,nlev) ! vertical tendency (derived)       
     real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! Mean vertical flux from dynamics
! Semi-Implicit diagnostic variables, computed during
! explict component, and then reused in Helmholtz component.
     real (kind=real_kind) :: grad_lnps(np,np,2)  ! gradient of log surface pressure               
     real (kind=real_kind) :: zeta(np,np,nlev)    ! relative vorticity                             
     real (kind=real_kind) :: div(np,np,nlev,timelevels)     ! divergence                          

!these fields are used internally during tracer advection for consistency and limiters
!field dp (same for tracers and dynamics at the beginning of each time step)
!is to construct dp_tracers at physics timestep
     real (kind=real_kind) :: dp(np,np,nlev)       
!divdp_proj is DSSed divdp
     real (kind=real_kind) :: divdp(np,np,nlev) 
     real (kind=real_kind) :: divdp_proj(np,np,nlev) 
!-------------------------------------------------------- 

     ! Primitive equations forcings
#ifdef CAM
     real (kind=real_kind) :: FQ(np,np,nlev,qsize_d, 1) ! F-Tracers   
     real (kind=real_kind) :: FM(np,np,2,nlev, 1)       ! F-momentum  
     real (kind=real_kind) :: FT(np,np,nlev, 1)         ! F-Temperature
     real (kind=real_kind) :: omega_prescribed(np,np,nlev) ! vertical tendency (prescribed)       
#else
     ! when does forcing ever need more than 1 timelevel?  
     real (kind=real_kind) :: FQ(np,np,nlev,qsize_d, timelevels) ! F-Tracers  
     real (kind=real_kind) :: FM(np,np,2,nlev, timelevels) ! F-momentum       
     real (kind=real_kind) :: FT(np,np,nlev, timelevels)    ! F-Temperature   
#endif
!   ^
!   |
! endif CAM 
     real (kind=real_kind) :: pecnd(np,np,nlev)         ! pressure pert. from condensate
     real (kind=real_kind) :: FQps(np,np,timelevels)         ! implied forcing of FQ on ps_v
                                                            ! used when conserving dry mass in the presence of precipitation 
  end type derived_state_t

  type, public :: elem_accum_t
!     real (kind=real_kind) :: u(np,np,nlev)       ! zonal velocity on sphere
!     real (kind=real_kind) :: T(np,np,nlev)       ! temperature
!     real (kind=real_kind) :: Q(np,np,nlev) ! tracers
!     real (kind=real_kind) :: ke(np,np,nlev)      ! kinetic energy

!
! **** ENERGY DIAGNOSTICS ****    
!
! Energy equation:   
! KE_t  = T1 + T2  + D1   + Err   +  vertical & horizontal advection terms
! IE_t  = S1 + D2                 +  vertical & horizontal advection terms
! PE_t  = S2        
!
!  KEvert*  =  KE net vertical advection    (should be zero) 
!  KEhoriz* =  KE net horizonatl advection  (should be zero)
!  IEvert*  =  IE net vertical advection    (should be zero) 
!  IEhoriz* =  IE net horizonatl advection  (should be zero)
!
! With leapfrog, energy equations are all exact except KE 
! (has an Err term that goes to zero as dt**2)
!
! Transfer terms:
! T1   = -< dp/dn u, RT_v/p grad_p >     KE<->IE:   T1 + T2-T2_s = S1
! T2   = -< dp/dn u, grad_phi >          KE<->PE:   T2_s         = S2
! T2_s = -< dp/dn u, grad_phis >
! S1   = < Cp_star dp/dn , RT omega_p/Cp_star >  
! S2   = -< div (u dp/dn), phis >                
!
#ifdef ENERGY_DIAGNOSTICS
     real (kind=real_kind) :: KEvert1(np,np)      ! term from continuity equ
     real (kind=real_kind) :: KEvert2(np,np)      ! term from momentum equ
     real (kind=real_kind) :: IEvert1(np,np)      ! term from continuity equ
     real (kind=real_kind) :: IEvert2(np,np)      ! term from T equ
     real (kind=real_kind) :: IEvert1_wet(np,np)  ! wet term from continuity equ
     real (kind=real_kind) :: IEvert2_wet(np,np)  ! wet term from T equ

     real (kind=real_kind) :: KEhorz1(np,np)      ! at time t
     real (kind=real_kind) :: KEhorz2(np,np)      ! after calling time_advance, these will be at time t-1
     real (kind=real_kind) :: IEhorz1(np,np)    
     real (kind=real_kind) :: IEhorz2(np,np)    
     real (kind=real_kind) :: IEhorz1_wet(np,np)    
     real (kind=real_kind) :: IEhorz2_wet(np,np)    

     real (kind=real_kind) :: T1(np,np)    
     real (kind=real_kind) :: T2(np,np) 
     real (kind=real_kind) :: T2_s(np,np)    
     real (kind=real_kind) :: S1(np,np)    
     real (kind=real_kind) :: S1_wet(np,np)    
     real (kind=real_kind) :: S2(np,np)    

     ! the KE conversion term and diffusion term
     real (kind=real_kind) :: DIFF(np,np,2,nlev) ! net hypervis term
     real (kind=real_kind) :: DIFFT(np,np,nlev) ! net hypervis term
     real (kind=real_kind) :: CONV(np,np,2,nlev) ! dpdn u dot CONV = T1 + T2
#endif
!   ^
!   |
! endif  ENERGY_DIAGNOSTICS

     ! the "4" timelevels represents data computed at:
     !  1  t-.5   
     !  2  t+.5   after dynamics
     !  3  t+.5   after forcing
     !  4  t+.5   after Robert
     ! after calling TimeLevelUpdate, all time above decrease by 1.0
     real (kind=real_kind) :: KEner(np,np,4)       
     real (kind=real_kind) :: PEner(np,np,4)       
     real (kind=real_kind) :: IEner(np,np,4)    
     real (kind=real_kind) :: IEner_wet(np,np,4)    
     real (kind=real_kind) :: Qvar(np,np,qsize_d,4)  ! Q variance at half time levels   
     real (kind=real_kind) :: Qmass(np,np,qsize_d,4) ! Q mass at half time levels
     real (kind=real_kind) :: Q1mass(np,np,qsize_d)  ! Q mass at full time levels

  end type elem_accum_t

#elif defined _PRIMDG
  ! definition of elem_state_t for DG Primitive Equations version of model
  type, public :: elem_state_t
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: phis(np,np)                  ! surface geopotential
     real (kind=real_kind) :: gradps(np,np,2)              ! gradient of surface geopotential
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)   ! contravarient comp
     !=======================================================================================================!
     real (kind=real_kind) :: couv(np,np,2,nlev) !covariant velocities
     real (kind=real_kind) :: uv(np,np,2,nlev)   !redundant copy of v (eliminate?)
     real (kind=real_kind) :: uv0(np,np,2,nlev)  !temp variable for velocities (eliminate?)
     real (kind=real_kind) :: pgrads(np,np,2,nlev)
     !=======================================================================================================!
     real (kind=real_kind) :: psi(np,np,nlev)
     real (kind=real_kind) :: phi(np,np,nlev)
     real (kind=real_kind) :: ht(np,np,nlev)
     real (kind=real_kind) :: T(np,np,nlev,timelevels) ! Temperature
     real (kind=real_kind) :: Q(np,np,nlev,timelevels) ! Moist
     real (kind=real_kind) :: pt3d(np,np,nlev)         ! Potential Temperature
     real (kind=real_kind) :: qt3d(np,np,nlev)
     real (kind=real_kind) :: peta(np,np,nlev)
     real (kind=real_kind) :: dp3d(np,np,nlev)
     !=======================================================================================================!
     real (kind=real_kind) :: zeta(np,np,nlev)
     real (kind=real_kind) :: pr3d(np,np,nlev+1)
     real (kind=real_kind) :: pr3d_ref(np,np,nlev+1)
     real (kind=real_kind) :: gp3d(np,np,nlev+1)
     !=======================================================================================================!
     real (kind=real_kind) :: ptop(np,np)
     real (kind=real_kind) :: sgp(np,np)
     real (kind=real_kind) :: tbar(nlev)
     !=======================================================================================================!
  end type elem_state_t

  ! this type is just a place holder at the moment.  In the future some of the fields in elem_state_t
  ! should be moved into derived_state_t to mirror the _PRIM verion
  type, public :: derived_state_t
     real (kind=real_kind) :: dummmy
     real (kind=real_kind) :: vstar(np,np,2,nlev) ! velocity on Lagrangian surfaces
  end type derived_state_t


  type, public :: elem_accum_t
    real (kind=real_kind) :: u(np,np,nlev)       ! zonal velocity on sphere
    real (kind=real_kind) :: T(np,np,nlev)       ! temperature
    real (kind=real_kind) :: ke(np,np,nlev)      ! kinetic energy
  end type

! neither _PRIM nor _PRIMDG has been defined, 
! we are in shallow water world
#else
  !
  ! SHALLOW WATER STRUCT
  !
  type, public :: elem_state_t
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: ps(np,np)           ! surface geopotential
     real (kind=real_kind) :: gradps(np,np,2)     ! gradient of surface geopotential     
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)   ! contravarient comp
     !=======================================================================================================!
     ! This are all fields added for DG (Shallow-Water)							!
     !=======================================================================================================!
#ifdef _SWDG     
     real (kind=real_kind) :: couv(np,np,2,nlev) ! covarient velocity
     real (kind=real_kind) :: uv(np,np,2,nlev)   ! spherical (u,v) velocity  vector        
     real (kind=real_kind) :: psi(np,np,nlev)    ! (ht+hs)*metdet(G)
     real (kind=real_kind) :: ht(np,np,nlev)     ! height variable     
     real (kind=real_kind) :: hs(np,np,nlev)     ! mountain height
#endif
! endif  _SWDG 
     !=======================================================================================================!
  end type elem_state_t

  type, public :: derived_state_t
     !=======================================================================================================!
     ! all shallow water diagnostic variables should be moved here.
     ! %state is for prognostic variables
     !=======================================================================================================!
     real (kind=real_kind) :: vstar(np,np,2,nlev) ! velocity on Lagrangian surfaces
  end type derived_state_t

!=======================================================================================================! 
#endif
!   ^
!   |
! endif  _PRIM, elif  _PRIMDR, else shallow water.

  type, public :: index_t
     integer(kind=int_kind) :: ia(npsq),ja(npsq)
     integer(kind=int_kind) :: is,ie
     integer(kind=int_kind) :: NumUniquePts
     integer(kind=int_kind) :: UniquePtOffset
  end type index_t

  type, public :: element_t
     integer(kind=int_kind) :: LocalId
     integer(kind=int_kind) :: GlobalId

     !=====================================
     !Add the link list hooks
     !=====================================
     type (element_t), pointer :: prev
     type (element_t), pointer :: next

     ! Coordinate values of element points
     type (spherical_polar_t) :: spherep(np,np)           ! Spherical coordinates of GLL points

     ! equ-angular gnomonic projection coordinates
     type (cartesian2D_t)     :: cartp(np,np)  ! gnomonic coordinates of GLL points 
     type (cartesian2D_t)     :: corners(4)    ! gnomonic coordinates of element corners
     real (kind=real_kind)    :: u2qmap(4,2)   ! bilinear map from ref element to quad in cubedsphere coordinates
                                               ! SHOULD BE REMOVED 
     ! 3D cartesian coordinates
     type (cartesian3D_t)     :: corners3D(4)  

     ! element diagnostics
     real (kind=real_kind)    :: area               ! Area of element
     real (kind=real_kind)    :: max_eig        ! max singular value of metinv
     real (kind=real_kind)    :: min_eig        ! min singular value of metinv
     real (kind=real_kind)    :: max_eig_ratio  ! max ratio of singular values
     real (kind=real_kind)    :: dx_short       ! short length scale
     real (kind=real_kind)    :: dx_long        ! long length scale

     real (kind=real_kind)    :: variable_hyperviscosity(np,np)  ! hyperviscosity based on above
     real (kind=real_kind)    :: hv_courant ! hyperviscosity courant number
     real (kind=real_kind)    :: tensorVisc(2,2,np,np)    !og, matrix V for tensor viscosity

     ! Edge connectivity information
     integer(kind=int_kind) :: node_numbers(4)
     integer(kind=int_kind) :: node_multiplicity(4)      ! number of elements sharing corner node

     type (GridVertex_t)       :: vertex                 ! Element grid vertex information
     type (EdgeDescriptor_t) :: desc

     type (elem_state_t)       :: state

     type (derived_state_t)    :: derived
#if defined _PRIM || defined _PRIMDG
     type (elem_accum_t)       :: accum
#endif
     ! Metric terms 
     real (kind=real_kind)    :: met(2,2,np,np)      ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metinv(2,2,np,np)   ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metdet(np,np)       ! g = SQRT(det(g_ij)) on velocity and pressure grid
     real (kind=real_kind)    :: rmetdet(np,np)      ! 1/metdet on velocity pressure grid
     real (kind=real_kind)    :: D(2,2,np,np)        ! Map covariant field on cube to vector field on the sphere
     real (kind=real_kind)    :: Dinv(2,2,np,np)     ! Map vector field on the sphere to covariant v on cube


     ! Convert vector fields from spherical to rectangular components
     ! The transpose of this operation is its pseudoinverse.
     real (kind=real_kind)    :: vec_sphere2cart(np,np,3,2)

     ! Mass matrix terms for an element on a cube face
     real (kind=real_kind)    :: mp(np,np)          ! mass matrix on  velocity and pressure grid
     real (kind=real_kind)    :: rmp(np,np)         ! inverse mass matrix on velocity and pressure grid

     ! Mass matrix terms for an element on the sphere
     ! This mass matrix is used when solving the equations in weak form
     ! with the natural (surface area of the sphere) inner product
     real (kind=real_kind)    :: spheremp(np,np)          ! mass matrix on velocity and pressure grid
     real (kind=real_kind)    :: rspheremp(np,np)         ! inverse mass matrix on velocity and pressure grid

     integer(kind=long_kind)  :: gdofP(np,np)        ! Global degree of freedom (P-grid)

     ! Coreolis term
     real (kind=real_kind)    :: fcor(np,np)        ! coreolis term

     ! Solver weights (used only for non-staggered grid
!     real (kind=real_kind)    :: solver_wts(np,np)  ! REMOVE. 
                                   ! non-conservative operatores: use spheremp*rspheremp
                                   ! conservative operators: use spheremp

     type (index_t) :: idxP
     type (index_t),pointer :: idxV
     integer :: FaceNum

     ! force element_t to be a multiple of 8 bytes.  
     ! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
     ! check core file for:
     ! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c ESR=0x01800000 CCR0=0x4800a002)
     integer :: dummy
  end type element_t


  public :: element_coordinates
  public :: element_var_coordinates
  public :: element_var_coordinates3D
  public :: GetColumnIdP,GetColumnIdV
  public :: allocate_element_desc
contains

  ! =======================================
  !  GetColumnIdP:
  !  
  !  Gives a unique identifier for a Physics 
  !  column on the P-grid
  ! =======================================
  function GetColumnIdP(elem,i,j) result(col_id)
    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id

    col_id = elem%gdofP(i,j)

  end function GetColumnIdP
  ! =======================================
  !  GetColumnIdV:
  !  
  !  Gives a unique identifier for a Physics 
  !  column on the V-grid
  ! =======================================
  function GetColumnIdV(elem,i,j) result(col_id)
    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id

    col_id = elem%gdofP(i,j)

  end function GetColumnIdV

  ! =======================================
  ! element_coordinates:
  !
  ! Initialize 2D rectilinear element 
  ! colocation points
  !
  ! =======================================

  function element_coordinates(start,end,points) result(cart)
    use kinds, only : longdouble_kind
    type (cartesian2D_t), intent(in) :: start
    type (cartesian2D_t), intent(in) :: end
    real (kind=longdouble_kind), intent(in) :: points(:)

    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))
    type (cartesian2D_t) :: length, centroid
    real (kind=longdouble_kind) :: y 
    integer i,j

    length%x   = 0.50D0*(end%x-start%x)
    length%y   = 0.50D0*(end%y-start%y)
    centroid%x = 0.50D0*(end%x+start%x) 
    centroid%y = 0.50D0*(end%y+start%y) 
    do j=1,SIZE(points)
       y = centroid%y + length%y*points(j)
       do i=1,SIZE(points)
          cart(i,j)%x = centroid%x + length%x*points(i)
          cart(i,j)%y = y
       end do
    end do

  end function element_coordinates

  function element_var_coordinates(c,points) result(cart)
    use kinds, only : longdouble_kind
    type (cartesian2D_t), intent(in) :: c(4)
    real (kind=longdouble_kind), intent(in) :: points(:)
    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))

    real (kind=longdouble_kind) :: p(size(points))
    real (kind=longdouble_kind) :: q(size(points))
    integer i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x 
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y 
       end do
    end do

  end function element_var_coordinates


  function element_var_coordinates3d(c,points) result(cart)
    use kinds, only : longdouble_kind
    type (cartesian3D_t), intent(in) :: c(4)
    real (kind=longdouble_kind), intent(in) :: points(:)
    type (cartesian3D_t) :: cart(SIZE(points),SIZE(points))

    real (kind=longdouble_kind) :: p(size(points))
    real (kind=longdouble_kind) :: q(size(points)),r
    integer i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x 
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y 
          cart(i,j)%z = p(i)*p(j)*c(1)%z &
                      + q(i)*p(j)*c(2)%z &
                      + q(i)*q(j)*c(3)%z &
                      + p(i)*q(j)*c(4)%z 

          ! project back to sphere:
          r = distance(cart(i,j))
          cart(i,j)%x = cart(i,j)%x/r
          cart(i,j)%y = cart(i,j)%y/r
          cart(i,j)%z = cart(i,j)%z/r
       end do
    end do

  end function element_var_coordinates3d


  subroutine allocate_element_desc(elem)
    type (element_t), intent(inout)   :: elem(:)
    integer                           :: num, j

    num      = SIZE(elem)
    
    do j=1,num
       allocate(elem(j)%desc%putmapP(max_neigh_edges))
       allocate(elem(j)%desc%getmapP(max_neigh_edges))
       allocate(elem(j)%desc%putmapP_ghost(max_neigh_edges))
       allocate(elem(j)%desc%getmapP_ghost(max_neigh_edges))
       allocate(elem(j)%desc%reverse(max_neigh_edges))
    end do
  end subroutine allocate_element_desc


end module element_mod
