!--------------------------------------------------
!ice3d_lib.F90
!This file is adapted from Frank Pattyn's Ice3D model
!Code updated to Fortran 90 and integrated with Glimmer
!by Tim Bocek and Jesse Johnson
!See F. Pattyn, "A new three-dimensional higher-order thermomechanical
!        ice sheet model: basic sensitivity, ice stream development,
!        and ice flow across subglacial lakes", Journal of Geophysical
!        Research, Volume 108, no. B8, 2003.
! and
! Tim Bocek's thesis
!----------------------------------------------------

#include "glide_nan.inc"
#include "glide_mask.inc"

!The following defines contain debug options that aid in examining the sparse
!matrix.  They should be disabled unless the higher-order code needs debugging.

!Define to output a NetCDF file of the partial iterations
#define OUTPUT_PARTIAL_ITERATIONS
!#define VERY_VERBOSE
!#define NOSHELF
!#define ENFORCE_PBC

!If defined, a text field contianing debug output will be written to the disk in the
!directory in which this file was run.
!#define DEBUG_FIELDS

module ice3d_lib
    use glimmer_global
    use glimmer_physcon, only: pi, grav, rhoi, rhoo, scyr
    use glide_deriv
    use glide_types
    use glide_grids
    use glimmer_sparse_type
    use glimmer_sparse
    use glimmer_log
    use xls
    use glide_nonlin

    implicit none
    real(dp), parameter :: SMALL=1.D-10, ZIP=1.D-30, BIG=1.D10
    integer, parameter :: MAXITER = 1500              ! For the non-linear iteration
    real(dp), parameter :: toler_adjust_factor=1.0

    real(dp) :: plastic_bed_regularization = 1e-2

#ifdef DEBUG_FIELDS
    !This is a debugging field to output what computation is done at what point.
    integer, dimension(:,:), allocatable :: pointtype
#endif

#ifdef VERY_VERBOSE
    logical, parameter :: sparverbose = .true.
#else
    logical, parameter :: sparverbose = .false.
#endif

!------------------------------------------------------
!   lookup coordinates in sparse matrix
!
!   Reference:
!     1   ->   i - 1, j - 1, k     ->   IM1_JM1_K
!x    2   ->   i - 1, j, k - 1
!x    3   ->   i - 1, j, k
!x    4   ->   i - 1, j, k + 1
!     5   ->   i - 1, j + 1, k
!x    6   ->   i, j - 1, k - 1
!x    7   ->   i, j - 1, k
!x    8   ->   i, j - 1, k + 1
!x    9   ->   i, j, k - 2
!x   10   ->   i, j, k - 1
!x   11   ->   i, j, k
!x   12   ->   i, j, k + 1
!x   13   ->   i, j, k + 2
!x   14   ->   i, j + 1, k - 1
!x   15   ->   i, j + 1, k
!x   16   ->   i, j + 1, k + 1
!    17   ->   i + 1, j - 1, k
!x   18   ->   i + 1, j, k - 1
!x   19   ->   i + 1, j, k
!x   20   ->   i + 1, j, k + 1
!    21   ->   i + 1, j + 1, k
!    22   ->   i - 2, j , k
!    23   ->   i + 2, j , k
!    24   ->   i , j - 2, k
!    25   ->   i , j + 2, k
!------------------------------------------------------
    integer, parameter :: STENCIL_SIZE = 25
    integer, parameter :: IM1_JM1_K = 1
    integer, parameter :: IM1_J_KM1 = 2
    integer, parameter :: IM1_J_K = 3
    integer, parameter :: IM1_J_KP1 = 4
    integer, parameter :: IM1_JP1_K = 5
    integer, parameter :: I_JM1_KM1 = 6
    integer, parameter :: I_JM1_K = 7
    integer, parameter :: I_JM1_KP1 = 8
    integer, parameter :: I_J_KM2 = 9
    integer, parameter :: I_J_KM1 = 10
    integer, parameter :: I_J_K = 11
    integer, parameter :: I_J_KP1 = 12
    integer, parameter :: I_J_KP2 = 13
    integer, parameter :: I_JP1_KM1 = 14
    integer, parameter :: I_JP1_K = 15
    integer, parameter :: I_JP1_KP1 = 16
    integer, parameter :: IP1_JM1_K = 17
    integer, parameter :: IP1_J_KM1 = 18
    integer, parameter :: IP1_J_K = 19
    integer, parameter :: IP1_J_KP1 = 20
    integer, parameter :: IP1_JP1_K = 21
    ! These are for second order up or downwinding
    integer, parameter :: IM2_J_K = 22
    integer, parameter :: IP2_J_K = 23
    integer, parameter :: I_JM2_K = 24
    integer, parameter :: I_JP2_K = 25

contains
    !Initialize rescaled coordinate coefficients
    subroutine init_rescaled_coordinates(dhdx,dhbdx,dhdy,dhbdy,surf,h,hb,&
               dzdx,dzdy,d2zdx2, d2zdy2, d2hdx2, d2hdy2, zeta,ax,ay,bx,by,cxy,dx,dy, &
               direction_x, direction_y)
!
        INTEGER MAXX,MAXY,NZETA
        real(dp), dimension(:,:), intent(in) :: dhdx
        real(dp), dimension(:,:), intent(in) :: dhbdx
        real(dp), dimension(:,:), intent(in) :: dhdy
        real(dp), dimension(:,:), intent(in) :: dhbdy
        real(dp), dimension(:,:), intent(in) :: surf
        real(dp), dimension(:,:), intent(in) :: h
        real(dp), dimension(:,:), intent(in) :: hb
        real(dp), dimension(:,:), intent(in) :: dzdx
        real(dp), dimension(:,:), intent(in) :: dzdy
        real(dp), dimension(:),   intent(in) :: zeta
        real(dp), dimension(:,:,:), intent(inout) :: ax
        real(dp), dimension(:,:,:), intent(inout) :: ay
        real(dp), dimension(:,:,:), intent(inout) :: bx
        real(dp), dimension(:,:,:), intent(inout) :: by
        real(dp), dimension(:,:,:), intent(inout) :: cxy
        
        real(dp), intent(in) :: dx
        real(dp), intent(in) :: dy
        
        real(dp), dimension(:,:), intent(in) :: d2zdx2,d2zdy2,d2hdx2,d2hdy2
        
        real(dp), dimension(:,:), intent(in) :: direction_x, direction_y

        real(dp) :: d2z_dxy, d2h_dxy

        INTEGER :: i,j,k

        !jvj checked for consistency. All calls to this function are made from
        !the same grids


        !The Pattyn model currently does not support CISM's rescaling of
        !variables, see glimmer_paramets.F90
#ifndef NO_RESCALE
        write(*,*) "Pattyn-Bocek model does not work with scaling yet.  Please re-compile with the flag -DNO_RESCALE."
        stop
#endif

#ifdef DEBUG_FIELDS
        call write_xls("h.txt",h)
        call write_xls("hb.txt",hb)
        call write_xls("surf.txt",surf)
        call write_xls("dhdx.txt",dhdx)
        call write_xls("dhdy.txt",dhdy)
        call write_xls("dhbdx.txt",dhbdx)
        call write_xls("dhbdy.txt",dhbdy)
        call write_xls("dzdx.txt",dzdx)
        call write_xls("dzdy.txt",dzdy)
        call write_xls("d2zdx2.txt",d2zdx2)
        call write_xls("d2zdy2.txt",d2zdy2)
        call write_xls("d2hdx2.txt",d2hdx2)
        call write_xls("d2hdy2.txt",d2hdy2)
#endif

        MAXY = size(h, 2)
        MAXX = size(h, 1)
        NZETA = size(zeta, 1)

        !Get a field of 2nd derivatives of the surface onto a staggered
        !grid
        !Determination of scaling factors ax, ay, bx, by, cxy
        !These are given by (40) through (42) in Pattyn 2002.
        !At boundaries or when ice thickness (h) is 0, these
        !factors are set to 0.
        do i = 1, MAXX
            do j = 1, MAXY
                if ( (i > 1) .and. (i < MAXX) .and. (j > 1) .and. (j < MAXY)&
                     .and. (h(i,j) > 0.)) then
                    do k = 1, NZETA
                        ax(k,i,j) = (dzdx(i, j) - zeta(k) * dhdx(i, j)) / h(i, j)
                        ay(k,i,j) = (dzdy(i, j) - zeta(k) * dhdy(i, j)) / h(i, j)
              
                        bx(k,i,j) = ( d2zdx2(i, j) - &
                                      zeta(k) * d2hdx2(i, j) - & 
                                      2. * ax(k,i,j) * dhdx(i, j) ) / h(i,j)
                         
                        by(k,i,j) = ( d2zdy2(i, j) - &
                                      zeta(k) * d2hdy2(i, j) - & 
                                      2. * ay(k,i,j) * dhdy(i, j) ) / h(i,j)
                        
                        !cxy requires cross derivatives of z and h!
                        !We compute these here by differentiating dzdy and dhdy w.r.t. x,
                        !as this lets us take upwinding needs into account
                        if (direction_x(i,j) > 0) then
                            d2z_dxy = dfdx_2d_downwind(dzdy, i, j, dx)
                            d2h_dxy = dfdx_2d_downwind(dhdy, i, j, dx)
                        else if (direction_x(i,j) < 0) then
                            d2z_dxy = dfdx_2d_upwind(dzdy, i, j, dx)
                            d2h_dxy = dfdx_2d_upwind(dhdy, i, j, dx)
                        else
                            d2z_dxy = dfdx_2d(dzdy, i, j, dx)
                            d2h_dxy = dfdx_2d(dhdy, i, j, dx)
                        end if 
                        cxy(k,i,j) = ( d2z_dxy - & 
                                      zeta(k) * d2h_dxy - &
                                      ax(k,i,j)*dhdy(i, j) - &
                                      ay(k,i,j)*dhdx(i, j) ) / h(i,j)
                         
                    end do
                else
                    ax(:,i,j)=0.
                    ay(:,i,j)=0.
                    bx(:,i,j)=0.
                    by(:,i,j)=0.
                    cxy(:,i,j)=0.
                endif
            end do
        end do
#ifdef DEBUG_FIELDS
        call write_xls_3d("ax.txt",ax)
        call write_xls_3d("ay.txt",ay)
        call write_xls_3d("bx.txt",bx)
        call write_xls_3d("by.txt",by)
        call write_xls_3d("cxy.txt",cxy)
#endif
    end subroutine
!
!
!------------------------------------------------------
!   Velocity estimation according to 0th order model (SIA)
!------------------------------------------------------
!
      SUBROUTINE veloc1(dzdx,dzdy,h,flwa,zeta,uvel,vvel,u,v,&
        MAXX,MAXY,NZETA,FLOWN,PERIODIC_X,PERIODIC_Y)
!
        INTEGER MAXX,MAXY,NZETA
        real(dp) :: dzdx(:,:)
        real(dp) :: dzdy(:,:),h(:,:),zeta(:)
        real(dp) :: flwa(:,:,:),uvel(:,:,:)
        real(dp) :: vvel(:,:,:),u(:,:)
        real(dp) :: v(:,:),FLOWN
        logical :: PERIODIC_X, PERIODIC_Y
!
      INTEGER i,j,k
      real(dp) :: diff1(NZETA),d,grad,z,diffus,us,vs
!
      diff1 = 0
      
      do 10 i=2,MAXX-1
        do 20 j=2,MAXY-1
          grad=(dzdx(i,j)**2)+(dzdy(i,j)**2)
          us=0.
          vs=0.
          d = (RHOI*GRAV*h(i,j))**FLOWN
          !Accumulated vertical integration using trapezoid rule??
          do 30 k=(NZETA-1),1,-1
            z = (0.5*(zeta(k+1)+zeta(k)))**FLOWN
            diff1(k)=d*grad*h(i,j)*(flwa(k,i,j+1)+flwa(k,i,j))*z*&
              (zeta(k)-zeta(k+1))+diff1(k+1)
   30     CONTINUE
          diffus=0.
          !Vertical averaging of diffusivity (this works because the
          !differences between successive vertical layers must sum to 1!)
          do 40 k=2,NZETA
            diffus=diffus+0.5*(diff1(k)+diff1(k-1))*&
              (zeta(k)-zeta(k-1))
   40     CONTINUE
          !Estimation of velocity from SIA diffusivity
          do 50 k=1,NZETA
            uvel(k,i,j)=diff1(k)*dzdx(i,j)+us
            vvel(k,i,j)=diff1(k)*dzdy(i,j)+vs
   50     CONTINUE
        !Estimation of vertical average velocity from vertical average diff.
        u(i,j)=diffus*dzdx(i,j)
        v(i,j)=diffus*dzdy(i,j)
   20   CONTINUE
   10 CONTINUE

    END subroutine
     !
!
!------------------------------------------------------
!   Velocity estimation according to higher order model
!------------------------------------------------------
!
    SUBROUTINE veloc2(efvs,uvel,vvel,flwa,dzdx,dzdy,h,ax,ay,&
                zeta,bx,by,cxy,beta,&
                dhbdx,dhbdy,FLOWN,ZIP,VEL2ERR,&
                TOLER, options, STAGGERED, delta_x, delta_y, &
                point_mask, active_points, geometry_mask, kinematic_bc_u, kinematic_bc_v, &
                marine_bc_normal)
                

        INTEGER MAXY,MAXX,NZETA,MANIFOLD

        real(dp), dimension(:,:,:), intent(inout) :: efvs
        real(dp), dimension(:,:,:), intent(inout) :: uvel
        real(dp), dimension(:,:,:), intent(inout) :: vvel
        real(dp), dimension(:,:,:), intent(in) :: flwa
        real(dp), dimension(:,:), intent(in) :: dzdx
        real(dp), dimension(:,:), intent(in) :: dzdy
        real(dp), dimension(:,:), intent(in) :: h
        real(dp), dimension(:,:,:), intent(in) :: ax
        real(dp), dimension(:,:,:), intent(in) :: ay
        real(dp), dimension(:,:,:), intent(in) :: bx
        real(dp), dimension(:,:,:), intent(in) :: by
        real(dp), dimension(:,:,:), intent(in) :: cxy
        real(dp), dimension(:), intent(in) :: zeta
        real(dp), dimension(:,:), intent(in) :: dhbdx
        real(dp), dimension(:,:), intent(in) :: dhbdy
        real(dp), dimension(:,:), intent(in) :: beta
        type(glide_options), intent(in) :: options 
        real(dp), intent(in) :: FLOWN,ZIP,VEL2ERR,TOLER,delta_x, delta_y
        logical, intent(in) :: STAGGERED !Whether the model is run on a staggered grid or a colocated grid.
                                         !This is passed so that corrections arising from averaging can be made.
        
        integer, dimension(:,:), intent(in) :: point_mask 
        !*FD Contains a unique nonzero number on each point of ice that should be computed or that is
        !*FD ajdacent to a point that should be computed. Other numbers contain
        !*FD zeros
        integer, intent(in) :: active_points !*FD The number of points that should be computed.  This is equivalent to the number of nonzero entries in point_mask.

        integer, intent(in), dimension(:,:) :: geometry_mask !*FD point_mask field as described in glide_mask.inc

        !Contains NaN everywhere except where a kinematic boundary is to be
        !applied, in which case contains the value at the boundary
        real(dp), dimension(:,:,:), intent(inout) :: kinematic_bc_u, kinematic_bc_v
        
        !Contains NaN everywhere except on the marine ice edge of an ice shelf,
        !where this should contain the angle of the normal to the marine edge,
        !in radians, with 0=12 o'clock, pi/2=3 o'clock, etc.
        real(dp), dimension(:,:), intent(in)   :: marine_bc_normal

        INTEGER l,lacc,m,iter,DU1,DU2,DV1,DV2,ijktot
        real(dp) :: error,tot,teta

        PARAMETER (DU1=1,DU2=2,DV1=3,DV2=4)

        real(dp), dimension(4,size(efvs,1)*size(efvs,2)*size(efvs,3)) :: em

        real(dp), dimension(2*size(efvs,1)*size(efvs,2)*size(efvs,3)) :: correction_vec
        !Velocity estimates computed for the *current* iteration.  uvel and
        !vvel, comparitively, hold the velocity estimates for the *last*
        !iteration.
        real(dp), dimension(size(efvs,1),size(efvs,2),size(efvs,3)) :: ustar, vstar
        real(dp), dimension(size(efvs,2),size(efvs,3)) :: tau !Basal traction, to be computed from the provided beta parameter
       
        real(dp), dimension(size(efvs,1),size(efvs,2),size(efvs,3)) :: dudx, dudy, dudz
        real(dp), dimension(size(efvs,1),size(efvs,2),size(efvs,3)) :: dvdx, dvdy, dvdz
        real(dp), dimension(size(efvs,2),size(efvs,3)) :: direction_x, direction_y

        logical :: cont = .true.

        type(sparse_matrix_type) :: matrix_u
        type(sparse_solver_workspace) :: matrix_workspace_u
        type(sparse_solver_options) :: matrix_options_u

        type(sparse_matrix_type) :: matrix_v
        type(sparse_solver_workspace) :: matrix_workspace_v
        type(sparse_solver_options) :: matrix_options_v


#ifdef OUTPUT_PARTIAL_ITERATIONS 
        integer :: ncid_debug
#endif
        !For timing the algorithm
        real(dp) :: solve_start_time, solve_end_time, iter_start_time, iter_end_time
        real(dp) :: max_vel 

        call cpu_time(solve_start_time)

#ifdef DEBUG_FIELDS
        allocate(pointtype(size(efvs, 2), size(efvs, 3)))
#endif

        !Get the sizes from the efvs field.  Calling code must make sure that
        !these all agree.
        maxx = size(efvs,2)
        maxy = size(efvs,3)
        nzeta = size(efvs,1)
        ijktot = active_points*nzeta

        error=VEL2ERR
        m=1
        l=1
        lacc=0
        em = 0
        correction_vec = 0
        tot = 0
        !Set up sparse matrix options
        call sparse_solver_default_options(options%which_ho_sparse, matrix_options_u)
        call sparse_solver_default_options(options%which_ho_sparse, matrix_options_v)
        matrix_options_u%base%tolerance=TOLER
        matrix_options_u%base%maxiters  = 500
        matrix_options_v%base%tolerance=TOLER
        matrix_options_v%base%maxiters  = 500

        !Create the sparse matrix
        call new_sparse_matrix(ijktot, ijktot*STENCIL_SIZE, matrix_u)
        call sparse_allocate_workspace(matrix_u, matrix_options_u, matrix_workspace_u, ijktot*STENCIL_SIZE)

        call new_sparse_matrix(ijktot, ijktot*STENCIL_SIZE, matrix_v)
        call sparse_allocate_workspace(matrix_v, matrix_options_v, matrix_workspace_v, ijktot*STENCIL_SIZE)


        call write_xls("h.txt",h)
        call write_xls_3d("flwa.txt",flwa)

#ifdef DEBUG_FIELDS
        call write_xls("dzdx.txt",dzdx)
        call write_xls("dzdy.txt",dzdy)
        call write_xls_3d("ax.txt",ax)
        call write_xls_3d("ay.txt",ay)
        call write_xls_3d("bx.txt",bx)
        call write_xls_3d("by.txt",by)
        call write_xls_3d("cxy.txt",cxy)
        call write_xls_3d("uvel_sia.txt",uvel)
        call write_xls_3d("vvel_sia.txt",vvel)
        call write_xls("beta.txt",beta)
        call write_xls("dhbdx.txt",dhbdx)
        call write_xls("dhbdy.txt",dhbdy)
        call write_xls("latbc.txt",marine_bc_normal)
        call write_xls_int("geometry_mask.txt",geometry_mask)
        call write_xls_int("point_mask.txt",point_mask)
        call write_xls_3d("kinematic_bc_u.txt",kinematic_bc_u)
        call write_xls_3d("kinematic_bc_v.txt",kinematic_bc_v)
        call write_xls("marine_bc_norms.txt",marine_bc_normal)
        call write_xls_direction_guide("direction_guide.txt",3,3)
        call write_xls("normal_x.txt", sin(marine_bc_normal))
        call write_xls("normal_y.txt", -cos(marine_bc_normal))
        write(*,*) "ZETA=",zeta
#endif

        !Copy the velocity estimate from the previous iteration as the current velocity
        ustar=uvel
        vstar=vvel

        !Deal with periodic boundary conditions
        call periodic_boundaries_3d(uvel,options%periodic_ew,options%periodic_ns)
        call periodic_boundaries_3d(vvel,options%periodic_ew,options%periodic_ns)
        !Compute basal traction

#ifdef OUTPUT_PARTIAL_ITERATIONS        
        ncid_debug = begin_iteration_debug(maxx, maxy, nzeta)
        call iteration_debug_step(ncid_debug, 0, efvs, uvel, vvel, geometry_mask, marine_bc_normal)
#endif


        !  =========================================
        !  Non-linear iteration on velocities
        !  =========================================

        print *, "Entering the non-linear iteration on velocities (1st order Pattyn)."
        print "(a,es10.3)","Error Tolerance is:", error
        write (*,*) "================================================================"
        write (*,*) "Nonlin.      Lin. Solv.     Error            Max        CPU time"
        write (*,*) " Iter.         Iter.                       Velocity             "
        write (*,*) "================================================================"

        nonlinear_iteration: do while (cont .and. l <= MAXITER)
            call cpu_time(iter_start_time)
            
            lacc=lacc+1
           
            !Every 10 iterations, we relax the error requirement somewhat (1/2
            !an order of magnitude) under the assumption that we won't converge 
            !given the current error tolerance.
            !(TODO: This leads to exponential increase in the tolerance - can we
            !get by without it?  Is this too liberal of a relaxation?)
            if (l == m*10) then
                error = error*toler_adjust_factor
#if DEBUG
                write(*,*) "Error tolerance is now", error
#endif
                m = m+1
            endif
             
            !Compute the basal traction.  This is done by either passing through
            !the beta parameter, or by using Shoof's plastic bed model (Schoof
            !2004).
            if (options%which_ho_bstress == HO_BSTRESS_LINEAR) then
                tau = beta
            else
                call plastic_bed(tau, beta, uvel(nzeta,:,:), vvel(nzeta,:,:))
            end if

            !Compute velocity derivatives
            call velderivs(uvel, vvel, dzdx, dzdy, geometry_mask, &
                           delta_x, delta_y, zeta, .false., &
                           direction_x, direction_y, &
                           dudx, dudy, dudz, dvdx, dvdy, dvdz)
#ifdef DEBUG_FIELDS
            call write_xls_3d("dudx.txt",dudx)
            call write_xls_3d("dudy.txt",dudy)
            call write_xls_3d("dudz.txt",dudz)
            call write_xls_3d("dvdx.txt",dvdx)
            call write_xls_3d("dvdy.txt",dvdy)
            call write_xls_3d("dvdz.txt",dvdz)
#endif


            !Compute viscosity
            if (options%which_ho_efvs == HO_EFVS_FULL) then
                call calcefvs(efvs,flwa,h,ax,ay,FLOWN,ZIP,dudx, dudy, dudz, dvdx, dvdy, dvdz)
            else if (options%which_ho_efvs == HO_EFVS_CONSTANT) then
                efvs = 1d6
            end if


            !Apply periodic boundary conditions to the viscosity
            call periodic_boundaries_3d(efvs,options%periodic_ew,options%periodic_ns)
            
            !Sparse matrix routine for determining velocities.  The new
            !velocities will get spit into ustar and vstar, while uvel and vvel
            !will still hold the old velocities.

            iter=sparuv(efvs,dzdx,dzdy,ax,ay,bx,by,cxy,h,&
                uvel,vvel,dudx,dudy,dudz,dvdx,dvdy,dvdz,&
                ustar,vstar,tau,dhbdx,dhbdy,ijktot,MAXY,&
                MAXX,NZETA,TOLER, options, delta_x, delta_y, zeta, point_mask, &
                geometry_mask,matrix_u, matrix_workspace_u, matrix_options_u, &
                matrix_v, matrix_workspace_v, matrix_options_v, &
                kinematic_bc_u, kinematic_bc_v, &
                marine_bc_normal,direction_x,direction_y, STAGGERED,l)

            !Apply periodic boundary conditions to the computed velocity
            call write_xls("ustar_beforebc.txt", ustar(1,:,:))
            
            call periodic_boundaries_3d(ustar,options%periodic_ew,options%periodic_ns)
            call periodic_boundaries_3d(vstar,options%periodic_ew,options%periodic_ns)

            call write_xls("ustar_afterbc.txt", ustar(1,:,:))

#ifdef OUTPUT_PARTIAL_ITERATIONS
            call iteration_debug_step(ncid_debug, l, efvs, ustar, vstar, geometry_mask, marine_bc_normal)
#endif

    !Apply unstable manifold correction.  This function returns
    !true if we need to keep iterating, false if we reached convergence
    cont = umc_correct_vels(ustar, uvel, &
            vstar, vvel, correction_vec, &
            maxy, maxx, nzeta, error, &
            tot, teta)    

call cpu_time(iter_end_time)

max_vel = maxval(sqrt(ustar**2 + vstar**2))

    print "(i5,10x,i5,7x,es10.3,5x,f10.2,5x,f8.3)", l, iter, tot, max_vel, iter_end_time - iter_start_time

    l = l + 1

    end do nonlinear_iteration

    if ( l >= MAXITER) then
    !call write_log("Maximum iterations exceeded in Pattyn velocity solve", GM_FATAL)
    write(*,*) "Maximum iterations exceeded. Moving on"
    end if

#ifdef OUTPUT_PARTIAL_ITERATIONS
      call end_debug_iteration(ncid_debug)
#endif
      
      call sparse_destroy_workspace(matrix_u, matrix_options_u, matrix_workspace_u)
      call sparse_destroy_workspace(matrix_v, matrix_options_v, matrix_workspace_v)
      call del_sparse_matrix(matrix_u) 
      call del_sparse_matrix(matrix_v) 

      call cpu_time(solve_end_time)

#ifdef DEBUG_FIELDS
      deallocate(pointtype)
#endif

      write(*,*) "Pattyn higher-order solve took",solve_end_time - solve_start_time
      return
      END subroutine

!----------------------------------------------------------------------------------------
! BOP
!
! !IROUTINE: unstable_manifold_correction
!
! !INTERFACE:
    function umc_correct_vels(u_new, u_old, &
                                          v_new, v_old, vec_correction, &
                                          maxx, maxy, nzeta, toler, &
                                          tot_out, theta_out)
        ! !RETURN VALUE:
        logical :: umc_correct_vels !whether another iteration step is needed
 
! !PARAMETERS:
        real(dp), dimension(:,:,:), intent(in) :: u_new  !Computed u component from this iteration
        real(dp), dimension(:,:,:), intent(inout) :: u_old !Computed u component from last iteration
        real(dp), dimension(:,:,:), intent(in) :: v_new !Computed v component from this iteration
        real(dp), dimension(:,:,:), intent(inout) :: v_old !Computed v component from last iteration
        real(dp), dimension(:), intent(inout) :: vec_correction !Old correction vector for v
        integer, intent(in) :: maxy, maxx, nzeta !Grid size
        real(dp), intent(in) :: toler !Error tolerance for the iteration
        real(dp), intent(out):: tot_out !Optional output of error
        real(dp), intent(out):: theta_out !Optional output of angle

        real(dp), dimension(maxy*maxx*nzeta*2) :: vec_new, vec_old
        integer :: linearize_idx
        
        linearize_idx = 1
        call linearize_3d(vec_new, linearize_idx, u_new)
        call linearize_3d(vec_new, linearize_idx, v_new)
        
        linearize_idx = 1
        call linearize_3d(vec_old, linearize_idx, u_old)
        call linearize_3d(vec_old, linearize_idx, v_old)

        umc_correct_vels = unstable_manifold_correction(vec_new, &
                vec_old, vec_correction, maxy*maxx*nzeta*2, toler, tot_out, theta_out)
        
        linearize_idx = 1
        call delinearize_3d(vec_old, linearize_idx, u_old)
        call delinearize_3d(vec_old, linearize_idx, v_old)
        
    end function umc_correct_vels


    !Reduces a 3d higher order velocity estimate to a 2d velocity field.
    subroutine vel_2d_from_3d(uvel, vvel, u, v, zeta)
        real(dp), dimension(:,:,:), intent(in) :: uvel
        real(dp), dimension(:,:,:), intent(in) :: vvel
        real(dp), dimension(:,:),  intent(out) :: u
        real(dp), dimension(:,:),  intent(out) :: v
        real(dp), dimension(:), intent(in) :: zeta
    
        integer :: maxx, maxy, nzeta, i, j, k

        maxx = size(uvel,2)
        maxy = size(uvel,1)
        nzeta = size(uvel,3)

        u = 0
        v = 0
        do i=1,MAXY
            do j=1,MAXX
                do k=2,NZETA
                    !Average velocity for layer (k-1...k) * thickness of layer (k-1..k)
                    !Does this lead to a vertical velocity average???
                    u(i,j) = u(i,j) + (uvel(i,j,k)+uvel(i,j,k-1))*(zeta(k)-zeta(k-1))*0.5
                    v(i,j) = v(i,j) + (vvel(i,j,k)+vvel(i,j,k-1))*(zeta(k)-zeta(k-1))*0.5
                end do
            end do
        end do
    end subroutine vel_2d_from_3d


!-----------------------------------------------------
!   plastic bed computation
!-----------------------------------------------------
    subroutine plastic_bed(tau, tau0, ubas, vbas)
        real(dp), dimension(:,:), intent(out) :: tau
        real(dp), dimension(:,:), intent(in)  :: tau0
        real(dp), dimension(:,:), intent(in)  :: ubas
        real(dp), dimension(:,:), intent(in)  :: vbas
        
        integer :: maxy, maxx, i, j

        maxy = size(ubas,1)
        maxx = size(ubas,2)
        !TODO: Vectorize
        do i = 1,maxx
            do j = 1,maxy
                tau(i,j) = tau0(i,j) / (sqrt(ubas(i,j)**2 + vbas(i,j)**2 + plastic_bed_regularization ** 2))
            end do
        end do
    end subroutine plastic_bed


    !Compute x,y,z derivative fields of u and v, upwinding if necessary at the
    !ice shelf front
    subroutine velderivs(uvel, vvel, dzdx, dzdy, geometry_mask, & 
                         dx, dy, levels, UPSTREAM, &
                         direction_x, direction_y, &
                         dudx, dudy, dudz, dvdx, dvdy, dvdz)
        use glide_mask, only: upwind_from_mask

        real(dp), dimension(:,:,:), intent(in) :: uvel
        real(dp), dimension(:,:,:), intent(in) :: vvel
        real(dp), dimension(:,:), intent(in) :: dzdx
        real(dp), dimension(:,:), intent(in) :: dzdy
        integer, dimension(:,:), intent(in) :: geometry_mask
        real(dp), intent(in) :: dx
        real(dp), intent(in) :: dy
        real(dp), dimension(:), intent(in) :: levels
        logical, intent(in) :: UPSTREAM
        real(dp), dimension(:,:), intent(out) :: direction_x
        real(dp), dimension(:,:), intent(out) :: direction_y

        real(dp), dimension(:,:,:), intent(out) :: dudx
        real(dp), dimension(:,:,:), intent(out) :: dudy
        real(dp), dimension(:,:,:), intent(out) :: dudz
        real(dp), dimension(:,:,:), intent(out) :: dvdx        
        real(dp), dimension(:,:,:), intent(out) :: dvdy        
        real(dp), dimension(:,:,:), intent(out) :: dvdz
 
        real(dp), dimension(size(uvel,1),size(uvel,2),size(uvel,3)) :: uvel_test, vvel_test

        !TEST of upwinding shelf boundary derivatives:
        !Everywhere where there's no ice, we set the velocity
        !to NaN.  If a location with no ice is ever used, this should
        !propegate.
        uvel_test = uvel
        vvel_test = vvel

#if 0
        do k=1,size(uvel,3)
            where(.not. GLIDE_HAS_ICE(geometry_mask))
                uvel_test(:,:,k) = NaN
                vvel_test(:,:,k) = NaN
            endwhere
        end do
#endif

        !TODO: Implement option for upstream differencing
        if (UPSTREAM) then
            !Upstream the number 
            call df_field_3d(uvel, dx, dy, levels, dudx, dudy, dudz, dzdy, dzdx)
            call df_field_3d(vvel, dx, dy, levels, dvdx, dvdy, dvdz, dzdy, dzdx)
        else
            direction_x = 0
            direction_y = 0
            call upwind_from_mask(geometry_mask, direction_x, direction_y)
            !call write_xls("direction_x.txt",direction_x)
            !call write_xls("direction_y.txt",direction_y)
            
            call df_field_3d(uvel_test, dx, dy, levels, dudx, dudy, dudz, direction_x, direction_y)
            call df_field_3d(vvel_test, dx, dy, levels, dvdx, dvdy, dvdz, direction_x, direction_y)
         end if

    end subroutine

!
!
!------------------------------------------------------
!   nonlinear viscosity term
!------------------------------------------------------
!
      subroutine calcefvs(efvs, flwa, h, ax, ay, FLOWN, ZIP, &
                        dudx, dudy, dudz, dvdx, dvdy, dvdz)
!
        real(dp),                   intent(in)  :: FLOWN,ZIP
        real(dp), dimension(:,:,:), intent(out) :: efvs
        real(dp), dimension(:,:,:), intent(in)  :: flwa
        real(dp), dimension(:,:),   intent(in)  :: h
        real(dp), dimension(:,:,:), intent(in)  :: ax
        real(dp), dimension(:,:,:), intent(in)  :: ay
        
        real(dp), dimension(:,:,:), intent(in) :: dudx, dudy, dudz
        real(dp), dimension(:,:,:), intent(in) :: dvdx, dvdy, dvdz
       !
        INTEGER :: i,j,k, MAXX, MAXY, NZETA
        real(dp) :: macht
        real(dp) :: exx,eyy,exy,exz,eyz,eeff
!
        MAXY = size(efvs, 3)
        MAXX = size(efvs, 2)
        NZETA = size(efvs, 1)

        macht=(1.-FLOWN)/(2.*FLOWN)
     
        !Compute efvs term from the acceleration fields
        do i=1,MAXX
            do j=1,MAXY
                do k=1,NZETA
          
                    if (h(i,j) > 0.) then
                        exx=dudx(k,i,j)+ax(k,i,j)*dudz(k,i,j)
                        eyy=dvdy(k,i,j)+ay(k,i,j)*dvdz(k,i,j)
                        exy=0.5*((dudy(k,i,j)+ay(k,i,j)*dudz(k,i,j))+(dvdx(k,i,j)+ax(k,i,j)*dvdz(k,i,j)))
                        exz=-0.5*dudz(k,i,j)/h(i,j)
                        eyz=-0.5*dvdz(k,i,j)/h(i,j)
                    else
                        exx=dudx(k,i,j)
                        eyy=dvdy(k,i,j)
                        exy=0.5*(dudy(k,i,j)+dvdx(k,i,j))
                        exz=0.
                        eyz=0.
                    endif
            
                    eeff=(exx*exx)+(eyy*eyy)+(exx*eyy)+(exy*exy)+(exz*exz)+(eyz*eyz)
                    efvs(k,i,j)=0.5 * (flwa(k,i,j)**(-1./FLOWN)) * ((eeff+ZIP)**macht)
                end do
            end do
        end do
        return
    end subroutine

    function sparuv(efvs,dzdx,dzdy,ax,ay,bx,by,cxy,h,uvel,vvel,dudx,dudy,dudz,&
                    dvdx,dvdy,dvdz,ustar,vstar,beta,dhbdx,dhbdy,IJKTOT,MAXY,  &
                    MAXX,NZETA,TOLER,options,GRIDX,GRIDY,zeta, point_mask,    & 
                    geometry_mask, matrix_u, matrix_workspace_u,              &
                    matrix_options_u, matrix_v, matrix_workspace_v,           &
                    matrix_options_v,kinematic_bc_u, kinematic_bc_v,          &
                    latbc_normal, direction_x, direction_y, STAGGERED,l)

        INTEGER IJKTOT,MAXY,MAXX,NZETA
        real(dp), dimension(:,:,:), intent(in) :: efvs
        real(dp), dimension(:,:), intent(in) :: dzdx
        real(dp), dimension(:,:), intent(in) :: dzdy
        real(dp), dimension(:,:,:), intent(in) :: ax
        real(dp), dimension(:,:,:), intent(in) :: ay
        real(dp), dimension(:,:,:), intent(in) :: bx
        real(dp), dimension(:,:,:), intent(in) :: by
        real(dp), dimension(:,:,:), intent(in) :: cxy
        real(dp), dimension(:,:), intent(in) :: h
        real(dp), dimension(:,:,:), intent(in), target :: uvel
        real(dp), dimension(:,:,:), intent(in), target :: vvel
        
        real(dp), dimension(:,:,:), intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz
        
        real(dp), dimension(:,:,:), intent(out), target :: ustar
        real(dp), dimension(:,:,:), intent(out), target :: vstar
        real(dp), dimension(:,:), intent(in) :: dhbdx
        real(dp), dimension(:,:), intent(in) :: dhbdy
        real(dp), dimension(:,:), intent(in) :: beta
        real(dp), dimension(:), intent(in) :: zeta
        integer, dimension(:,:), intent(in) :: geometry_mask
        real(dp), dimension(:,:,:), intent(in), target :: kinematic_bc_u
        real(dp), dimension(:,:,:), intent(in), target :: kinematic_bc_v
        real(dp), dimension(:,:), intent(in)   :: latbc_normal !On the marine ice front, this is the angle of the normal to the ice front
        real(dp), intent(in) :: toler
        real(dp), intent(in) :: gridx
        real(dp), intent(in) :: gridy
        real(dp), dimension(:,:), intent(in) :: direction_x,direction_y
        integer, intent(in) :: l

        !Sparse matrix variables.  These are passed in so that allocation can be
        !done once per velocity solve instead of once per iteration
        type(sparse_matrix_type), target, intent(inout) :: matrix_u,matrix_v

        type(sparse_solver_workspace), target,intent(in):: matrix_workspace_u,&
                                                           matrix_workspace_v

        type(sparse_solver_options), target,intent(in)  :: matrix_options_u,  & 
                                                           matrix_options_v

        type(glide_options), intent(in) :: options
        logical, intent(in)::STAGGERED
        real(dp) :: rhs
        
        INTEGER i,j,k,m,sparuv,iter, ierr
        real(dp) :: d(IJKTOT),x(IJKTOT),coef(STENCIL_SIZE),err
      
        integer, dimension(:,:) :: point_mask
        
        integer :: stencil_center_idx
        integer :: si, sj, sk
        real(dp), dimension(:,:,:), pointer :: velpara, velperp, velpara_star, kinematic_bc_para
        type(sparse_matrix_type), pointer :: matrix
        type(sparse_solver_workspace), pointer :: matrix_workspace
        type(sparse_solver_options), pointer :: matrix_options

        integer :: whichcomponent
        character(1) :: componentstr
#ifdef DEBUG_FIELDS
        pointtype = 0
#endif
        sparuv = 0
        do whichcomponent = 1,2
            if (whichcomponent == 1) then !Set up to compute u component
                velpara => uvel
                velperp => vvel
                velpara_star => ustar
                componentstr = "u"
                kinematic_bc_para => kinematic_bc_u
                matrix => matrix_u
                matrix_workspace => matrix_workspace_u
                matrix_options => matrix_options_u
            else !Compute v component
                velpara => vvel
                velperp => uvel
                velpara_star => vstar
                componentstr = "v"
                kinematic_bc_para => kinematic_bc_v
                matrix => matrix_v
                matrix_workspace => matrix_workspace_v
                matrix_options => matrix_options_v
            end if
            !Initialize sparse matrix & vectors
            d=0
            x=0
            call sparse_clear(matrix_u)
            call sparse_clear(matrix_v)
#ifdef VERY_VERBOSE
            write(*,*)"Begin Matrix Assembly"
#endif

            do i=1,MAXX
                do j=1,MAXY
                    if (point_mask(i,j) /= 0) then
                        do k=1,NZETA
                            coef = 0
                            stencil_center_idx = csp_masked(I_J_K,i,j,k,point_mask,NZETA) 
                            if (.not. GLIDE_HAS_ICE( geometry_mask(i,j) ) .or. &
                                      GLIDE_IS_LAND_MARGIN( geometry_mask(i,j) ) .or. &
                                      GLIDE_IS_THIN( geometry_mask(i,j) ) &
                                         .and. .not. options%ho_include_thinice) then
                                !No ice - "pass through"
                                coef(I_J_K)=1.
                                !Normally, we use uvel(i,j,k) as our initial guess.
                                !However, in this case, we know the velocity
                                !should be 0, so we'll replace our uvel with that
                                velpara(k,i,j) = 0
                                rhs = 0
#ifdef DEBUG_FIELDS
                                pointtype(i,j) = 1
#endif
                            else if (.not. IS_NAN(kinematic_bc_para(k,i,j))) then
                                !If a kinematic boundary condition was specified
                                !for this location, hold the location at the
                                !specified value.
                                coef(I_J_K)=1.
                                rhs = kinematic_bc_para(k,i,j)
#ifdef DEBUG_FIELDS
                                pointtype(i,j) = 2
#endif
                            else if ((i.eq.1).or.(i.eq.MAXX).or.(j.eq.1).or.(j.eq.MAXY)) then
                                ! Boundary condition at the edges of the domain.
                                ! If we don't have a kinematic boundary already
                                ! specified, we just "pass through" the initial
                                ! guess.  If this lies on a periodic boundary, 
                                ! enforce it in the sparse matrix 
                                ! (which means we need to do some insertions 
                                ! ourselves...)
                                rhs=velpara(k,i,j)
 
#ifdef ENFORCE_PBC
                                if (i .eq. 1 .and. options%periodic_ew) then
                                    rhs = 0
                                    call sparse_insert_val(matrix, &
                                               stencil_center_idx, &
                                               csp_stenciled(j, MAXX - 1, k,point_mask,NZETA), &
                                               -1d0)
                                else if (i .eq. MAXX .and. options%periodic_ew) then
                                    rhs = 0
                                    call sparse_insert_val(matrix, &
                                               stencil_center_idx, &
                                               csp_stenciled(j, 2, k,point_mask,NZETA), &
                                               -1d0)                                
                                else if (j .eq. 1 .and. options%periodic_ns) then
                                    rhs = 0
                                    call sparse_insert_val(matrix, &
                                               stencil_center_idx, &
                                               csp_stenciled(MAXY - 1, i, k,point_mask,NZETA), &
                                               -1d0) 
                                else if (j .eq. MAXY .and. options%periodic_ns) then
                                    rhs = 0
                                    call sparse_insert_val(matrix, &
                                               stencil_center_idx, &
                                               csp_stenciled(2, i, k,point_mask,NZETA), &
                                               -1d0)                                                            
                                end if
#endif
                                coef(I_J_K)=1.

#ifdef DEBUG_FIELDS
                                pointtype(i,j) = 3
#endif
                            else
#ifdef DEBUG_FIELDS
                                pointtype(i,j) = 4
#endif
                                call sparse_setup(componentstr, i, j, k, efvs, dzdx, dzdy, ax, ay, bx, by, cxy, &
                                     h, gridx, gridy, zeta, uvel, vvel, dudx, dudy, dudx, dvdx, dvdy, dvdz, &
                                     dhbdx, dhbdy, beta, geometry_mask, &
                                     latbc_normal, maxx, maxy, Nzeta, coef, rhs, direction_x, direction_y, STAGGERED, &
                                     options%which_ho_source)
                            endif
                            d(stencil_center_idx)=rhs
                            !Preliminary benchmarks indicate the we actually reach
                            !convergance faster if we use a 0 initial guess rather than
                            !the previous velocity.  More testing is needed.
                            !To use the current velocity as the guess, uncomment this
                            !line here and and in the Y velocity section of this
                            !function.
                            !x=csp(I_J_K,i,j,k,MAXX,NZETA))=uvel(i,j,k) !Use current velocity as guess of solution
                            if (abs(coef(I_J_K)) < SMALL) then
                                write(*,*) "WARNING: 0 on diagonal at position",i,j,k
                                write(*,*) "component:",componentstr
                                write(*,*) "Thickness at this point:", h(i,j)
                                write(*,*) "Mask at this point:", geometry_mask(i,j)
                                stop
                            end if

                            do m=1,STENCIL_SIZE
                                if (abs(coef(m)) > SMALL) then
                                    !Because of the transposition of Pattyn's original code compared to
                                    !Glimmer/CISM, we need to interpret indices slightly differently:
                                    !i = y index (would normally be z index)
                                    !k = z index (would normally be y index)
                                    !j = x as usual
                                    si = stencil_y(m,j)
                                    sj = stencil_x(m,i)
                                    sk = stencil_z(m,k)

                                    if (si > 0 .and. si <= maxy .and. &
                                        sj > 0 .and. sj <= maxx .and. &
                                        sk > 0 .and. sk <= nzeta) then
                                            if (point_mask(sj,si) == 0) then
                                                write(*,*) "ERROR: point is off mask."
                                                write(*,*) "component:",componentstr
                                                write(*,*) "location:",i,j,k
                                                write(*,*) "stencil:",si,sj,sk
                                                write(*,*) "position and coefficient in stencil:", m, coef(m)
                                                write(*,*) "h for stencil center and bad point:",h(i,j),h(si,sj)
                                                write(*,*) "Mask for stencil center and bad point:",&
                                                           geometry_mask(i,j), geometry_mask(si,sj)
                                                write(*,*) "Coefficient:", coef(m)
                                            end if
                                            call sparse_insert_val(matrix, &
                                                    stencil_center_idx, &
                                                    csp_stenciled(si,sj,sk,point_mask,NZETA), &
                                                    coef(m)) 
                                    end if
                                end if
                            end do !End loop through stencil
                        end do !END k loop
                    end if !End mask check
                end do !End j loop
            end do !End i loop
#ifdef VERY_VERBOSE
            write(*,*)"End Matrix Assembly"
            write(*,*)"Begin Matrix Solve"
#endif

#ifdef DEBUG_FIELDS
            call write_xls_int("point_type.txt", pointtype)

#endif
            call sparse_solver_preprocess(matrix, matrix_options, matrix_workspace)             
            ierr = sparse_solve(matrix, d, x, matrix_options, matrix_workspace,  err, iter, verbose=sparverbose)

            sparuv = sparuv + iter
       
            if (ierr /= 0) then
                !The sparse solve failed.  If a fallback was specified, use it
                if (options%which_ho_sparse_fallback >= 0 .and. &
                    options%which_ho_sparse_fallback /= options%which_ho_sparse) then
                    write(*,*)"Sparse solve failed, falling back on alternate method"
                    call sparse_easy_solve(matrix, d, x, err, iter, options%which_ho_sparse_fallback, &
                                           __FILE__, __LINE__)
                else
                    call handle_sparse_error(matrix, matrix_options, ierr, __FILE__, __LINE__)      
                end if
            end if
            call sparse_solver_postprocess(matrix, matrix_options, matrix_workspace)
#ifdef VERY_VERBOSE
            write(*,*)"End Matrix Solve"
#endif
            !Delinearize the solution
            do i=1,MAXX
                do j=1,MAXY
                    do k=1,NZETA
                        if (point_mask(i,j) /= 0) then
                            !write(*,*)csp_masked(11,i,j,k,point_mask,nzeta)
                            velpara_star(k,i,j)=x(csp_masked(11,i,j,k,point_mask,nzeta))
                        else
                            velpara_star(k,i,j)=0
                        end if
                    end do
                end do
            end do
        end do !END whichcomponent loop
    end function sparuv
!
!
!
!------------------------------------------------------
!   lookup coordinates in sparse matrix
!
!   Reference:
!     1   ->   i - 1, j - 1, k
!     2   ->   i - 1, j, k - 1
!     3   ->   i - 1, j, k
!     4   ->   i - 1, j, k + 1
!     5   ->   i - 1, j + 1, k
!     6   ->   i, j - 1, k - 1
!     7   ->   i, j - 1, k
!     8   ->   i, j - 1, k + 1
!     9   ->   i, j, k - 2
!    10   ->   i, j, k - 1
!    11   ->   i, j, k
!    12   ->   i, j, k + 1
!    13   ->   i, j, k + 2
!    14   ->   i, j + 1, k - 1
!    15   ->   i, j + 1, k
!    16   ->   i, j + 1, k + 1
!    17   ->   i + 1, j - 1, k
!    18   ->   i + 1, j, k - 1
!    19   ->   i + 1, j, k
!    20   ->   i + 1, j, k + 1
!    21   ->   i + 1, j + 1, k
!------------------------------------------------------
!
      !Given a point and a stencil location,
      !the following three functions return the equivalent location
      
      !NOTE: When referring to stencils, i refers to the y coordinate
      function stencil_y(pos, i)
        integer, intent(in) :: pos, i
        integer :: stencil_y
        if (pos <= 5) then
            stencil_y = i - 1
        else if (pos <= 16 .or. pos == 24 .or. pos == 25) then
            stencil_y = i
        else if (pos == 22) then
            stencil_y = i - 2
        else if (pos == 23) then
            stencil_y = i + 2
        else
            stencil_y = i + 1
        end if
      end function

      !NOTE: When referring to stencils, j refers to the x coordinate
      function stencil_x(pos, j)
        integer, intent(in) :: pos, j
        integer :: stencil_x
        select case(pos)
            case(1, 6, 7, 8, 17)
                stencil_x = j - 1
            case(5, 14, 15, 16, 21)
                stencil_x = j + 1
            case(24)
                stencil_x = j - 2
            case(25)
                stencil_x = j + 2
            case default
                stencil_x = j
        end select
      end function

      function stencil_z(pos, k)
        integer, intent(in) :: pos, k
        integer :: stencil_z
        select case(pos)
            case(2, 6, 10, 14, 18)
                stencil_z = k - 1
            case(4, 8, 12, 16, 20)
                stencil_z = k + 1
            case(9)
                stencil_z = k - 2
            case(13)
                stencil_z = k + 2
            case default
                stencil_z = k
        end select
      end function
 !
      function csp_masked(pos, i, j, k, point_mask, nzeta)
         integer, intent(in) :: pos
         integer, intent(in) :: i !X coordinate
         integer, intent(in) :: j !Y coordinate
         integer, intent(in) :: k !Z coordinate
         integer, dimension(:,:), intent(in) :: point_mask
         integer, intent(in) :: nzeta

         integer :: csp_masked
         !TODO: Should i and j be transposed here?
         csp_masked = (point_mask(stencil_x(pos,i), stencil_y(pos,j)) - 1) * nzeta &
                      + stencil_z(pos, k)
      end function

      !Returns the coordinate without performing stencil lookups (this assumes
      !those have already been done)
      !NOTE: si should refer to the y coordinate, and sj should refer to the x coordinate!!
      function csp_stenciled(si, sj, sk, point_mask, nzeta)
            integer, intent(in) :: si, sj, sk
            integer, dimension(:,:), intent(in) :: point_mask
            integer, intent(in) :: nzeta

            integer :: csp_stenciled

            csp_stenciled = (point_mask(sj, si) - 1)*nzeta + sk
      end function

!
!------------------------------------------------------
!   central difference calculation for uvel
!------------------------------------------------------
!

    !This differencing scheme function uses the following nomenclature.
    !This reduces the code duplication that would have otherwise happened!!!
    !u - Refers to velocity in the x direction
    !v - Refers to velocity in the y direction
    !para - Refers to a vector that is parallel to the component of the
    !       velocity being computed
    !perp - Refers to a vector that is perpendicular to the component of
    !       the velocity being computed
    !Under this nomenclature, d_para_d_para means either "dudx" or "dvdy".
    !d_perp_d_perp means "dvdy" or "dudx" (in each case, for computing u and v)
    !Note: x and y might be referred to directly when they don't change
    !      depending on the differencing scheme being used
    !
    !The "component" argument should be a single character indicating
    !whether the "u" or "v" component is being requested.
    !DON'T PANIC, this works!  It's been tested!  Really!  The bug is somewhere
    !else!
    subroutine sparse_setup(component, i,j,k,efvs,dzdx,dzdy,ax,ay,bx,by,cxy,&
        h, dx, dy, dz, uvel, vvel, dudx_field, dudy_field, dudz_field, &
        dvdx_field, dvdy_field, dvdz_field, &
        dhbdx,dhbdy,beta,mask,latbc_normal,MAXX,MAXY,Ndz,&
        coef, rhs, direction_x, direction_y, STAGGERED, WHICH_SOURCE)
!
        integer, intent(in) :: i,j,k, MAXY, MAXX, Ndz

        !Output array
        real(dp), dimension(STENCIL_SIZE), intent(out) :: coef

        !Output RHS value
        real(dp), intent(out) :: rhs

        !Viscosity
        real(dp), dimension(:,:,:), intent(in) :: efvs

        real(dp), dimension(:,:,:), intent(in) :: dudx_field, dudy_field, dudz_field 
        real(dp), dimension(:,:,:), intent(in) :: dvdx_field, dvdy_field, dvdz_field

        !Surface Gradients
        real(dp), dimension(:,:), target, intent(in) :: dzdx, dzdy

        !Rescaling Factors
        real(dp), dimension(:,:,:), target, intent(in) :: ax, ay, bx, by, cxy
        
        !Ice Thickness
        real(dp), dimension(:,:), intent(in) :: h

        !Grid Spacing (Z is an irregular grid)
        real(dp), intent(in) :: dx, dy
        real(dp), dimension(:), intent(in) :: dz

        !Current velocities
        real(dp), dimension(:,:,:), target, intent(in) :: uvel, vvel

        !Bedrock Gradients
        real(dp), dimension(:,:), target, intent(in) :: dhbdx, dhbdy

        !Basal Traction
        real(dp), dimension(:,:), intent(in) :: beta

        integer, dimension(:,:), intent(in) :: mask

        !Contains the angle of the normal to the marine margine, NaN
        !everywhere not on the margin
        real(dp), dimension(:,:), intent(in) :: latbc_normal

        logical, intent(in) :: STAGGERED

        integer, intent(in) :: WHICH_SOURCE

        real(dp), dimension(:,:),intent(in) :: direction_x, direction_y

        character(1), intent(in) :: component
!
        !Temporary calculations done in the ABSOLUTE (using u/x, v/y) coordinate
        !system
        real(dp), target :: defvsdx, defvsdy, defvsdz, defvsdx2, defvsdy2
        real(dp), target :: dudx, dudy, dudz, dudx2, dudy2, dudz2, dudxz, dudyz, dudxy
        real(dp), target :: dvdx, dvdy, dvdz, dvdx2, dvdy2, dvdz2, dvdxz, dvdyz, dvdxy

        !Temporary calculations done in the RELATIVE (using para, perp)
        !coordinate system
        real(dp), pointer :: dpara_dpara, dpara_dperp, dpara_dx, dpara_dy, dpara_dz
        real(dp), pointer :: dpara_dpara2, dpara_dparaz, dpara_dz2, dpara_dperp2, dpara_dperpz
        real(dp), pointer :: dpara_dyz, dpara_dxz, dpara_dx2, dpara_dy2
        real(dp), pointer :: dperp_dpara, dperp_dperp, dperp_dx, dperp_dy, dperp_dz
        real(dp), pointer :: dperp_dxy, dperp_dxz, dperp_dyz, dperp_dz2
        real(dp), pointer :: defvs_dpara2, defvs_dperp2

        !Fields done in the RELATIVE coordinate system
        real(dp), dimension(:,:,:), pointer :: a_para, a_perp, b_para, b_perp
        real(dp), dimension(:,:,:), pointer :: vel_perp, vel_para
        real(dp), dimension(:,:), pointer :: dz_dpara, dz_dperp, dhb_dpara, dhb_dperp

        !Cached difference calculations for the nonstaggered Z grid
        real(dp) dz_down1, dz_down2, dz_down3, dz_up1, dz_up2, dz_up3
        real(dp) dz_cen1, dz_cen2, dz_cen3, dz_sec1, dz_sec2, dz_sec3
        
        !Set up a system of pointers to encapsulate the nomenclature described
        !above
        !TODO: Go through the code and figure out which of these bindings
        !      can just become the variable that gets assigned to
        if (component == "u") then
            dpara_dpara => dudx
            dpara_dperp => dudy
            dpara_dx    => dudx
            dpara_dy    => dudy
            dpara_dz    => dudz
            
            dpara_dpara2 => dudx2
            dpara_dparaz => dudxz
            dpara_dz2    => dudz2
            dpara_dx2    => dudx2
            dpara_dy2    => dudy2
            dpara_dperp2 => dudy2
            dpara_dperpz => dudyz
            dpara_dyz    => dudyz
            dpara_dxz    => dudxz

            dperp_dpara => dvdx
            dperp_dperp => dvdy
            dperp_dx    => dvdx
            dperp_dy    => dvdy
            dperp_dz    => dvdz

            dperp_dxy => dvdxy
            dperp_dxz => dvdxz
            dperp_dyz => dvdyz
            dperp_dz2 => dvdz2

            defvs_dpara2 => defvsdx2
            defvs_dperp2 => defvsdy2

            a_para => ax
            a_perp => ay
            b_para => bx
            b_perp => by

            dz_dpara  => dzdx
            dz_dperp  => dzdy
            dhb_dpara => dhbdx
            dhb_dperp => dhbdy

            vel_perp => vvel
            vel_para => uvel
        else
            dpara_dpara => dvdy
            dpara_dperp => dvdx
            dpara_dx    => dvdx
            dpara_dy    => dvdy
            dpara_dz    => dvdz
            
            dpara_dpara2 => dvdy2
            dpara_dparaz => dvdyz
            dpara_dx2    => dvdx2
            dpara_dy2    => dvdy2
            dpara_dz2    => dvdz2
            dpara_dperp2 => dvdx2
            dpara_dperpz => dvdxz
            dpara_dyz    => dvdyz
            dpara_dxz    => dvdxz

            dperp_dpara => dudy
            dperp_dperp => dudx
            dperp_dx    => dudx
            dperp_dy    => dudy
            dperp_dz    => dudz

            dperp_dxy => dudxy
            dperp_dxz => dudxz
            dperp_dyz => dudyz
            dperp_dz2 => dudz2

            defvs_dpara2 => defvsdy2
            defvs_dperp2 => defvsdx2

            a_para => ay
            a_perp => ax
            b_para => by
            b_perp => bx

            dz_dpara  => dzdy
            dz_dperp  => dzdx
            dhb_dpara => dhbdy
            dhb_dperp => dhbdx

            vel_perp => uvel
            vel_para => vvel
        end if
#ifndef NOSHELF        
        if (GLIDE_IS_CALVING(mask(i,j)) .and. &
            WHICH_SOURCE /= HO_SOURCE_DISABLED) then !Marine margin dynamic (Neumann) boundary condition
            call sparse_marine_margin(component,i,j,k,h,latbc_normal, dzdx, dzdy, dhbdx, dhbdy, &
                                      vel_perp, efvs, dx, dy, ax, ay, dz,coef, rhs, &
                                      dudx_field,dudy_field,dudz_field,dvdx_field,dvdy_field,dvdz_field,&
                                      direction_x, direction_y, mask, STAGGERED, WHICH_SOURCE)
        else  &
#endif
        if (k.eq.1) then !Upper boundary condition (stress-free surface)
            !Finite difference coefficients for an irregular Z grid, downwinded
            dz_down1=(2.*dz(k)-dz(k+1)-dz(k+2))/(dz(k+1)-dz(k))/(dz(k+2)-dz(k))
            dz_down2=(dz(k+2)-dz(k))/(dz(k+2)-dz(k+1))/(dz(k+1)-dz(k))
            dz_down3=(dz(k)-dz(k+1))/(dz(k+2)-dz(k+1))/(dz(k+2)-dz(k))
 
            dpara_dpara = 4*dz_dpara(i,j)
            dpara_dz = 4*a_para(k,i,j)*dz_dpara(i,j) + a_perp(k,i,j)*dz_dperp(i,j) + 1./h(i,j)
            dpara_dperp = dz_dperp(i,j)

            dperp_dpara = dz_dperp(i,j)
            dperp_dperp = 2*dz_dpara(i,j)
            dperp_dz = 2.*a_perp(k,i,j)*dz_dpara(i,j)+a_para(k,i,j)*dz_dperp(i,j)

            coef(IM1_J_K) = -.5*dpara_dy/dy
            coef(I_JM1_K) = -.5*dpara_dx/dx
            coef(I_J_K) = dpara_dz*dz_down1
            coef(I_J_KP1) = dpara_dz*dz_down2
            coef(I_J_KP2) = dpara_dz*dz_down3
            coef(I_JP1_K) = dpara_dx*.5/dx
            coef(IP1_J_K) = dpara_dy*.5/dy
            rhs = -dperp_dx * dfdx_3d(vel_perp,i,j,k,dx) &
                -  dperp_dy * dfdy_3d(vel_perp,i,j,k,dy) &
                -  dperp_dz * dfdz_3d_downwind_irregular(vel_perp,i,j,k,dz)
                
        else if (k.eq.Ndz) then !Lower boundary condition (Basal sliding)
            !Finite difference coefficients for an irregular Z grid, upwinded
            dz_up1=(dz(k)-dz(k-1))/(dz(k-1)-dz(k-2))/(dz(k)-dz(k-2))
            dz_up2=(dz(k-2)-dz(k))/(dz(k)-dz(k-1))/(dz(k-1)-dz(k-2))
            dz_up3=(2.*dz(k)-dz(k-1)-dz(k-2))/(dz(k)-dz(k-1))/(dz(k)-dz(k-2))

            if (beta(i,j) >= BIG) then !No sliding at the base
                coef(I_J_K)=1.
                rhs=0.
            else !Compute sliding 
                dpara_dpara = 4.*dhb_dpara(i,j)
                dpara_dz = 4.*a_para(k,i,j)*dhb_dpara(i,j) + a_perp(k,i,j)*dhb_dperp(i,j) + 1./h(i,j)
                dpara_dperp = dhb_dperp(i,j)

                dperp_dpara = dhb_dperp(i,j)
                dperp_dperp = 2.*dhb_dpara(i,j)
                dperp_dz = 2.*a_perp(k,i,j)*dhb_dpara(i,j) + a_para(k,i,j)*dhb_dperp(i,j)

                !If we have a non-zero friction, then we need to scale the
                !strain rates by the viscosity.  If there is no friction, then
                !we have a stress-free base and these coefficients disappear.
                if (beta(i,j) >= SMALL) then
                   dpara_dpara = dpara_dpara*efvs(k,i,j)
                   dpara_dperp = dpara_dperp*efvs(k,i,j)
                   dpara_dz = dpara_dz*efvs(k,i,j)
                   dperp_dpara = dperp_dpara*efvs(k,i,j)
                   dperp_dperp = dperp_dperp*efvs(k,i,j)
                   dperp_dz = dperp_dz*efvs(k,i,j)
                end if

                coef(IM1_J_K) = -.5*dpara_dy/dy
                coef(I_JM1_K) = -.5*dpara_dx/dx
                coef(I_J_KM2) = dpara_dz*dz_up1
                coef(I_J_KM1) = dpara_dz*dz_up2
                coef(I_J_K)   = dpara_dz*dz_up3
                coef(I_JP1_K) = .5*dpara_dx/dx
                coef(IP1_J_K) = .5*dpara_dy/dy
                rhs = -dperp_dx * dfdx_3d(vel_perp, i, j, k, dx) &
                    -  dperp_dy * dfdy_3d(vel_perp, i, j, k, dy) &
                    -  dperp_dz * dfdz_3d_upwind_irregular(vel_perp, i, j, k, dz)

                !Adjust center of stencil for the basal friction parameter (this
                !occurs in all cases b/c the ice shelf case necessarily has beta = 0
                !CORRECTION: Adding in normalization term that Pattyn has
                !neglected.
                coef(11) = coef(11) + beta(i,j)! * sqrt(1 + dhbdx(i,j)**2 + dhbdy(i,j)**2)
                
            endif
        else !Interior of the ice (e.g. not a boundary condition)
            !Finite difference coefficients for an irregular Z grid
            dz_cen1 = (dz(k)-dz(k+1))/(dz(k)-dz(k-1))/(dz(k+1)-dz(k-1))
            dz_cen2 = (dz(k+1)-2.*dz(k)+dz(k-1))/(dz(k)-dz(k-1))/(dz(k+1)-dz(k))
            dz_cen3 = (dz(k)-dz(k-1))/(dz(k+1)-dz(k))/(dz(k+1)-dz(k-1))
            dz_sec1 = 2./(dz(k)-dz(k-1))/(dz(k+1)-dz(k-1))
            dz_sec2 = 2./(dz(k)-dz(k+1))/(dz(k)-dz(k-1))
            dz_sec3 = 2./(dz(k+1)-dz(k))/(dz(k+1)-dz(k-1))

            defvsdx=dfdx_3d(efvs,i,j,k,dx)
            defvsdy=dfdy_3d(efvs,i,j,k,dy)
            defvsdz=dfdz_3d_irregular(efvs,i,j,k,dz)
            defvsdx2=defvsdx+ax(k,i,j)*defvsdz
            defvsdy2=defvsdy+ay(k,i,j)*defvsdz

            dpara_dpara = 4.*defvs_dpara2
            dpara_dz = 4.*a_para(k,i,j)*defvs_dpara2 + 4.*b_para(k,i,j)*efvs(k,i,j) &
                     +    a_perp(k,i,j)*defvs_dperp2 +    b_perp(k,i,j)*efvs(k,i,j) &
                     +    defvsdz/(h(i,j)*h(i,j))
            dpara_dpara2 = 4.*efvs(k,i,j)
            dpara_dz2 = efvs(k,i,j)*(4.*a_para(k,i,j)**2 + a_perp(k,i,j)**2 + 1/h(i,j)**2)
            dpara_dparaz = 8.*a_para(k,i,j)*efvs(k,i,j)
            dpara_dperp=defvs_dperp2
            dpara_dperp2=efvs(k,i,j)
            dpara_dperpz=2.*efvs(k,i,j)*a_perp(k,i,j)

            dperp_dpara = defvs_dperp2
            dperp_dperp = 2 * defvs_dpara2
            dperp_dz  = 2 * a_perp(k,i,j) * defvs_dpara2  + a_para(k,i,j) * defvs_dperp2 + 3 * efvs(k,i,j) * cxy(k,i,j) 
            dperp_dxy = 3 * efvs(k,i,j)
            dperp_dxz = 3 * efvs(k,i,j) * ay(k,i,j)
            dperp_dyz = 3 * efvs(k,i,j) * ax(k,i,j)
            dperp_dz2 = 3 * efvs(k,i,j) * ax(k,i,j)*ay(k,i,j)

            coef(IM1_J_KM1)  = dpara_dyz*dz_cen1*(-.5/dy)
            coef(IM1_J_K)  = (dpara_dyz*dz_cen2 + dpara_dy)*(-.5/dy) + dpara_dy2/(dy**2)
            coef(IM1_J_KP1)  = dpara_dyz*dz_cen3*(-.5/dy)
          
            coef(I_JM1_KM1)  = dpara_dxz*dz_cen1*(-.5/dx)
            coef(I_JM1_K)  = (dpara_dx + dpara_dxz*dz_cen2)*(-.5/dx) + dpara_dx2 / dx**2
            coef(I_JM1_KP1)  = dpara_dxz*dz_cen3*(-.5/dx)
          
            coef(I_J_KM1) = dpara_dz*dz_cen1 + dpara_dz2*dz_sec1
            coef(I_J_K) = dpara_dz*dz_cen2 + dpara_dz2*dz_sec2 + dpara_dx2*(-2/dx**2) + dpara_dy2*(-2/dy**2)
            coef(I_J_KP1) = dpara_dz*dz_cen3 + dpara_dz2*dz_sec3
          
            coef(I_JP1_KM1) = dpara_dxz * dz_cen1 * (.5/dx)
            coef(I_JP1_K) = (dpara_dx + dpara_dxz*dz_cen2) * (.5/dx) + dpara_dx2 / dx**2
            coef(I_JP1_KP1) = dpara_dxz * dz_cen3 * (.5/dx)
          
            coef(IP1_J_KM1) = dpara_dyz * dz_cen1 * (.5/dy)
            coef(IP1_J_K) = (dpara_dyz*dz_cen2 + dpara_dy) * (.5/dy) + dpara_dy2 / dy**2
            coef(IP1_J_KP1) = dpara_dyz * dz_cen3 * (.5/dy)
          
            rhs =RHOI * GRAV * dz_dpara(i,j) &
                - dperp_dx  * dfdx_3d(vel_perp,i,j,k,dx)&
                - dperp_dy  * dfdy_3d(vel_perp,i,j,k,dy)&
                - dperp_dz  * dfdz_3d_irregular(vel_perp,i,j,k,dz)&
                - dperp_dxy * d2fdxy_3d(vel_perp,i,j,k,dx,dy)&
                - dperp_dxz * d2fdxz_3d(vel_perp,i,j,k,dx,dz)&
                - dperp_dyz * d2fdyz_3d(vel_perp,i,j,k,dy,dz)&
                - dperp_dz2 * d2fdz2_3d_irregular(vel_perp,i,j,k,dz)
        
        end if

    end subroutine sparse_setup

    !Computes finite differences for the marine margin
    subroutine sparse_marine_margin(component,i,j,k,h,normals, dzdx, dzdy, dhbdx, dhbdy, &
                                    vel_perp, efvs, dx, dy, ax, ay, zeta,coef, rhs, &
                                    dudx,dudy,dudz,dvdx,dvdy,dvdz,direction_x,direction_y, geometry_mask, STAGGERED, &
                                    WHICH_SOURCE)
        character(*), intent(in) :: component !*FD Either "u" or "v"
        integer, intent(in) :: i,j,k !*FD Point that the boundary condition is computed for
        real(dp), dimension(:,:), intent(in) :: h !*FD Ice thickness field
        real(dp), dimension(:,:), intent(in) :: normals !*FD Pre-computed angles that are normal to the marine margin
        real(dp), dimension(:,:), intent(in) :: dzdx, dzdy, dhbdx, dhbdy
        real(dp), dimension(:,:,:), intent(in) :: vel_perp !*FD Velocity component that is perpendicular to the one being computed
        real(dp), dimension(:,:,:), intent(in) :: efvs !*FD effective viscosity

        real(dp), intent(in) :: dx, dy !*FD grid spacing in x and y directions
        
        real(dp), intent(in), dimension(:,:,:), target :: ax, ay !*FD coefficients introduced during vertical rescaling
        
        real(dp), dimension(:), intent(in) :: zeta !*FD Irregular grid levels in vertical coordinate
        
        real(dp), dimension(STENCIL_SIZE), intent(inout) :: coef !*FD Output array of stencil coefficients
        real(dp),                intent(out)   :: rhs !*FD Element of the right-hand side vector corresponding to this point
        
        real(dp), dimension(:,:,:), intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz
        
        integer, dimension(:,:), intent(in) :: geometry_mask

        real(dp),dimension(:,:) :: direction_x, direction_y

        logical, intent(in) :: STAGGERED
        integer, intent(in) :: WHICH_SOURCE

        !Whether or not we need to upwind or downwind lateral derivatives
        !Upwind is for when there's no ice on positive side,
        !downwind is for when there's no ice on negative side
        logical :: upwind_x, upwind_y, downwind_x, downwind_y

        !x and y components of the normal vector with unit length
        real(dp), target :: n_x, n_y

        !Balance of hydrostatic and cryostatic pressures at this level
        real(dp) :: source
        
        !Derivative of the perpendicular component (v if computing u and vice
        !versa) with respect to the parallel direction (x if computing u, y if
        !computing v)
        real(dp), target :: dperp_dx
        real(dp), target :: dperp_dy
        real(dp), target :: dperp_dz

        !The derivatives of the velocity component being computed with respect
        !to each dimension.  These end up as coefficients in the stencil
        !definition
        real(dp), target :: dx_coeff, dy_coeff, dz_coeff

        !Cached difference calculations for the nonstaggered Z grid
        real(dp) dz_down1, dz_down2, dz_down3, dz_up1, dz_up2, dz_up3
        real(dp) dz_cen1, dz_cen2, dz_cen3 
        real(dp), dimension(:,:,:), pointer :: a_para, a_perp

        real(dp), pointer :: n_para, n_perp, dperp_dpara, dperp_dperp
        real(dp), pointer :: para_coeff, perp_coeff

        !Surface (upper or lower) derivatives.  These are only nonzero at k=1 or k=nzeta
        real(dp) :: ds_dpara = 0, ds_dperp = 0

        !Get the normal unit vector
        !The angles are defined in radians, with 0 = 12 o' clock, pi/2 = 3 o' clock,
        !etc.
        !We will need to convert this into an angle system with 0 = 3 o' clock,
        !pi/2 = 12 o' clock, etc.
        !At first I thought:
        !theta = normals(i,j) + 2*(PI/4 - normals(i,j))
        !But really this reduces to swapping sine and cosine (easy proof to work
        !through).
        n_x = sin(normals(i,j))
        n_y = -cos(normals(i,j))

       !Clamp n_x and n_y so that small values register as zero
        !(The need to do this arises from the fact that, sin(k*pi)
        !cannot be exactly 0)
        if (abs(n_x) < 1d-10) then
            n_x = 0
        end if

        if (abs(n_y) < 1d-10) then
            n_y = 0
        end if
   
        !Divide the source term by viscosity, less error-prone to do it here
        !rather than multiply everything else by it (though it may give
        !iterative solvers more fits maybe?)
        source = source_term(i, j, k, H, zeta, WHICH_SOURCE, STAGGERED)/efvs(k,i,j)
        ! jvj this is just a trick to see if there is an effective viscosity
        ! instability at the front.
        ! initial tests indicate that there is, ie this is much more stable.
        !source = source_term(i, j, k, H, zeta, WHICH_SOURCE, STAGGERED)**3 * 1.e-18

        !Determine whether to use upwinded or downwinded derivatives for the
        !horizontal directions.
        upwind_x = .false.
        upwind_y = .false.
        downwind_x = .false.
        downwind_y = .false.

#ifdef DEBUG_FIELDS
        pointtype(i,j) = 5
#endif

        if (component=="u") then
            n_para => n_x
            n_perp => n_y

            dperp_dpara => dperp_dx
            dperp_dperp => dperp_dy
            
            a_para => ax
            a_perp => ay

            para_coeff => dx_coeff
            perp_coeff => dy_coeff
        
            dperp_dx = dvdx(k,i,j)
            dperp_dy = dvdy(k,i,j)
            dperp_dz = dvdz(k,i,j)

            !Set the ds terms
            !Not sure if this is actually needed yet
            !if (k == 1) then !Surface
            !    ds_dpara = dzdx(i,j)
            !    ds_dperp = dzdy(i,j)
            !else if (k == size(zeta)) then !Bed
            !    ds_dpara = dhbdx(i,j)
            !    ds_dperp = dhbdy(i,j)
            !end if
        else 
            n_para => n_y
            n_perp => n_x
            
            dperp_dpara => dperp_dy
            dperp_dperp => dperp_dx
            
            a_para => ay
            a_perp => ax

            para_coeff => dy_coeff
            perp_coeff => dx_coeff

            dperp_dx = dudx(k,i,j)
            dperp_dy = dudy(k,i,j)
            dperp_dz = dudz(k,i,j)

            !Set the ds terms
            !if (k == 1) then !Surface
            !    ds_dpara = dzdy(i,j)
            !    ds_dperp = dzdx(i,j)
            !else if (k == size(zeta)) then !Bed
            !    ds_dpara = dhbdy(i,j)
            !    ds_dperp = dhbdx(i,j)
            !end if
        end if
        
        if (direction_y(i,j) > 0) then
           downwind_y = .true.
        else if (direction_y(i,j) < 0) then
            upwind_y = .true.
        end if

        if (direction_x(i,j) > 0) then
            downwind_x = .true.
        else if (direction_x(i,j) < 0) then
            upwind_x = .true.
        end if
       
        para_coeff = 4*(n_para + ds_dpara)
        perp_coeff =   (n_perp + ds_dperp)
        
        dz_coeff   =   (4*a_para(k,i,j)*(n_para + ds_dpara) + a_perp(k,i,j)*(n_perp + ds_dperp))

        rhs = source * n_para &
                - 2*(n_para+ds_dpara)*dperp_dperp & 
                -   (n_perp+ds_dperp)*dperp_dpara &
                -   (2*a_para(k,i,j)*(n_perp+ds_dperp) + a_perp(k,i,j)*(n_para+ds_dpara))*dperp_dz



        !Enter the z component into the finite difference scheme.
        !If we are on the top of bottom of the ice shelf, we will need
        !to upwind or downwind respectively
        !These terms arise from the vertical rescaling of the system;
        !any vertical explicitness of the system its self has been
        !thrown out based on the assumption of plug flow.
        if (k==1) then !Top of the ice shelf
            !Finite difference coefficients for an irregular Z grid, downwinded
            dz_down1=(2.*zeta(k)-zeta(k+1)-zeta(k+2))/(zeta(k+1)-zeta(k))/(zeta(k+2)-zeta(k))
            dz_down2=(zeta(k+2)-zeta(k))/(zeta(k+2)-zeta(k+1))/(zeta(k+1)-zeta(k))
            dz_down3=(zeta(k)-zeta(k+1))/(zeta(k+2)-zeta(k+1))/(zeta(k+2)-zeta(k))
        
            coef(I_J_KP2) = coef(I_J_KP2) + dz_coeff*dz_down3
            coef(I_J_KP1) = coef(I_J_KP1) + dz_coeff*dz_down2
            coef(I_J_K)   = coef(I_J_K)   + dz_coeff*dz_down1
        else if (k==size(zeta)) then !Bottom of the ice shelf
            !Finite difference coefficients for an irregular Z grid, upwinded
            dz_up1=(zeta(k)-zeta(k-1))/(zeta(k-1)-zeta(k-2))/(zeta(k)-zeta(k-2))
            dz_up2=(zeta(k-2)-zeta(k))/(zeta(k)-zeta(k-1))/(zeta(k-1)-zeta(k-2))
            dz_up3=(2.*zeta(k)-zeta(k-1)-zeta(k-2))/(zeta(k)-zeta(k-1))/(zeta(k)-zeta(k-2))

            coef(I_J_KM2) = coef(I_J_KM2) + dz_coeff*dz_up1
            coef(I_J_KM1) = coef(I_J_KM1) + dz_coeff*dz_up2
            coef(I_J_K)   = coef(I_J_K)   + dz_coeff*dz_up3
        else !Middle of the ice shelf, use centered difference
            dz_cen1 = (zeta(k)-zeta(k+1))/(zeta(k)-zeta(k-1))/(zeta(k+1)-zeta(k-1))
            dz_cen2 = (zeta(k+1)-2.*zeta(k)+zeta(k-1))/(zeta(k)-zeta(k-1))/(zeta(k+1)-zeta(k))
            dz_cen3 = (zeta(k)-zeta(k-1))/(zeta(k+1)-zeta(k))/(zeta(k+1)-zeta(k-1))

            coef(I_J_KM1) = coef(I_J_KM1) + dz_coeff*dz_cen1
            coef(I_J_K)   = coef(I_J_K)   + dz_coeff*dz_cen2
            coef(I_J_KP1) = coef(I_J_KP1) + dz_coeff*dz_cen3
        end if
        
        !Apply finite difference approximations of the x and y derivatives.  The
        !coefficients have already been computed above, so we just need to make
        !sure to use the correct difference form (upwind or downwind, etc.)

        if (downwind_y) then
            coef(IP2_J_K) = coef(IP2_J_K) - 0.5 * dy_coeff/dy
            coef(IP1_J_K) = coef(IP1_J_K) + 2.0 * dy_coeff/dy
            coef(I_J_K)   = coef(I_J_K)   - 1.5 * dy_coeff/dy
        else if (upwind_y) then
            coef(I_J_K)   = coef(I_J_K)   + 1.5 * dy_coeff/dy
            coef(IM1_J_K) = coef(IM1_J_K) - 2.0 * dy_coeff/dy
            coef(IM2_J_K) = coef(IM2_J_K) + 0.5 * dy_coeff/dy
        else
            coef(IP1_J_K) = coef(IP1_J_K) + .5*dy_coeff/dy
            coef(IM1_J_K) = coef(IM1_J_K) - .5*dy_coeff/dy  
        end if

        if (downwind_x) then
            coef(I_JP2_K) = coef(I_JP2_K) - 0.5 * dx_coeff/dx
            coef(I_JP1_K) = coef(I_JP1_K) + 2.0 * dx_coeff/dx
            coef(I_J_K)   = coef(I_J_K)   - 1.5 * dx_coeff/dx 
        else if (upwind_x) then
            coef(I_J_K)   = coef(I_J_K)   + 1.5 * dx_coeff/dx
            coef(I_JM1_K) = coef(I_JM1_K) - 2.0 * dx_coeff/dx
            coef(I_JM2_K) = coef(I_JM2_K) + 0.5 * dx_coeff/dx
        else
            coef(I_JP1_K) = coef(I_JP1_K) + 0.5 * dx_coeff/dx
            coef(I_JM1_K) = coef(I_JM1_K) - 0.5 * dx_coeff/dx  
        end if

    end subroutine

    !Computes the source term for an ice shelf from hydrostatic and cryostatic pressure
    function source_term(i, j, k, H, zeta, which_source, staggered)
        integer, intent(in) :: i,j,k
        real(dp),dimension(:,:), intent(in) :: H
        real(dp), dimension(:), intent(in) :: zeta
        
        integer, intent(in) :: which_source
        logical, intent(in) :: staggered

        real(dp) :: source_term
        
        real(dp) :: integrate_from

        !Determine the hydrostatic pressure at the current location
        !If we are at the part of the ice shelf above the water, 
        !then we are at 1ATM pressure which we just assume to be 0 
        !throughout the model anyway

        source_term = 0

        if (which_source == HO_SOURCE_AVERAGED) then
            !Vertically averaged model
            if (k == 1) then
                source_term = 0
            else
                source_term = .5 * rhoi * grav * h(i,j) * (1 - rhoi/rhoo)
            end if

        else if (which_source == HO_SOURCE_EXPLICIT) then
            !In order to get the correct source term, I integrate from the level above this
            !one to the current level so that all of the pressure is accounted for.
            if (k == 1) then
                source_term = 0
            else
                source_term = .5 * rhoi * grav * H(i,j) * (zeta(k)**2 - zeta(k-1)**2)

                !When we consider the contribution of hydrostatic pressure to the 
                !source term, we only want to integrate that part of the current
                !layer that is below sea level.
                integrate_from = max(1-rhoi/rhoo, zeta(k-1))

                if (zeta(k) > (1 - rhoi/rhoo)) then
                    source_term = source_term + rhoo * grav * H(i,j) * &
                        ( ( 1-rhoi/rhoo) * (zeta(k) - integrate_from) - &
                        .5*(zeta(k)**2 - integrate_from**2) )
                end if
                source_term = source_term / (zeta(k) - zeta(k-1))
            end if
        end if  

    end function

!
!
!------------------------------------------------------
!   Calculation of stress field based on velocities
!------------------------------------------------------
!
      SUBROUTINE stressf(efvs,uvel,vvel,flwa,h,ax,ay,dx,dy,&
        dz,txz,tyz,txx,tyy,txy,FLOWN,ZIP,PERIODIC_X, PERIODIC_Y)
!
        INTEGER MAXY,MAXX,NZETA
        real(dp) :: FLOWN,ZIP
        real(dp), dimension(:,:,:), intent(inout) :: efvs
        real(dp), dimension(:,:,:), intent(in) :: uvel
        real(dp), dimension(:,:,:), intent(in) :: vvel
        real(dp), dimension(:,:,:), intent(in) :: flwa
        real(dp), dimension(:,:),   intent(in) :: h
        real(dp), dimension(:,:,:), intent(in) :: ax
        real(dp), dimension(:,:,:), intent(in) :: ay
        real(dp), dimension(:,:,:), intent(out) :: txz 
        real(dp), dimension(:,:,:), intent(out) :: tyz
        real(dp), dimension(:,:,:), intent(out) :: txx
        real(dp), dimension(:,:,:), intent(out) :: tyy
        real(dp), dimension(:,:,:), intent(out) :: txy
        real(dp),                   intent(in) :: dx
        real(dp),                   intent(in) :: dy
        real(dp), dimension(:),     intent(in) :: dz
        
        logical :: PERIODIC_X, PERIODIC_Y

      INTEGER i,j,k
      real(dp), dimension(:,:,:), allocatable :: dudz, dvdz, dudx, dvdx, dudy, dvdy

      real(dp) :: exx,eyy,exy,exz,eyz,eeff
      real(dp) :: macht
!
      
      macht=(1.-FLOWN)/(2.*FLOWN)
      
      MAXY = size(uvel, 3)
      MAXX = size(uvel, 2)
      NZETA = size(dz)
      
            !Allocate temporary derivative data
      allocate(dudx(NZETA, MAXX, MAXY))
      allocate(dudy(NZETA, MAXX, MAXY))
      allocate(dudz(NZETA, MAXX, MAXY))
      allocate(dvdx(NZETA, MAXX, MAXY))
      allocate(dvdy(NZETA, MAXX, MAXY))
      allocate(dvdz(NZETA, MAXX, MAXY))
      
      
      !TODO: Figure out what the UPSTREAM flag does!!
      !Because the values being derived (uvel and vvel) already
      !exist on the staggered grid, we use the unstaggered version
      !of the field derivative function
      call df_field_3d(uvel, dx, dy, dz, dudy, dudx, dudz)
      call df_field_3d(vvel, dx, dy, dz, dvdy, dvdx, dvdz)
      !TODO: Can we somehow avoid the code duplication with calcefvs?
      do  i=1,MAXX
        do  j=1,MAXY
          do  k=1,NZETA
            if (h(i,j).gt.0.) then
              exx=dudx(k,i,j)+ax(k,i,j)*dudz(k,i,j)
              eyy=dvdy(k,i,j)+ay(k,i,j)*dvdz(k,i,j)
              exy=0.5*((dudy(k,i,j)+ay(k,i,j)*dudz(k,i,j))+(dvdx(k,i,j)+ax(k,i,j)*dvdz(k,i,j)))
              exz=-0.5*dudz(k,i,j)/h(i,j)
              eyz=-0.5*dvdz(k,i,j)/h(i,j)
            else
              exx=dudx(k,i,j)
              eyy=dvdy(k,i,j)
              exy=0.5*(dudy(k,i,j)+dvdx(k,i,j))
              exz=0.
              eyz=0.
            endif
            eeff=(exx*exx)+(eyy*eyy)+(exx*eyy)+(exy*exy)+(exz*exz)+&
              (eyz*eyz)
            if (eeff < ZIP) eeff=ZIP
            efvs(k,i,j)=0.5*(flwa(k,i,j)**(-1./FLOWN))*(eeff**macht)
            txz(k,i,j)=2.0*efvs(k,i,j)*exz
            tyz(k,i,j)=2.0*efvs(k,i,j)*eyz
            txx(k,i,j)=2.0*efvs(k,i,j)*exx
            tyy(k,i,j)=2.0*efvs(k,i,j)*eyy
            txy(k,i,j)=2.0*efvs(k,i,j)*exy
        end do
      end do
    end do


      !Adjust boundaries if using periodic boundary conditions.
      !REFACTORED CODE BELOW
      call periodic_boundaries_3d(txz,PERIODIC_X,PERIODIC_Y)
      call periodic_boundaries_3d(tyz,PERIODIC_X,PERIODIC_Y)
      call periodic_boundaries_3d(txx,PERIODIC_X,PERIODIC_Y)
      call periodic_boundaries_3d(tyy,PERIODIC_X,PERIODIC_Y)
      call periodic_boundaries_3d(txy,PERIODIC_X,PERIODIC_Y)
      
            !Free temporary derivative data
      deallocate(dudx)
      deallocate(dudy)
      deallocate(dudz)
      deallocate(dvdx)
      deallocate(dvdy)
      deallocate(dvdz)
      
      return
      END subroutine

#ifdef OUTPUT_PARTIAL_ITERATIONS
!------------------------------------------------------------------------------
!The following are methods to write iterations to a NetCDF file.  
!The NetCDF file used is not managed through Glimmer's NetCDF interface
!I know this is code duplication, it is for DEBUG PURPOSES ONLY
!------------------------------------------------------------------------------
function begin_iteration_debug(nx, ny, nz)
    use netcdf
    use glimmer_ncdf, only: nc_errorhandle
    integer :: nx, ny, nz
    integer :: err
    integer :: ncid

    integer :: xdim, ydim, zdim, iterdim, dims(4), dims2d(2), varid
    integer :: begin_iteration_debug

    ncid = 0
    err = nf90_create("iterdebug.nc", NF90_CLOBBER, ncid)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_dim(ncid, "x", nx, xdim)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_dim(ncid, "y", ny, ydim)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_dim(ncid, "z", nz, zdim)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_dim(ncid, "iter", NF90_UNLIMITED, iterdim)
    call nc_errorhandle(__FILE__, __LINE__, err)

    dims = (/xdim, ydim, zdim, iterdim/)
    dims2d = (/xdim, ydim/)

    err = nf90_def_var(ncid, "efvs", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_var(ncid, "uvel", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_var(ncid, "vvel", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
    err = nf90_def_var(ncid, "velnorm", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 

    err = nf90_def_var(ncid, "mask", NF90_INT, dims2d, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
    err = nf90_def_var(ncid, "bc_normal", NF90_DOUBLE, dims2d, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
  
    err = nf90_enddef(ncid)
    call nc_errorhandle(__FILE__, __LINE__, err)

    begin_iteration_debug=ncid


end function

subroutine iteration_debug_step(ncid, iter, efvs, uvel, vvel, geometry_mask, marine_bc_normal)
    use netcdf
    use glimmer_ncdf, only: nc_errorhandle

    integer :: ncid, iter
    real(dp), dimension(:,:,:) :: efvs, uvel, vvel
    integer, dimension(:,:) :: geometry_mask
    real(dp), dimension(:,:) :: marine_bc_normal

    integer :: varid, err

    integer :: nx, ny, nz
    integer :: i,j,k, start(4), count(4), start2d(2)
    
    nx = size(efvs, 2)
    ny = size(efvs, 3)
    nz = size(efvs, 1)

    count = (/1,1,1,1/)
    do i = 1,nx
        do j = 1,ny
            do k = 1,nz
                start=(/i,j,k,iter+1/)
       err = nf90_inq_varid(ncid, "efvs", varid)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_put_var(ncid, varid, efvs(k,i,j), start)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_inq_varid(ncid, "uvel", varid)
       call nc_errorhandle(__FILE__, __LINE__, err) 
       
       err = nf90_put_var(ncid, varid, uvel(k,i,j), start)
       
       err = nf90_inq_varid(ncid, "vvel", varid)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_put_var(ncid, varid, vvel(k,i,j), start)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_inq_varid(ncid, "velnorm", varid)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_put_var(ncid, varid, sqrt(uvel(k,i,j)**2 + vvel(k,i,j)**2), start)
       call nc_errorhandle(__FILE__, __LINE__, err)

        end do
       if (iter .EQ. 1) then
           start2d=(/i,j/)
           err = nf90_inq_varid(ncid, "mask", varid)
           call nc_errorhandle(__FILE__, __LINE__, err)
           err = nf90_put_var(ncid, varid, geometry_mask(i,j), start2d)
           call nc_errorhandle(__FILE__, __LINE__, err)

           err = nf90_inq_varid(ncid, "bc_normal", varid)
           call nc_errorhandle(__FILE__, __LINE__, err)
           err = nf90_put_var(ncid, varid, marine_bc_normal(i,j), start2d)
           call nc_errorhandle(__FILE__, __LINE__, err)
       end if
        end do
    end do
   
    !Close and open the dataset so that changes get written to it
    !if (mod(iter+1, 10) == 0) then
        call end_debug_iteration(ncid)

        err = nf90_open("iterdebug.nc", ior(NF90_WRITE,NF90_SHARE), ncid)
        call nc_errorhandle(__FILE__, __LINE__, err)
    !end if

end subroutine

subroutine end_debug_iteration(ncid)
    use netcdf
    use glimmer_ncdf, only: nc_errorhandle

    integer :: ncid

    integer :: err

    err = nf90_close(ncid)
    call nc_errorhandle(__FILE__, __LINE__, err)
end subroutine
#endif

subroutine write_xls_direction_guide(filename, nx, ny)
    integer :: nx, ny
    character(*) :: filename
    
    real(dp), dimension(nx, ny) :: field

    field = 0
    field(1,1) = 1
    field(2,1) = 1

    call write_xls(filename, field)

end subroutine


!
!---------------------------------------------------
!
!                 END OF SUBROUTINES
!
!---------------------------------------------------
end module ice3d_lib

