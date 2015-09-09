
#ifdef HAVE_CONFIG_H 
#include "config.inc" 
#endif
#include "glide_nan.inc"
#include "glide_mask.inc"

#define shapedbg(x) write(*,*) "x", shape(x)

module glide_velo_higher
    !Includes for the higher-order velocity computations that this calls out to
    use ice3d_lib
    use glam_strs2, only: glam_velo_fordsiapstr, JFNK, umask

    !globals
    use glimmer_global, only : dp
    use glimmer_paramets, only : vis0, vis0_glam 
    use glimmer_physcon, only: gn

    !Other modules that this needs to call out to
    use glide_types
    use glide_vertint
    use glide_grids, only: stagvarb
    use glide_mask
    use glide_grids
    implicit none
    
    !TODO: Parameterize the following globals
    real(dp), parameter :: VEL2ERR  = 1e-4
    real(dp), parameter :: TOLER    = 1e-6
    integer,  parameter :: CONVT    = 4
    real(dp), parameter :: SHTUNE   = 1.D-16
    integer,  parameter :: UPSTREAM = 0
    integer,  parameter :: MANIFOLD = 1
contains
        
    subroutine init_velo_hom_pattyn(model)
        type(glide_global_type) :: model

        !Init the beta field based on the selected option
        !If we selected to use 1/btrc, this needs to be computed later
        !If we selected to read the beta field, it's been done already and no
        !action is needed
        select case (model%options%which_ho_beta_in)
            case(HO_BETA_ALL_NAN)
                model%velocity_hom%beta = NAN
            case(HO_BETA_USE_SOFT)
                where (model%velocity%bed_softness /= 0)
                    model%velocity_hom%beta = 1 / model%velocity%bed_softness
                elsewhere
                    model%velocity_hom%beta = NAN
                end where
        end select

    end subroutine

    subroutine calc_slip_ratio(flwa, thick, slip_ratio_coef, beta)
        real(dp), dimension(:,:), intent(in) :: flwa !*FD Glen's A at the base
        real(dp), dimension(:,:), intent(in) :: thick !*FD Ice thickness
        real(dp),                 intent(in) :: slip_ratio_coef
        real(dp), dimension(:,:), intent(out):: beta

        beta = 1.0D0/(slip_ratio_coef*flwa*thick)
    end subroutine  

    !This is a temporary wrapper function to get all the HO setup and calling
    !code out of glide_thck so it keeps the clutter down.  I'm just passing the
    !model right now because the analysis on what needs to be passed has yet to
    !be done.
    subroutine run_ho_diagnostic(model)
        use glide_thckmask
        use glide_mask
        
        type(glide_global_type),intent(inout) :: model
        !For HO masking
        logical :: empty
        integer :: totpts
        integer, save :: tstep ! JFL to be removed
        real(sp), dimension(model%general%ewn-1, model%general%nsn-1) :: stagmassb

        !TEMPORARY arrays, these should at some point be placed in Model
        !probably
        integer, dimension(model%general%ewn-1, model%general%nsn-1)  :: geom_mask_stag
        real(dp), dimension(model%general%ewn-1, model%general%nsn-1) :: latbc_norms_stag

        tstep = tstep + 1 ! JFL to be removed
        ! NOTE that choice of non-linear solver, "NL_solver" flag (1 = Picard; 2 = JFNK),  
        ! is now input as run time option (model%options%which_ho_nonlinear)        

        !Beta field computations that change in time
        if (model%options%which_ho_beta_in == HO_BETA_USE_BTRC) then
           where (model%velocity%btrc /= 0)
               model%velocity_hom%beta = 1/model%velocity%btrc
            elsewhere
                model%velocity_hom%beta = NAN
            end where
        else if (model%options%which_ho_beta_in == HO_BETA_SLIP_RATIO) then
            call calc_slip_ratio(model%temper%flwa(model%general%upn,:,:), model%geomderv%stagthck, &
                                 model%paramets%slip_ratio, model%velocity_hom%beta)
        end if

        !Compute the "geometry mask" (type of square) for the staggered grid

        call glide_set_mask(model%numerics, model%geomderv%stagthck, model%geomderv%stagtopg, &
                                model%general%ewn-1, model%general%nsn-1, model%climate%eus, &
                                geom_mask_stag) 

        !Augment masks with kinematic boundary condition info
        call augment_kinbc_mask(model%geometry%thkmask, model%velocity_hom%kinbcmask)
        call augment_kinbc_mask(geom_mask_stag, model%velocity_hom%kinbcmask)

        !Compute the normal vectors to the marine margin for the staggered grid
        call glide_marine_margin_normal(model%geomderv%stagthck, geom_mask_stag, latbc_norms_stag)

        ! save the final mask to 'dynbcmask' for exporting to netCDF output file
        model%velocity_hom%dynbcmask = geom_mask_stag

        if (model%options%which_ho_diagnostic == HO_DIAG_PATTYN_STAGGERED) then
#ifdef VERY_VERBOSE
            write(*,*)"Running Pattyn staggered"
#endif
            !Compute the "point mask" (mask that assigns a unique value to each
            !grid point on which ice dynamics are enabled) for the
            !staggered grid.  
            stagmassb = 0

            
            call glide_maskthck(model%geomderv%stagthck, stagmassb, .true., model%numerics%thklim,&
                                model%geometry%dom, model%velocity_hom%velmask, totpts, empty)
                 
               
            call velo_hom_pattyn(model%general%ewn, model%general%nsn, model%general%upn, &
                                 model%numerics%dew, model%numerics%dns, model%numerics%sigma, &
                                 model%geomderv%stagthck, model%geomderv%stagusrf, model%geomderv%staglsrf, &
                                 model%geomderv%dthckdew, model%geomderv%dthckdns, &
                                 model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                                 model%geomderv%dlsrfdew, model%geomderv%dlsrfdns, & 
                                 model%geomderv%d2usrfdew2, model%geomderv%d2usrfdns2, &
                                 model%geomderv%d2thckdew2, model%geomderv%d2thckdns2, &
                                 model%velocity_hom%velmask, totpts, &
                                 geom_mask_stag, &
                                 model%temper%flwa, real(gn, dp), model%velocity_hom%beta, &
                                 model%options, &
                                 latbc_norms_stag,&
                                 model%velocity_hom%uvel, model%velocity_hom%vvel, &
                                 model%velocity_hom%is_velocity_valid, &
                                 model%velocity_hom%uflx, model%velocity_hom%vflx, &
                                 model%velocity_hom%efvs, model%velocity_hom%tau, &
                                 model%velocity_hom%gdsx, model%velocity_hom%gdsy)

            call hom_diffusion_pattyn(model%velocity_hom%uvel, model%velocity_hom%vvel, model%geomderv%stagthck, &
                                      model%geomderv%dusrfdew, model%geomderv%dusrfdns, model%numerics%sigma, &
                                      model%velocity_hom%diffu_x, model%velocity_hom%diffu_y, &
                                      model%velocity_hom%uflx,    model%velocity_hom%vflx) 

        else if (model%options%which_ho_diagnostic == HO_DIAG_PATTYN_UNSTAGGERED) then
#ifdef VERY_VERBOSE
            write(*,*)"Running Pattyn unstaggered"
#endif
             call velo_hom_pattyn_nonstag(model%general%ewn, model%general%nsn, model%general%upn, &
                                         model%numerics%dew, model%numerics%dns, model%numerics%sigma, &
                                         model%geometry%thck, model%geometry%usrf, model%geometry%lsrf, &
                                         model%geomderv%dusrfdew_unstag, model%geomderv%dusrfdns_unstag, &
                                         model%geomderv%dthckdew_unstag, model%geomderv%dthckdns_unstag, &
                                         model%geomderv%dlsrfdew_unstag, model%geomderv%dlsrfdns_unstag, &
                                         model%geomderv%d2usrfdew2_unstag, model%geomderv%d2usrfdns2_unstag, &
                                         model%geomderv%d2thckdew2_unstag, model%geomderv%d2thckdns2_unstag, &
                                         model%geometry%mask, model%geometry%totpts, &
                                         model%geometry%thkmask, &
                                         model%temper%flwa, real(gn, dp), &
                                         model%velocity_hom%beta, model%options, &
                                         model%geometry%marine_bc_normal, &
                                         model%velocity_hom%uvel, model%velocity_hom%vvel, &
                                         model%velocity_hom%is_velocity_valid)

            call hom_diffusion_pattyn(model%velocity_hom%uvel, model%velocity_hom%vvel, model%geomderv%stagthck, &
                                      model%geomderv%dusrfdew, model%geomderv%dusrfdns, model%numerics%sigma, &
                                      model%velocity_hom%diffu_x, model%velocity_hom%diffu_y, &
                                      model%velocity_hom%uflx,    model%velocity_hom%vflx) 
            
        else if (model%options%which_ho_diagnostic == HO_DIAG_PP) then

          if ( model%options%which_ho_nonlinear == HO_NONLIN_PICARD ) then ! Picard (standard solver)

            call glam_velo_fordsiapstr( model%general%ewn,       model%general%nsn,                 &
                                        model%general%upn,                                          &
                                        model%numerics%dew,      model%numerics%dns,                &
                                        model%numerics%sigma,    model%numerics%stagsigma,          &
                                        model%geometry%thck,     model%geometry%usrf,               &
                                        model%geometry%lsrf,     model%geometry%topg,               &
                                        model%geomderv%dthckdew, model%geomderv%dthckdns,           &
                                        model%geomderv%dusrfdew, model%geomderv%dusrfdns,           &
                                        model%geomderv%dlsrfdew, model%geomderv%dlsrfdns,           & 
                                        model%geomderv%stagthck, model%temper%flwa*vis0/vis0_glam,  &
                                        model%basalproc%minTauf,                                    & 
                                        model%velocity_hom%btraction,                               & 
                                        geom_mask_stag,                                             &
                                        model%options%which_ho_babc,                                &
                                        model%options%which_ho_efvs,                                &
                                        model%options%which_ho_resid,                               &
                                        model%options%which_ho_nonlinear,                           &
                                        model%options%which_ho_sparse,                              &
                                        model%options%periodic_ew,                                  &
                                        model%options%periodic_ns,                                  &
                                        model%velocity_hom%beta,                                    & 
                                        model%velocity_hom%uvel, model%velocity_hom%vvel,           &
                                        model%velocity_hom%uflx, model%velocity_hom%vflx,           &
                                        model%velocity_hom%efvs, tstep)

          else if ( model%options%which_ho_nonlinear == HO_NONLIN_JFNK ) then ! JFNK (solver in development...)

            call JFNK                  ( model%general%ewn,       model%general%nsn,                 &
                                        model%general%upn,                                          &
                                        model%numerics%dew,      model%numerics%dns,                &
                                        model%numerics%sigma,    model%numerics%stagsigma,          &
                                        model%geometry%thck,     model%geometry%usrf,               &
                                        model%geometry%lsrf,     model%geometry%topg,               &
                                        model%geomderv%dthckdew, model%geomderv%dthckdns,           &
                                        model%geomderv%dusrfdew, model%geomderv%dusrfdns,           &
                                        model%geomderv%dlsrfdew, model%geomderv%dlsrfdns,           & 
                                        model%geomderv%stagthck, model%temper%flwa*vis0/vis0_glam,  &
                                        model%basalproc%minTauf,                                    & 
                                        model%velocity_hom%btraction,                               & 
                                        geom_mask_stag,                                             &
                                        model%options%which_ho_babc,                                &
                                        model%options%which_ho_efvs,                                &
                                        model%options%which_ho_resid,                               &
                                        model%options%which_ho_nonlinear,                           &
                                        model%options%which_ho_sparse,                              &
                                        model%options%periodic_ew,                                  &
                                        model%options%periodic_ns,                                  &
                                        model%velocity_hom%beta,                                    & 
                                        model%velocity_hom%uvel, model%velocity_hom%vvel,           &
                                        model%velocity_hom%uflx, model%velocity_hom%vflx,           &
                                        model%velocity_hom%efvs, tstep)
           else
              call write_log('Invalid which_ho_nonlinear option.',GM_FATAL)
            end if

        end if
        !Compute the velocity norm - this is independant of the methods used to compute the u and v components so
        !we put it out here
        model%velocity_hom%velnorm = sqrt(model%velocity_hom%uvel**2 + model%velocity_hom%vvel**2)
        model%velocity_hom%is_velocity_valid = .true.
        
    end subroutine

    subroutine velo_hom_pattyn(ewn, nsn, upn, dew, dns, sigma, &
                               thck, usrf, lsrf, dthckdew, dthckdns, dusrfdew, dusrfdns, &
                               dlsrfdew, dlsrfdns, &
                               d2zdx2, d2zdy2, d2hdx2, d2hdy2, &
                               point_mask, totpts, geometry_mask, flwa, flwn, btrc, &
                               options, &
                               marine_bc_normal, &
                               uvel, vvel, valid_initial_guess, uflx, vflx, efvs, tau, gdsx, gdsy)
                           
        integer, intent(in)  :: ewn !*FD Number of cells X
        integer, intent(in)  :: nsn !*FD Number of cells Y
        integer, intent(in)  :: upn !*FD Number of cells Z
        real(dp), intent(in) :: dew
        real(dp), intent(in) :: dns
        real(dp), dimension(:), intent(in) :: sigma !*FD Sigma coord for rescaled Z dimension
        real(dp), dimension(:,:), intent(in) :: thck !*FD Thickness, on staggered grid
        real(dp), dimension(:,:), intent(in) :: usrf !*FD Upper surface profile
        real(dp), dimension(:,:), intent(in) :: lsrf !*FD Lower surface profile
        real(dp), dimension(:,:), intent(in) :: dthckdew !*FD X thickness gradient
        real(dp), dimension(:,:), intent(in) :: dthckdns !*FD Y thickness gradient
        real(dp), dimension(:,:), intent(in) :: dusrfdew !*FD X surface gradient
        real(dp), dimension(:,:), intent(in) :: dusrfdns !*FD Y surface gradient
        real(dp), dimension(:,:), intent(in) :: dlsrfdew !*FD X bed gradient
        real(dp), dimension(:,:), intent(in) :: dlsrfdns !*FD Y bed gradient
        real(dp), dimension(:,:), intent(in) :: d2zdx2, d2zdy2, d2hdx2, d2hdy2
        integer,  dimension(:,:), intent(in) :: point_mask     !*FD Numbers points in the staggered grid that are included in computation
        integer, intent(in) :: totpts
        integer, dimension(:,:), intent(in) :: geometry_mask
        real(dp), dimension(:,:,:), intent(in) :: flwa !*FD Glen's A (rate factor) - Used for thermomechanical coupling
        real(dp), dimension(:,:), intent(in)   :: btrc !*FD Basal Traction, either betasquared or tau0
        real(dp), dimension(:,:), intent(in)   :: marine_bc_normal
        real(dp), intent(in) :: flwn !*FD Exponent in Glenn power law
        type(glide_options), intent(in) :: options
        real(dp), dimension(:,:,:), intent(inout) :: uvel 
        real(dp), dimension(:,:,:), intent(inout) :: vvel
        logical, intent(in) :: valid_initial_guess !*Whether or not the given uvel or vvel are appropriate initial guesses.  If not we'll have to roll our own.
        real(dp), dimension(:,:), intent(out) :: uflx
        real(dp), dimension(:,:), intent(out) :: vflx
        real(dp), dimension(:,:,:), intent(out) :: efvs !*FD Effective viscosity
        type(glide_tensor), intent(inout)         :: tau
        real(dp), dimension(:,:,:), intent(out) :: gdsx !*FD X driving stress
        real(dp), dimension(:,:,:), intent(out) :: gdsy !*FD Y driving stress

        !Second derivative of surface
                !Arrays for rescaled coordinate parameters
        real(dp), dimension(upn, ewn-1, nsn-1) :: ax, ay, bx, by, cxy
        
        real(dp), dimension(ewn, nsn)::u,v

        real(dp), dimension(ewn-1, nsn-1) :: direction_x, direction_y
        real(dp), dimension(upn, ewn-1, nsn-1) :: flwa_stag

        real(dp), dimension(upn, ewn-1, nsn-1) :: kinematic_bc_u, kinematic_bc_v
    
        integer :: k


        direction_x = 0
        direction_y = 0


        !Construct fields that contain kinematic boundary condition velocities where they are specified,
        !and NaN where kinematic boundaries should be computed
        !TODO: Integrate this more elegantly (using the mask in Pattyn's model rather than construct the old-style
        !fields outside of it)
        do k = 1, upn
            where (GLIDE_IS_DIRICHLET_BOUNDARY(geometry_mask))
                kinematic_bc_u(k,:,:) = uvel(k,:,:)
                kinematic_bc_v(k,:,:) = vvel(k,:,:)
            elsewhere
                kinematic_bc_u(k,:,:) = NaN
                kinematic_bc_v(k,:,:) = NaN
            endwhere
        end do

        call stagvarb_3d(flwa, flwa_stag, ewn, nsn, upn)

        !Compute a new geometry mask based on staggered geometry
       
        !Compute rescaled coordinate parameters (needed because Pattyn uses an
        !irregular Z grid and scales so that 0 is the surface, 1 is the bed)
        call init_rescaled_coordinates(dthckdew,dlsrfdew,dthckdns,dlsrfdns,usrf,thck,lsrf,&
                                               dusrfdew,dusrfdns,d2zdx2,d2zdy2,d2hdx2,d2hdy2,&
                                               sigma,ax,ay,bx,by,cxy,dew,dns,direction_x,direction_y)
        !"Spin up" estimate with Pattyn's SIA model runs if we don't already
        !have a good initial guess
        if (.not. valid_initial_guess) then
            call veloc1(dusrfdew, dusrfdns, thck, flwa_stag, sigma, uvel, vvel, u, v, ewn-1, nsn-1, upn, &
                        FLWN, options%periodic_ew, options%periodic_ns)
            !If we are performing the plastic bed iteration, the SIA is not
            !enough and we need to spin up a better estimate by shoehorning the
            !tau0 values into a linear bed estimate
            if (options%which_ho_bstress == HO_BSTRESS_PLASTIC) then
                call veloc2(efvs, uvel, vvel, flwa_stag, dusrfdew, dusrfdns, thck, ax, ay, &
                        sigma, bx, by, cxy, btrc/100, dlsrfdew, dlsrfdns, FLWN, ZIP, VEL2ERR, &
                        TOLER, options, .true., dew, dns,point_mask, &
                        totpts,geometry_mask, &
                        kinematic_bc_u, kinematic_bc_v, marine_bc_normal)
            end if
        end if

        !Higher order velocity estimation
        !I am assuming that efvs (effective viscosity) is the same as mu
        !A NOTE ON COORDINATE TRANSPOSITION:
        !Because of the transposition, ewn=maxy and nsn=maxx.  However, veloc2
        !passes maxy in *first*, so these really get passed in the same order
        !that they normally would.
        call veloc2(efvs, uvel, vvel, flwa_stag, dusrfdew, dusrfdns, thck, ax, ay, &
                    sigma, bx, by, cxy, btrc, dlsrfdew, dlsrfdns, FLWN, ZIP, VEL2ERR, &
                    TOLER, options, .true., dew, dns, &
                    point_mask,totpts, geometry_mask, kinematic_bc_u, kinematic_bc_v, marine_bc_normal)
       
        !Final computation of stress field for output
        !call stressf(mu_t, uvel_t, vvel_t, flwa_t, stagthck_t, ax, ay, dew, dns, sigma, & 
        !             tau_xz_t, tau_yz_t, tau_xx_t, tau_yy_t, tau_xy_t, flwn, zip, options%periodic_ew, options%periodic_ns) 
    
    end subroutine velo_hom_pattyn

    subroutine hom_diffusion_pattyn(uvel_hom, vvel_hom, stagthck, dusrfdew, dusrfdns, sigma, diffu_x, diffu_y, uflx, vflx)
        !*FD Estimate of higher-order diffusivity vector given Pattyn's approach
        !*FD Implements (54) in Pattyn 2003, including vertical integration
        !*FD Also returns fluxes computed from the velocities, as these are part of the
        !*FD HO diffusivity computation anyway
        real(dp), dimension(:,:,:), intent(in) :: uvel_hom
        real(dp), dimension(:,:,:), intent(in) :: vvel_hom
        real(dp), dimension(:,:),   intent(in) :: stagthck
        real(dp), dimension(:,:),   intent(in) :: dusrfdew
        real(dp), dimension(:,:),   intent(in) :: dusrfdns
        real(dp), dimension(:),     intent(in) :: sigma
        real(dp), dimension(:,:),   intent(out):: diffu_x
        real(dp), dimension(:,:),   intent(out):: diffu_y
        real(dp), dimension(:,:),   intent(out):: uflx
        real(dp), dimension(:,:),   intent(out):: vflx
        !Local Variables
        real(dp), dimension(size(uvel_hom, 2), size(uvel_hom, 3)) :: u_int !*FD Vertically integrated u velocity
        real(dp), dimension(size(uvel_hom, 2), size(uvel_hom, 3)) :: v_int !*FD Vertically integrated v velocity

        call vertint_output2d(uvel_hom, u_int, sigma)
        call vertint_output2d(vvel_hom, v_int, sigma)

        uflx = u_int*stagthck
        vflx = v_int*stagthck

        diffu_x = uflx*dusrfdew
        diffu_y = vflx*dusrfdns

    end subroutine hom_diffusion_pattyn


    !*FD Version of higher order velo call that computes velocities on ice grid
    !*FD initially and then staggers them.  This might introduce less error than
    !computing velocities on a staggered grid; it is thought that the averaging
    !there causes problems.
    subroutine velo_hom_pattyn_nonstag(ewn, nsn, upn, dew, dns, sigma, thck, usrf, lsrf, &
                              dusrfdew, dusrfdns, dthckdew, dthckdns, dlsrfdew, dlsrfdns, &
                              d2zdx2, d2zdy2, d2hdx2, d2hdy2, point_mask, totpts, geometry_mask, &
                              flwa, flwn, btrc, options, marine_bc_normal,&
                              uvel, vvel, valid_initial_guess)
                            
        integer, intent(in)  :: ewn !*FD Number of cells X
        integer, intent(in)  :: nsn !*FD Number of cells Y
        integer, intent(in)  :: upn !*FD Number of cells Z
        real(dp), intent(in) :: dew !*FD Grid spacing X
        real(dp), intent(in) :: dns !*FD Grid spacing Y
        real(dp), dimension(:), intent(in) :: sigma !*FD Sigma coord for rescaled Z dimension
        real(dp), dimension(:,:), intent(in) :: thck !*FD Thickness, on non-staggered grid
        real(dp), dimension(:,:), intent(in) :: usrf !*FD Upper surface profile
        real(dp), dimension(:,:), intent(in) :: lsrf !*FD Lower surface profile
        real(dp), dimension(:,:), intent(in) :: dthckdew !*FD X thickness gradient
        real(dp), dimension(:,:), intent(in) :: dthckdns !*FD Y thickness gradient
        real(dp), dimension(:,:), intent(in) :: dusrfdew !*FD X surface gradient
        real(dp), dimension(:,:), intent(in) :: dusrfdns !*FD Y surface gradient
        real(dp), dimension(:,:), intent(in) :: dlsrfdew !*FD X bed gradient
        real(dp), dimension(:,:), intent(in) :: dlsrfdns !*FD Y bed gradient
        real(dp), dimension(:,:), intent(in) :: d2zdx2, d2zdy2, d2hdx2, d2hdy2
        integer,  dimension(:,:), intent(in) :: point_mask     !*FD Numbers points in the staggered grid that are included in computation
        integer, intent(in) :: totpts
        integer, dimension(:,:), intent(in) :: geometry_mask
        real(dp), dimension(:,:,:), intent(in) :: flwa !*FD Glen's A (rate factor) - Used for thermomechanical coupling
        real(dp), dimension(:,:), intent(in)   :: btrc !*FD Basal Traction, either betasquared or tau0
        real(dp), dimension(:,:), intent(in)   :: marine_bc_normal
        real(dp), intent(in) :: flwn !*FD Exponent in Glenn power law
        type(glide_options), intent(in) :: options 
        real(dp), dimension(:,:,:), intent(out) :: uvel 
        real(dp), dimension(:,:,:), intent(out) :: vvel
        logical, intent(in) :: valid_initial_guess !*Whether or not the given uvel or vvel are appropriate initial guesses.  If not we'll have to roll our own.
        
        !Arrays for rescaled coordinate parameters
        real(dp), dimension(upn, ewn, nsn) :: ax, ay, bx, by, cxy
        
        real(dp), dimension(ewn, nsn)::u,v

        !whether to upwind or downwind derivatives.
        real(dp), dimension(ewn, nsn) :: direction_x, direction_y
        
        !Arrays to hold unstaggered data
        real(dp), dimension(ewn,nsn) :: btrc_unstag
        
        real(dp), dimension(upn, ewn, nsn) :: kinematic_bc_u_unstag, kinematic_bc_v_unstag
        
        real(dp), dimension(upn, ewn, nsn) :: efvs, uvel_unstag, vvel_unstag
        
        integer :: k 

        !Determine whether to upwind or downwind derivatives at points on the
        !interior of the model domain (this is mainly important for the marine
        !margin
        !Since this is a call to Pattyn's model, derivatives need to be
        !transposed.
        direction_y = 0
        direction_x = 0
        call upwind_from_mask(geometry_mask, direction_x, direction_y)

        call unstagger_field_3d(uvel, uvel_unstag, options%periodic_ew, options%periodic_ns)
        call unstagger_field_3d(vvel, vvel_unstag, options%periodic_ew, options%periodic_ns)
       
        call unstagger_field_2d(btrc, btrc_unstag, options%periodic_ew, options%periodic_ns)

        !Construct fields that contain kinematic boundary condition velocities where they are specified,
        !and NaN where kinematic boundaries should be computed
        do k = 1, upn
            where (GLIDE_IS_DIRICHLET_BOUNDARY(geometry_mask))
                kinematic_bc_u_unstag(k,:,:) = uvel_unstag(k,:,:)
                kinematic_bc_v_unstag(k,:,:) = vvel_unstag(k,:,:)
            elsewhere
                kinematic_bc_u_unstag(k,:,:) = NaN
                kinematic_bc_v_unstag(k,:,:) = NaN
            endwhere
        end do

        !Compute rescaled coordinate parameters (needed because Pattyn uses an
        !irregular Z grid and scales so that 0 is the surface, 1 is the bed)
        call init_rescaled_coordinates(dthckdew,dlsrfdew,dthckdns,dlsrfdns,usrf,thck,lsrf,&
                                               dusrfdew,dusrfdns,d2zdx2,d2zdy2,d2hdx2,d2hdy2,&
                                               sigma,ax,ay,bx,by,cxy,dew,dns,direction_x,direction_y)
       
        !"Spin up" estimate with Pattyn's SIA model runs if we don't already
        !have a good initial guess
        if (.not. valid_initial_guess) then
            call veloc1(dusrfdew, dusrfdns, thck, flwa, sigma, uvel_unstag, vvel_unstag, u, v, ewn, nsn, upn, &
                        FLWN, options%periodic_ew, options%periodic_ns)
            !If we are performing the plastic bed iteration, the SIA is not
            !enough and we need to spin up a better estimate by shoehorning the
            !tau0 values into a linear bed estimate
            if (options%which_ho_bstress == HO_BSTRESS_PLASTIC) then
                call veloc2(efvs, uvel, vvel, flwa, dusrfdew, dusrfdns, thck, ax, ay, &
                        sigma, bx, by, cxy, btrc_unstag, dlsrfdew, dlsrfdns, FLWN, ZIP, VEL2ERR, &
                        TOLER, options, .false., dew, dns,point_mask,totpts,geometry_mask,&
                        kinematic_bc_u_unstag, kinematic_bc_v_unstag, marine_bc_normal)
            end if
        end if
        
        !Higher order velocity estimation
        !I am assuming that efvs (effective viscosity) is the same as mu
        !A NOTE ON COORDINATE TRANSPOSITION:
        !Because of the transposition, ewn=maxy and nsn=maxx.  However, veloc2
        !passes maxy in *first*, so these really get passed in the same order
        !that they normally would.
        call veloc2(efvs, uvel_unstag, vvel_unstag, flwa, dusrfdew, dusrfdns, thck, ax, ay, &
                    sigma, bx, by, cxy, btrc_unstag, dlsrfdew, dlsrfdns, FLWN, ZIP, VEL2ERR, &
                    TOLER, options, .false., dew, dns, &
                    point_mask,totpts,geometry_mask,kinematic_bc_u_unstag, kinematic_bc_v_unstag, marine_bc_normal)

        !Final computation of stress field for output
        !call stressf(mu_t, uvel_t, vvel_t, flwa_t, stagthck_t, ax, ay, dew, dns, sigma, & 
        !             tau_xz_t, tau_yz_t, tau_xx_t, tau_yy_t, tau_xy_t, flwn, zip, periodic_ew, periodic_ns) 
        
        
        call stagvarb_3d_mask(uvel_unstag, uvel,ewn,nsn,upn,geometry_mask)
        call stagvarb_3d_mask(vvel_unstag, vvel,ewn,nsn,upn,geometry_mask)
        
    end subroutine velo_hom_pattyn_nonstag


end module glide_velo_higher
