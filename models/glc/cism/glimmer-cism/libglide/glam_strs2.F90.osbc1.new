!CLEANUP - glam_strs2.F90
! Changed 'real (kind = dp)' to 'real(dp)'.  (Did this in other modules too.)
! Changed glam_velo_fordsiapstr to glam_velo_solver
! Change JFNK to JFNK_velo_solver
!
! "glam_strs2.F90"
!
! 3d velocity calculation based on Blatter/Pattyn, 1st-order equations, by Tony Payne (Univ.
! of Bristol) and Steve Price (Univ. of Bristol / Los Alamos Nat. Lab.). Boundary conditions
! available include periodic (lateral), free surface, zero slip at bed, specified basal 
! traction at bed, and specified basal yield stress at bed (all three of which are implemented
! through various verions of the specified traction b.c.)
! include macros for glide mask definitions
#include "glide_mask.inc"
#include "config.inc"

!GlobalIDs are for distributed TRILINOS variable IDs
#ifdef TRILINOS
#define globalIDs
#endif

!TODO - This module is complex and hard to understand.
!       In particular, there are chunks of code that are used more than once, for Picard as well as JFNK.
!       It would be better to combine these chunks of code into subroutines that can be called
!        from multiple places in the code--or even better, to remove the extra chunks of code
!        if they are no longer needed.

!***********************************************************************
module glam_strs2
!***********************************************************************

use iso_c_binding
use glimmer_paramets, only : dp
use glimmer_physcon,  only : gn, rhoi, rhoo, grav, pi, scyr

!SCALING - What used to be called vis0_glam is now called vis0 and is used by both dycores.
!          Note: if thk0 = 1, then tau0 = rhoi*grav
!          It would not be hard to remove tau0, vis0, and evs0 from this module and elsewhere in the code.
          
use glimmer_paramets, only : thk0, len0, vel0, vis0, tim0, evs0, tau0

use glimmer_log,      only : write_log
use glide_mask
use glimmer_sparse_type
use glimmer_sparse
use glide_types

implicit none

  logical, save :: lateralboundry = .false.
  integer, dimension(6), save :: loc_latbc

  real(dp), allocatable, dimension(:,:,:),     save :: flwafact
  real(dp), allocatable, dimension(:),         save :: dups
  real(dp), allocatable, dimension(:,:,:,:,:), save :: corr
  real(dp), allocatable, dimension(:,:,:,:),   save :: usav
  real(dp), dimension(2),                      save :: usav_avg
  real(dp), allocatable, dimension(:,:,:),     save :: tvel
  real(dp), allocatable, dimension(:),         save :: dup, dupm

  integer, dimension(:,:), allocatable :: uindx

  ! regularization constant for eff. strain rate to avoid infinite visc.
  ! NOTE: would be good to explore how small this really needs to be, as 
  ! code converges much better when this value is made larger.

!SCALING - Is this constant OK?
  real(dp), parameter :: effstrminsq = (1.0d-20 * tim0)**2
  real(dp) :: homotopy = 0.0

  real(dp) :: p1, p2, p3    ! variants of Glen's "n" (e.g. n, (1-n)/n)
  real(dp) :: dew2, dns2, dew4, dns4

  ! combinations of coeffs. used in momentum balance calcs
  real(dp) :: cdxdy
  real(dp), dimension(2) :: cdxdx
  real(dp), dimension(:),   allocatable :: cdsds, cds
  real(dp), dimension(:), allocatable :: cvert, fvert
  real(dp), dimension(:,:), allocatable :: cdsdx

  real(dp), dimension(:), allocatable :: dsigmadew, dsigmadns
  real(dp), dimension(:), allocatable :: d2sigmadew2, d2sigmadns2, d2sigmadewdns
  real(dp) :: d2sigmadewdsigma, d2sigmadnsdsigma

  ! vectors of coeffs. used for switching symmetric solution subroutines between calc.
  ! of x-comp of vel or y-comp of vel
  real(dp), dimension(2), parameter ::   &
           oneorfour = (/ 1.d0, 4.d0 /),     &
           fourorone = (/ 4.d0, 1.d0 /),     &
           oneortwo  = (/ 1.d0, 2.d0 /),     &
           twoorone  = (/ 2.d0, 1.d0 /)

  real(dp), allocatable, dimension(:,:,:), save  :: ughost 
  real(dp), allocatable, dimension(:,:,:), save  :: vghost

  ! coeff. for forward differencing template, used for stress bcs at lateral boundaries
  real(dp), dimension(3), parameter ::   &
           onesideddiff = (/ -3.d0, 4.d0, -1.d0 /)

  ! geometric 2nd and cross-derivs
  real(dp), dimension(:,:), allocatable :: &
              d2thckdew2, d2usrfdew2, d2thckdns2, d2usrfdns2, d2thckdewdns, d2usrfdewdns

  ! variables for plastic-till basal BC iteration using Newton method
  real(dp), dimension(:,:,:), allocatable :: velbcvect, plastic_coeff_lhs, plastic_coeff_rhs, &
                                                     plastic_rhs, plastic_resid
  real(dp), dimension(:,:,:,:), allocatable :: ghostbvel

  ! variables for use in sparse matrix calculation
  real(dp), dimension(:), allocatable :: pcgval, rhsd, rhsx
  integer, dimension(:), allocatable :: pcgcol, pcgrow
  integer, dimension(2) :: pcgsize
  ! additional storage needed for off diagonal blocks when using JFNK for nonlinear iteration 
  real(dp), dimension(:), allocatable :: pcgvaluv, pcgvalvu
  integer, dimension(:), allocatable :: pcgcoluv, pcgrowuv, pcgcolvu, pcgrowvu
  integer :: ct, ct2

!*sfp* NOTE: these redefined here so that they are "in scope" and can avoid being passed as args
  integer :: whatsparse ! needed for putpgcg()
  integer :: nonlinear  ! flag for indicating type of nonlinar iteration (Picard vs. JFNK)

  logical, save :: storeoffdiag = .false. ! true only if using JFNK solver and block, off diag coeffs needed
  logical, save :: calcoffdiag = .false. 
  logical, save :: inisoln = .false.      ! true only if a converged solution (velocity fields) exists

  real(dp) :: linearSolveTime = 0
  real(dp) :: totalLinearSolveTime = 0 ! total linear solve time

  ! AGS: partition information for distributed solves
  ! JEFF: Moved to module-level scope for globalIDs
  integer, allocatable, dimension(:) :: myIndices
  real(dp), allocatable, dimension(:) :: myX, myY, myZ
  integer, allocatable, dimension(:,:,:) :: loc2_array
  integer :: mySize = -1

  ! JEFF: Debugging Output Variables
  integer :: overallloop = 1

!***********************************************************************

contains

!***********************************************************************

!TODO - Is this subroutine still needed?

subroutine dumpvels(name, uvel, vvel)
    !JEFF routine to track the uvel and vvel calculations in Picard Iteration for debugging
    !3/28/11
    use parallel
    implicit none

    character(*) :: name
    real(dp), dimension(:,:,:), intent(inout) :: uvel, vvel  ! horiz vel components: u(z), v(z)

    if (distributed_execution()) then
       if (this_rank == 0) then
           write(*,*) name, "Proc 0 uvel & vvel (1,7:8,16:17)", uvel(1,7:8,16:17), vvel(1,7:8,16:17)
       else
           write(*,*) name, "Proc 1 uvel & vvel (1,7:8,0:1)", uvel(1,7:8,0:1), vvel(1,7:8,0:1)
       endif
    else
       write(*,*) name, "Parallel uvel & vvel (1,5:6,15:16)", uvel(1,5:6,15:16), vvel(1,5:6,15:16)
    endif 
end subroutine dumpvels


subroutine glam_velo_init( ewn,   nsn,   upn,    &
                           dew,   dns,           &
                           sigma)

    ! Allocate arrays and initialize variables.
    implicit none

    integer, intent(in) :: ewn, nsn, upn
    real(dp), intent(in) :: dew, dns

    real(dp), dimension(:), intent(in)  :: sigma

    integer :: up

    allocate( dup(upn) )
    allocate( dupm(upn) )
    allocate( cvert(upn) )
    allocate( cdsdx(upn,2) )
    allocate( cdsds(upn) )
    allocate( cds(upn) )
    allocate( fvert(upn) )
    allocate(ughost(2,ewn-1,nsn-1))
    allocate(vghost(2,ewn-1,nsn-1))

    ! NOTE: "dup", the sigma coordinate spacing is defined as a vector to allow it to 
    ! be read in from file for use with non-constant vertical grid spacing. Currently, this
    ! is not working, so the code will not give accurate results if the sigma coordinate is
    ! not regularly spaced. 
    dup = (/ ( (sigma(2)-sigma(1)), up = 1, upn) /)
    dupm = - 0.25d0 / dup

!whl - Moved stagsigma calculation to glide_setup module
!!    stagsigma(1:upn-1) = (sigma(1:upn-1) + sigma(2:upn)) / 2.d0

    ! p1 = -1/n   - used with rate factor in eff. visc. def.
    ! p2 = (1-n)/2n   - used with eff. strain rate in eff. visc. def. 
    ! p3 = (1-n)/n   !TODO - Remove p3?  It is never used.

    p1 = -1.d0 / real(gn,dp)
    p2 = (1.d0 - real(gn,dp)) / (2.d0 * real(gn,dp))
    p3 = (1.d0 - real(gn,dp)) / real(gn,dp)

    dew2 = 2.d0 * dew; dns2 = 2.d0 * dns        ! 2x the standard grid spacing
    dew4 = 4.d0 * dew; dns4 = 4.d0 * dns        ! 4x the standard grid spacing

    allocate(dsigmadew(upn),  dsigmadns(upn))
    allocate(d2sigmadew2(upn),d2sigmadns2(upn),d2sigmadewdns(upn))

    allocate (d2thckdew2(ewn-1,nsn-1),d2thckdns2(ewn-1,nsn-1),d2thckdewdns(ewn-1,nsn-1), &
              d2usrfdew2(ewn-1,nsn-1),d2usrfdns2(ewn-1,nsn-1),d2usrfdewdns(ewn-1,nsn-1))

    allocate(flwafact(1:upn-1,ewn,nsn))  ! NOTE: the vert dim here must agree w/ that of 'efvs'

    allocate(dups(upn))

  ! allocate/initialize variables for plastic-till basal BC iteration using Newton method
    allocate(velbcvect(2,ewn-1,nsn-1),plastic_coeff_rhs(2,ewn-1,nsn-1),plastic_coeff_lhs(2,ewn-1,nsn-1), &
            plastic_rhs(2,ewn-1,nsn-1), plastic_resid(1,ewn-1,nsn-1) )
    allocate(ghostbvel(2,3,ewn-1,nsn-1))        !! for saving the fictious basal vels at the bed !!

    plastic_coeff_rhs(:,:,:) = 0.d0
    plastic_coeff_lhs(:,:,:) = 0.d0
    plastic_rhs(:,:,:) = 0.d0
    plastic_resid(:,:,:) = 0.d0
    ghostbvel(:,:,:,:) = 0.d0
    velbcvect(:,:,:) = 0.d0

    flwafact = 0.d0

     ! define constants used in various FD calculations associated with the 
     ! subroutine 'findcoefst'   
     call calccoeffsinit(upn, dew, dns)

    dups = (/ (sigma(up+1) - sigma(up), up=1,upn-1), 0.d0 /)

end subroutine glam_velo_init


!***********************************************************************

! Note that this is the driver subroutine, called from 'run_ho_diagnostic' in
! 'glide_velo_higher.F90'. In turn, 'run_ho_model' is called from 'inc_remap_driver' in
! 'glam.F90', and 'inc_remap_driver' is called from 'glide_tstep_ps' in 'glide.F90'.

subroutine glam_velo_solver(ewn,      nsn,    upn,  &
                            dew,      dns,          &
                            sigma,    stagsigma,    &
                            thck,     usrf,         &
                            lsrf,     topg,         &
                            dthckdew, dthckdns,     &
                            dusrfdew, dusrfdns,     &
                            dlsrfdew, dlsrfdns,     &
                            stagthck, flwa,         &
                            mintauf,                &
                            btraction,              &
                            umask,                  &
                            whichbabc,              &
                            whichefvs,              &
                            whichresid,             &
                            whichnonlinear,         &
                            whichsparse,            &
                            periodic_ew,periodic_ns,&
                            beta,                   &
                            uvel,     vvel,         &
                            uflx,     vflx,         &
                            efvs )

  use parallel

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:),   intent(inout)  :: umask

!TODO - Make umask intent in?
  ! NOTE: 'inout' status to 'umask' should be changed to 'in' at some point, 
  ! but for now this allows for some minor internal hacks to CISM-defined mask  

  real(dp), intent(in) :: dew, dns

  real(dp), dimension(:),     intent(in)  :: sigma, stagsigma       ! sigma coords
  real(dp), dimension(:,:),   intent(in)  :: thck, usrf, lsrf, topg ! geom vars
  real(dp), dimension(:,:),   intent(in)  :: dthckdew, dthckdns     ! thick grads
  real(dp), dimension(:,:),   intent(in)  :: dusrfdew, dusrfdns     ! upper surf grads
  real(dp), dimension(:,:),   intent(in)  :: dlsrfdew, dlsrfdns     ! basal surf grads
  real(dp), dimension(:,:),   intent(in)  :: stagthck               ! staggered thickness
  real(dp), dimension(:,:),   intent(in)  :: minTauf                ! till yield stress
  real(dp), dimension(:,:,:), intent(inout) :: btraction            ! consistent basal traction array
  real(dp), dimension(:,:,:), intent(in)  :: flwa                   ! flow law rate factor

  ! This is the betasquared field from CISM (externally specified), and should eventually
  ! take the place of the subroutine 'calcbetasquared' below. For now, there is simply an option
  ! in the subroutine 'calcbetasquared' (case 9) to use this external, CISM specified value for
  ! the betasquared field as opposed to one of the values calculated internally.
  real(dp), dimension(:,:),   intent(in)  :: beta

  integer, intent(in) :: whichbabc    ! options for betasquared field to use
  integer, intent(in) :: whichefvs    ! options for efvs calculation (calculate it or make it uniform)
  integer, intent(in) :: whichresid   ! options for method to use when calculating vel residul
  integer, intent(in) :: whichnonlinear  ! options for which method for doing elliptic solve
  integer, intent(in) :: whichsparse  ! options for which method for doing elliptic solve
  logical, intent(in) :: periodic_ew, periodic_ns  ! options for applying periodic bcs or not

  real(dp), dimension(:,:,:), intent(inout) :: uvel, vvel  ! horiz vel components: u(z), v(z)
  real(dp), dimension(:,:),   intent(out) :: uflx, vflx  ! horiz fluxs: u_bar*H, v_bar*H
  real(dp), dimension(:,:,:), intent(out) :: efvs        ! effective viscosity

  integer :: ew, ns, up     ! counters for horiz and vert do loops

  real(dp), parameter :: minres = 1.0d-4    ! assume vel fields converged below this resid 
  real(dp), parameter :: NL_tol = 1.0d-06   ! to have same criterion than with JFNK
  real(dp), save, dimension(2) :: resid     ! vector for storing u resid and v resid 
  real(dp) :: plastic_resid_norm = 0.d0    ! norm of residual used in Newton-based plastic bed iteration

  integer, parameter :: cmax = 300                  ! max no. of iterations
  integer :: counter, linit                         ! iteration counter, ???
  character(len=100) :: message                     ! error message

  ! variables used for incorporating generic wrapper to sparse solver
  type(sparse_matrix_type) :: matrix
  real(dp), dimension(:), allocatable :: answer, uk_1, vk_1, F
  real(dp) :: err, L2norm, L2square, NL_target
  integer :: iter, pic
  integer , dimension(:), allocatable :: g_flag ! jfl flag for ghost cells

  ! variables for when to stop outer loop when using Picard for nonlinear iteration 
  real(dp) :: outer_it_criterion, outer_it_target

  ! variables for debugging output JEFF
  character(3) :: loopnum
  character(3) :: looptime
  real(dp) :: multiplier

 call t_startf("PICARD_pre")
  ! RN_20100125: assigning value for whatsparse, which is needed for putpcgc()
  whatsparse = whichsparse

  ! assign value for nonlinear iteration flag
  nonlinear = whichnonlinear

!TODO - Note: d2usrfdew2 and d2usrfdns2 are needed at all locally owned velocity points.
!       I am not sure where and why the upwind 2nd derivatives are computed.

  ! calc geometric 2nd deriv. for generic input variable 'ipvr', returns 'opvr'
  call geom2ders(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdew2, d2usrfdns2)
  call geom2ders(ewn, nsn, dew, dns, thck, stagthck, d2thckdew2, d2thckdns2)

  ! calc geometric 2nd cross-deriv. for generic input variable 'ipvr', returns 'opvr'
  call geom2derscros(dew, dns, thck, stagthck, d2thckdewdns)
  call geom2derscros(dew, dns, usrf, stagthck, d2usrfdewdns)

  allocate(uindx(ewn-1,nsn-1))

  ! If a point from the 2d array 'mask' is associated with a non-zero ice thickness
  ! assign it a unique number. If not assign a zero.             
  uindx = indxvelostr(ewn, nsn, upn, umask,pcgsize(1))

!!!!!!!!!! Boundary conditions HACKS section !!!!!!!!!!!!!

!TODO - Replace these hacks with something more robust.

!! A hack of the boundary condition mask needed for the Ross Ice Shelf exp.
!! The quick check of whether or not this is the Ross experiment is to look
!! at the domain size.
 if( ewn == 151 .and. nsn == 115 )then
    call not_parallel(__FILE__, __LINE__)
    do ns=1,nsn-1; do ew=1,ewn-1
        if( umask(ew,ns) == 21 .or. umask(ew,ns) == 5 )then
            umask(ew,ns) = 73
        endif
    end do; end do
 end if

!! hack for basal processes submodel test case, to avoid floatation at downstream
!! end yet still allow for application of a floating ice bc there
!  do ns=1,nsn-1; do ew=1,ewn-1
!      if( umask(ew,ns) == 37 )then
!          umask(ew,ns) = 41
!      endif
!  end do; end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! allocate space for storing temporary across-flow comp of velocity
  allocate(tvel(upn,ewn-1,nsn-1))
  tvel = 0.d0

  ! allocate space for variables used by 'mindcrash' function (unstable manifold correction)
  allocate(corr(upn,ewn-1,nsn-1,2,2),usav(upn,ewn-1,nsn-1,2))
  ! and initialize them
  corr = 0.d0
  usav = 0.d0

  ! make an initial guess at the size of the sparse matrix
  pcgsize(2) = pcgsize(1) * 20

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

#ifdef TRILINOS
  if (whatsparse == STANDALONE_TRILINOS_SOLVER) then
#ifdef globalIDs
     if (main_task) write(*,*) "Using GlobalIDs..."
	 ! JEFF: Define myIndices in terms of globalIDs
     allocate(myIndices(pcgsize(1)))  ! myIndices is an integer vector with a unique ID for each layer for ice grid points
     allocate(myX(pcgsize(1))) ! Coordinates of nodes, used by ML preconditioner
     allocate(myY(pcgsize(1))) 
     allocate(myZ(pcgsize(1))) 
     call distributed_create_partition(ewn, nsn, (upn + 2) , uindx, pcgsize(1), myIndices, myX, myY, myZ)  ! Uses uindx mask to determine ice grid points.
     mySize = pcgsize(1)  ! Set variable for inittrilinos

     !write(*,*) "GlobalIDs myIndices..."
     !write(*,*) "pcgsize = ", pcgsize(1)
     !write(*,*) "myIndices = ", myIndices
     !call parallel_stop(__FILE__, __LINE__)
#else
     ! AGS: Get partition -- later this will be known by distributed glimmer
     call dopartition(pcgsize(1), mySize) 
     allocate(myIndices(mySize))
     allocate(myX(1)) ! Coordinates (only set for globalIDs)
     allocate(myY(1)) 
     allocate(myZ(1)) 
     call getpartition(mySize, myIndices) 

     if (distributed_execution()) then
         if (main_task) write(*,*) "Distributed Version cannot be run without globalIDs.  Stopping."
         call not_parallel(__FILE__, __LINE__)  ! Fatal if running without GlobalIDs in MPI
     endif

     !write(*,*) "Trilinos Generated Partition Map myIndices..."
     !write(*,*) "pcgsize = ", pcgsize(1), " mySize = ", mySize
     !write(*,*) "myIndices = ", myIndices
     !call parallel_stop(__FILE__, __LINE__)
#endif

     ! Now send this partition to Trilinos initialization routines
     call inittrilinos(20, mySize, myIndices, myX, myY, myZ) 

     ! Set if need full solution vector returned or just owned portion
#ifdef globalIDs
     ! We default Trilinos to return full solution vector.
     ! Per AGS, for both distributed and serial globalIDs cases, we must return just owned portion in order to prevent a permutation of the results.
     ! This was built into parallel_set_trilinos_return_vect that handled the difference between _single and _mpi versions.
     ! AGS: no longer needed, now the default
     ! call returnownedvector()
!     call parallel_set_trilinos_return_vect
#endif

     !No Triad matrix needed in this case -- save on memory alloc
     pcgsize(2) = 1

     ! JEFF: deallocate myIndices after the solve loop, because used in translation between globalIDs and local indices
     ! deallocate(myIndices)
  endif
#else
  if (whatsparse == STANDALONE_TRILINOS_SOLVER) then
         write(*,*) 'Error: Trilinos sparse solver requires Trilinos build'
         stop
  endif
#endif

!==============================================================================
! RN_20100126: End of the block
!==============================================================================

  ! allocate sparse matrix variables
  allocate (pcgrow(pcgsize(2)),pcgcol(pcgsize(2)),rhsd(pcgsize(1)), &
            pcgval(pcgsize(2)))

  allocate(matrix%row(pcgsize(2)), matrix%col(pcgsize(2)), &
            matrix%val(pcgsize(2)), answer(pcgsize(1)))

  allocate( uk_1(pcgsize(1)), vk_1(pcgsize(1)), &
            F(2*pcgsize(1)), g_flag(pcgsize(1)) ) ! jfl for res calc.

  ! set residual and iteration counter to initial values
  resid = 1.d0
  counter = 1
  L2norm = 1.0d20

  ! intialize outer loop test vars
  outer_it_criterion = 1.0
  outer_it_target = 0.0

  if (main_task) then
     ! print some info to the screen to update on iteration progress
     print *, ' '
     print *, 'Running Payne/Price higher-order dynamics solver'
     print *, ' '
     if( whichresid == 3 )then
       print *, 'iter #     resid (L2 norm)       target resid'
     else
       print *, 'iter #     uvel resid         vvel resid       target resid'
     end if
     print *, ' '
  endif

 call t_stopf("PICARD_pre")
  ! ****************************************************************************************
  ! START of Picard iteration
  ! ****************************************************************************************
 call t_startf("PICARD_iter")

  call ghost_preprocess( ewn, nsn, upn, uindx, ughost, vghost, &
                         uk_1, vk_1, uvel, vvel, g_flag) ! jfl_20100430

  ! Picard iteration; continue iterating until resid falls below specified tolerance
  ! or the max no. of iterations is exceeded

  !JEFF Guarantees at least one loop
  outer_it_criterion = 1.0
  outer_it_target = 0.0

  do while ( outer_it_criterion >= outer_it_target .and. counter < cmax)    ! use L2 norm for resid calculation
 call t_startf("PICARD_in_iter")

  ! choose outer loop stopping criterion
  if( counter > 1 )then
    if( whichresid == 3 )then
      outer_it_criterion = L2norm
      outer_it_target = NL_target
    else
      outer_it_criterion = maxval(resid)
      outer_it_target = minres   
    end if
  else
    outer_it_criterion = 1.0d10
    outer_it_target = 1.0d-12
  end if

#ifdef JEFFTEST
    !JEFF Debugging Output to see what differences in final vvel and tvel.
    write(loopnum,'(i3.3)') counter
    write(Looptime, '(i3.3)') overallloop
    loopnum = trim(loopnum)  ! Trying to get rid of spaces in name.
    Looptime = trim(Looptime)
    call distributed_print("uvela_ov"//Looptime//"_pic"//loopnum//"_tsk", uvel)

    call distributed_print("vvela_ov"//Looptime//"_pic"//loopnum//"_tsk", vvel)

    ! call dumpvels("Before findefvsstr", uvel, vvel)

    ! call distributed_print("preefvs_ov"//Looptime//"_pic"//loopnum//"_tsk", efvs)
#endif

!HALO - To avoid parallel halo calls for efvs within glam_strs2, we need to compute efvs in one layer of halo cells

    ! calc effective viscosity using previously calc vel. field
    call findefvsstr(ewn,  nsn,  upn,      &
                     stagsigma,  counter,  &
                     whichefvs,  efvs,     &
                     uvel,       vvel,     &
                     flwa,       thck,     &
                     dusrfdew,   dthckdew, &
                     dusrfdns,   dthckdns, &
                     umask)

    ! calculate coeff. for stress balance in y-direction 
    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     2,           efvs,           &
                     vvel,        uvel,           &
                     thck,        dusrfdns,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta,        btraction,      &
                     counter, 0 )

    ! put vels and coeffs from 3d arrays into sparse vector format
    call solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, vvel )

!==============================================================================
! jfl 20100412: residual for v comp: Fv= A(u^k-1,v^k-1)v^k-1 - b(u^k-1,v^k-1)  
!==============================================================================

!TODO - Is L2square summed correctly in res_vect?
    !JEFF - The multiplication Ax is done across all nodes, but Ax - b is only 
    !       computed locally, so L2square needs to be summed.
 call t_startf("PICARD_res_vect")
    call res_vect( matrix, vk_1, rhsd, size(rhsd), g_flag, L2square, whichsparse ) 
 call t_stopf("PICARD_res_vect")

    L2norm  = L2square
    F(1:pcgsize(1)) = vk_1(:)
      
!   call output_res(ewn,nsn,upn,uindx,counter,size(vk_1),vk_1, 2) ! JFL

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

 call t_startf("PICARD_solvea")
  if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
     call sparse_easy_solve(matrix, rhsd, answer, err, iter, whichsparse)
#ifdef TRILINOS
  else
     call solvewithtrilinos(rhsd, answer, linearSolveTime) 
     totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
     ! write(*,*) 'Total linear solve time so far', totalLinearSolveTime
#endif
  endif
 call t_stopf("PICARD_solvea")

!==============================================================================
! RN_20100129: End of the block
!==============================================================================

    vk_1 = answer ! jfl for residual calculation

    ! put vels and coeffs from sparse vector format (soln) back into 3d arrays
    call solver_postprocess( ewn, nsn, upn, 2, uindx, answer, tvel, ghostbvel )

    ! NOTE: y-component of velocity that comes out is called "tvel", to differentiate it
    ! from the y-vel solution from the previous iteration, which is maintained as "vvel". 
    ! This is necessary since we have not yet solved for the x-comp of vel, which needs the
    ! old prev. guess as an input (NOT the new guess).

!TODO - Can we eliminate the periodic option here?  Periodic BC should be handled automatically.

! implement periodic boundary conditions in y (if flagged)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( periodic_ns )then
        call not_parallel(__FILE__, __LINE__)

        tvel(:,:,nsn-1) = tvel(:,:,2)
        tvel(:,:,1) = tvel(:,:,nsn-2)
    end if
    if( periodic_ew )then
        call not_parallel(__FILE__, __LINE__)

        tvel(:,ewn-1,:) = tvel(:,2,:)
        tvel(:,1,:) = tvel(:,ewn-2,:)
    end if
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! calculate coeff. for stress balance calc. in x-direction 
    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     1,           efvs,           &
                     uvel,        vvel,           &
                     thck,        dusrfdew,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta,        btraction,      &
                     counter, 0 )

    ! put vels and coeffs from 3d arrays into sparse vector format
    call solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, uvel )

!==============================================================================
! jfl 20100412: residual for u comp: Fu= C(u^k-1,v^k-1)u^k-1 - d(u^k-1,v^k-1)  
!==============================================================================

 call t_startf("PICARD_res_vect")
    call res_vect( matrix, uk_1, rhsd, size(rhsd), g_flag, L2square, whichsparse ) 
 call t_stopf("PICARD_res_vect")

    L2norm = sqrt(L2norm + L2square)
    F(pcgsize(1)+1:2*pcgsize(1)) = uk_1(:) ! F = [ Fv, Fu ]

!    print *, 'L2 with/without ghost (k)= ', counter, &
!              sqrt(DOT_PRODUCT(F,F)), L2norm
!    if (counter <= 2) NL_target = NL_tol * L2norm
!    if (counter == 1) NL_target = NL_tol * L2norm
    if (counter == 1) NL_target = 1.0d-4 

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

 call t_startf("PICARD_solveb")
  if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
     call sparse_easy_solve(matrix, rhsd, answer, err, iter, whichsparse)
#ifdef TRILINOS
  else
     call solvewithtrilinos(rhsd, answer, linearSolveTime) 
     totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
     ! write(*,*) 'Total linear solve time so far', totalLinearSolveTime
#endif
  endif
 call t_stopf("PICARD_solveb")

!==============================================================================
! RN_20100129: End of the block
!==============================================================================

    uk_1 = answer ! jfl for residual calculation

    ! put vels and coeffs from sparse vector format (soln) back into 3d arrays
    call solver_postprocess( ewn, nsn, upn, 1, uindx, answer, uvel, ghostbvel )

    ! call fraction of assembly routines, passing current vel estimates (w/o manifold
    ! correction!) to calculate consistent basal tractions

    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     2,           efvs,           &
                     tvel,        uvel,           &
                     thck,        dusrfdns,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta,        btraction,      &
                     counter, 1 )

   call findcoefstr(ewn,  nsn,   upn,             &
                     dew,  dns,   sigma,          &
                     1,           efvs,           &
                     uvel,        tvel,           &
                     thck,        dusrfdew,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta,        btraction,      &
                     counter, 1 )

    !JEFF Commented out plasticbediteration() per Steve Price. December 2010
    ! call plasticbediteration( ewn, nsn, uvel(upn,:,:), tvel(upn,:,:), btraction, minTauf, &
    !                          plastic_coeff_lhs, plastic_coeff_rhs, plastic_rhs, plastic_resid )

    ! apply unstable manifold correction to converged velocities

 call t_startf("PICARD_mindcrsh")
    call mindcrshstr(1,whichresid,uvel,counter,resid(1))
    vvel = tvel
    call mindcrshstr(2,whichresid,vvel,counter,resid(2))
 call t_stopf("PICARD_mindcrsh")

!TODO - I'm pretty sure these updates *are* needed.
!       
    ! coordinate halos for updated uvel and vvel
    call staggered_parallel_halo(uvel)
    call staggered_parallel_halo(vvel)

    !call dumpvels("After mindcrsh", uvel, vvel)

!TODO - Remove periodic BC stuff?

! implement periodic boundary conditions in x (if flagged)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( periodic_ns )then
        call not_parallel(__FILE__, __LINE__)

        uvel(:,:,nsn-1) = uvel(:,:,2)
        uvel(:,:,1) = uvel(:,:,nsn-2)
    end if
    if( periodic_ew )then
        call not_parallel(__FILE__, __LINE__)

        uvel(:,ewn-1,:) = uvel(:,2,:)
        uvel(:,1,:) = uvel(:,ewn-2,:)
    end if
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (this_rank == 0) then
        ! Can't use main_task flag because main_task is true for all processors in case of parallel_single
        ! output the iteration status: iteration number, max residual, and location of max residual
        ! (send output to the screen or to the log file, per whichever line is commented out) 
        if( whichresid == 3 )then
            print '(i4,3g20.6)', counter, L2norm, NL_target    ! Output when using L2norm for convergence
            print '(a,i4,3g20.6)', "sup-norm uvel, vvel=", counter, resid(1), resid(2), minres
            !write(message,'(i4,3g20.6)') counter, L2norm, NL_target
            !call write_log (message)
        else
            print '(i4,3g20.6)', counter, resid(1), resid(2), minres
            !write(message,'(" * strs ",i3,3g20.6)') counter, resid(1), resid(2), minres
            !call write_log (message)
        end if
    endif

    counter = counter + 1   ! advance the iteration counter
 call t_stopf("PICARD_in_iter")

  end do  ! while ( outer_it_criterion >= outer_it_target .and. counter < cmax)

  inisoln = .true.

  ! ****************************************************************************************
  ! END of Picard iteration
  ! ****************************************************************************************
 call t_stopf("PICARD_iter")

 call t_startf("PICARD_post")
  call ghost_postprocess( ewn, nsn, upn, uindx, uk_1, vk_1, &
                          ughost, vghost )

!TODO - If needed, these loops should be over locally owned velocity points: (ilo-1:ihi, jlo-1:jhi).
!       However, I don't think uflx and vflx are needed; they are not used by remapping subroutine.

  do ns = 1+staggered_shalo,size(umask,2)-staggered_nhalo
      do ew = 1+staggered_whalo,size(umask,1)-staggered_ehalo
      ! calc. fluxes from converged vel. fields (needed for input to thickness evolution subroutine)
         if (umask(ew,ns) > 0) then
             uflx(ew,ns) = vertintg(upn, sigma, uvel(:,ew,ns)) * stagthck(ew,ns)
             vflx(ew,ns) = vertintg(upn, sigma, vvel(:,ew,ns)) * stagthck(ew,ns)
         end if
      end do
  end do

  !JEFF: Coordinate halos
  !JEFF: umask is marked as INOUT and is updated for the Ross Ice Shelf experiment, but for no other, so don't update halos
  !JEFF: uvel, vvel, uflx, and vflx are calculated in this routine, but only for "owned" grid cells, so update halos to get neighboring values.

  !call staggered_parallel_halo(uvel) (called earlier)
  !call staggered_parallel_halo(vvel) (called earlier)

!HALO - Do we need halo updates for btraction and efvs?
!       I think we don't need an update for efvs, because it is already computed in a layer of halo cells.
!       And I think we don't need an update for btraction, because it is computed in bodyset for all
!        locally owned velocity points.

  call parallel_halo(efvs)
  call parallel_halo(btraction)

!HALO - Pretty sure we don't need these updates; uflx and vflx are not used.
  call staggered_parallel_halo(uflx)
  call staggered_parallel_halo(vflx)

#ifdef JEFFTEST    
  !JEFF Debugging Output to see what differences in final vvel and tvel.
    write(CurrTimeLoopStr, '(i3.3)') CurrTimeLoop
    call distributed_print("uvel_post_ov"//CurrTimeLoopStr//"_tsk", uvel)

    call distributed_print("vvel_post_ov"//CurrTimeLoopStr//"_tsk", vvel)
#endif
  ! JEFF: Deallocate myIndices which is used to intialize Trilinos
  if (whatsparse == STANDALONE_TRILINOS_SOLVER) then
     deallocate(myIndices)
     deallocate(myX)
     deallocate(myY)
     deallocate(myZ)
  endif

  ! de-allocation sparse matrix solution variables 
  deallocate(tvel)
  deallocate(uindx,corr,usav)
  deallocate(pcgval,pcgrow,pcgcol,rhsd)
  deallocate(matrix%row, matrix%col, matrix%val)
  deallocate(answer)
  deallocate(uk_1, vk_1, F, g_flag) ! jfl 

  !JEFF Debugging output
  overallloop = overallloop + 1
 call t_stopf("PICARD_post")

  return

end subroutine glam_velo_solver

!***********************************************************************

!TODO - Can we pass arguments explicitly instead of passing 'model'?

subroutine JFNK_velo_solver  (model,umask)

  use parallel

  use iso_c_binding 
  use glide_types, only : glide_global_type, pass_through

  implicit none

  type(glide_global_type) ,intent(inout) :: model

!TODO - Can we make the mask intent in?

  integer, dimension(:,:),   intent(inout)  :: umask  !*sfp* replaces the prev., internally calc. mask
                                                      ! ... 'inout' status allows for a minor alteration
                                                      ! to cism defined mask, which don't necessarily 
                                                      ! associate all/any boundaries as a unique mask value.

! new glide_global_type variables for everything needed to pass thru trilinos NOX to calc_F as a pointer.
  type(pass_through) ,target  :: resid_object
  type(pass_through) ,pointer :: fptr=>NULL()
  type(c_ptr)                 :: c_ptr_to_object

!KJE for NOX
  integer(c_int) :: xk_size
  real(dp), dimension(:), allocatable :: xk_1, xk_1_plus
  real(dp), dimension(:), allocatable :: vectx
  integer ,dimension(:) ,allocatable :: gx_flag, g_flag

! split off of derived types

!TODO - Should the following be passed in explicitly?

! intent(in)
  integer :: ewn, nsn, upn
  real(dp) :: dew, dns

  real(dp), dimension(:)     ,pointer :: sigma, stagsigma
  real(dp), dimension(:,:)   ,pointer :: thck, usrf, lsrf, topg
  real(dp), dimension(:,:)   ,pointer :: dthckdew, dthckdns
  real(dp), dimension(:,:)   ,pointer :: dusrfdew, dusrfdns
  real(dp), dimension(:,:)   ,pointer :: dlsrfdew, dlsrfdns
  real(dp), dimension(:,:)   ,pointer :: stagthck
  real(dp), dimension(:,:,:) ,pointer :: flwa
  real(dp), dimension(:,:)   ,pointer :: minTauf
  real(dp), dimension(:,:,:) ,pointer :: btraction            ! consistent basal traction array
  
!TODO - Anything to update here?
  !*sfp* This is the betasquared field from CISM (externally specified), and should eventually
  ! take the place of the subroutine 'calcbetasquared' below (for now, using this value instead
  ! will simply be included as another option within that subroutine) 
  real(dp), dimension(:,:)  ,pointer :: beta 

  integer :: whichbabc
  integer :: whichefvs
  integer :: whichresid
  integer :: whichsparse
  integer :: whichnonlinear
  logical :: periodic_ew, periodic_ns

!TODO - Should the following be passed out explicitly?
! intent(out)
  real(dp), dimension(:,:,:) ,pointer :: uvel, vvel
  real(dp), dimension(:,:)   ,pointer :: uflx, vflx
  real(dp), dimension(:,:,:) ,pointer :: efvs

  integer :: ew, ns, up, nele, k

  real(dp), parameter :: NL_tol = 1.0d-06

  integer, parameter :: img = 20, img1 = img+1
  integer :: kmax = 1000

  character(len=100) :: message

!*sfp* needed to incorporate generic wrapper to solver
  type(sparse_matrix_type) :: matrixA, matrixC, matrixtp, matrixAuv, matrixAvu
  real(dp), dimension(:), allocatable :: answer, uk_1, vk_1
  real(dp), dimension(:), allocatable :: vectp, uk_1_plus, vk_1_plus
  real(dp), dimension(:), allocatable :: dx, F, F_plus
  real(dp), dimension(:), allocatable :: wk1, wk2, rhs
  real(dp), dimension(:,:), allocatable :: vv, wk
  real(dp) :: L2norm, L2norm_wig, tol, gamma_l, epsilon,NL_target
  real(dp) :: crap
  integer :: tot_its, itenb, maxiteGMRES, iout, icode

!  interface
!    subroutine noxsolve(vectorSize,vector,v_container) bind(C,name='noxsolve')
!      use iso_c_binding
!          integer(c_int)                :: vectorSize
!          real(c_double)  ,dimension(*) :: vector
!          type(c_ptr)                   :: v_container
!      end subroutine noxsolve
!
!    subroutine noxinit(vectorSize,vector,comm,v_container) bind(C,name='noxinit')
!      use iso_c_binding
!          integer(c_int)                :: vectorSize,comm
!          real(c_double)  ,dimension(*) :: vector
!          type(c_ptr)                   :: v_container
!      end subroutine noxinit
!
!      subroutine noxfinish() bind(C,name='noxfinish')
!       use iso_c_binding
!      end subroutine noxfinish
!  end interface

!TODO - Could eliminate pointers if arguments are passed in explicitly.
 call t_startf("JFNK_pre")
  ewn = model%general%ewn
  nsn = model%general%nsn
  upn = model%general%upn
  dew = model%numerics%dew
  dns = model%numerics%dns
  sigma => model%numerics%sigma(:)
  stagsigma => model%numerics%stagsigma(:)
  thck => model%geometry%thck(:,:)
  usrf => model%geometry%usrf(:,:)
  lsrf => model%geometry%lsrf(:,:)
  topg => model%geometry%topg(:,:)
  dthckdew => model%geomderv%dthckdew(:,:)
  dthckdns => model%geomderv%dthckdns(:,:)
  dusrfdew => model%geomderv%dusrfdew(:,:)
  dusrfdns => model%geomderv%dusrfdns(:,:)
  dlsrfdew => model%geomderv%dlsrfdew(:,:)
  dlsrfdns => model%geomderv%dlsrfdns(:,:)
  stagthck => model%geomderv%stagthck(:,:)
  flwa => model%temper%flwa(:,:,:)
  mintauf => model%basalproc%minTauf(:,:)
  btraction => model%velocity%btraction(:,:,:)
  whichbabc = model%options%which_ho_babc
  whichefvs = model%options%which_ho_efvs
  whichresid = model%options%which_ho_resid
  whichsparse = model%options%which_ho_sparse
  whichnonlinear = model%options%which_ho_nonlinear
  periodic_ew = model%options%periodic_ew
  periodic_ns = model%options%periodic_ns
  beta => model%velocity%beta(:,:)

  uvel => model%velocity%uvel(:,:,:)
  vvel => model%velocity%vvel(:,:,:)
  uflx => model%velocity%uflx(:,:)
  vflx => model%velocity%vflx(:,:)
  efvs => model%stress%efvs(:,:,:)

  ! RN_20100125: assigning value for whatsparse, which is needed for putpcgc()
!TODO - Can we use just one variable for each of these options?
  whatsparse = whichsparse
  nonlinear = whichnonlinear

!TODO - Much of the following code is a copy of code above.  
!       Can we get by with a single copy?  I'm thinking of operations that are done once, before the iterations begin.

  ! *sfp** geometric 1st deriv. for generic input variable 'ipvr',
  !      output as 'opvr' (includes 'upwinding' for boundary values)
  call geom2ders(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdew2, d2usrfdns2)
  call geom2ders(ewn, nsn, dew, dns, thck, stagthck, d2thckdew2, d2thckdns2)

  ! *sfp** geometric (2nd) cross-deriv. for generic input variable 'ipvr', output as 'opvr'
  call geom2derscros(dew, dns, thck, stagthck, d2thckdewdns)
  call geom2derscros(dew, dns, usrf, stagthck, d2usrfdewdns)

!TODO - Do these derivatives have to go in the model derived type and the residual object?
  ! put d2's into model derived type structure to eventually go into resid_object
  model%geomderv%d2thckdew2 = d2thckdew2
  model%geomderv%d2thckdns2 = d2thckdns2
  model%geomderv%d2usrfdew2 = d2usrfdew2
  model%geomderv%d2usrfdns2 = d2usrfdns2

  ! *sfp** make a 2d array identifying if the associated point has zero thickness,
  !      has non-zero thickness and is interior, or has non-zero thickness
  !      and is along a boundary

  !*sfp* This subroutine has been altered from its original form (was a function, still included
  ! below w/ subroutine but commented out) to allow for a tweak to the CISM calculated mask (adds
  ! in an unique number for ANY arbitrary boundary, be it land, water, or simply at the edge of
  ! the calculation domain). 

  allocate(uindx(ewn-1,nsn-1))

  ! *sfp** if a point from the 2d array 'mask' is associated with non-zero ice thickness,
  !      either a boundary or interior point, give it a unique number. If not, give it a zero			 
  uindx = indxvelostr(ewn, nsn, upn, umask, pcgsize(1))

  L2norm = 1.0d20
 
  ! *sfp** an initial guess at the size of the sparse matrix
  pcgsize(2) = pcgsize(1) * 20

  ! Structure to become NOX implementation for JFNK solve
  xk_size=2*pcgsize(1)

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

#ifdef TRILINOS
  if (whatsparse == STANDALONE_TRILINOS_SOLVER) then
#ifdef globalIDs
     if (main_task) write(*,*) "Using GlobalIDs..."
	 ! JEFF: Define myIndices in terms of globalIDs
     allocate(myIndices(pcgsize(1)))  ! myIndices is an integer vector with a unique ID for each layer for ice grid points
     allocate(myX(pcgsize(1))) ! Coordinates of nodes, used by ML preconditioner
     allocate(myY(pcgsize(1))) 
     allocate(myZ(pcgsize(1))) 
     call distributed_create_partition(ewn, nsn, (upn + 2) , uindx, pcgsize(1), myIndices, myX, myY, myZ)  ! Uses uindx mask to determine ice grid points.
     mySize = pcgsize(1)  ! Set variable for inittrilinos

     !write(*,*) "GlobalIDs myIndices..."
     !write(*,*) "pcgsize = ", pcgsize(1)
     !write(*,*) "myIndices = ", myIndices
     !call parallel_stop(__FILE__, __LINE__)
#else
     ! AGS: Get partition -- later this will be known by distributed glimmer
     call dopartition(pcgsize(1), mySize) 
     allocate(myIndices(mySize))
     allocate(myX(1)) ! Coordinates (only set for globalIDs)
     allocate(myY(1)) 
     allocate(myZ(1)) 
     call getpartition(mySize, myIndices) 

     if (distributed_execution()) then
         if (main_task) write(*,*) "Distributed Version cannot be run without globalIDs.  Stopping."
         call not_parallel(__FILE__, __LINE__)  ! Fatal if running without GlobalIDs in MPI
     endif

     !write(*,*) "Trilinos Generated Partition Map myIndices..."
     !write(*,*) "pcgsize = ", pcgsize(1), " mySize = ", mySize
     !write(*,*) "myIndices = ", myIndices
     !call parallel_stop(__FILE__, __LINE__)
#endif

     call inittrilinos(25, mySize, myIndices, myX, myY, myZ)   !TODO - Why 25 instead of 20 as above?

     ! Set if need full solution vector returned or just owned portion
#ifdef globalIDs
     ! We default Trilinos to return full solution vector.
     ! Per AGS, for both distributed and serial globalIDs cases, we must return just owned portion in order to prevent a permutation of the results.
     ! This was built into parallel_set_trilinos_return_vect that handled the difference between _single and _mpi versions.
     ! AGS: no longer needed, now the default
     ! call returnownedvector()
!     call parallel_set_trilinos_return_vect
#endif

     ! Triad sparse matrix not used in this case, so save on memory
     pcgsize(2) = 1

     ! JEFF: deallocate myIndices after the solve loop, because used in translation between globalIDs and local indices
     ! deallocate(myIndices)
  endif
#endif

!TODO This is the end of the block of code that is (mostly) cut and pasted from above.

!==============================================================================
! RN_20100126: End of the block
!==============================================================================

  ! For NOX JFNK
  allocate( vectx(2*pcgsize(1)), xk_1(2*pcgsize(1)), xk_1_plus(2*pcgsize(1)), gx_flag(2*pcgsize(1)) )
  ! allocate( uk_1(pcgsize(1)), vk_1(pcgsize(1)), g_flag(pcgsize(1)))

  ! *sfp** allocate space matrix variables
  allocate (pcgrow(pcgsize(2)),pcgcol(pcgsize(2)), rhsd(pcgsize(1)), rhsx(2*pcgsize(1)), &
            pcgval(pcgsize(2)))
  allocate(matrixA%row(pcgsize(2)), matrixA%col(pcgsize(2)), &
            matrixA%val(pcgsize(2)), answer(pcgsize(1)))
  allocate(matrixC%row(pcgsize(2)), matrixC%col(pcgsize(2)), &
            matrixC%val(pcgsize(2)))
  allocate(matrixtp%row(pcgsize(2)), matrixtp%col(pcgsize(2)), &
            matrixtp%val(pcgsize(2)))

  !*sfp* allocation for storage of (block matrix) off diagonal terms in coeff sparse matrix
  ! (these terms are usually sent to the RHS and treated as a source term in the operator splitting
  ! done in the standard Picard iteration)
  allocate(pcgrowuv(pcgsize(2)),pcgcoluv(pcgsize(2)),pcgvaluv(pcgsize(2)))
  allocate(pcgrowvu(pcgsize(2)),pcgcolvu(pcgsize(2)),pcgvalvu(pcgsize(2)))
  allocate(matrixAuv%row(pcgsize(2)),matrixAuv%col(pcgsize(2)),matrixAuv%val(pcgsize(2)))
  allocate(matrixAvu%row(pcgsize(2)),matrixAvu%col(pcgsize(2)),matrixAvu%val(pcgsize(2)))

  allocate( uk_1(pcgsize(1)), vk_1(pcgsize(1)),g_flag(pcgsize(1)) )
  allocate( vectp(pcgsize(1)), uk_1_plus(pcgsize(1)), vk_1_plus(pcgsize(1)) )
  allocate( F(2*pcgsize(1)), F_plus(2*pcgsize(1)), dx(2*pcgsize(1)) )
  allocate( wk1(2*pcgsize(1)), wk2(2*pcgsize(1)), rhs(2*pcgsize(1)) )
  allocate( vv(2*pcgsize(1),img1), wk(2*pcgsize(1), img) )

  call init_resid_type(resid_object, model, uindx, umask, d2thckdewdns, d2usrfdewdns, &
                       pcgsize, gx_flag, matrixA, matrixC, L2norm, ewn, nsn)
  fptr => resid_object
  c_ptr_to_object = c_loc(fptr)

  call ghost_preprocess_jfnk( ewn, nsn, upn, uindx, ughost, vghost, &
                         xk_1, uvel, vvel, gx_flag, pcgsize(1)) ! jfl_20100430

! SLAP incompatibility with main_task
if (main_task) then
  print *, ' '
  print *, 'Running Payne/Price higher-order dynamics with JFNK solver' 
end if

  calcoffdiag = .false.    ! save off diag matrix components
 call t_stopf("JFNK_pre")

!==============================================================================
! Beginning of Newton loop. Solves F(x) = 0 for x where x = [v, u] and
!                                                       F = [Fv(u,v), Fu(u,v)] 
!==============================================================================

!TODO - Anything to do here?  Is NOX now standard?
!
!       Should we eliminate the standalone JFNK solver?
!       We want to maintain a serial SLAP solver for glam, but does it have to support JFNK?

! UNCOMMENT these lines to switch to NOX's JFNK
! AGS: To Do:  send in distributed xk_1, or myIndices array, for distributed nox
#ifdef TRILINOS 

 call t_startf("JFNK_noxinit")
  call noxinit(xk_size, xk_1, 1, c_ptr_to_object)
 call t_stopf("JFNK_noxinit")

 call t_startf("JFNK_noxsolve")
  call noxsolve(xk_size, xk_1, c_ptr_to_object)
 call t_stopf("JFNK_noxsolve")

 call t_startf("JFNK_noxfinish")
  call noxfinish()
 call t_stopf("JFNK_noxfinish")

  kmax = 0     ! turn off native JFNK below
#endif

!==============================================================================
! JFNK loop: calculate F(u^k-1,v^k-1)
!==============================================================================

!TODO - Remove not_parallel calls if this will always be run with SLAP?

  ! This do loop is only used for SLAP.  Do not parallelize.
  do k = 1, kmax

   call t_startf("JFNK_SLAP")

!WHL - commenting out this call for now
!!!    call not_parallel(__FILE__, __LINE__)

!    calcoffdiag = .true.    ! save off diag matrix components
!    calcoffdiag = .false.    ! save off diag matrix components

    call calc_F (xk_1, F, xk_size, c_ptr_to_object, 0)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    L2norm = fptr%L2norm
    matrixA = fptr%matrixA
    matrixC = fptr%matrixC

!   calcoffdiag = .false.    ! next time calling calc_F, DO NOT save off diag matrix components

    L2norm_wig = sqrt(DOT_PRODUCT(F,F)) ! with ghost

!==============================================================================
! -define nonlinear target (if k=1)
! -check at all k if target is reached
!==============================================================================

    if (k == 1) NL_target = NL_tol * (L2norm_wig + 1.0e-2)

    print *, 'L2 w/ghost (k)= ',k,L2norm_wig,L2norm

    if (L2norm_wig < NL_target) exit ! nonlinear convergence criterion

!==============================================================================
! solve J(u^k-1,v^k-1)dx = -F(u^k-1,v^k-1) with fgmres, dx = [dv, du]  
!==============================================================================

    rhs = -1.d0*F

    dx  = 0.d0 ! initial guess

    call forcing_term (k, L2norm_wig, gamma_l)

    tol = gamma_l * L2norm_wig ! setting the tolerance for fgmres

    epsilon = 1.d-07 ! for J*vector approximation

    maxiteGMRES = 300
      
    iout   = 0    ! set  higher than 0 to have res(ite)

    icode = 0

!TODO - Can we avoid GOTO and CONTINUE?  Very old-style Fortran.

 10 CONTINUE
      
      call fgmres (2*pcgsize(1),img,rhs,dx,itenb,vv,wk,wk1,wk2, &
                   tol,maxiteGMRES,iout,icode,tot_its)

      IF ( icode == 1 ) THEN   ! precond step: use of Picard linear solver
                               ! wk2 = P^-1*wk1
        call apply_precond_nox( wk2, wk1, xk_size, c_ptr_to_object )
        GOTO 10
      ELSEIF ( icode >= 2 ) THEN  ! matvec step: Jacobian free approach
                                  ! J*wk1 ~ wk2 = (F_plus - F)/epsilon
         
! form  v^k-1_plus = v^k-1 + epsilon*wk1v. We use solver_postprocess to 
! transform vk_1_plus from a vector to a 3D field. (same idea for u^k-1_plus)
        vectx(:) = wk1(1:2*pcgsize(1)) ! for v and u
        xk_1_plus = xk_1 + epsilon*vectx

! form F(x + epsilon*wk1) = F(u^k-1 + epsilon*wk1u, v^k-1 + epsilon*wk1v)
        call calc_F (xk_1_plus, F_plus, xk_size, c_ptr_to_object, 1)

! put approximation of J*wk1 in wk2

        wk2 =  ( F_plus - F ) / epsilon

        GOTO 10
      ENDIF

!------------------------------------------------------------------------
! End of FGMRES method    
!------------------------------------------------------------------------
      if (tot_its == maxiteGMRES) then
        print *,'WARNING: FGMRES has not converged'
        stop
      endif

! icode = 0 means that fgmres has finished and sol contains the app. solution

!------------------------------------------------------------------------
! Update solution vectors (x^k = x^k-1 + dx) and 3D fields 
!------------------------------------------------------------------------
      xk_1 = xk_1 + dx(1:2*pcgsize(1))

  call t_stopf("JFNK_SLAP")
 end do   ! k = 1, kmax 

 call t_startf("JFNK_post")
! (need to update these values from fptr%uvel,vvel,stagthck etc)
  call solver_postprocess_jfnk( ewn, nsn, upn, uindx, xk_1, vvel, uvel, ghostbvel, pcgsize(1) )
  call ghost_postprocess_jfnk( ewn, nsn, upn, uindx, xk_1, ughost, vghost, pcgsize(1) )

    ! call fraction of assembly routines, passing current vel estimates (w/o manifold
    ! correction!) to calculate consistent basal tractions
    !
    ! *SFP* NOTE that if wanting to use basal tractions for the Newton method of converging on a
    ! coulomb-friction basasl BC, must update basal tractions estimate at EACH nonlinear iteration.
    ! In this case, the following two calls need to sit INSIDE of the do loop above. They are left
    ! out here because the current implementation of NOX skips to the end of this do loop, in order
    ! to skip JFs original implementation of JFNK (and jumping out of the do loop means these calls
    ! are skipped if they are inside of the do loop).
    !
    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     2,           efvs,           &
                     vvel,        uvel,           &
                     thck,        dusrfdns,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta,        btraction,      &
                     k, 1 )

   call findcoefstr(ewn,  nsn,   upn,             &
                     dew,  dns,   sigma,          &
                     1,           efvs,           &
                     uvel,        vvel,           &
                     thck,        dusrfdew,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta,        btraction,      &
                     k, 1 )

  inisoln = .true.

  print*,"Solution vector norm after JFNK = " ,sqrt(DOT_PRODUCT(xk_1,xk_1))

!TODO - The remaining code in this subroutine is cut and pasted from above.
!       Can we encapsulate this repeated code in a subroutine?

!TODO - I don't think uflx and vflux are needed.

  do ns = 1+staggered_shalo,size(umask,2)-staggered_nhalo
      do ew = 1+staggered_whalo,size(umask,1)-staggered_ehalo
      ! *sfp** calc. fluxes from converged vel. fields (for input to thickness evolution subroutine)
         if (umask(ew,ns) > 0) then
             uflx(ew,ns) = vertintg(upn, sigma, uvel(:,ew,ns)) * stagthck(ew,ns)
             vflx(ew,ns) = vertintg(upn, sigma, vvel(:,ew,ns)) * stagthck(ew,ns)
         end if
      end do
  end do

  ! JEFF: Deallocate myIndices which is used to intialize Trilinos
  if (whatsparse == STANDALONE_TRILINOS_SOLVER) then
     deallocate(myIndices)
     deallocate(myX)
     deallocate(myY)
     deallocate(myZ)
  endif

  ! *sfp* de-allocation of sparse matrix solution variables 
  deallocate(uindx)
  deallocate(pcgval,pcgrow,pcgcol,rhsd, rhsx)
  deallocate(pcgvaluv,pcgrowuv,pcgcoluv)
  deallocate(pcgvalvu,pcgrowvu,pcgcolvu)
  deallocate(matrixA%row, matrixA%col, matrixA%val)
  deallocate(matrixAuv%row, matrixAuv%col, matrixAuv%val)
  deallocate(matrixAvu%row, matrixAvu%col, matrixAvu%val)
  deallocate(matrixC%row, matrixC%col, matrixC%val)
  deallocate(matrixtp%row, matrixtp%col, matrixtp%val)
!  deallocate(uk_1, vk_1, xk_1, g_flag, gx_flag)
  deallocate(answer, dx, vectp, uk_1_plus, vk_1_plus, xk_1_plus, gx_flag )
  deallocate(F, F_plus)
  deallocate(wk1, wk2)
  deallocate(vv, wk)

 !model%velocity%uvel = uvel
 !model%velocity%vvel = vvel
 !model%velocity%uflx = uflx
 !model%velocity%vflx = vflx
 !model%stress%efvs = efvs

!HALO - Not sure whether these are needed.  Where does JFNK do its parallel halo updates for uvel, vvel?

 !PW following are needed for glam_velo_fordsiapstr - putting here until can be convinced
 !   that they are not needed (or that they should be delayed until later)
  call staggered_parallel_halo(uvel)
  call staggered_parallel_halo(vvel)

!HALO - Not sure we need these two updates
!       I think we do not need an update for efvs, because it is already computed in a layer of halo cells.
!       And I think we don't need an update for btraction, because it is computed in bodyset for all
!        locally owned velocity points.
  call parallel_halo(efvs)
  call parallel_halo(btraction)

!HALO - Probably do not need these two updates
  call staggered_parallel_halo(uflx)
  call staggered_parallel_halo(vflx)

 call t_stopf("JFNK_post")

  return

end subroutine JFNK_velo_solver

!***********************************************************************

function indxvelostr(ewn,  nsn,  upn,  &
                     mask, pointno)
!if a point from the 2d array 'mask' is associated with non-zero ice thickness, 
! (either a boundary or interior point) give it a unique number. If not, give it a zero.
  use parallel
  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, intent(in), dimension(:,:) :: mask
  integer, intent(out) :: pointno

  integer :: ew, ns
  integer, dimension(size(mask,1),size(mask,2)) :: indxvelostr

  pointno = 1

!TODO - I think these should be over locally owned velocity points: (ilo-1:ihi, jlo-1:jhi).

  do ns = 1+staggered_shalo,size(mask,2)-staggered_nhalo
     do ew = 1+staggered_whalo,size(mask,1)-staggered_ehalo
        if ( GLIDE_HAS_ICE( mask(ew,ns) ) ) then
          indxvelostr(ew,ns) = pointno
          pointno = pointno + 1
        else
          indxvelostr(ew,ns) = 0
        end if
      end do
  end do

  ! add two ghost points at upper and lower boundaries (needed for sfc and basal bcs)
  pointno = (pointno - 1) * (upn + 2)

  return

end function indxvelostr

!***********************************************************************

subroutine findefvsstr(ewn,  nsn, upn,       &
                       stagsigma, counter,   &
                       whichefvs, efvs,      &
                       uvel,      vvel,      &
                       flwa,      thck,      &
                       dusrfdew,  dthckdew,  &
                       dusrfdns,  dthckdns,  &
                       mask)

  ! calculate the eff. visc.    
  use parallel
  implicit none

  integer, intent(in) :: ewn, nsn, upn
  real(dp), intent(in), dimension(:)     :: stagsigma
  real(dp), intent(in), dimension(:,:,:) :: uvel, vvel, flwa
  real(dp), intent(inout), dimension(:,:,:) :: efvs
  real(dp), intent(in), dimension(:,:) :: thck, dthckdew, dusrfdew, &
                                                  dusrfdns, dthckdns
  integer, intent(in), dimension(:,:) :: mask
  integer, intent(in) :: whichefvs, counter

  integer :: ew, ns, up

  real(dp), dimension(size(efvs,1)) :: effstr, ugradup, vgradup, &
                                               ugradew, ugradns, vgradew, vgradns

  integer, dimension(2) :: mew, mns

  ! This is the factor 1/4(X0/H0)^2 in front of the term ((dv/dz)^2+(du/dz)^2) 
  real(dp), parameter :: f1 = 0.25d0 * (len0 / thk0)**2

  if (counter == 1) then

!  if (main_task) then
!    print *, 'nsn=', nsn
!    print *, 'ewn=', ewn
!    print *, 'uvel shape =', shape(uvel)
!    print *, 'vvel shape =', shape(vvel)
!    print *, 'thck shape =', shape(thck)
!    print *, 'efvs shape =', shape(efvs)
!    print *, 'flwafact shape =', shape(flwafact)
!  endif

!TODO - If we are not supporting glam_strs2 together with the old Glimmer temperature routines,
!       then we can assume that temp and flwa live on the staggered vertical grid.

!whl - If temp and flwa live on the staggered vertical grid (like the effective viscosity),
!      then the size of flwa is (upn-1), and vertical averaging of flwa is not needed here.

!HALO - These loops should be over locally owned cells, plus one halo layer (ilo-1:ihi+1, jlo-1:jhi+1).

     if (size(flwa,1)==upn-1) then   ! temperature and flwa live on staggered vertical grid

        do ns = 2,nsn-1
        do ew = 2,ewn-1
           if (thck(ew,ns) > 0.d0) then
              ! This is the rate factor term in the expression for the eff. visc: 1/2*A^(-1/n).
              ! If both temperature and eff. visc. live on a staggered grid in the vertical, then
              !  no vertical averaging is needed.
              flwafact(1:upn-1,ew,ns) = 0.5d0 * flwa(1:upn-1,ew,ns)**p1
           end if
        end do
        end do

     else  ! size(flwa,1)=upn; temperature and flwa live on unstaggered vertical grid

       do ns = 2,nsn-1
       do ew = 2,ewn-1
          if (thck(ew,ns) > 0.d0) then
             ! this is the rate factor term in the expression for the eff. visc: 1/2*A^(-1/n),
             ! which is averaged to midpoints in the vertical (i.e. it lives on a staggered 
             ! grid in the vertical, which is the case for "efvs" as well).
             forall (up = 1:upn-1) flwafact(up,ew,ns) = 0.5d0 * (sum(flwa(up:up+1,ew,ns)) / 2.d0)**p1
          end if
       end do
       end do

     end if   ! present(flwa_vstag)
  endif       ! counter

  select case(whichefvs)

  case(0)       ! calculate eff. visc. using eff. strain rate

!HALO - These loops should be over (ilo-1:ihi+1, jlo-1:jhi+1).
!       In other words, efvs is required in locally owned cells, plus one halo layer.
!       This is equivalent to the loop below, provided nhalo = 2.
!       As is, this code should *not* be run with nhalo = 1.
!       To run with nhalo = 1, we would need to compute efvs in locally owned cells and then do a halo update. 

  do ns = 2,nsn-1
      do ew = 2,ewn-1
        if (thck(ew,ns) > 0.d0) then
    ! The hsum() is on the unstaggered grid picking up the four points.  
    ! Then there is a derivative in the vertical direction.  
            ugradup = vertideriv(upn, hsum(uvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns))
            vgradup = vertideriv(upn, hsum(vvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns))

            ugradew = horizderiv(upn,  stagsigma,        &
                         sum(uvel(:,ew-1:ew,ns-1:ns),3), &
                         dew4, ugradup,                  &
                         sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                         sum(dthckdew(ew-1:ew,ns-1:ns)))

            vgradew = horizderiv(upn,  stagsigma,        &
                         sum(vvel(:,ew-1:ew,ns-1:ns),3), &
                         dew4, vgradup,                  &
                         sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                         sum(dthckdew(ew-1:ew,ns-1:ns)))

            ugradns = horizderiv(upn,  stagsigma,        &
                         sum(uvel(:,ew-1:ew,ns-1:ns),2), &
                         dns4, ugradup,                  &
                         sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                         sum(dthckdns(ew-1:ew,ns-1:ns)))

            vgradns = horizderiv(upn,  stagsigma,        &
                         sum(vvel(:,ew-1:ew,ns-1:ns),2), &
                         dns4, vgradup,                  &
                         sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                         sum(dthckdns(ew-1:ew,ns-1:ns)))

            ! "effstr" = eff. strain rate squared
            effstr = ugradew**2 + vgradns**2 + ugradew*vgradns + &
                         0.25d0 * (vgradew + ugradns)**2 + &
!                         f1 * (ugradup**2 + vgradup**2)      ! make line ACTIVE for "capping" version (see note below)   
                         f1 * (ugradup**2 + vgradup**2) + effstrminsq ! make line ACTIVE for new version

    ! -----------------------------------------------------------------------------------
    ! NOTES on capping vs. non-capping version of eff. strain rate calc.
    ! -----------------------------------------------------------------------------------
    !
    ! Set eff. strain rate (squared) to some min value where it falls below some 
    ! threshold value, 'effstrminsq'. Commented out the old version below, which "caps" 
    ! the min eff strain rate (and thus the max eff visc) in favor of a version that 
    ! leads to a "smooth" description of eff strain rate (and eff visc). The change for 
    ! new version is that the value of 'effstrminsq' simply gets added in with the others
    ! (e.g. how it is done in the Pattyn model). The issues w/ the capping approach are 
    ! discussed (w.r.t. sea ice model) in: Lemieux and Tremblay, JGR, VOL. 114, C05009, 
    ! doi:10.1029/2008JC005017, 2009). Long term, the capping version should probably be
    !  available as a config file option or possibly removed altogether.   

    ! Old "capping" scheme       ! these lines must be active to use the "capping" scheme for the efvs calc
!            where (effstr < effstrminsq)
!                   effstr = effstrminsq
!            end where

    ! Note that the vert dims are explicit here, since glide_types defines this 
    ! field as having dims 1:upn. This is something that we'll have to decide on long-term;
    ! should efvs live at cell centroids in the vert (as is assumed in this code)
    ! or should we be doing some one-sided diffs at the sfc/bed boundaries so that it has vert dims 
    ! of upn? For now, we populate ONLY the first 1:upn-1 values of the efvs vector and leave the one
    ! at upn empty (the Pattyn/Bocek/Johnson core would fill all values, 1:upn).

    ! NOTE also that efvs lives on the non-staggered grid in the horizontal. That is, in all of the 
    ! discretizations conducted below, efvs is explicitly averaged from the normal horiz grid onto the 
    ! staggered horiz grid (Thus, in the calculations, efvs is treated as if it lived on the staggered 
    ! horiz grid, even though it does not). 

            ! Below, p2=(1-n)/2n. The 1/2 is from taking the sqr root of the squared eff. strain rate
            efvs(1:upn-1,ew,ns) = flwafact(1:upn-1,ew,ns) * effstr**p2 + homotopy
!            efvs(:,ew,ns) = flwafact(:,ew,ns) * effstr**p2

        else
           efvs(:,ew,ns) = effstrminsq ! if the point is associated w/ no ice, set to min value
        end if

       end do   ! end ew
   end do       ! end ns

  case(1)       ! set the eff visc to some const value 

!HALO - Make this loop consistent with loop above.

!   *sfp* changed default setting for linear viscosity so that the value of the rate
!   factor is taken into account
  do ns = 2,nsn-1
      do ew = 2,ewn-1
       if (thck(ew,ns) > 0.d0) then
! KJE code used to have this
!       efvs(1:upn-1,ew,ns) = 0.5d0 * flwa(1:upn-1,ew,ns)**(-1.d0)
        efvs(1:upn-1,ew,ns) = flwafact(1:upn-1,ew,ns)
        else
           efvs(:,ew,ns) = effstrminsq ! if the point is associated w/ no ice, set to min value
       end if
      end do
  end do

  end select

! JEFF Halo does NOT verify, because used a staggered array to hold unstaggered data.  The unstaggered data fits because remove two rows and columns of data.
! The current parallel_halo(efvs) routine won't update a staggered array.
! I think it is OK, so I'm passing for now.
!  if (.NOT. parallel_halo_verify(efvs)) then
!      ! efvs is an unstaggered grid in an staggered declaration.  Steve Price said he reused the variable from the other core.
!      write(*,*) "Halo Verify failed for findefvstr()"
!      call parallel_stop(__FILE__, __LINE__)
!  endif

  return
end subroutine findefvsstr

!***********************************************************************

function vertideriv(upn, varb, thck)

  implicit none

  integer, intent(in) :: upn
  real(dp), intent(in), dimension(:) :: varb
  real(dp), intent(in) :: thck

  real(dp), dimension(size(varb)-1) :: vertideriv
  !'dupm' is defined as -1/(2*del_sigma), in which case it seems like 
  ! there should be a '-' in front of this expression ... but note that
  ! the negative sign is implicit in the fact that the vertical index 
  ! increases moving downward in the ice column (up=1 is the sfc, 
  ! up=upn is the bed).

  vertideriv(1:upn-1) = dupm * (varb(2:upn) - varb(1:upn-1)) / thck

  return

end function vertideriv

!***********************************************************************

function horizderiv(upn,     stagsigma,   &
                    varb,    grid,        &
                    dvarbdz, dusrfdx, dthckdx)

  implicit none

  integer, intent(in) :: upn
  real(dp), dimension(:), intent(in) :: stagsigma
  real(dp), dimension(:,:), intent(in) :: varb
  real(dp), dimension(:), intent(in) :: dvarbdz
  real(dp), intent(in) :: dusrfdx, dthckdx, grid

  real(dp) :: horizderiv(size(varb,1)-1)

  horizderiv = (varb(1:upn-1,2) + varb(2:upn,2) - varb(1:upn-1,1) - varb(2:upn,1)) / grid - &
                dvarbdz * (dusrfdx - stagsigma * dthckdx) / 4.d0

  return

end function horizderiv

!***********************************************************************

function getlocrange(upn, indx)

  implicit none

  integer, intent(in) :: upn
  integer, intent(in) :: indx
  integer, dimension(2) :: getlocrange

  getlocrange = (indx - 1) * (upn + 2) + 1 + (/ 1, upn /)

  return

end function getlocrange

!***********************************************************************

function getlocationarray(ewn, nsn, upn, mask, indxmask)

  use parallel

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:), intent(in) :: mask
  integer, dimension(:,:), intent(in) :: indxmask
  integer, dimension(ewn,nsn,2) :: getlocationarray
  integer :: ew, ns

#ifdef globalIDs
  ! Returns in (:,:,1) the global ID bases for each grid point, including 
  ! halos and those without ice.
  ! Since the code checks elsewhere whether ice occurs at a given grid point, 
  ! this information is not encoded here. For the local indices (see below)
  ! the mask information is used since ice-free grid points are not indexed
  ! locally

!TODO - Not sure if these loops are correct because I don't understand what this subroutine is doing.
!       Is the input mask on the scalar (ice) grid? 
  do ns=1,nsn
    do ew=1,ewn
      getlocationarray(ew,ns,1) = parallel_globalID(ns, ew, upn + 2)  ! Extra two layers for ghost layers
    end do
  end do

  ! Returns in (:,:,2) the local index base for each ice grid point 
  !  (same indices as those used in myIndices)
  ! indxmask is ice mask with non-zero values for cells with ice.
  ! If a point (ew,ns) doesn't have ice, then value is set to 0.
  ! If a point (ew,ns) is in the halo, value is also set to 0.
  ! upn+2 is the total number of vertical layers including any ghosts
  ! (logic modelled after distributed_create_partition)

  ! initialize to zero (in order to set halo and ice-free cells to zero)
  getlocationarray(:,:,2) = 0

!TODO - Should this be over locally owned velocity points?

  ! Step through indxmask, but exclude halo
  do ns = 1+staggered_shalo,size(indxmask,2)-staggered_nhalo
    do ew = 1+staggered_whalo,size(indxmask,1)-staggered_ehalo
      if ( indxmask(ew,ns) /= 0 ) then
        getlocationarray(ew,ns,2) = (indxmask(ew,ns) - 1) * (upn+2) + 1
      endif
    end do
  end do

#else
  integer, dimension(ewn,nsn) :: temparray
  integer :: cumsum

  ! initialize to zero
  cumsum = 0
  temparray = 0
  getlocationarray = 0

!TODO - Should this be over locally owned velocity points?

  do ns=1+staggered_shalo,size(mask,2)-staggered_nhalo
    do ew=1+staggered_whalo,size(mask,1)-staggered_ehalo
      if ( GLIDE_HAS_ICE( mask(ew,ns) ) ) then
        cumsum = cumsum + ( upn + 2 )
        getlocationarray(ew,ns,1) = cumsum
        temparray(ew,ns) = upn + 2
      else
        getlocationarray(ew,ns,1) = 0
        temparray(ew,ns) = 1
      end if
    end do
  end do

  getlocationarray(:,:,1) = ( getlocationarray(:,:,1) + 1 ) - temparray
  getlocationarray(:,:,2) = getlocationarray(:,:,1)
#endif

    return

end function getlocationarray

!***********************************************************************

function slapsolvstr(ewn, nsn, upn, &
                     vel, uindx, its, answer )

! *sp* routine to solve Ax=b sparse matrix problem 

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  real(dp), dimension(:,:,:), intent(in) :: vel
  integer, dimension(:,:), intent(in) :: uindx

  real(dp), dimension(:), intent(out) :: answer

  real(dp), dimension(size(vel,1),size(vel,2),size(vel,3)) :: slapsolvstr
  integer, intent(inout) :: its

  integer :: ew, ns

  real(dp), dimension(:), allocatable :: rwork
  integer, dimension(:), allocatable :: iwork

  real(dp), parameter :: tol = 1.0d-12
  real(dp) :: err

  integer, parameter :: isym = 0, itol = 2, itmax = 100
  integer, dimension(2) :: loc
  integer :: iter, ierr, mxnelt

! ** move to values subr   

  pcgsize(2) = ct - 1

  call ds2y(pcgsize(1),pcgsize(2),pcgrow,pcgcol,pcgval,isym)

!** plot the matrix to check that it has the correct form
!call dcpplt(pcgsize(1),pcgsize(2),pcgrow,pcgcol,pcgval,isym,ulog)      

  mxnelt = 60 * pcgsize(1); allocate(rwork(mxnelt),iwork(mxnelt))

!**     solve the problem using the SLAP package routines     
!**     -------------------------------------------------
!**     n ... order of matrix a (in)
!**     b ... right hand side vector (in)                        
!**     x ... initial quess/final solution vector (in/out)                        
!**     nelt ... number of non-zeroes in A (in)
!**     ia, ja ... sparse matrix format of A (in)
!**     a ... matrix held in SLAT column format (in)
!**     isym ... storage method (0 is complete) (in)
!**     itol ... convergence criteria (2 recommended) (in)                     
!**     tol ... criteria for convergence (in)
!**     itmax ... maximum number of iterations (in)
!**     iter ... returned number of iterations (out)
!**     err ... error estimate of solution (out)
!**     ierr ... returned error message (0 is ok) (out)
!**     iunit ... unit for error writes during iteration (0 no write) (in)
!**     rwork ... workspace for SLAP routines (in)
!**     mxnelt ... maximum array and vector sizes (in)
!**     iwork ... workspace for SLAP routines (in)

!TODO - Are loop bounds OK? Since this is for the serial SLAP solver, I think so.

! *sp* initial estimate for vel. field?
  do ns = 1,nsn-1
  do ew = 1,ewn-1
   if (uindx(ew,ns) /= 0) then
    loc = getlocrange(upn, uindx(ew,ns))
    answer(loc(1):loc(2)) = vel(:,ew,ns)
    answer(loc(1)-1) = vel(1,ew,ns)
    answer(loc(2)+1) = vel(upn,ew,ns)
   end if
  end do
  end do

  call dslucs(pcgsize(1),rhsd,answer,pcgsize(2),pcgrow,pcgcol,pcgval, &
              isym,itol,tol,itmax,iter,err,ierr,0,rwork,mxnelt,iwork,mxnelt)

  if (ierr /= 0) then
    print *, 'pcg error ', ierr, itmax, iter, tol, err
    ! stop
  end if

  deallocate(rwork,iwork)

!TODO - Are loop bounds OK? Since this is for the serial SLAP solver, I think so.

  do ns = 1,nsn-1
  do ew = 1,ewn-1
     if (uindx(ew,ns) /= 0) then
       loc = getlocrange(upn, uindx(ew,ns))
       slapsolvstr(:,ew,ns) = answer(loc(1):loc(2))
     else
       slapsolvstr(:,ew,ns) = 0.d0
     end if
  end do
  end do

  its = its + iter

  return

end function slapsolvstr

! *****************************************************************************

subroutine solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, vel )

  ! Puts sparse matrix variables in SLAP triad format into "matrix" derived type, 
  ! so that it can be passed to the generic solver wrapper, "sparse_easy_solve". 
  ! Takes place of the old, explicit solver interface to SLAP linear solver.
  use parallel

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  real(dp), dimension(:,:,:), intent(in) :: vel
  integer, dimension(:,:), intent(in) :: uindx
  type(sparse_matrix_type), intent(inout) :: matrix
  real(dp), dimension(:), intent(out) :: answer

  integer :: ew, ns
  integer, dimension(2) :: loc

  pcgsize(2) = ct - 1

  matrix%order = pcgsize(1)
  matrix%nonzeros = pcgsize(2)
  matrix%symmetric = .false.

  matrix%row = pcgrow
  matrix%col = pcgcol
  matrix%val = pcgval

!HALO - Not sure about the loops here.  Should be over locally owned velocity points?

  ! Initial estimate for vel. field; take from 3d array and put into
  ! the format of a solution vector.
  do ns = 1+staggered_shalo,size(uindx,2)-staggered_nhalo
   do ew = 1+staggered_whalo,size(uindx,1)-staggered_ehalo
        if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            answer(loc(1):loc(2)) = vel(:,ew,ns)
            answer(loc(1)-1) = vel(1,ew,ns)
            answer(loc(2)+1) = vel(upn,ew,ns)

            !JEFF Verifying Trilinos Input
            ! write(*,*) "Initial answer at (", ew, ", ", ns, ") = ", answer(loc(1)-1:loc(2)+1)
        end if
    end do
  end do

end subroutine solver_preprocess

!***********************************************************************

subroutine solver_postprocess( ewn, nsn, upn, pt, uindx, answrapped, ansunwrapped, ghostbvel )

  ! Unwrap the vels from the solution vector and place into a 3d array.
  use parallel

  implicit none

  integer, intent(in) :: ewn, nsn, upn, pt
  integer, dimension(:,:), intent(in) :: uindx
  real(dp), dimension(:), intent(in) :: answrapped
  real(dp), dimension(upn,ewn-1,nsn-1), intent(out) :: ansunwrapped
  real(dp), dimension(:,:,:,:), intent(inout) :: ghostbvel   

  integer, dimension(2) :: loc
  integer :: ew, ns

!HALO - Not sure about the loops here.  Should be over locally owned velocity points?

  do ns = 1+staggered_shalo,size(uindx,2)-staggered_nhalo
      do ew = 1+staggered_whalo,size(uindx,1)-staggered_ehalo
          if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            ansunwrapped(:,ew,ns) = answrapped(loc(1):loc(2))
            !! save the fictitious basal velocities for basal traction calculation !!
            ghostbvel(pt,:,ew,ns) = answrapped( loc(2)-1:loc(2)+1 )  
          else
            ansunwrapped(:,ew,ns) = 0.d0
          end if
      end do
  end do

end subroutine solver_postprocess

!***********************************************************************

subroutine solver_postprocess_jfnk( ewn, nsn, upn, uindx, answrapped, ansunwrappedv, &
                                    ansunwrappedu, ghostbvel, pcg1 )

   ! Unwrap the vels from the solution vector and place into a 3d array.
   use parallel

   implicit none

   integer :: pcg1
   integer, intent(in) :: ewn, nsn, upn
   integer, dimension(:,:), intent(in) :: uindx
   real(dp), dimension(:), intent(in) :: answrapped
   real(dp), dimension(upn,ewn-1,nsn-1), intent(out) :: ansunwrappedv, ansunwrappedu
   real(dp), dimension(:,:,:,:), intent(inout) :: ghostbvel

   integer, dimension(2) :: loc
   integer :: ew, ns

!HALO - Not sure about the loops here.  Should be over locally owned velocity points?

   do ns = 1+staggered_shalo,size(uindx,2)-staggered_nhalo
       do ew = 1+staggered_whalo,size(uindx,1)-staggered_ehalo
           if (uindx(ew,ns) /= 0) then
             loc = getlocrange(upn, uindx(ew,ns))
             ansunwrappedv(:,ew,ns) = answrapped(loc(1):loc(2))
             ansunwrappedu(:,ew,ns) = answrapped(pcg1+loc(1):pcg1+loc(2))
             !! save the fictitious basal velocities for basal traction calculation !!
             ghostbvel(2,:,ew,ns) = answrapped( loc(2)-1:loc(2)+1 )
             ghostbvel(1,:,ew,ns) = answrapped( pcg1+loc(2)-1:pcg1+loc(2)+1 )
           else
             ansunwrappedv(:,ew,ns) = 0.d0
             ansunwrappedu(:,ew,ns) = 0.d0
           end if
       end do
   end do
end subroutine solver_postprocess_jfnk


!***********************************************************************

subroutine form_matrix( matrix ) ! for JFNK solver

  ! Puts sparse matrix variables in SLAP triad format into "matrix" 
  ! derived type. Similar to solver_preprocess but does not form answer vector

  implicit none

!  integer, intent(in) :: ewn, nsn, upn
  type(sparse_matrix_type), intent(inout) :: matrix

  pcgsize(2) = ct - 1

  matrix%order = pcgsize(1)
  matrix%nonzeros = pcgsize(2)
  matrix%symmetric = .false.

  matrix%row = pcgrow
  matrix%col = pcgcol
  matrix%val = pcgval

end subroutine form_matrix

!***********************************************************************

subroutine forcing_term ( k, L2normk_1, gamma_l ) 

  ! Calculates the forcing term (i.e. the factor that multiplies the initial
  ! L2 norm to determine the tolerance for the linear solve in the JFNK solver)
  ! at iteration k given the L2norm at k-1 and k-2.
  ! jfl, 10 Sept 2010

  ! See eq 2.6 in S.C. Eisenstat, H.F. Walker, Choosing the forcing terms in
  ! an inexact Newton method, SIAM J. Sci. Comput. 17 (1996) 16-32.

  implicit none
      
  integer, intent(in) :: k
  real(dp), intent(in) :: L2normk_1 ! L2 norm at k-1
  real(dp), intent(out):: gamma_l
  real(dp) :: gamma_ini, gamma_min, expo
  real(dp), save :: L2normk_2      ! L2 norm at k-2

      gamma_ini = 0.9d0
      gamma_min = 0.01d0
      expo      = 2d0

      if (k == 1) then
         gamma_l = gamma_ini
      else
         gamma_l = (L2normk_1 / L2normk_2)**expo
      endif

      if (gamma_l > gamma_ini) gamma_l = gamma_ini
      if (gamma_l < gamma_min) gamma_l = gamma_min
      
      L2normk_2 = L2normk_1

end subroutine forcing_term

!***********************************************************************

subroutine apply_precond( matrixA, matrixC, nu1, nu2, wk1, wk2, whichsparse ) 

  ! Apply preconditioner operator for JFNK solver: wk2 = P^-1 *wk1 
  ! The preconditioner operator is in fact taken from the Picard solver
  ! There is a splitting of the v (A matrix) and u (C matrix) equations
  ! Each component is solved to a loose tolerance (as opposed to Picard)

  implicit none

  integer, intent(in) :: nu1, nu2, whichsparse
  integer :: iter
  type(sparse_matrix_type), intent(in) :: matrixA, matrixC
  real(dp), dimension(nu2), intent(in) :: wk1
  real(dp), dimension(nu2), intent(out):: wk2
  real(dp), dimension(nu1) :: answer, vectp
  real(dp) :: err

! precondition v component 
       
!TODO - Add decimal points to real variables below

!      answer = 0d0 ! initial guess
      answer = 0.d0 ! initial guess
      vectp(:) = wk1(1:nu1) ! rhs for precond v
      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
         call sparse_easy_solve(matrixA, vectp, answer, err, iter, whichsparse, nonlinear_solver = nonlinear)
#ifdef TRILINOS
      else
         call restoretrilinosmatrix(0); 
         call solvewithtrilinos(vectp, answer, linearSolveTime) 
         totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
         write(*,*) 'Total linear solve time so far', totalLinearSolveTime
#endif
      endif
      wk2(1:nu1) = answer(:)

! precondition u component 
       
!      answer = 0d0 ! initial guess
      answer = 0.d0 ! initial guess
      vectp(:) = wk1(nu1+1:nu2) ! rhs for precond u
      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
         call sparse_easy_solve(matrixC, vectp, answer, err, iter, whichsparse, nonlinear_solver = nonlinear)
#ifdef TRILINOS
      else
         call restoretrilinosmatrix(1); 
         call solvewithtrilinos(vectp, answer, linearSolveTime) 
         totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
         write(*,*) 'Total linear solve time so far', totalLinearSolveTime
#endif
      endif
      wk2(nu1+1:nu2) = answer(:)

end subroutine apply_precond

!***********************************************************************

subroutine apply_precond_nox( wk2_nox, wk1_nox, xk_size, c_ptr_to_object )  bind(C, name='apply_precond_nox')

  ! Apply preconditioner operator for JFNK solver through NOX: wk2 = P^-1 *wk1 
  ! The preconditioner operator is in fact taken from the Picard solver
  ! There is a splitting of the v (A matrix) and u (C matrix) equations
  ! Each component is solved to a loose tolerance (as opposed to Picard)

  implicit none

! variables coming through from NOX
  integer(c_int) ,intent(in) ,value  :: xk_size
  real (c_double) ,intent(in)        :: wk1_nox(xk_size)
  real (c_double) ,intent(out)       :: wk2_nox(xk_size)
  type(pass_through) ,pointer        :: fptr=>NULL()
  type(c_ptr) ,intent(inout)         :: c_ptr_to_object

  integer :: nu1, nu2, whichsparse
  integer :: iter
  type(sparse_matrix_type) :: matrixA, matrixC
  real(dp), dimension(xk_size) :: wk1
  real(dp), dimension(xk_size) :: wk2
  real(dp), allocatable, dimension(:) :: answer, vectp
  real(dp) :: err

  call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr= resid_object

  matrixA = fptr%matrixA
  matrixC = fptr%matrixC
  whichsparse = fptr%model%options%which_ho_sparse
  pcgsize = fptr%pcgsize
           
  nu1 = pcgsize(1)
  nu2 = 2*pcgsize(1)
  allocate ( answer(nu1) )
  allocate ( vectp(nu1) )
  wk1  = wk1_nox

! ID as a test
!  wk2_nox  = wk1

! precondition v component 
       
!TODO - Add decimal points to real variables below

!      answer = 0d0 ! initial guess
      answer = 0.d0 ! initial guess
      vectp(:) = wk1(1:nu1) ! rhs for precond v
      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
         call sparse_easy_solve(matrixA, vectp, answer, err, iter, whichsparse, nonlinear_solver = nonlinear)
#ifdef TRILINOS
      else
         call restoretrilinosmatrix(0);
         call solvewithtrilinos(vectp, answer, linearSolveTime)
         totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
!         write(*,*) 'Total linear solve time so far', totalLinearSolveTime
#endif
      endif
      wk2(1:nu1) = answer(:)

! precondition u component 
       
      answer = 0.d0 ! initial guess
      vectp(:) = wk1(nu1+1:nu2) ! rhs for precond u
      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
         call sparse_easy_solve(matrixC, vectp, answer, err, iter, whichsparse, nonlinear_solver = nonlinear)
#ifdef TRILINOS
      else
         call restoretrilinosmatrix(1);
         call solvewithtrilinos(vectp, answer, linearSolveTime)
         totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
!         write(*,*) 'Total linear solve time so far', totalLinearSolveTime
#endif
      endif
      wk2(nu1+1:nu2) = answer(:)

  wk2_nox  = wk2

end subroutine apply_precond_nox

!***********************************************************************

subroutine reset_effstrmin (esm_factor) bind(C, name='reset_effstrmin')
  use iso_c_binding  
  real (c_double), intent(in):: esm_factor
 
  ! esm_factor of 0 leads to desired target. Valid values are [0,10]
!  effstrminsq = effstrminsq_target * 10.0**(2.0 * esm_factor)
  
  ! Homotopy parameter needs to be zero when esm_factor hits zero
  if (esm_factor > 1.0e-10) then
    homotopy = 10**( esm_factor - 9.0 )
  else
    homotopy = 0.0;
  endif

end subroutine reset_effstrmin

!***********************************************************************

!TODO - There is more repeated code here, making code maintenance difficult.
!       Would it be possible to package the repeated code into a subroutine called
!        from multiple places?

 subroutine calc_F (xtp, F, xk_size, c_ptr_to_object, ispert) bind(C, name='calc_F')

  ! Calculates either F(x) or F(x+epsilon*vect) for the JFNK method
  ! Recall that x=[v,u]
  ! xtp is both vtp and utp in one vector

  use iso_c_binding  
  use glide_types ,only : glide_global_type, pass_through
  use parallel

  implicit none

   integer(c_int) ,intent(in) ,value  :: xk_size
! ispert is 0 for base calculations, 1 for perturbed calculations
   integer(c_int) ,intent(in) ,value  :: ispert 
   real(c_double)  ,intent(in)        :: xtp(xk_size)
   real(c_double)  ,intent(out)       :: F(xk_size)
   type(pass_through) ,pointer        :: fptr=>NULL()
   type(c_ptr) ,intent(inout)         :: c_ptr_to_object

  integer :: ewn, nsn, upn, counter, whichbabc, whichefvs, i
  integer  ,dimension(2)   :: pcgsize
  integer  ,dimension(:) ,allocatable :: gxf ! 0 :reg cell
  integer  ,dimension(:,:) ,allocatable :: ui, um
  real(dp) :: dew, dns
  real(dp), dimension(:)  ,pointer :: sigma, stagsigma
  real(dp), dimension(:,:) ,pointer :: thck, dusrfdew, dthckdew, dusrfdns, dthckdns, &
                                         dlsrfdew, dlsrfdns, stagthck, lsrf, topg, minTauf, beta
  real(dp), dimension(:,:) ,pointer ::  d2usrfdew2, d2thckdew2, d2usrfdns2, d2thckdns2
  real(dp), dimension(:,:,:) ,pointer :: efvs, btraction
  real(dp), dimension(:,:,:) ,pointer :: uvel, vvel, flwa
  type(sparse_matrix_type) :: matrixA, matrixC
  real(dp), dimension(:) ,allocatable :: vectx
  real(dp), dimension(:) ,allocatable :: vectp

  real(dp) :: L2square
!  real(dp), intent(inout):: L2norm
  real(dp) :: L2norm

 call t_startf("Calc_F")
  call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr= resid_object
           
  ewn = fptr%model%general%ewn
  nsn = fptr%model%general%nsn
  upn = fptr%model%general%upn
  whichbabc = fptr%model%options%which_ho_babc
  whichefvs = fptr%model%options%which_ho_efvs
  dew = fptr%model%numerics%dew
  dns = fptr%model%numerics%dns
  sigma => fptr%model%numerics%sigma(:)
  stagsigma => fptr%model%numerics%stagsigma(:)
  thck => fptr%model%geometry%thck(:,:)
  lsrf => fptr%model%geometry%lsrf(:,:)
  topg => fptr%model%geometry%topg (:,:)
  stagthck => fptr%model%geomderv%stagthck(:,:)
  dthckdew => fptr%model%geomderv%dthckdew(:,:)
  dthckdns => fptr%model%geomderv%dthckdns(:,:)
  dusrfdew => fptr%model%geomderv%dusrfdew(:,:)
  dusrfdns => fptr%model%geomderv%dusrfdns(:,:)
  dlsrfdew => fptr%model%geomderv%dlsrfdew(:,:)
  dlsrfdns => fptr%model%geomderv%dlsrfdns(:,:)
  d2thckdew2 => fptr%model%geomderv%d2thckdew2(:,:)
  d2thckdns2 => fptr%model%geomderv%d2thckdns2(:,:)
  d2usrfdew2 => fptr%model%geomderv%d2usrfdew2(:,:)
  d2usrfdns2 => fptr%model%geomderv%d2usrfdns2(:,:)
  minTauf => fptr%model%basalproc%minTauf(:,:)
  beta => fptr%model%velocity%beta(:,:)
!intent (inout) terms
  btraction => fptr%model%velocity%btraction(:,:,:)
  flwa => fptr%model%temper%flwa(:,:,:)
  efvs => fptr%model%stress%efvs(:,:,:)
  uvel => fptr%model%velocity%uvel(:,:,:)
  vvel => fptr%model%velocity%vvel(:,:,:)
  L2norm = fptr%L2norm

  allocate( ui(ewn-1,nsn-1), um(ewn-1,nsn-1) )
  ui= fptr%ui
  um = fptr%um

  pcgsize = fptr%pcgsize
  allocate( gxf(2*pcgsize(1)) )

  gxf = fptr%gxf
! temporary to test JFNK -  need to take out
  counter = 1

  d2usrfdewdns = fptr%d2usrfcross
  d2thckdewdns = fptr%d2thckcross

  matrixA = fptr%matrixA
  matrixC = fptr%matrixC
  allocate( vectp( pcgsize(1)) )
  allocate( vectx(2*pcgsize(1)) )

  call solver_postprocess_jfnk( ewn, nsn, upn, ui, &
                                xtp, vvel, uvel, ghostbvel, pcgsize(1) )

    ! coordinate halos for updated uvel and vvel
    call staggered_parallel_halo(uvel)
    call staggered_parallel_halo(vvel)

    call findefvsstr(ewn,  nsn,  upn,       &
                     stagsigma,  counter,  &
                     whichefvs,  efvs,     &
                     uvel,       vvel,     &
                     flwa,       thck,     &
                     dusrfdew,   dthckdew, &
                     dusrfdns,   dthckdns, &
                     um)

!==============================================================================
! jfl 20100412: residual for v comp: Fv= A(utp,vtp)vtp - b(utp,vtp)  
!==============================================================================

    ! *sfp** calculation of coeff. for stress balance calc. 
    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     2,           efvs,           &
                     vvel,        uvel,           &
                     thck,        dusrfdns,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     ui,       um,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta,        btraction,      &
                     counter, 0 )

    rhsx(1:pcgsize(1)) = rhsd ! Fv

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
      call form_matrix ( matrixA ) ! to get A(utp,vtp)
#ifdef TRILINOS
    else
      if (ispert == 0) then
        call savetrilinosmatrix(0); 
      endif
#endif
    end if
    
    vectp = xtp(1:pcgsize(1))

 call t_startf("Calc_F_res_vect")
    call res_vect(matrixA, vectp, rhsd, pcgsize(1),  gxf, L2square, whatsparse)
 call t_stopf("Calc_F_res_vect")
    L2norm=L2square

    F(1:pcgsize(1)) = vectp(1:pcgsize(1)) 

!==============================================================================
! jfl 20100412: residual for u comp: Fu= C(utp,vtp)utp - d(utp,vtp)  
!==============================================================================

    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     1,           efvs,           &
                     uvel,        vvel,           &
                     thck,        dusrfdew,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     ui,       um,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta,        btraction,      &
                     counter, 0 )

    rhsx(pcgsize(1)+1:2*pcgsize(1)) = rhsd ! Fv

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
      call form_matrix ( matrixC ) ! to get C(utp,vtp)
#ifdef TRILINOS
    else
      if (ispert == 0) then
        call savetrilinosmatrix(1); 
      endif
#endif
    end if
    
    vectp(1:pcgsize(1)) = xtp(pcgsize(1)+1:2*pcgsize(1))

 call t_startf("Calc_F_res_vect")
    call res_vect(matrixC, vectp, rhsd, pcgsize(1), gxf, L2square, whatsparse)
 call t_stopf("Calc_F_res_vect")
    L2norm = sqrt(L2norm + L2square)

    F(pcgsize(1)+1:2*pcgsize(1)) = vectp(1:pcgsize(1)) 

!      vectx = xtp 
!
!   call res_vect_jfnk(matrixA, matrixC, vectx, rhsx, pcgsize(1), 2*pcgsize(1), gxf, L2square, whatsparse)
!    L2norm = L2square
!    F = vectx 

    call solver_postprocess_jfnk( ewn, nsn, upn, ui, xtp, vvel, uvel, ghostbvel, pcgsize(1) )

  fptr%model%velocity%btraction => btraction(:,:,:)
  fptr%model%temper%flwa => flwa(:,:,:)
  fptr%model%stress%efvs => efvs(:,:,:)
  fptr%model%velocity%uvel => uvel(:,:,:)
  fptr%model%velocity%vvel => vvel(:,:,:)

  fptr%L2norm = L2norm
  fptr%matrixA = matrixA
  fptr%matrixC = matrixC
 call t_stopf("Calc_F")

end subroutine calc_F

!***********************************************************************

subroutine ghost_preprocess( ewn, nsn, upn, uindx, ughost, vghost, & 
                             uk_1, vk_1, uvel, vvel, g_flag)

! puts vel values in  uk_1, vk_1 (including ghost values) and creates the
! ghost flag vector. uk_1, vk_1 and the ghost flag vector are used for 
! the residual calculation (jfl 20100430)

  use parallel

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:), intent(in) :: uindx
  integer, dimension(:), intent(out) :: g_flag 
  real(dp), dimension(2,ewn-1,nsn-1), intent(in) ::ughost,vghost 
  real(dp), dimension(:,:,:), intent(in) :: uvel, vvel
  real(dp), dimension(:), intent(out) :: uk_1, vk_1 

  integer :: ew, ns
  integer, dimension(2) :: loc

  g_flag = 0

!TODO - Loops should be over locally owned velocity points?

  do ns = 1+staggered_shalo,size(uindx,2)-staggered_nhalo
   do ew = 1+staggered_whalo,size(uindx,1)-staggered_ehalo
        if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            uk_1(loc(1):loc(2)) = uvel(:,ew,ns)
            uk_1(loc(1)-1)      = ughost(1,ew,ns) ! ghost at top
            uk_1(loc(2)+1)      = ughost(2,ew,ns) ! ghost at base

            vk_1(loc(1):loc(2)) = vvel(:,ew,ns)
            vk_1(loc(1)-1)      = vghost(1,ew,ns) ! ghost at top
            vk_1(loc(2)+1)      = vghost(2,ew,ns) ! ghost at base

            g_flag(loc(1)-1) = 1 ! ghost at top
            g_flag(loc(2)+1) = 2 ! ghost at base
        end if
    end do
  end do

end subroutine ghost_preprocess

!***********************************************************************

 subroutine ghost_preprocess_jfnk( ewn, nsn, upn, uindx, ughost, vghost, &
                              xk_1, uvel, vvel, gx_flag, pcg1)

 ! puts vel values in  xk_1 (including ghost values) and creates the
 ! ghost flag vector. xk_1 and the ghost flag vector are used for
 ! the residual calculation (jfl 20100430), adapted to combine uk, vk (kje 20101002)
   use parallel

   implicit none

   integer, intent(in) :: ewn, nsn, upn
   integer, dimension(:,:), intent(in) :: uindx
   integer, dimension(:), intent(out) :: gx_flag
   real(dp), dimension(2,ewn-1,nsn-1), intent(in) ::ughost,vghost
   real(dp), dimension(:,:,:), intent(in) :: uvel, vvel
   real(dp), dimension(:), intent(out) :: xk_1

   integer :: ew, ns, pcg1
   integer, dimension(2) :: loc
   
   gx_flag = 0

!TODO - Loops should be over locally owned velocity points?
   
   do ns = 1+staggered_shalo,size(uindx,2)-staggered_nhalo
    do ew = 1+staggered_whalo,size(uindx,1)-staggered_ehalo
         if (uindx(ew,ns) /= 0) then
             loc = getlocrange(upn, uindx(ew,ns))
             xk_1(pcg1+loc(1):pcg1+loc(2)) = uvel(:,ew,ns)
             xk_1(pcg1+loc(1)-1)      = ughost(1,ew,ns) ! ghost at top
             xk_1(pcg1+loc(2)+1)      = ughost(2,ew,ns) ! ghost at base

             xk_1(loc(1):loc(2)) = vvel(:,ew,ns)
             xk_1(loc(1)-1)      = vghost(1,ew,ns) ! ghost at top
             xk_1(loc(2)+1)      = vghost(2,ew,ns) ! ghost at base

! independent of u and v
             gx_flag(loc(1)-1) = 1 ! ghost at top
             gx_flag(loc(2)+1) = 2 ! ghost at base 
         end if
     end do
   end do
 
 end subroutine ghost_preprocess_jfnk

!***********************************************************************

subroutine ghost_postprocess( ewn, nsn, upn, uindx, uk_1, vk_1, &
                              ughost, vghost )

! puts ghost values (which are now in  uk_1 and vk_1) into ughost and 
! vghost so that they can be used fro the next time step (jfl 20100430)
  use parallel

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:), intent(in) :: uindx
  real(dp), dimension(:), intent(in) :: uk_1, vk_1
  real(dp), dimension(2,ewn-1,nsn-1), intent(out) :: ughost,vghost

  integer :: ew, ns
  integer, dimension(2) :: loc

!TODO - Loops should be over locally owned velocity points?

  do ns = 1+staggered_shalo,size(uindx,2)-staggered_nhalo
      do ew = 1+staggered_whalo,size(uindx,1)-staggered_ehalo
          if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            ughost(1,ew,ns) = uk_1(loc(1)-1) ! ghost at top
            ughost(2,ew,ns) = uk_1(loc(2)+1) ! ghost at base
            vghost(1,ew,ns) = vk_1(loc(1)-1) ! ghost at top
            vghost(2,ew,ns) = vk_1(loc(2)+1) ! ghost at base
          else 
            ughost(1,ew,ns) = 0.d0
            ughost(2,ew,ns) = 0.d0
            vghost(1,ew,ns) = 0.d0
            vghost(2,ew,ns) = 0.d0
          end if
      end do
  end do
end subroutine ghost_postprocess

!***********************************************************************

 subroutine ghost_postprocess_jfnk( ewn, nsn, upn, uindx, xk_1, &
                               ughost, vghost, pcg1 )

 ! puts ghost values (which are now in  uk_1 and vk_1) into ughost and
 ! vghost so that they can be used fro the next time step (jfl 20100430)
 ! update to use combined uk and vk = xk (kje 20101003)
   use parallel

   implicit none

   integer, intent(in) :: ewn, nsn, upn, pcg1
   integer, dimension(:,:), intent(in) :: uindx
   real(dp), dimension(:), intent(in) :: xk_1
   real(dp), dimension(2,ewn-1,nsn-1), intent(out) :: ughost,vghost
   
   integer :: ew, ns
   integer, dimension(2) :: loc

!TODO - Loops should be over locally owned velocity points?
   
   do ns = 1+staggered_shalo,size(uindx,2)-staggered_nhalo
       do ew = 1+staggered_whalo,size(uindx,1)-staggered_ehalo
           if (uindx(ew,ns) /= 0) then
             loc = getlocrange(upn, uindx(ew,ns))
             ughost(1,ew,ns) = xk_1(pcg1+loc(1)-1) ! ghost at top
             ughost(2,ew,ns) = xk_1(pcg1+loc(2)+1) ! ghost at base
             vghost(1,ew,ns) = xk_1(loc(1)-1) ! ghost at top
             vghost(2,ew,ns) = xk_1(loc(2)+1) ! ghost at base
           else 
             ughost(1,ew,ns) = 0.d0
             ughost(2,ew,ns) = 0.d0
             vghost(1,ew,ns) = 0.d0
             vghost(2,ew,ns) = 0.d0
           end if
       end do
   end do
 end subroutine ghost_postprocess_jfnk

!***********************************************************************

subroutine mindcrshstr(pt,whichresid,vel,counter,resid)

  ! Function to perform 'unstable manifold correction' (see Hindmarsh and Payne, 1996,
  ! "Time-step limits for stable solutions of the ice-sheet equation", Annals of
  ! Glaciology, 23, p.74-85)
  use parallel

  implicit none

  real(dp), intent(inout), dimension(:,:,:) :: vel
  integer, intent(in) :: counter, pt, whichresid

  real(dp), intent(out) :: resid

!TODO - critlimit is never used
  real(dp), parameter :: ssthres = 5.d0 * pi / 6.d0, &
                         critlimit = 10.d0 / (scyr * vel0), &
                         small = 1.0d-16

  real(dp), intrinsic :: abs, acos

  real(dp) :: temp_vel

  integer, dimension(2), save :: new = 1, old = 2
  !JEFF integer :: locat(3)
  integer ew, ns, nr

  integer, dimension(size(vel,1),size(vel,2),size(vel,3)) :: vel_ne_0
  real(dp) :: sum_vel_ne_0

! Note: usav and corr initialized to zero upon allocation; following probably
! not necessary, but occurs only once (per nonlinear solve)
  if (counter == 1) then
    usav(:,:,:,pt) = 0.d0
    corr(:,:,:,old(pt),pt) = 0.d0
  end if

  ! RESIDUAL CALCULATION

  select case (whichresid)
  ! options for residual calculation method, as specified in configuration file
  ! (see additional notes in "higher-order options" section of documentation)
  ! case(0): use max of abs( vel_old - vel ) / vel )
  ! case(1): use max of abs( vel_old - vel ) / vel ) but ignore basal vels
  ! case(2): use mean of abs( vel_old - vel ) / vel )
  ! case(3): use max of abs( vel_old - vel ) / vel ) (in addition to L2 norm calculated externally)

!TODO - All loops in this subroutine should be over locally owned velocity points?

   case(0)
    ! resid = maxval( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel /= 0.d0)
    resid = 0.d0
    do ns = 1 + staggered_shalo, size(vel, 3) - staggered_nhalo
      do ew = 1 + staggered_whalo, size(vel, 2) - staggered_ehalo
        do nr = 1, size(vel, 1)
          if (vel(nr,ew,ns) /= 0.d0) then
            resid = max(resid, abs(usav(nr,ew,ns,pt) - vel(nr,ew,ns)) / vel(nr,ew,ns))
          endif
        enddo
      enddo
    enddo

    resid = parallel_reduce_max(resid)
    !locat is only used in diagnostic print statement below.
    !locat = maxloc( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel /= 0.d0)

   case(1)
    ! nr = size( vel, dim=1 ) ! number of grid points in vertical ...
    ! resid = maxval( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ), MASK = vel /= 0.d0)
    resid = 0.d0
    do ns = 1 + staggered_shalo, size(vel, 3) - staggered_nhalo
      do ew = 1 + staggered_whalo, size(vel, 2) - staggered_ehalo
        do nr = 1, size(vel, 1) - 1
          if (vel(nr,ew,ns) /= 0.d0) then
            resid = max(resid, abs(usav(nr,ew,ns,pt) - vel(nr,ew,ns)) / vel(nr,ew,ns))
          endif
        enddo
      enddo
    enddo

    resid = parallel_reduce_max(resid)
    !locat = maxloc( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
    !        MASK = vel /= 0.d0)

   case(2)
    call not_parallel(__FILE__, __LINE__)
    !JEFF This has not been translated to parallel.
    resid = 0.d0
    nr = size( vel, dim=1 )
    vel_ne_0 = 0
    where ( vel /= 0.d0 ) vel_ne_0 = 1

    ! include basal velocities in resid. calculation when using MEAN
    ! JEFF Compute sums across nodes in order to compute mean.
    resid = sum( abs((usav(:,:,:,pt) - vel ) / vel ), &
            MASK = vel /= 0.d0)

    resid = parallel_reduce_sum(resid)
    sum_vel_ne_0 = sum( vel_ne_0 )
    sum_vel_ne_0 = parallel_reduce_sum(sum_vel_ne_0)

    resid = resid / sum_vel_ne_0

    ! ignore basal velocities in resid. calculation when using MEAN
    ! resid = sum( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),   &
    !           MASK = vel /= 0.d0) / sum( vel_ne_0(1:nr-1,:,:) )

    ! NOTE that the location of the max residual is somewhat irrelevent here
    !      since we are using the mean resid for convergence testing
    ! locat = maxloc( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel /= 0.d0)

   case(3)
    ! resid = maxval( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel /= 0.d0)
    resid = 0.d0
    do ns = 1 + staggered_shalo, size(vel, 3) - staggered_nhalo
      do ew = 1 + staggered_whalo, size(vel, 2) - staggered_ehalo
        do nr = 1, size(vel, 1)
          if (vel(nr,ew,ns) /= 0.d0) then
            resid = max(resid, abs(usav(nr,ew,ns,pt) - vel(nr,ew,ns)) / vel(nr,ew,ns))
          endif
        enddo
      enddo
    enddo

    resid = parallel_reduce_max(resid)
    !locat = maxloc( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel /= 0.d0)

  end select

  ! Additional debugging line, useful when trying to determine if convergence is being consistently
  ! help up by the residual at one or a few particular locations in the domain.
  !print '("* ",i3,g20.6,3i6,g20.6)', counter, resid, locat, vel(locat(1),locat(2),locat(3))*vel0

  ! SAVE VELOCITY AND CALCULATE CORRECTION

  corr(:,:,:,new(pt),pt) = vel(:,:,:) - usav(:,:,:,pt)  ! changed

!  if (counter > 1) then
!    where (acos((corr(:,:,:,new(pt),pt) * corr(:,:,:,old(pt),pt)) / &
!          (abs(corr(:,:,:,new(pt),pt)) * abs(corr(:,:,:,old(pt),pt)) + small)) > &
!           ssthres .and. corr(:,:,:,new(pt),pt) - corr(:,:,:,old(pt),pt) /= 0.d0 )
!      mindcrshstr = usav(:,:,:,pt) + &
!                    corr(:,:,:,new(pt),pt) * abs(corr(:,:,:,old(pt),pt)) / &
!                    abs(corr(:,:,:,new(pt),pt) - corr(:,:,:,old(pt),pt))
!!      mindcrshstr = vel; ! jfl uncomment this and comment out line above
!!                         ! to avoid the unstable manifold correction
!    elsewhere
!      mindcrshstr = vel;
!    end where
!  else
!    mindcrshstr = vel;
!  end if
!  usav(:,:,:,pt) = vel
!  vel = mindcrshstr

  if (counter > 1) then

    ! Replace where clause with explicit, owned variables for each processor.
    do ns = 1 + staggered_shalo, size(vel, 3) - staggered_nhalo
      do ew = 1 + staggered_whalo, size(vel, 2) - staggered_ehalo
        do nr = 1, size(vel, 1)
          temp_vel = vel(nr,ew,ns)
          if (acos((corr(nr,ew,ns,new(pt),pt) * corr(nr,ew,ns,old(pt),pt)) / &
               (abs(corr(nr,ew,ns,new(pt),pt)) * abs(corr(nr,ew,ns,old(pt),pt)) + small)) > &
              ssthres .and. corr(nr,ew,ns,new(pt),pt) - corr(nr,ew,ns,old(pt),pt) /= 0.d0) then
            vel(nr,ew,ns) = usav(nr,ew,ns,pt) + &
                corr(nr,ew,ns,new(pt),pt) * abs(corr(nr,ew,ns,old(pt),pt)) / &
                    abs(corr(nr,ew,ns,new(pt),pt) - corr(nr,ew,ns,old(pt),pt))
          endif
          usav(nr,ew,ns,pt) = temp_vel
        enddo
      enddo
    enddo
  else

    usav(:,:,:,pt) = vel

  end if

  ! UPDATE POINTERS

  !*sfp* Old version
  ! if (new(pt) == 1) then; old(pt) = 1; new(pt) = 2; else; old(pt) = 1; new(pt) = 2; end if

  !*sfp* correction from Carl Gladdish
  if (new(pt) == 1) then; old(pt) = 1; new(pt) = 2; else; old(pt) = 2; new(pt) = 1; end if

  return

end subroutine mindcrshstr

!***********************************************************************

!TODO - Why are there two of these subroutines?  Can we remove one of them?

function mindcrshstr2(pt,whichresid,vel,counter,resid)

  ! Function to perform 'unstable manifold correction' (see Hindmarsch and Payne, 1996,
  ! "Time-step limits for stable solutions of the ice-sheet equation", Annals of  
  ! Glaciology, 23, p.74-85)
  
  ! Alternate unstable manifold scheme, based on DeSmedt, Pattyn, and De Goen, J. Glaciology 2010
  ! Written by Carl Gladdish
  
!TODO - something to do here?
  use parallel  ! Use of WHERE statements is causing inconsistencies on the halos in parallel.  Rewrite like mindcrshstr()
  implicit none
  
  real(dp), intent(in), dimension(:,:,:) :: vel
  integer, intent(in) :: counter, pt, whichresid 
  real(dp), intent(out) :: resid
  
  real(dp), dimension(size(vel,1),size(vel,2),size(vel,3)) :: mindcrshstr2
  
  integer, parameter :: start_umc = 3
  real(dp), parameter :: cvg_accel = 2.d0
  real(dp), parameter :: small = 1.0d-16
  
  real(dp) in_prod, len_new, len_old, mean_rel_diff, sig_rel_diff
  real(dp) :: theta 
  real(dp), intrinsic :: abs, acos
  
  integer, dimension(2), save :: new = 1, old = 2
  integer :: locat(3)
  
  integer :: nr
  integer,      dimension(size(vel,1),size(vel,2),size(vel,3)) :: vel_ne_0
  real(dp),dimension(size(vel,1),size(vel,2),size(vel,3)) :: rel_diff
  
  call not_parallel(__FILE__, __LINE__)

  if (counter == 1) then
    usav(:,:,:,pt) = 0.d0
    corr(:,:,:,:,:) = 0.d0
  end if
  
  corr(:,:,:,new(pt),pt) = vel - usav(:,:,:,pt)
  
  if (counter >= start_umc) then
  
  in_prod = sum( corr(:,:,:,new(pt),pt) * corr(:,:,:,old(pt),pt) )
  len_new = sqrt(sum( corr(:,:,:,new(pt),pt) * corr(:,:,:,new(pt),pt) ))
  len_old = sqrt(sum( corr(:,:,:,old(pt),pt) * corr(:,:,:,old(pt),pt) ))
  
  theta = acos( in_prod / (len_new * len_old + small) )
    
   if (theta  < (1.0/8.0)*pi) then
        mindcrshstr2 = usav(:,:,:,pt) + cvg_accel * corr(:,:,:,new(pt),pt)
!        print *, theta/pi, 'increased correction'
   else if(theta < (19.0/20.0)*pi) then
        mindcrshstr2 = vel
!        print *, theta/pi, 'standard correction'
   else
        mindcrshstr2 = usav(:,:,:,pt) + (1.0/cvg_accel) * corr(:,:,:,new(pt),pt)
!        print *, theta/pi, 'decreasing correction'
   end if
   
  else
  
    mindcrshstr2 = vel;
 !   print *, 'Not attempting adjustment to correction'
     
  end if
  

  ! now swap slots for storing the previous correction
  if (new(pt) == 1) then
      old(pt) = 1; new(pt) = 2
  else
      old(pt) = 2; new(pt) = 1
  end if

  if (counter == 1) then
        usav_avg = 1.d0
  else
        usav_avg(1) = sum( abs(usav(:,:,:,1)) ) / size(vel)  ! a x-dir transport velocity scale
        usav_avg(2) = sum( abs(usav(:,:,:,2)) ) / size(vel)  ! a y-dir transport velocity scale
  end if

!  print *, 'usav_avg(1)',usav_avg(1),'usav_avg(2)',usav_avg(2)

  select case (whichresid)

  ! options for residual calculation method, as specified in configuration file
  ! (see additional notes in "higher-order options" section of documentation)
  ! case(0): use max of abs( vel_old - vel ) / vel )
  ! case(1): use max of abs( vel_old - vel ) / vel ) but ignore basal vels
  ! case(2): use mean of abs( vel_old - vel ) / vel )

   case(0)
    rel_diff = 0.d0
    vel_ne_0 = 0
    where ( mindcrshstr2 /= 0.d0 )
        vel_ne_0 = 1
        rel_diff = abs((usav(:,:,:,pt) - mindcrshstr2) / mindcrshstr2) &
                           * usav_avg(pt)/sqrt(sum(usav_avg ** 2.0))
    end where

    resid = maxval( rel_diff, MASK = mindcrshstr2 /= 0.d0 )
    locat = maxloc( rel_diff, MASK = mindcrshstr2 /= 0.d0 )

!    mean_rel_diff = sum(rel_diff) / sum(vel_ne_0)
!    sig_rel_diff = sqrt( sum((rel_diff - mean_rel_diff) ** 2.0 )/ sum(vel_ne_0) )
!    print *, 'mean', mean_rel_diff, 'sig', sig_rel_diff

    !write(*,*) 'locat', locat
    !call write_xls('resid1.txt',abs((usav(1,:,:,pt) - mindcrshstr2(1,:,:)) / (mindcrshstr2(1,:,:) + 1e-20)))

   case(1)
    !**cvg*** should replace vel by mindcrshstr2 in the following lines, I belive
    nr = size( vel, dim=1 ) ! number of grid points in vertical ...
    resid = maxval( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
                        MASK = vel /= 0.d0)
    locat = maxloc( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
            MASK = vel /= 0.d0)

   case(2)
    !**cvg*** should replace vel by mindcrshstr2 in the following lines, I believe
    nr = size( vel, dim=1 )
    vel_ne_0 = 0
    where ( vel /= 0.d0 ) vel_ne_0 = 1

    ! include basal velocities in resid. calculation when using MEAN
    resid = sum( abs((usav(:,:,:,pt) - vel ) / vel ), &
            MASK = vel /= 0.d0) / sum( vel_ne_0 )

    ! ignore basal velocities in resid. calculation when using MEAN
    ! resid = sum( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),   &
    !           MASK = vel /= 0.d0) / sum( vel_ne_0(1:nr-1,:,:) )

    ! NOTE that the location of the max residual is somewhat irrelevent here
    !      since we are using the mean resid for convergence testing
    locat = maxloc( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel /= 0.d0)

  end select

  usav(:,:,:,pt) = mindcrshstr2

    ! Additional debugging line, useful when trying to determine if convergence is being consistently
    ! held up by the residual at one or a few particular locations in the domain.
!    print '("* ",i3,g20.6,3i6,g20.6)', counter, resid, locat, vel(locat(1),locat(2),locat(3))*vel0

  return

end function mindcrshstr2

!***********************************************************************

subroutine findcoefstr(ewn,  nsn,   upn,            &
                       dew,  dns,   sigma,          &
                       pt,          efvs,           &
                       thisvel,     othervel,       &
                       thck,        thisdusrfdx,    &
                       dusrfdew,    dthckdew,       &
                       d2usrfdew2,  d2thckdew2,     &
                       dusrfdns,    dthckdns,       &
                       d2usrfdns2,  d2thckdns2,     &
                       d2usrfdewdns,d2thckdewdns,   &
                       dlsrfdew,    dlsrfdns,       &
                       stagthck,    whichbabc,      &
                       uindx,       mask,           &
                       lsrf,        topg,           &
                       minTauf,     flwa,           &
                       beta,        btraction,      &
                       count, assembly )

  ! Main subroutine for determining coefficients that go into the LHS matrix A 
  ! in the expression Au = b. Calls numerous other subroutines, including boundary
  ! condition subroutines, which determine "b".

  use parallel

  implicit none

  integer, intent(in) :: ewn, nsn, upn, count, assembly
  real(dp), intent(in) :: dew, dns
  real(dp), dimension(:), intent(in) :: sigma

  real(dp), dimension(:,:,:), intent(in) :: efvs, thisvel, &
                                                    othervel
  real(dp), dimension(:,:), intent(in) :: stagthck, thisdusrfdx,     &
                                                  dusrfdew,   dthckdew,      &
                                                  d2usrfdew2, d2thckdew2,    &
                                                  dusrfdns,   dthckdns,      &
                                                  d2usrfdns2, d2thckdns2,    &
                                                  d2usrfdewdns,d2thckdewdns, &
                                                  dlsrfdew,   dlsrfdns,      &
                                                  thck, lsrf, topg

  real(dp), dimension(:,:), intent(in) :: minTauf
  real(dp), dimension(:,:), intent(in) :: beta
  real(dp), dimension(:,:,:), intent(in) :: flwa
  real(dp), dimension(:,:,:), intent(inout) :: btraction

  integer, dimension(:,:), intent(in) :: mask, uindx
  integer, intent(in) :: pt, whichbabc

  real(dp), dimension(ewn-1,nsn-1) :: betasquared
  real(dp), dimension(2,2,2) :: localefvs
  real(dp), dimension(3,3,3) :: localothervel
  real(dp), dimension(upn) :: boundaryvel
  real(dp) :: flwabar

  integer, dimension(6,2) :: loc2
  integer, dimension(2) :: loc2plusup
  integer, dimension(3) :: shift
  integer :: ew, ns, up, up_start

  ct = 1        ! index to count the number of non-zero entries in the sparse matrix
  ct2 = 1

  if( assembly == 1 )then   ! for normal assembly (assembly=0), start vert index at sfc and go to bed
      up_start = upn        ! for boundary traction calc (assembly=1), do matrix assembly on for equations at bed
  else
      up_start = 1
  end if

  ! calculate/specify the map of 'betasquared', for use in the basal boundary condition. 
  ! Input to the subroutine 'bodyset' (below) ... 

  call calcbetasquared (whichbabc,              &
                        dew,        dns,        &
                        ewn,        nsn,        &
                        lsrf,       topg,       &
                        thck,                   &
                        thisvel(upn,:,:),       &
                        othervel(upn,:,:),      &
                        minTauf, beta,          &
                        betasquared )
! intent(out) betasquared

  ! Note loc2_array is defined only for non-halo ice grid points.
  ! JEFFLOC returns an array with starting indices into solution vector for each ice grid point.
  allocate(loc2_array(ewn,nsn,2))
  loc2_array = getlocationarray(ewn, nsn, upn, mask, uindx)
!  !!!!!!!!! useful for debugging !!!!!!!!!!!!!!
!    print *, 'loc2_array = '
!    print *, loc2_array
!    pause

!HALO - This loop should be over locally owned velocity points: (ilo-1:ihi,jlo-1:jhi)
!       Note: efvs has been computed in a layer of halo cells, so we have its value in all
!             neighbors of locally owned velocity points.

  ! JEFFLOC Do I need to restrict to non-halo grid points?
  do ns = 1+staggered_shalo,size(mask,2)-staggered_nhalo
    do ew = 1+staggered_whalo,size(mask,1)-staggered_ehalo

     ! Calculate the depth-averaged value of the rate factor, needed below when applying an ice shelf
     ! boundary condition (complicated code so as not to include funny values at boundaries ...
     ! ... kind of a mess and could be redone or made into a function or subroutine).
     ! SUM has the definition SUM(ARRAY, DIM, MASK) where MASK is either scalar or the same shape as ARRAY
     ! JEFFLOC Concerned about the edges at (ew+1, ns), (ew, ns+1), and (ew+1,ns+1)

!SCALING - The following is OK because flwa*vis0 is equal to the dimensional flow factor.
!          The product will still equal the dimensional flow factor when vis0 = 1.

     flwabar = ( sum( flwa(:,ew,ns),     1, flwa(1,ew,ns)*vis0     < 1.0d-10 )/real(upn) + &
                 sum( flwa(:,ew,ns+1),   1, flwa(1,ew,ns+1)*vis0   < 1.0d-10 )/real(upn)  + &
                 sum( flwa(:,ew+1,ns),   1, flwa(1,ew+1,ns)*vis0   < 1.0d-10 )/real(upn)  + &
                 sum( flwa(:,ew+1,ns+1), 1, flwa(1,ew+1,ns+1)*vis0 < 1.0d-10 )/real(upn) ) / &
               ( sum( flwa(:,ew,ns)/flwa(:,ew,ns),         1, flwa(1,ew,ns)*vis0     < 1.0d-10 )/real(upn) + &
                 sum( flwa(:,ew,ns+1)/flwa(:,ew,ns+1),     1, flwa(1,ew,ns+1)*vis0   < 1.0d-10 )/real(upn) + &
                 sum( flwa(:,ew+1,ns)/flwa(:,ew+1,ns),     1, flwa(1,ew+1,ns)*vis0   < 1.0d-10 )/real(upn) + &
                 sum( flwa(:,ew+1,ns+1)/flwa(:,ew+1,ns+1), 1, flwa(1,ew+1,ns+1)*vis0 < 1.0d-10 )/real(upn) )

    loc2(1,:) = loc2_array(ew,ns,:)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!TODO - Sometimes we may want to solve for the velocity at the domain boundary,
!       or at a land margin.  Can we allow this?
!       
    if ( GLIDE_HAS_ICE(mask(ew,ns)) .and. .not. &
         GLIDE_IS_COMP_DOMAIN_BND(mask(ew,ns)) .and. .not. &
         GLIDE_IS_MARGIN(mask(ew,ns)) .and. .not. &
         GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) .and. .not. &
         GLIDE_IS_CALVING(mask(ew,ns) ) .and. .not. &
         GLIDE_IS_THIN(mask(ew,ns) ) ) &
    then
!    print *, 'In main body ... ew, ns = ', ew, ns
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        call calccoeffs( upn,               sigma,              &
                         stagthck(ew,ns),                       &
                         dusrfdew(ew,ns),   dusrfdns(ew,ns),    &
                         dthckdew(ew,ns),   dthckdns(ew,ns),    &
                         d2usrfdew2(ew,ns), d2usrfdns2(ew,ns),  &
                         d2usrfdewdns(ew,ns),                   &
                         d2thckdew2(ew,ns), d2thckdns2(ew,ns),  &
                         d2thckdewdns(ew,ns))

        ! get index of cardinal neighbours
        loc2(2,:) = loc2_array(ew+1,ns,:)
        loc2(3,:) = loc2_array(ew-1,ns,:)
        loc2(4,:) = loc2_array(ew,ns+1,:)
        loc2(5,:) = loc2_array(ew,ns-1,:)

        ! this loop fills coeff. for all vertical layers at index ew,ns (including sfc. and bed bcs)
        do up = up_start, upn

            ! Function to adjust indices at sfc and bed so that most correct values of 'efvs' and 'othervel'
            ! are passed to function. Because of the fact that efvs goes from 1:upn-1 rather than 1:upn
            ! we simply use the closest values. This could probably be improved upon at some point
            ! by extrapolating values for efvs at the sfc and bed using one-sided diffs, and it is not clear
            ! how important this simplfication is.
            !JEFFLOC indshift() returns three-element shift index for up, ew, and ns respectively.
            !JEFFLOC It does get passed loc2_array, but it doesn't use it.  Further, the shifts can be at most 1 unit in any direction.
            shift = indshift( 0, ew, ns, up, ewn, nsn, upn, loc2_array(:,:,1), stagthck(ew-1:ew+1,ns-1:ns+1) )

!HALO - Note that ew and ns below are locally owned velocity points, i.e. in the range (ilo-1:ihi, jlo-1:jhi)
!HALO - This means we need efvs in one layer of halo cells.

!JEFFLOC As long as not accessing halo ice points, then won't shift off of halo of size at least 1.
!JEFFLOC Completed scan on 11/23.  Testing change of definition of loc2_array.

            call bodyset(ew,  ns,  up,        &
                         ewn, nsn, upn,       &
                         dew,      dns,       &
                         pt,       loc2_array,&
                         loc2,     stagthck,  &
                         thisdusrfdx,         &
                         dusrfdew, dusrfdns,  &
                         dlsrfdew, dlsrfdns,  &
                         efvs(up-1+shift(1):up+shift(1),ew:ew+1,ns:ns+1),  &
                         othervel(up-1+shift(1):up+1+shift(1),  &
                         ew-1+shift(2):ew+1+shift(2),  &
                         ns-1+shift(3):ns+1+shift(3)), &
                         thisvel(up-1+shift(1):up+1+shift(1),  &
                         ew-1+shift(2):ew+1+shift(2),  &
                         ns-1+shift(3):ns+1+shift(3)), &
                         betasquared(ew,ns),           &
                         btraction,                    &
                         whichbabc, assembly )

        end do  ! upn

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!TODO - Not sure COMP_DOMAIN_BND condition is needed

    elseif ( GLIDE_IS_CALVING( mask(ew,ns) ) .and. .not. &
             GLIDE_IS_COMP_DOMAIN_BND(mask(ew,ns) ) .and. .not. &
             GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) .and. .not. &
             GLIDE_IS_THIN(mask(ew,ns) ) ) &
    then
!    print *, 'At a SHELF boundary ... ew, ns = ', ew, ns
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        call calccoeffs( upn,               sigma,              &
                         stagthck(ew,ns),                       &
                         dusrfdew(ew,ns),   dusrfdns(ew,ns),    &
                         dthckdew(ew,ns),   dthckdns(ew,ns),    &
                         d2usrfdew2(ew,ns), d2usrfdns2(ew,ns),  &
                         d2usrfdewdns(ew,ns),                   &
                         d2thckdew2(ew,ns), d2thckdns2(ew,ns),  &
                         d2thckdewdns(ew,ns))

        do up = up_start, upn
           lateralboundry = .true.
           shift = indshift(  1, ew, ns, up,                   &
                               ewn, nsn, upn,                  &
                               loc2_array(:,:,1),              &
                               stagthck(ew-1:ew+1,ns-1:ns+1) )

            call bodyset(ew,  ns,  up,        &
                         ewn, nsn, upn,       &
                         dew,      dns,       &
                         pt,       loc2_array,&
                         loc2,     stagthck,  &
                         thisdusrfdx,         &
                         dusrfdew, dusrfdns,  &
                         dlsrfdew, dlsrfdns,  &
                         efvs(up-1+shift(1):up+shift(1),ew:ew+1,ns:ns+1),  &
                         othervel(up-1+shift(1):up+1+shift(1),  &
                         ew-1+shift(2):ew+1+shift(2),  &
                         ns-1+shift(3):ns+1+shift(3)), &
                         thisvel(up-1+shift(1):up+1+shift(1),  &
                         ew-1+shift(2):ew+1+shift(2),  &
                         ns-1+shift(3):ns+1+shift(3)), &
                         betasquared(ew,ns),           &
                         btraction,                    &
                         whichbabc, assembly,          &              
                         abar=flwabar, cc=count )
        end do
        lateralboundry = .false.

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!TODO - Here we deal with cells on the computational domain boundary.
!       Currently the velocity is always set to a specified value on this boundary.
!       With open (non-Dirichlet) BCs, we might want to solve for these velocities,
!        using the code above to compute the matrix elements.

    elseif ( GLIDE_HAS_ICE(mask(ew,ns)) .and. ( GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) .or. &
             GLIDE_IS_COMP_DOMAIN_BND(mask(ew,ns)) ) .or. GLIDE_IS_LAND_MARGIN(mask(ew,ns)) .or. &
             GLIDE_IS_THIN(mask(ew,ns)) ) &
    then
!    print *, 'At a NON-SHELF boundary ... ew, ns = ', ew, ns

!WHLdebug
!    print*, ' '
!    print*, 'At a NON-SHELF boundary ... ew, ns = ', ew, ns
!    print*, 'LAND_MARGIN =', GLIDE_IS_LAND_MARGIN(mask(ew,ns))

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ! Put specified value for vel on rhs. NOTE that this is NOT zero by default 
        ! unless the initial guess is zero. It will be set to whatever the initial value 
        ! for the vel at location up,ew,ns is in the initial array!
        loc2plusup = loc2(1,:)
        call valueset(0.d0, loc2plusup)

        loc2plusup = loc2(1,:) + upn + 1
        call valueset(0.d0, loc2plusup)

        do up = up_start, upn
           loc2plusup = loc2(1,:) + up
           call valueset( thisvel(up,ew,ns), loc2plusup )     ! vel at margin set to initial value 
!           call valueset( 0.d0 )                ! vel at margin set to 0 
        end do

    end if

    end do;     ! ew 
  end do        ! ns

  deallocate(loc2_array)

  return

end subroutine findcoefstr

!***********************************************************************

subroutine bodyset(ew,  ns,  up,           &
                   ewn, nsn, upn,          &
                   dew,      dns,          &
                   pt,       loc2_array,   &
                   loc2,     stagthck,     &
                   thisdusrfdx,            &
                   dusrfdew, dusrfdns,     &
                   dlsrfdew, dlsrfdns,     &
                   local_efvs,             &
                   local_othervel,         &
                   local_thisvel,          &
                   betasquared,            &
                   btraction,              &
                   whichbabc, assembly,    &
                   abar, cc)

  ! This subroutine does the bulk of the work in calling the appropriate discretiztion routines,
  ! which determine the values for coefficients that will go into the sparse matrix, for points
  ! on and inside of the boundaries.

  use glimmer_paramets, only: evs0, evs_scale
  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, intent(in) :: ew, ns, up
  real(dp), intent(in) :: dew, dns
  integer, intent(in) :: pt, whichbabc, assembly
  integer, dimension(ewn,nsn,2), intent(in) :: loc2_array
  integer, dimension(6,2), intent(in) :: loc2

  real(dp), dimension(:,:), intent(in) :: stagthck
  real(dp), dimension(:,:), intent(in) :: dusrfdew, dusrfdns
  real(dp), dimension(:,:), intent(in) :: dlsrfdew, dlsrfdns
  real(dp), dimension(:,:), intent(in) :: thisdusrfdx
  real(dp), dimension(2,2,2), intent(in) :: local_efvs
  ! "local_othervel" is the other vel component (i.e. u when v is being calc and vice versa),
  ! which is taken as a known value (terms involving it are moved to the RHS and treated as sources)
  real(dp), dimension(3,3,3), intent(in) :: local_othervel, local_thisvel
  real(dp), intent(in) :: betasquared
  real(dp), dimension(:,:,:), intent(inout) :: btraction
  real(dp), intent(in), optional :: abar
  integer, intent(in), optional :: cc

  ! storage space for coefficients that go w/ the discretization at the local point up, ew, ns.
  ! Note that terms other than 'g' are used for storing particular parts needed for calculation
  ! of the basal traction vector.
  real(dp), dimension(3,3,3) :: g, h, g_cros, g_vert, g_norm, g_vel_lhs, g_vel_rhs

  ! source term for the rhs when using ice shelf lateral boundary condition,
  ! e.g. source = rho*g*H/(2*Neff) * ( 1 - rho_i / rho_w ) for ice shelf
  real(dp) :: source

  real(dp) :: slopex, slopey    ! local sfc (or bed) slope terms

  ! lateral boundary normal and vector to indicate use of forward
  ! or bacward one-sided diff. when including specified stress lateral bcs
  real(dp), dimension(2) :: fwdorbwd, normal

  real(dp) :: nz   ! z dir normal vector component at sfc or bed (takes diff value for each)

  integer, dimension(2) :: bcflag  ! indicates choice of sfc and basal bcs ...

  real(dp) :: scalebabc

  integer, dimension(2) :: loc2plusup

  loc2plusup = loc2(1,:) + up

  if( lateralboundry )then

  ! *********************************************************************************************
  ! lateral boundary conditions 

  ! if at sfc or bed, source due to seawater pressure is 0 and bc normal vector
  ! should contain sfc/bed slope components, e.g. (-ds/dx, -ds/dy, 1) or (db/dx, db/dy, -1)
     source = 0.d0

     call getlatboundinfo( ew,  ns,  up,                                 &
                           ewn, nsn, upn,                                &
                           stagthck(ew-1:ew+1, ns-1:ns+1),               &
                           loc2_array(:,:,1), fwdorbwd, normal, loc_latbc)

     if( up == 1 .or. up == upn )then

        if( up == 1 )then                ! specify necessary variables and flags for free sfc
           bcflag = (/1,0/)
           loc2plusup = loc2(1,:) + up - 1   ! reverse the sparse matrix / rhs vector row index by 1 ...
           slopex = -dusrfdew(ew,ns); slopey = -dusrfdns(ew,ns); nz = 1.d0
        else                             ! specify necessary variables and flags for basal bc
   
           if( whichbabc == 6 )then
                bcflag = (/0,0/)             ! flag for u=v=0 at bed; doesn't work well so commented out here...
                                             ! better to specify very large value for betasquared below
           elseif( whichbabc >=0 .and. whichbabc <= 5 )then
                bcflag = (/1,1/)              ! flag for specififed stress at bed: Tau_zx = betasquared * u_bed,
                                              ! where betasquared is MacAyeal-type traction parameter
           elseif( whichbabc == 7 )then

                write(*,*)"ERROR: This option is not supported in the current release of CISM."
                write(*,*)"       A future release will support use of Newton iteration on " 
                write(*,*)"       plastic till basal BC."
                stop

                bcflag = (/1,2/)              ! plastic bed iteration using Newton implementation
 
           end if   

           loc2plusup = loc2(1,:) + up + 1   ! advance the sparse matrix / rhs row vector index by 1 ...
           slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns); nz = -1.d0

        end if

!        !! Hack to avoid bad sfc and basal bc normal vectors !!
        slopex = 0.d0; slopey = 0.d0

        g = normhorizmainbc_lat(dew,           dns,             &
                                slopex,        slopey,          &
                                dsigmadew(up), dsigmadns(up),   &
                                pt,            2,               &
                                dup(up),       local_efvs,      &
                                oneorfour,     fourorone,       &
                                onesideddiff,                   &
                                normal,        fwdorbwd)

        ! add on coeff. associated w/ du/dsigma  
        g(:,3,3) = g(:,3,3) &
                 + vertimainbc( stagthck(ew,ns), bcflag, dup(up),     &
                                local_efvs,      betasquared,   g_vert,    nz, &
                                plastic_coeff=plastic_coeff_lhs(pt,ew,ns)  )

        !! scale basal bc coeffs when using JFNK solver 
        scalebabc = scalebasalbc( g, bcflag, lateralboundry, betasquared, local_efvs )
        g = g / scalebabc

        ! put the coeff. for the b.c. equation in the same place as the prev. equation
        ! (w.r.t. cols), on a new row ...
        call fillsprsebndy( g, loc2plusup(1), loc_latbc, up, normal, pt )

        ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
        ! which results from moving them from the LHS over to the RHS, has been moved
        ! inside of "croshorizmainbc_lat".
        ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
        ! which results from moving them from the LHS over to the RHS, has been moved
        ! inside of "croshorizmainbc_lat".
        rhsd(loc2plusup(2)) = sum( croshorizmainbc_lat(dew,           dns,           &
                                                   slopex,        slopey,        &
                                                   dsigmadew(up), dsigmadns(up), &
                                                   pt,            2,             &
                                                   dup(up),       local_othervel,&
                                                   local_efvs,                   &
                                                   oneortwo,      twoorone,      &
                                                   onesideddiff,                 &
                                                   normal,fwdorbwd)              &
                                                 * local_othervel ) /scalebabc
!OSBC: note that efvs is now passed in to subroutine above, since visc. is included in LHS matrix
! coeffs rather than divide through onto the RHS solution vector.

!TODO: These are not currently being used. They are for getting the off diag. blocks of coeffs. from the
!LHS matrix to improve preconditioning. However, we never got them working. Consider removing here and in
! related locations below? Can grep for this by looking for the related if construct below.
!        if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag )then
!            storeoffdiag = .true.
!            h = -croshorizmainbc_lat(dew,           dns,           & 
!                                slopex,        slopey,        &
!                                dsigmadew(up), dsigmadns(up), & 
!                                pt,            2,             & 
!                                dup(up),       local_othervel,& 
!                                oneortwo,      twoorone,      & 
!                                onesideddiff,                 &
!                                 normal,fwdorbwd) / scalebabc 
!            call fillsprsebndy( h, loc2plusup(1), loc_latbc, up, normal, pt ) 
!            storeoffdiag = .false.
!        end if     

    end if     ! up = 1 or up = upn (IF at lateral boundary and IF at surface or bed)

    ! If in main body and at ice/ocean boundary, calculate depth-averaged stress
    ! due to sea water, bc normal vector components should be boundary normal 
    loc2plusup = loc2(1,:) + up

    ! for this bc, the normal vector components are not the sfc/bed slopes but are taken
    ! from a normal to the shelf front in map view (x,y plane); slopex,slopey are simply renamed here 
    slopex = normal(1)
    slopey = normal(2)

    ! Two options here, (1) use the 1d solution that involves the rate factor (not accurate for 
    !                       3d domains, but generally more stable 
    !                   (2) use the more general solution that involves the eff. visc. and normal
    !                       vector orientation at lateral boundary
    ! ... only one of these options should be active at a time (comment the other lines out)
    ! The default setting is (2), the more general case that should also work in the 1d case.

    ! In some cases, the two options can be used together to improve performance, e.g. for the Ross
    ! ice shelf experiment, a number of early iterations use the more simple bc (option 1) and then
    ! when the solution has converged a bit, we switch to the more realistic implementation (option 2).
    ! That is achieved in the following if construct ...

!OSBC: NOTE that for the one-side bc implementation, the eff. visc. remains w/ the RHS matrix coeffs, so it
! is not necessary here to fiddle with how it is calculated for the source term in the RHS vector. For this
! reason, we instead go straight to the floating ice bc source term calc for the full 2d code (below).

!    if( cc < 2 .and. .not. inisoln )then  ! This should be the default option for the shelf BC source term.
                                          ! If no previous guess for the eff. visc. exists, this option
                                          ! uses the 1d version of the BC for one iteration in order to 
                                          ! "precondition" the soln for the next iteration. W/o this option
                                          ! active, the 2d version of the BC fails, presumably because of
                                          ! eff. visc. terms in the denom. of the source term which are either
                                          ! too large (inf) or to small (~0).

! These options are primarily for debugging the shelf BC source term
!    if( cc >= 0 )then        ! - use this to use only the 1d version
    if( cc > 1000000 )then   ! - use this to go straight to the full 2d version of the bc

    ! --------------------------------------------------------------------------------------
    ! (1) source term (strain rate at shelf/ocean boundary) from Weertman's analytical solution 
    ! --------------------------------------------------------------------------------------
    ! See eq. 2, Pattyn+, 2006, JGR v.111; eq. 8, Vieli&Payne, 2005, JGR v.110). Note that this 
    ! contains the 1d assumption that ice is not spreading lateraly !(assumes dv/dy = 0 for u along flow)

    source = abar * vis0 * ( 1.d0/4.d0 * rhoi * grav * stagthck(ew,ns)*thk0 * ( 1.d0 - rhoi/rhoo))**3.d0

    ! multiply by 4 so that case where v=0, du/dy = 0, LHS gives: du/dx = du/dx|_shelf 
    ! (i.e. LHS = 4*du/dx, requires 4*du/dx_shelf)
    source = source * 4.d0

    ! split source based on the boundary normal orientation and non-dimensinoalize
    ! Note that it is not really appropriate to apply option (1) to 2d flow, since terms other than du/dx in 
    ! eff. strain rate are ignored. For 2d flow, should use option (2) below. 
     source = source * normal(pt)
     source = source * tim0 ! make source term non-dim
    ! --------------------------------------------------------------------------------------

  else

    ! --------------------------------------------------------------------------------------
    ! (2) source term (strain rate at shelf/ocean boundary) from MacAyeal depth-ave solution. 
    ! --------------------------------------------------------------------------------------

    source = (rhoi*grav*stagthck(ew,ns)*thk0) / tau0 / 2.d0 * ( 1.d0 - rhoi / rhoo )

    ! terms after "/" below count number of non-zero efvs cells ... needed for averaging of the efvs at boundary 

!SCALING - Units of efvs are evs0 = tau0 / (vel0/len0)
!          Multiply efvs by evs0/evs_scale so we get the same result in these two cases:
!           (1) Old Glimmer with scaling:         evs0 = evs_scale = tau0/(vel0/len0), and efvs is non-dimensional
!           (2) New Glimmer-CISM without scaling: evs0 = 1, evs_scale = tau0/(vel0/len0), and efvs is in Pa s

!!!    source = source / ( sum(local_efvs, local_efvs > 1.0d-12) / &    ! OLD version
!!!             sum( (local_efvs/local_efvs), local_efvs > 1.0d-12 ) )

!    source = source / ( sum(local_efvs, local_efvs*evs0/evs_scale > 1.0d-12) / &     ! NEW version
!             sum( (local_efvs/local_efvs), local_efvs*evs0/evs_scale > 1.0d-12 ) )
!OSBC: This scaling is not needed in the one-sided bc version of the code, since the eff. visc. are
! kept w/ the LHS coeff. and not moved to the denom. of the RHS vector source terms.

    source = source * normal(pt) ! partition according to normal vector at lateral boundary
                                 ! NOTE that source term is already non-dim here 
    ! --------------------------------------------------------------------------------------

  end if

!OSBC: For the one-sided diff version of the code, need to pass in the efvs. since it is
! part of the LHS matrix coeffs. now.
    g = normhorizmainbc_lat(dew,           dns,        &
                            slopex,        slopey,     &
                            dsigmadew(up), dsigmadns(up),  &
                            pt,            1,          &
                            dup(up),       local_efvs, &
                            oneorfour,     fourorone,  &
                            onesideddiff,              &
                            normal,        fwdorbwd)

    ! put the coeff. for the b.c. equation in the same place as the prev. equation
    ! (w.r.t. cols), on a new row ...
    call fillsprsebndy( g, loc2plusup(1), loc_latbc, up, normal, pt )

    ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
    ! which results from moving them from the LHS over to the RHS, has been moved
    ! inside of "croshorizmainbc_lat".
    rhsd(loc2plusup(2)) = sum( croshorizmainbc_lat(dew,           dns,            &
                                               slopex,        slopey,         &
                                               dsigmadew(up), dsigmadns(up),  &
                                               pt,            1,              &
                                               dup(up),       local_othervel, &
                                               local_efvs,                    &
                                               oneortwo,      twoorone,       &
                                               onesideddiff,                  &
                                               normal,        fwdorbwd)       &
                                              * local_othervel ) + source
!OSBC: For the one-sided diff version of the code, need to pass in the efvs (above) since it is
! part of the LHS matrix coeffs. now.

!TODO: As noted above in line 3594, consider removing this?
!     if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag )then
!         storeoffdiag = .true.
!         h = -croshorizmainbc_lat(dew,           dns,            &
!                                 slopex,        slopey,         &
!                                 dsigmadew(up), dsigmadns(up),  &
!                                 pt,            1,              &
!                                 dup(up),       local_othervel, &
!                                 oneortwo,      twoorone,       &
!                                 onesideddiff,                  &
!                                 normal,        fwdorbwd)
!         call fillsprsebndy( h, loc2plusup(1), loc_latbc, up, normal, pt )
!         storeoffdiag = .false.
!     end if

  else   ! NOT at a lateral boundary 

! *********************************************************************************************
! normal discretization for points inside of lateral boundary and inside main body of ice sheet

!OSBC: Changes in this section are to replace ghost cells w/ one-sided diffs at sfc/basal indices.

!OSBC: This if construct skips the normal discretization for the RHS and LHS for the sfc and basal indices
     if( up /= upn .and. up /= 1 )then
         g = normhorizmain(pt,up,local_efvs)
         g(:,2,2) = g(:,2,2) + vertimain(hsum(local_efvs),up)
!OSBC: NOTE that version of 'fillspremain' for one-sided bcs needs additional index to specify column shift of
! coeffs. of rows in LHS matrix
         call fillsprsemain(g,loc2plusup(1),loc2(:,1),up,pt,0)
         ! NOTE that in the following expression, the "-" sign on the crosshoriz terms,
         ! which results from moving them from the LHS over to the RHS, is explicit and
         ! hast NOT been moved inside of "croshorizmin" (as is the case for the analogous
         ! boundary condition routines).
         rhsd(loc2plusup(2)) = thisdusrfdx(ew,ns) - &
                                            sum(croshorizmain(pt,up,local_efvs) * local_othervel)
     end if

!OSBC: Replace ghost cells w/ one-sided diffs at sfc/basal indices.
! The follow two if constructs set the ghost cells to have ones on the diag and zeros on the rhs,
! enforcing a zero vel bc for the ghost cells.
     if( up == upn  )then
        loc2plusup = loc2(1,:) + upn + 1    ! basal ghost cells
        call valueset(0.d0, loc2plusup)
     endif
     if( up == 1  )then
        loc2plusup = loc2(1,:)              ! sfc ghost cells
        call valueset(0.d0, loc2plusup)
     endif

!TODO: As noted above in line 3594, consider removing this?
!     if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag )then
!         storeoffdiag = .true.
!         h = croshorizmain(pt,up,local_efvs)   
!         call fillsprsemain(h,loc2plusup(1),loc2(:,1),up,pt)
!         storeoffdiag = .false.
!     end if     

  end if

! *********************************************************************************************
! higher-order sfc and bed boundary conditions in main body of ice sheet (NOT at lat. boundry)

!TODO - whichbabc option numbers should not be hardwired.
!       This code will break if we change the numbering.

  if(  ( up == upn  .or. up == 1 ) .and. .not. lateralboundry) then

        if( up == 1 )then                ! specify necessary variables and flags for free sfc
           bcflag = (/1,0/)
           loc2plusup = loc2(1,:) + up - 1   ! reverse the sparse matrix / rhs vector row index by 1 ...
           slopex = -dusrfdew(ew,ns); slopey = -dusrfdns(ew,ns); nz = 1.d0
        else                             ! specify necessary variables and flags for basal bc

           if( whichbabc == 6 )then
                bcflag = (/0,0/)             ! flag for u=v=0 at bed; doesn't work well so commented out here...
                                             ! better to specify very large value for betasquared below
           elseif( whichbabc >=0 .and. whichbabc <= 5 )then
                bcflag = (/1,1/)              ! flag for specififed stress at bed: Tau_zx = betasquared * u_bed,
                                              ! where betasquared is MacAyeal-type traction parameter
           elseif( whichbabc == 7 )then

                write(*,*)"ERROR: This option is not supported in the current release of CISM."
                write(*,*)"       A future release will support use of Newton iteration on " 
                write(*,*)"       plastic till basal BC."
                stop

                bcflag = (/1,2/)              ! plastic bed iteration using Newton implementation

           end if
           
           loc2plusup = loc2(1,:) + up + 1   ! advance the sparse matrix / rhs row vector index by 1 ...
           slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns); nz = -1.d0
        
        end if

!OSBC: OLD bc method, using ghost cells; this subroutine has been moved to bottom module.
!     g = normhorizmainbc(dew,           dns,     &
!                         slopex,        slopey,  &
!                         dsigmadew(up), dsigmadns(up),  &
!                         pt,            bcflag,  &
!                         dup(up),                &
!                         oneorfour,     fourorone)

!OSBC: NEW bc method, using one-sided diffs. Note that this subroutine
! has a new, slightly diff name than the original.  
     g = normhorizmainbcos(dew,           dns,          &
                         slopex,        slopey,         &
                         dsigmadew(up), dsigmadns(up),  &
                         pt,            bcflag,         &
                         dup(up),       local_efvs,     &
                         oneorfour,     fourorone)

     g_norm = g              ! save for basal traction calculation

     ! Now add on coeff. associated w/ du/dsigma

!OSBC: OLD bc method, using ghost cells; this subroutine has been moved however, because it
! it is still used in the case of sfc/basal bcs when at a lateral boundary.
!     g(:,2,2) = g(:,2,2)   &
!              + vertimainbc( stagthck(ew,ns),bcflag,dup(up),local_efvs,betasquared, &
!                            g_vert, nz, plastic_coeff=plastic_coeff_lhs(pt,ew,ns) )

!OSBC: NEW bc method, using one-sided diffs. Note that this subroutine
! has a new, slightly diff name than the original.  
     g(:,2,2) = g(:,2,2)   &
              + vertimainbcos( stagthck(ew,ns),bcflag,dup(up),local_efvs,betasquared, &
                            g_vert, nz, plastic_coeff=plastic_coeff_lhs(pt,ew,ns) )

     !! scale basal bc coeffs when using JFNK solver 
     scalebabc = scalebasalbc( g, bcflag, lateralboundry, betasquared, local_efvs )
     g = g / scalebabc

     ! put the coeff. for the b.c. equation in the same place as the prev. equation
     ! (w.r.t. cols), on a new row ...
     !call fillsprsemain(g,loc2plusup(1),loc2(:,1),up,pt)

!OSBC: The call above goes w/ the old ghost cell bc implementation. The new set of calls is below.
!
! Replace ghost cells w/ one-sided diffs at sfc/basal indices. This section shifts the LHS matrix coeffs for the sfc
! and basal bcs back on to the main diagonal, as opposed to staggered off the diag, which was necessary for the ghost 
! cell implementation.
     if( up == 1 .or. up == upn )then
       loc2plusup = loc2(1,:) + up  ! Need to reset this index since we want the bc on the actual row
                                    ! coinciding with the boundary at up=1
        if( up == 1 )then
          call fillsprsemain(g,loc2plusup(1),loc2(:,1),up,pt,1)
        else if( up == upn )then
          call fillsprsemain(g,loc2plusup(1),loc2(:,1),up,pt,-1)
        end if
     end if


     ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
     ! which results from moving them from the LHS over to the RHS, has been moved
     ! inside of "croshorizmainbc".
     if( bcflag(2) == 2 )then

! OSBC: NOTE that the one-sided basal bc implementation has not yet been updated to work with the 
! plastic bed basal boundary conditions using Newton methods. However, the code above will already 
! report an error and exit the code if those basal bcs are asked for in the configuration file.
!
! Since this is the only place in the one-side basal bc version of the code where 'crosshorizmainbc' is 
! called, I've (1) moved this subroutine to the end of the module, along w/ other subroutines/functions that
! are redundant for the time being, and (2) commented out the call here.
!          rhsd(loc2plusup(2)) = sum( croshorizmainbc(dew,           dns,            &
!                                             slopex,        slopey,         &
!                                             dsigmadew(up), dsigmadns(up),  &
!                                             pt,            bcflag,         &
!                                             dup(up),       local_othervel, &
!                                             local_efvs,                    &
!                                             oneortwo,      twoorone,       &
!                                             g_cros,                       &
!                                             velbc=velbcvect(pt,ew,ns),     &
!                                             plastic_coeff=plastic_coeff_rhs(pt,ew,ns) ) &
!                                              * local_othervel )            &
!                                             + plastic_rhs(pt,ew,ns) / (sum(local_efvs(2,:,:))/4.d0) &
!                                             *(len0/thk0) / scalebabc

!TODO: As noted above in line 3594, consider removing this?
!         if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag )then
!             storeoffdiag = .true.
!             h = -croshorizmainbc(dew,           dns,            &
!                                 slopex,        slopey,         &
!                                 dsigmadew(up), dsigmadns(up),  &
!                                 pt,            bcflag,         &
!                                 dup(up),       local_othervel, &
!                                 local_efvs,                    &
!                                 oneortwo, twoorone, g_cros ) / scalebabc
!             call fillsprsemain(h,loc2plusup(1),loc2(:,1),up,pt)
!             storeoffdiag = .false.
!         end if

     else if( bcflag(2) /= 2 )then

! OSBC: Need to reset this index since we want the bc coeffs. on the actual row coinciding 
! with the boundary at up=1.
     loc2plusup = loc2(1,:) + up    

! OSBC: for one-side bc implementation, call slightly diff. subroutine here.
!          rhsd(loc2plusup(2)) = sum( croshorizmainbc(dew,           dns,            &
         rhsd(loc2plusup(2)) = sum( croshorizmainbcos(dew,          dns,            &
                                             slopex,        slopey,         &
                                             dsigmadew(up), dsigmadns(up),  &
                                             pt,            bcflag,         &
                                             dup(up),       local_othervel, &
                                             local_efvs,                    &
                                             oneortwo, twoorone, g_cros )  &
                                              * local_othervel ) / scalebabc 

!TODO: As noted above in line 3594, consider removing this?
!         if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag)then
!             storeoffdiag = .true.
!             h = -croshorizmainbc(dew,           dns,            &
!                                 slopex,        slopey,         &
!                                 dsigmadew(up), dsigmadns(up),  &
!                                 pt,            bcflag,         &
!                                 dup(up),       local_othervel, &
!                                 local_efvs,                    &
!                                 oneortwo, twoorone, g_cros ) / scalebabc
!             call fillsprsemain(h,loc2plusup(1),loc2(:,1),up,pt)
!             storeoffdiag = .false.
!         end if

      end if

      ! The following calculates the basal traction AFTER an updated solution is obtain by passing the new
      ! values of uvel, vvel back to the matrix assembly routines, and thus obtaining updated values of the 
      ! relevant coefficients. The if construct allows the assembly routines to be called for only the vert
      ! layers that are needed to cacluate the basal traction (as opposed to all vert levels 1:upn).
      if( assembly == 1 )then

!OSBC: one-sided bc implementation does not need ghost vels for traction calculation (below, old calc.
! is commented out, followed by new calculation)
!      select case( pt )
!         case(1)
!           g_vel_lhs(:,:,:) = ghostbvel(1,:,ew-1:ew+1,ns-1:ns+1)
!           g_vel_rhs(:,:,:) = ghostbvel(2,:,ew-1:ew+1,ns-1:ns+1)
!         case(2)
!           g_vel_lhs(:,:,:) = ghostbvel(2,:,ew-1:ew+1,ns-1:ns+1)
!           g_vel_rhs(:,:,:) = ghostbvel(1,:,ew-1:ew+1,ns-1:ns+1)
!       end select

           g_vel_lhs = local_thisvel
           g_vel_rhs = local_othervel

!HALO - Since ew and ns are locally owned velocity points, we will have btraction at all such points.

!       btraction(pt,ew,ns) = sum( (g_norm+g_vert)*g_vel_lhs*thk0/len0*sum(local_efvs(2,:,:))/4.d0 ) &
!                           - sum( g_cros*g_vel_rhs*thk0/len0*sum(local_efvs(2,:,:))/4.d0 )
! NOTE that one-sided bc implementation also uses diff. scales, since efvs is now explicitly attached to the bc coeffs
! in the LHS matrix, as opposed to being divided through and placed in the RHS vector. 
       btraction(pt,ew,ns) = sum( (g_norm+g_vert)*g_vel_lhs*thk0/len0 ) &
                           - sum( g_cros*g_vel_rhs*thk0/len0 )

     end if

  end if   ! (up = 1 or up = upn) and lateralboundry = F

! *********************************************************************************************

  return

end subroutine bodyset

!***********************************************************************

subroutine valueset(local_value, loc2plusup)

  ! plugs given value into the right location in the rhs vector of matrix equation Ax=rhs

  implicit none

  real(dp), intent(in) :: local_value
  integer, dimension(2), intent(in) :: loc2plusup

  call putpcgc(1.d0,loc2plusup(1),loc2plusup(1))
  rhsd(loc2plusup(2)) = local_value

  return

end subroutine valueset

!***********************************************************************

subroutine calccoeffsinit (upn, dew, dns)

  ! determines constants used in various FD calculations associated with 'findcoefst'
  ! In general, the constants contain (1) grid spacing info, (2) numeric constants 
  ! used for averaging of eff. visc. from normal grid in horiz onto stag grid in horiz. 
  implicit none

  integer, intent(in) :: upn
  real(dp), intent(in) :: dew, dns

  ! this coefficient used in finite differences of vertical terms.
  cvert(:) = (len0**2) / (4.d0 * thk0**2 * dup**2)

  ! these coefficients used in finite differences of horizontal terms
  ! for d/dx(fdu/dx), d/dx(fdu/dy), d/dsigma(fdu/dx), d/dx(fdu/dsigma) and
  ! du/dsigma. 
  cdxdx = (/ 0.25d0 / dew**2, 0.25d0 / dns**2 /)
  cdsdx(:,1) = 0.0625d0 / (dew * dup); cdsdx(:,2) = 0.0625d0 / (dns * dup);
  cdsds = 0.25d0 / (dup * dup)
  cds = 0.0625d0 / dup
  cdxdy = 0.0625d0 / (dew * dns)

  return

end subroutine calccoeffsinit

!***********************************************************************

subroutine calccoeffs(upn,        sigma,                    &
                      stagthck,                             &
                      dusrfdew,   dusrfdns,                 &
                      dthckdew,   dthckdns,                 &
                      d2usrfdew2, d2usrfdns2, d2usrfdewdns, &
                      d2thckdew2, d2thckdns2, d2thckdewdns)

  ! Called from 'findcoefst' to find coefficients in stress balance equations
  ! Detemines coeficients needed for finite differencing.
  ! This is a column-based operation. In general these coefficients refer 
  ! to grid transformations and averaging of efvs to half grid points.

  implicit none

  integer, intent(in) :: upn
  real(dp), dimension(:), intent(in) :: sigma
  real(dp), intent(in) :: stagthck, dusrfdew, dusrfdns, dthckdew, dthckdns, &
                                  d2usrfdew2, d2usrfdns2, d2usrfdewdns, &
                                  d2thckdew2, d2thckdns2, d2thckdewdns

  fvert(:) = cvert(:) / stagthck**2

  dsigmadew = calcdsigmadx(upn, sigma, dusrfdew, dthckdew, stagthck)
  dsigmadns = calcdsigmadx(upn, sigma, dusrfdns, dthckdns, stagthck)

  d2sigmadew2 = calcd2sigmadxdy(upn,        sigma,      &
                                d2usrfdew2, d2thckdew2, &
                                dusrfdew,   dusrfdew,   &
                                dthckdew,   dthckdew,   &
                                stagthck)

  d2sigmadns2 = calcd2sigmadxdy(upn,        sigma,      &
                                d2usrfdns2, d2thckdns2, &
                                dusrfdns,   dusrfdns,   &
                                dthckdns,   dthckdns,   &
                                stagthck)

  d2sigmadewdns = calcd2sigmadxdy(upn,          sigma,         &
                                  d2usrfdewdns, d2thckdewdns,  &
                                  dusrfdew,     dusrfdns,      &
                                  dthckdew,     dthckdns,      &
                                  stagthck)

  d2sigmadewdsigma = calcd2sigmadxdsigma(dthckdew,stagthck)
  d2sigmadnsdsigma = calcd2sigmadxdsigma(dthckdns,stagthck)

  return

end subroutine calccoeffs

!***********************************************************************

function calcdsigmadx(upn,     sigma,    &
                      dusrfdx, dthckdx,  &
                      stagthck)

  implicit none

  integer, intent(in) :: upn
  real(dp), dimension(:), intent(in) :: sigma
  real(dp), intent(in) :: stagthck, dusrfdx, dthckdx
  real(dp), dimension(upn) :: calcdsigmadx

  calcdsigmadx = (dusrfdx - sigma * dthckdx) / stagthck

  return

end function calcdsigmadx

!***********************************************************************

function calcd2sigmadxdy(upn,        sigma,       &
                         d2usrfdxdy, d2thckdxdy,  &
                         dusrfdx,    dusrfdy,     &
                         dthckdx,    dthckdy,     &
                         stagthck)

  implicit none

  integer, intent(in) :: upn
  real(dp), dimension(:), intent(in) :: sigma
  real(dp), intent(in) :: d2usrfdxdy, d2thckdxdy, dusrfdx, dusrfdy, &
                                  dthckdx, dthckdy, stagthck
  real(dp), dimension(upn) :: calcd2sigmadxdy

  calcd2sigmadxdy = (stagthck * d2usrfdxdy - &
                     dusrfdx * dthckdy - dusrfdy * dthckdx + &
                     sigma * (2.d0 * dthckdx * dthckdy - &
                     stagthck * d2thckdxdy)) / stagthck**2

  return

end function calcd2sigmadxdy

!***********************************************************************

function calcd2sigmadxdsigma(dthckdx,stagthck)

  implicit none

  real(dp), intent(in) :: dthckdx, stagthck
  real(dp) :: calcd2sigmadxdsigma

  calcd2sigmadxdsigma = - dthckdx / stagthck

  return

end function calcd2sigmadxdsigma

!***********************************************************************

function vertimain(efvs,up)

  implicit none

  real(dp), dimension(2), intent(in) :: efvs

  real(dp), dimension(3) :: vertimain

  integer, intent(in) :: up

  vertimain(3) = fvert(up) * efvs(2)
  vertimain(1) = fvert(up) * efvs(1)
  vertimain(2) = - vertimain(3) - vertimain(1)

  return

end function vertimain

!***********************************************************************

function normhorizmain(which,up,efvs)

  ! Called from 'findcoefst' to calculate normal-stress grad terms 
  !      like: d/dx(f(du/dx)), d/dy(f(dv/dy)), etc.  
  ! ... calls FUNCTIONS: horiztermdxdx, horiztermdsdx, horiztermdxds,
  !                      horiztermdsds, horiztermds 
  ! determines coefficients from d/dx(fdu/dx) and d/dy(fdu/dy)

  implicit none

  integer, intent(in) :: which, up
  real(dp), dimension(:,:,:), intent(in) :: efvs

  real(dp), dimension(3,3,3) :: normhorizmain
  real(dp), dimension(3,3,3) :: g, h
  real(dp), dimension(2) :: sumefvsup, sumefvsew, sumefvsns
  real(dp) :: sumefvs

  g = 0.d0
  h = 0.d0

  sumefvsup = hsum(efvs)
  sumefvsew = sum(sum(efvs,3),1)
  sumefvsns = sum(sum(efvs,2),1)
  sumefvs = sum(efvs)

! for d(f.du/dx)/dx

  g(2,:,2) = horiztermdxdx(sumefvsew,cdxdx(1))
  g(:,1:3:2,2) = g(:,1:3:2,2) + horiztermdsdx(dsigmadew(up),sumefvsup,cdsdx(up,1))
  g(1:3:2,:,2) = g(1:3:2,:,2) + horiztermdxds(dsigmadew(up),sumefvsew,cdsdx(up,1))
  g(:,2,2) = g(:,2,2) + horiztermdsds(dsigmadew(up)**2,sumefvsup,cdsds(up))
  g(1:3:2,2,2) = g(1:3:2,2,2) + horiztermds(d2sigmadew2(up)+d2sigmadewdsigma*dsigmadew(up),sumefvs,cds(up))

! for d(f.du/dy)/dy 

  h(2,2,:) = horiztermdxdx(sumefvsns,cdxdx(2))
  h(:,2,1:3:2) = h(:,2,1:3:2) + horiztermdsdx(dsigmadns(up),sumefvsup,cdsdx(up,2))
  h(1:3:2,2,:) = h(1:3:2,2,:) + horiztermdxds(dsigmadns(up),sumefvsns,cdsdx(up,2))
  h(:,2,2) = h(:,2,2) + horiztermdsds(dsigmadns(up)**2,sumefvsup,cdsds(up))
  h(1:3:2,2,2) = h(1:3:2,2,2) + horiztermds(d2sigmadns2(up)+d2sigmadnsdsigma*dsigmadns(up),sumefvs,cds(up))

  normhorizmain = g * fourorone(which) + h * oneorfour(which)

  return

end function normhorizmain

!***********************************************************************

function croshorizmain(which,up,efvs)

  ! Called from 'findcoefst' to calculate cross-stress grad terms 
  !      like: d/dx(f(du/dy)), d/dy(f(dv/dx)), etc.  
  ! ... calls FUNCTIONS: horiztermdxdy, horiztermdsdx, horiztermdxds, 
  !                      horiztermdsds, horiztermds 
  ! determines coefficients from d/dx(fdu/dy) and d/dy(fdu/dx)

  implicit none

  integer, intent(in) :: which, up
  real(dp), dimension(:,:,:), intent(in) :: efvs

  real(dp), dimension(3,3,3) :: croshorizmain
  real(dp), dimension(3,3,3) :: g = 0.d0, h = 0.d0
  real(dp), dimension(2) :: sumefvsup, sumefvsew, sumefvsns
  real(dp) :: sumefvs

  g = 0.d0
  h = 0.d0

  sumefvsup = hsum(efvs)
  sumefvsew = sum(sum(efvs,3),1)
  sumefvsns = sum(sum(efvs,2),1)
  sumefvs = sum(efvs)

! for d(f.du/dy)/dx

  g(2,:,1:3:2) = horiztermdxdy(sumefvsew,cdxdy)
  g(:,2,1:3:2) = g(:,2,1:3:2) + horiztermdsdx(dsigmadew(up),sumefvsup,cdsdx(up,2))
  g(1:3:2,:,2) = g(1:3:2,:,2) + horiztermdxds(dsigmadns(up),sumefvsew,cdsdx(up,1))
  g(:,2,2) = g(:,2,2) + horiztermdsds(dsigmadew(up)*dsigmadns(up),sumefvsup,cdsds(up))
  g(1:3:2,2,2) = g(1:3:2,2,2) + horiztermds(d2sigmadewdns(up)+d2sigmadnsdsigma*dsigmadew(up),sumefvs,cds(up))

! for d(f.du/dx)/dy 

  h(2,1:3:2,:) = transpose(horiztermdxdy(sumefvsns,cdxdy))
  h(:,1:3:2,2) = h(:,1:3:2,2) + horiztermdsdx(dsigmadns(up),sumefvsup,cdsdx(up,1))
  h(1:3:2,2,:) = h(1:3:2,2,:) + horiztermdxds(dsigmadew(up),sumefvsns,cdsdx(up,2))
  h(:,2,2) = h(:,2,2) + horiztermdsds(dsigmadew(up)*dsigmadns(up),sumefvsup,cdsds(up))
  h(1:3:2,2,2) = h(1:3:2,2,2) + horiztermds(d2sigmadewdns(up)+d2sigmadewdsigma*dsigmadns(up),sumefvs,cds(up))

  croshorizmain = g * twoorone(which) + h * oneortwo(which)

  return

end function croshorizmain

!***********************************************************************

! ***************************************************************************
! start of functions to deal with higher-order boundary conditions at sfc and bed
! ***************************************************************************

function vertimainbc(thck, bcflag, dup, efvs, betasquared, g_vert, nz, plastic_coeff)

! altered form of 'vertimain' that calculates coefficients for higher-order
! b.c. that go with the 'normhorizmain' term: -(X/H)^2 * dsigma/dzhat * du/dsigma 

    implicit none

    real(dp), intent(in) :: dup, thck, betasquared 
    real(dp), intent(in) :: nz                      ! sfc normal vect comp in z-dir
    real(dp), intent(in), dimension(2,2,2) :: efvs
    real(dp), intent(out), dimension(3,3,3) :: g_vert
    real(dp), optional, intent(in) :: plastic_coeff
    integer, intent(in), dimension(2) :: bcflag

    real(dp) :: c
    real(dp), dimension(3) :: vertimainbc

    c = 0.d0
    g_vert = 0.d0

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    if( bcflag(1) == 1 )then

           c = nz / thck / (2*dup) * (len0**2 / thk0**2)   ! value of coefficient

           vertimainbc(:) = 0.d0
           vertimainbc(3) = -c
           vertimainbc(1) = c
           vertimainbc(2) = vertimainbc(3) + vertimainbc(1) ! should = 0

           ! this is the part of the vertimain coeff. block that we want to keep for calc
           ! of boundary tractions (note that it DOES NOT include terms from boundary forcing)
           g_vert(:,2,2) = vertimainbc

   ! for higher-order BASAL B.C. w/ specified basal traction, add on the necessary source term ...
    if( bcflag(2) == 1 )then

            ! last set of terms is mean visc. of ice nearest to the bed
            vertimainbc(2) = vertimainbc(2)   &
                           + ( betasquared / ( sum( efvs(2,:,:) ) / 4.d0 ) ) * (len0 / thk0)
    end if

    ! for higher-order BASAL B.C. w/ plastic yield stress iteration ...
    if( bcflag(2) == 2 )then

             ! last set of terms is mean visc. of ice nearest to the bed
            vertimainbc(2) = vertimainbc(2)   &
                           + ( plastic_coeff / ( sum( efvs(2,:,:) ) / 4.d0 ) ) * (len0 / thk0)
    end if

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this is not often implemented, as it is generally sufficient to implement 
    ! an "almost" no slip BC by just making the coeff. for betasquared very large (and the 
    ! the code converges more quickly/stably in this case than for actual no-slip).
    else if( bcflag(1) == 0 )then

           ! if u,v set to 0, there are no coeff. assoc. with du/digma terms ...
           vertimainbc(:) = 0.d0

    end if

    return

end function vertimainbc

!***********************************************************************

function vertimainbcos(thck, bcflag, dup, efvs, betasquared, g_vert, nz, plastic_coeff)

! altered form of 'vertimain' that calculates coefficients for higher-order
! b.c. that go with the 'normhorizmain' term: -(X/H)^2 * dsigma/dzhat * du/dsigma

    implicit none

    real (dp), intent(in) :: dup, thck, betasquared
    real (dp), intent(in) :: nz                      ! sfc normal vect comp in z-dir
    real (dp), intent(in), dimension(2,2,2) :: efvs
    real (dp), intent(out), dimension(3,3,3) :: g_vert
    real (dp), optional, intent(in) :: plastic_coeff
    integer, intent(in), dimension(2) :: bcflag

    real (dp) :: c
    real (dp), dimension(3) :: vertimainbcos
    real (dp) :: bar_sfc, bar_bed, efvsbar_bed, efvsbar_sfc

     ! averaging number for eff. visc. at domain edges
    bar_sfc = sum( (efvs(1,:,:)/efvs(1,:,:)), efvs(1,:,:) > effstrminsq )
    bar_bed = sum( (efvs(2,:,:)/efvs(2,:,:)), efvs(2,:,:) > effstrminsq )

    ! average visc. to use in coeff. calc.
    efvsbar_sfc = sum( efvs(1,:,:), efvs(1,:,:) > effstrminsq ) / bar_sfc
    efvsbar_bed = sum( efvs(2,:,:), efvs(2,:,:) > effstrminsq ) / bar_bed

    ! make the following lines active to turn OFF the visc. dependence in the LHS matrix coeffs.
    !efvsbar_sfc = 1.0d0; efvsbar_bed = 1.0d0

    c = 0.d0
    g_vert = 0.d0

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    if( bcflag(1) == 1 .and. bcflag(2) == 0 )then

           c = nz / thck / (2*dup) * (len0**2 / thk0**2) * efvsbar_sfc ! value of coefficient

           vertimainbcos(:) = 0.d0
           vertimainbcos(1) = 3.d0*c
           vertimainbcos(2) = -4.d0*c
           vertimainbcos(3) = c

           ! this is the part of the vertimain coeff. block that we want to keep for calc
           ! of boundary tractions (note that it DOES NOT include terms from boundary forcing)
           g_vert(:,2,2) = vertimainbcos

   end if

   ! for higher-order BASAL B.C. w/ specified basal traction, add on the necessary source term ...
   if( bcflag(1) == 1 .and. bcflag(2) == 1 )then

           c = nz / thck / (2*dup) * (len0**2 / thk0**2) * efvsbar_bed ! value of coefficient

           vertimainbcos(:) = 0.d0
           vertimainbcos(1) = -1.d0*c
           vertimainbcos(2) = 4.d0*c
           vertimainbcos(3) = -3.d0*c

            ! this is the part of the vertimain coeff. block that we want to keep for calc
            ! of boundary tractions (note that it DOES NOT include terms from boundary forcing)
            ! NOTE that here we do this BEFORE adding in the sliding coefficient, as in the standard
            ! expression for the BC, this term is on the RHS.
            g_vert(:,2,2) = vertimainbcos

           ! this is the part of the vertimain coeff. block that we want to keep for calc
           ! of boundary tractions (note that it DOES NOT include terms from boundary forcing)

            ! last set of terms is mean visc. of ice nearest to the bed
!            vertimainbcos(3) = vertimainbcos(3)   &
!                           + ( betasquared / efvsbar_bed ) * (len0 / thk0)
            vertimainbcos(3) = vertimainbcos(3)   &
                           + ( betasquared ) * (len0 / thk0)

    end if

    ! for higher-order BASAL B.C. w/ plastic yield stress iteration ...
    if( bcflag(2) == 2 )then

             ! last set of terms is mean visc. of ice nearest to the bed
!            vertimainbcos(2) = vertimainbcos(2)   &
!                           + ( plastic_coeff / efvsbar_bed ) * (len0 / thk0)
! ... but one-sided implementation doesn't need this, since efvs stays on LHS
            vertimainbcos(2) = vertimainbcos(2)   &
                           + ( plastic_coeff ) * (len0 / thk0)
    end if

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this is not often implemented, as it is generally sufficient to implement
    ! an "almost" no slip BC by just making the coeff. for betasquared very large (and the
    ! the code converges more quickly/stably in this case than for actual no-slip).
    if( bcflag(1) == 0 )then

           ! if u,v set to 0, there are no coeff. assoc. with du/digma terms ...
           vertimainbcos(:) = 0.d0

    end if

    return

end function vertimainbcos

!***********************************************************************

function normhorizmainbcos(dew,       dns,      &
                         dusrfdew,  dusrfdns,   &
                         dsigmadew, dsigmadns,  &
                         which,     bcflag,     &
                         dup,       efvs,       &
                         oneorfour, fourorone)

    ! Determines higher-order surface and basal boundary conditions for LHS of equation.
    ! Gives 3x3x3 coeff. array for either u or v component of velocity, depending on the
    ! value of the flag 'which'. Example of function call:
    !
    !  g = normhorizmainbc(dusrfew(ew,ns),dusrfnx(ew,ns),dsigmadew(up),dsigmadns(up),which,up,bcflag)
    !
    ! ... where g is a 3x3x3 array.
    !
    ! 'bcflag' is a 1 x 2 vector to indicate (1) which b.c. is being solved for (surface or bed) and
    ! (2), if solving for the bed b.c., which type of b.c. to use. For example, bcflag = [ 0, 0 ]
    ! denotes free sfc bc; bcflag = [ 1, 0 ] denotes basal bc w/ u=v=0, etc. (see also subroutine
    ! "bodyset"). "fourorone" and "oneorfour" are given by vectors: fourorone = [ 4 1 ]; oneorfour = [ 1 4 ].
    ! A single value is chosen from each vector and applied to the calculation of coefficients below.
    ! The "correct" value needed to satisfy the expression is chosen based on the "which" flag, which
    ! takes on a value of 1 for calculations in the x direction and a value of 2 for calculations in
    ! the y direction.

    implicit none

    real (kind = dp), intent(in) :: dew, dns
    real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real (kind = dp), intent(in), dimension(2) :: oneorfour, fourorone
    real (kind = dp), dimension(3,3,3) :: normhorizmainbcos
    real (kind = dp), dimension(3,3,3) :: g
    real (kind = dp) :: c

    integer, intent(in) :: which
    integer, intent(in), dimension(2) :: bcflag
    real (kind = dp), intent(in), dimension(2,2,2) :: efvs

    c = 0.d0
    g(:,:,:) = 0.d0

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this handles the case for specified stress at the bed as well, as we
    ! simply pass in a different value for the normal vector (slope) components (still
    ! called "dusrfdns", "dusrfdew" here, but args passed in are different).
    if( bcflag(1) == 1 .and. bcflag(2) == 0 )then


           ! first, coeff. that go with du/dsigma, and thus are associated
           ! with u(1,2,2) and u(3,2,2) ...
!           c = ( fourorone(which) * dusrfdew * dsigmadew   &
!               + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup)
           c = ( fourorone(which) * dusrfdew * dsigmadew   &
               + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup) * ( sum( efvs(1,:,:) ) / 4.d0 )

           g(1,2,2) = 3.d0*c
           g(2,2,2) = -4.d0*c
           g(3,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
!           c = fourorone(which) * dusrfdew / (2*dew)
           c = fourorone(which) * dusrfdew / (2*dew) * ( sum( efvs(1,:,:) ) / 4.d0 )
           g(1,3,2) = c
           g(1,1,2) = -c

!           c = oneorfour(which) * dusrfdns / (2*dns)
           c = oneorfour(which) * dusrfdns / (2*dns) * ( sum( efvs(1,:,:) ) / 4.d0 )
           g(1,2,3) = c
           g(1,2,1) = -c

    end if

    ! higher-order, specified traction basal bc, must use fwd rather than bwd one-sided
    ! diff in vertical direction
    if( bcflag(1) == 1 .and. bcflag(2) == 1 )then

           ! first, coeff. that go with du/dsigma, and thus are associated
           ! with u(1,2,2) and u(3,2,2) ...
!           c = ( fourorone(which) * dusrfdew * dsigmadew   &
!               + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup)
           c = ( fourorone(which) * dusrfdew * dsigmadew   &
               + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup) * ( sum( efvs(2,:,:) ) / 4.d0 )

           g(1,2,2) = -1.d0*c
           g(2,2,2) = 4.d0*c
           g(3,2,2) = -3.d0*c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
!           c = fourorone(which) * dusrfdew / (2*dew)
           c = fourorone(which) * dusrfdew / (2*dew) * ( sum( efvs(2,:,:) ) / 4.d0 )
           g(3,3,2) = c
           g(3,1,2) = -c

!           c = oneorfour(which) * dusrfdns / (2*dns)
           c = oneorfour(which) * dusrfdns / (2*dns) * ( sum( efvs(2,:,:) ) / 4.d0 )
           g(3,2,3) = c
           g(3,2,1) = -c

    end if

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! note that this requires that rhs(up) be set to 0 as well ...
    if( bcflag(1) == 0 )then

           g(:,:,:) = 0.d0
           g(2,2,2) = 1.d0;

    end if

    normhorizmainbcos = g

    return

end function normhorizmainbcos

!***********************************************************************

function croshorizmainbcos(dew,       dns,       &
                         dusrfdew,  dusrfdns,  &
                         dsigmadew, dsigmadns, &
                         which,     bcflag,    &
                         dup,       local_othervel,  &
                         efvs,                       &
                         oneortwo,  twoorone,        &
                         g_cros, velbc, plastic_coeff )

    ! As described for "normhorizmainbc" above. The vectors "twoorone" and
    ! "oneortwo" are given by: twoorone = [ 2 1 ]; oneortwo = [ 1 2 ];

    implicit none

    integer, intent(in) :: which
    integer, intent(in), dimension(:) :: bcflag

    real (kind = dp), intent(in) :: dew, dns
    real (kind = dp), intent(in), dimension(:) :: oneortwo, twoorone
    real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real (kind = dp), intent(in), dimension(:,:,:) :: local_othervel
    real (kind = dp), intent(in), dimension(:,:,:) :: efvs
    real (kind = dp), intent(in), optional :: velbc, plastic_coeff
    real (kind = dp), intent(out),dimension(:,:,:) :: g_cros


    real (kind = dp), dimension(3,3,3) :: g, croshorizmainbcos
    real (kind = dp) :: c
    integer :: nz

    c = 0.d0
    g(:,:,:) = 0.d0
    g_cros = g
    nz = 0

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this handles the case for specified stress at the bed as well, as we
    ! simply pass in a different value for the normal vector (slope) components (still
    ! called "dusrfdns", "dusrfdew" here, but args passed in are different).
    if( bcflag(1) == 1 .and. bcflag(2) == 0 )then

           ! first, coeff. that go with du/dsigma, and thus are associated
           ! with u(1,2,2) and u(3,2,2) ...
!           c = ( - twoorone(which) * dusrfdew * dsigmadns   &
!                 - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup)
           c = ( - twoorone(which) * dusrfdew * dsigmadns   &
                 - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup) * ( sum( efvs(1,:,:) ) / 4.d0 )

           g(1,2,2) = 3.d0*c
           g(2,2,2) = -4.d0*c
           g(3,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
!           c = - oneortwo(which) * dusrfdns / (2*dew)
           c = - oneortwo(which) * dusrfdns / (2*dew) * ( sum( efvs(1,:,:) ) / 4.d0 )
           g(1,3,2) = c
           g(1,1,2) = -c

!           c = - twoorone(which) * dusrfdew / (2*dns)
           c = - twoorone(which) * dusrfdew / (2*dns) * ( sum( efvs(1,:,:) ) / 4.d0 )
           g(1,2,3) = c
           g(1,2,1) = -c

    end if

    ! higher-order, specified traction basal bc, must use fwd rather than bwd one-sided
    ! diff in vertical direction
    if( bcflag(1) == 1 .and. bcflag(2) == 1 )then

           ! first, coeff. that go with du/dsigma, and thus are associated
           ! with u(1,2,2) and u(3,2,2) ...
!           c = ( - twoorone(which) * dusrfdew * dsigmadns   &
!                 - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup)
           c = ( - twoorone(which) * dusrfdew * dsigmadns   &
                 - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup) * ( sum( efvs(2,:,:) ) / 4.d0 )

           g(1,2,2) = -1.d0*c
           g(2,2,2) = 4.d0*c
           g(3,2,2) = -3.d0*c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
!           c = - oneortwo(which) * dusrfdns / (2*dew)
           c = - oneortwo(which) * dusrfdns / (2*dew) * ( sum( efvs(2,:,:) ) / 4.d0 )
           g(3,3,2) = c
           g(3,1,2) = -c


!           c = - twoorone(which) * dusrfdew / (2*dns)
           c = - twoorone(which) * dusrfdew / (2*dns) * ( sum( efvs(2,:,:) ) / 4.d0 )
           g(3,2,3) = c
           g(3,2,1) = -c

    end if

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! This forces the multiplication by 'local_otherval' in the main program
    ! to result in a value of 1, thus leaving the boundary vel. unchanged
    ! ... conditional makes sure there is no div by zero if the bc value IS also zero
    if( bcflag(1) == 0 )then

        g(:,:,:) = 0.d0

        where( local_othervel /= 0.d0 )
            g = 1
        elsewhere
            g = 0.d0
        endwhere

        nz = sum( g )
        g(:,:,:) = 0.d0

        where( local_othervel /= 0.d0 )
            g = ( velbc / nz ) / local_othervel
        elsewhere
            g = 0.d0
        endwhere

     end if

     ! NOTE: here we define 'g_cros' FIRST, because we want the value w/o the plastic
     ! bed coeff. included (needed for estimate of basal traction in plastic bed iteration)
     g_cros = g

     if( bcflag(2) == 2 )then        ! add on coeff. associated w/ plastic bed iteration

         g(2,2,2) = g(2,2,2) + plastic_coeff / ( sum( efvs(2,:,:) ) / 4.d0 ) * (len0 / thk0)

     end if

    croshorizmainbcos = g

    return

end function croshorizmainbcos

!***********************************************************************

function normhorizmainbc_lat(dew,       dns,   &
                             dusrfdew,  dusrfdns,  &
                             dsigmadew, dsigmadns, &
                             which,     what,      &
                             dup,       efvs,      &
                             oneorfour, fourorone, &
                             onesideddiff,         &
                             normal,    fwdorbwd)

    ! Analogous to "normhorizmainbc" but for the case of lateral stress (ice shelf)
    ! boundary conditions. Note that the basic form of the equations is the same. 
    ! What changes here is (1) the value of the normal vector that is passed in (at
    ! the sfc and bed we pass in the surface or basal slopes, while at the boundaries
    ! we use the normal vector orientation to the boundary in map view) and (2) we to
    ! to use one sided diffs at the lateral boundaries rather than centerd diffs.

    ! Note that we assume here that du/dz (and thus du/dsigma) is approx. 0 for an ice 
    ! shelf, and also that the sfc/basal slopes of an ice shelf are very flat at/near 
    ! the boundary. Thus, we assume flow is depth independent and we ignore gradients 
    ! in sigma. 

    implicit none

    real(dp), intent(in) :: dew, dns
    real(dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real(dp), intent(in), dimension(2) :: oneorfour, fourorone, normal, fwdorbwd
    real(dp), intent(in), dimension(3) :: onesideddiff
    real (kind = dp), intent(in), dimension(:,:,:) :: efvs

    integer, intent(in) :: which, what

    real(dp), dimension(3,3,3) :: normhorizmainbc_lat
    real(dp), dimension(3,3,3) :: g
    real(dp), dimension(2) :: whichbc
    real(dp) :: c
    real (kind = dp) :: bar, efvsbar

    c = 0.d0; g(:,:,:) = 0.d0; whichbc = (/ 0.d0, 1.d0 /)

    ! averaging number for eff. visc. at domain edges
    bar = sum( (efvs(:,:,:)/efvs(:,:,:)), efvs(:,:,:) > effstrminsq )

    ! average visc. to use in coeff. calc.
    efvsbar = sum( efvs(:,:,:), efvs(:,:,:) > effstrminsq ) / bar

    ! make the following lines active to turn OFF the visc. dependence in the LHS matrix coeffs.
!    efvsbar = 1.0d0; 

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    ! (also applies to basal stress bc) 

    ! first, coeff. that go with du/dsigma, and thus are associated with u(1,2,2) and u(3,2,2) ...
    ! ...note that these are stored in an empty column of 'g' (a corner column) so that we don't 
    ! overwrite these values in the case of fwd/bwd horiz. diffs., which require 3 spaces
    c = ( fourorone(which) * dusrfdew * dsigmadew    &
            + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup) * efvsbar
    g(3,3,3) = -c * whichbc(what)
    g(1,3,3) = c * whichbc(what)

    if( normal(1) == 0.d0 )then     ! centered in x ...

           c = fourorone(which) * dusrfdew / (2*dew) * efvsbar
           g(2,3,2) = c
           g(2,1,2) = -c

    elseif( normal(1) /= 0.d0 )then     ! forward/backward in x ...

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(1) * dusrfdew / (2*dew) * efvsbar
           g(2,2-int(fwdorbwd(1)),2) = c

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(2) * dusrfdew / (2*dew) * efvsbar
           g(2,2,2) = c

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(3) * dusrfdew / (2*dew) * efvsbar
           g(2,2+int(fwdorbwd(1)),2) = c

    end if

    if( normal(2) == 0.d0 ) then   ! centered in y ... 
                                       ! (NOTE that y coeff. are stored in g(1,:,:) )

           c = oneorfour(which) * dusrfdns / (2*dns) * efvsbar
           g(1,2,3) = c
           g(1,2,1) = -c

    elseif( normal(2) /= 0.d0) then ! forward/backward in y ...

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(1) * dusrfdns / (2*dns) * efvsbar
           g(1,2,2-int(fwdorbwd(2))) = c

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(2) * dusrfdns / (2*dns) * efvsbar
           g(1,2,2) = c

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(3) * dusrfdns / (2*dns) * efvsbar
           g(1,2,2+int(fwdorbwd(2))) = c

    end if

    normhorizmainbc_lat = g

    return

end function normhorizmainbc_lat

!***********************************************************************

function croshorizmainbc_lat (dew,       dns,       &
                              dusrfdew,  dusrfdns,  &
                              dsigmadew, dsigmadns, &
                              which,     what,      &
                              dup,       local_othervel,  &
                              efvs,                 &
                              oneortwo,  twoorone,  &
                              onesideddiff,         &
                              normal,    fwdorbwd)

    ! Analagous to "normhorizmainbc_lat" but for cross terms. See notes above.

    implicit none

    real(dp), intent(in) :: dew, dns
    real(dp), intent(in), dimension(2) :: oneortwo, twoorone, fwdorbwd, normal
    real(dp), intent(in), dimension(3) :: onesideddiff
    real(dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real(dp), intent(in), dimension(3,3,3) :: local_othervel
    real (kind = dp), intent(in), dimension(:,:,:) :: efvs

    integer, intent(in) :: which, what

    real(dp), dimension(3,3,3) :: g, croshorizmainbc_lat
    real(dp), dimension(3) :: gvert
    real(dp), dimension(2) :: whichbc
    real(dp) :: c

    integer, dimension(2) :: inormal

    real (kind = dp) :: bar, efvsbar

    ! averaging number for eff. visc. at domain edges
    bar = sum( (efvs(:,:,:)/efvs(:,:,:)), efvs(:,:,:) > effstrminsq )

    ! average visc. to use in coeff. calc.
    efvsbar = sum( efvs(:,:,:), efvs(:,:,:) > effstrminsq ) / bar

    ! make the following lines active to turn OFF the visc. dependence in the LHS matrix coeffs.
!    efvsbar = 1.0d0; 

    c = 0.d0
    g(:,:,:) = 0.d0
    gvert = 0.d0
    whichbc = (/ 0.d0, 1.d0 /)
    croshorizmainbc_lat = 0.d0

    ! first, coeff. that go with du/dsigma, and thus are associated with u(1,2,2) and u(3,2,2) 
    ! ... note that these are stored in a separate vector (to avoid being overwritten if stored in normal 'g')  
    c = ( - twoorone(which) * dusrfdew * dsigmadns   &
              - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup) * efvsbar
    gvert(3) = -c * whichbc(what)
    gvert(1) = c * whichbc(what)

    if( normal(1) == 0.d0 )then        ! centered in x ...

           c = -oneortwo(which) * dusrfdns / (2*dew) * efvsbar
           g(2,3,2) = c
           g(2,1,2) = -c

    elseif( normal(1) /= 0.d0 )then    ! forward/backward in x ...
                                           ! (NOTE that x coeff. are stored in g(2,:,:) )

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(1) * dusrfdns / (2*dew) * efvsbar
           g(2,2-int(fwdorbwd(1)),2) = c

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(2) * dusrfdns / (2*dew) * efvsbar
           g(2,2,2) = c

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(3) * dusrfdns / (2*dew) * efvsbar
           g(2,2+int(fwdorbwd(1)),2) = c

    end if

    if( normal(2) == 0.d0 )then    ! centered in y ...
                                       ! (NOTE that y coeff. are stored in g(1,:,:) )

           c = -twoorone(which) * dusrfdew / (2*dns) * efvsbar
           g(1,2,3) = c
           g(1,2,1) = -c

    elseif( normal(2) /= 0.d0 )then ! forward/backward in y ...

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(1) * dusrfdew / (2*dns) * efvsbar
           g(1,2,2-int(fwdorbwd(2))) = c

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(2) * dusrfdew / (2*dns) * efvsbar
           g(1,2,2) = c

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(3) * dusrfdew / (2*dns) * efvsbar
           g(1,2,2+int(fwdorbwd(2))) = c

    end if

    ! Now rearrange position of coefficients in structure 'g' so that they are multiplied by 
    ! the correct velocity component of 'local_othervel' in 'bodyset' ...
    ! ... this can be done by using the boundary normal vector to shift the indices of the rows/columns
    ! in 'g', in the appropriate direction. First, convert the boundary normal to an integer index ...
    inormal(1) = int( normal(1)/abs(normal(1)) )
    inormal(2) = int( normal(2)/abs(normal(2)) )
    if( abs( inormal(1) ) /= 1 )then; inormal(1) = 0; end if
    if( abs( inormal(2) ) /= 1 )then; inormal(2) = 0; end if

    croshorizmainbc_lat(2,:,2+inormal(2)) = g(2,:,2)    ! move x-coeffs. appropriate amount
    croshorizmainbc_lat(1,2+inormal(1),:) = g(1,2,:)    ! move y-coeffs. appropriate amount

    ! sum coeffs. that are in same column and flatten so that all coeff. are on level (2,:,:)   
    croshorizmainbc_lat(2,:,:) = croshorizmainbc_lat(2,:,:) + croshorizmainbc_lat(1,:,:)

    ! set remaining coeff. on this level to to 0 ...
    croshorizmainbc_lat(1,:,:) = 0.d0

    ! accounter for vertical terms stored seperately and temporarily in 'gvert'
    croshorizmainbc_lat(1,2+inormal(1),2+inormal(2)) = gvert(1) * whichbc(what)
    croshorizmainbc_lat(3,2+inormal(1),2+inormal(2)) = gvert(3) * whichbc(what)

    return

end function croshorizmainbc_lat

!***********************************************************************

! ---> the following routines are for derivatives in the main body 

function horiztermdxdx(efvs,fact)

  ! this is the d/dx(f.du/dx) and d/dy(f.du/dy) terms

  implicit none

  real(dp), dimension(2), intent(in) :: efvs
  real(dp), intent(in) :: fact

  real(dp), dimension(3) :: horiztermdxdx  

  horiztermdxdx(3) = efvs(2) * fact
  horiztermdxdx(1) = efvs(1) * fact
  horiztermdxdx(2) = - horiztermdxdx(3) - horiztermdxdx(1)

  return

end function horiztermdxdx

!***********************************************************************

function horiztermdxdy(efvs,fact)

  ! this is the d/dy(f.du/dx) and d/dx(f.du/dy) terms

  implicit none

  real(dp), dimension(2), intent(in) :: efvs
  real(dp), intent(in) :: fact

  real(dp), dimension(3,2) :: horiztermdxdy

  horiztermdxdy(3,2) = efvs(2) * fact 
  horiztermdxdy(2,2) = horiztermdxdy(3,2)
  horiztermdxdy(3,1) = - horiztermdxdy(3,2)
  horiztermdxdy(2,1) = - horiztermdxdy(3,2)

  horiztermdxdy(1,2) = - efvs(1) * fact
  horiztermdxdy(2,2) = horiztermdxdy(2,2) + horiztermdxdy(1,2)
  horiztermdxdy(2,1) = horiztermdxdy(2,1) - horiztermdxdy(1,2)
  horiztermdxdy(1,1) = - horiztermdxdy(1,2)

  return

end function horiztermdxdy

!***********************************************************************

function horiztermdsdx(dsigmadxy,efvs,fact)

  ! this is the d/ds(f.du/dx) and d/ds(f.du/dy) terms

  implicit none

  real(dp), dimension(2), intent(in) :: efvs
  real(dp), intent(in) :: dsigmadxy, fact

  real(dp), dimension(3,2) :: horiztermdsdx  

  horiztermdsdx(3,2) = dsigmadxy * efvs(2) * fact
  horiztermdsdx(2,2) = horiztermdsdx(3,2)
  horiztermdsdx(3,1) = - horiztermdsdx(3,2)
  horiztermdsdx(2,1) = - horiztermdsdx(3,2)

  horiztermdsdx(1,2) = - dsigmadxy * efvs(1) * fact
  horiztermdsdx(2,2) = horiztermdsdx(2,2) + horiztermdsdx(1,2)
  horiztermdsdx(2,1) = horiztermdsdx(2,1) - horiztermdsdx(1,2)
  horiztermdsdx(1,1) = - horiztermdsdx(1,2)

  return

end function horiztermdsdx

!***********************************************************************

function horiztermdxds(dsigmadxy,efvs,fact)

  ! this is the d/dx(f.du/ds) and d/dy(f.du/ds) terms

  implicit none

  real(dp), dimension(2), intent(in) :: efvs
  real(dp), intent(in) :: dsigmadxy, fact

  real(dp), dimension(2,3) :: horiztermdxds

  horiztermdxds(2,3) = dsigmadxy * efvs(2) * fact
  horiztermdxds(2,2) = horiztermdxds(2,3)
  horiztermdxds(1,3) = - horiztermdxds(2,3)
  horiztermdxds(1,2) = - horiztermdxds(2,3)

  horiztermdxds(2,1) = - dsigmadxy * efvs(1) * fact
  horiztermdxds(2,2) = horiztermdxds(2,2) + horiztermdxds(2,1)
  horiztermdxds(1,2) = horiztermdxds(1,2) - horiztermdxds(2,1)
  horiztermdxds(1,1) = - horiztermdxds(2,1)

  return

end function horiztermdxds

!***********************************************************************

function horiztermdsds(dsigmadxysq,efvs,fact)

  ! this is the d/ds(f.du/ds) term

  implicit none

  real(dp), dimension(2), intent(in) :: efvs
  real(dp), intent(in) :: dsigmadxysq, fact

  real(dp), dimension(3) :: horiztermdsds

  horiztermdsds(3) = dsigmadxysq * efvs(2) * fact
  horiztermdsds(1) = dsigmadxysq * efvs(1) * fact

  horiztermdsds(2) = - horiztermdsds(3) - horiztermdsds(1)

  return

end function horiztermdsds

!***********************************************************************

function horiztermds(d2sigmadxy2etc,efvs,fact)

  ! this is the f.du/ds term

  implicit none

  real(dp), intent(in) :: efvs, d2sigmadxy2etc, fact

  real(dp), dimension(2) :: horiztermds

  horiztermds(2) = d2sigmadxy2etc * efvs * fact
  horiztermds(1) = - horiztermds(2)

  return

end function horiztermds

! ---> end of routines for derivatives in the main body 

!***********************************************************************
 
subroutine fillsprsemain(inp,locplusup,ptindx,up,pt,osshift)

  ! scatter coefficients from 3x3x3 block "g" onto sparse matrix row
  implicit none

  real(dp), dimension(3,3,3), intent(in):: inp
  integer, intent(in) :: locplusup, up, pt
  integer, dimension(6), intent(in) :: ptindx
  integer, intent(in) :: osshift

  ! insert entries to "g" that are on same level
  call putpcgc(inp(2,2,2),ptindx(1)+up+osshift,locplusup,pt)
  call putpcgc(inp(2,3,2),ptindx(2)+up+osshift,locplusup,pt)
  call putpcgc(inp(2,1,2),ptindx(3)+up+osshift,locplusup,pt)
  call putpcgc(inp(2,2,3),ptindx(4)+up+osshift,locplusup,pt)
  call putpcgc(inp(2,2,1),ptindx(5)+up+osshift,locplusup,pt)

  ! add points for level above (that is, points in "g"  with a LARGER first index,
  ! which correspond to grid points that are CLOSER TO THE BED than at current level)
  call putpcgc(inp(3,2,2),ptindx(1)+up+1+osshift,locplusup,pt)
  call putpcgc(inp(3,3,2),ptindx(2)+up+1+osshift,locplusup,pt)
  call putpcgc(inp(3,1,2),ptindx(3)+up+1+osshift,locplusup,pt)
  call putpcgc(inp(3,2,3),ptindx(4)+up+1+osshift,locplusup,pt)
  call putpcgc(inp(3,2,1),ptindx(5)+up+1+osshift,locplusup,pt)

  ! add points for level below (that is, points in "g" with a SMALLER first index,
  ! which correspond to grid points that are CLOSER TO THE SURFACE than at current level)
  call putpcgc(inp(1,2,2),ptindx(1)+up-1+osshift,locplusup,pt)
  call putpcgc(inp(1,3,2),ptindx(2)+up-1+osshift,locplusup,pt)
  call putpcgc(inp(1,1,2),ptindx(3)+up-1+osshift,locplusup,pt)
  call putpcgc(inp(1,2,3),ptindx(4)+up-1+osshift,locplusup,pt)
  call putpcgc(inp(1,2,1),ptindx(5)+up-1+osshift,locplusup,pt)

  return

end subroutine fillsprsemain

!***********************************************************************

subroutine fillsprsebndy(inp,locplusup,ptindx,up,normal,pt)

  ! scatter coeff. from 3x3x3 block "g" onto sparse matrix row. This subroutine
  ! is specifically for the boundary conditions, which are handled differently
  ! than points in the "main" body of the domain (interior to boundaries).
  implicit none

  integer, intent(in) :: locplusup, up, pt
  integer, dimension(6), intent(in) :: ptindx
  real(dp), dimension(3,3,3), intent(in) :: inp
  real(dp), dimension(2), intent(in) :: normal

  ! at points where mixed centered and one-side diffs. would apply
  if( normal(1) == 0.d0 )then         ! at boundary normal to y, centered diffs in x 
    if( normal(2) == -1.d0 )then      ! at boundary w/ normal [0,-1]
           call putpcgc(inp(1,3,3),ptindx(5)+up-1,locplusup,pt)
           call putpcgc( inp(2,3,3)+inp(1,2,1),ptindx(5)+up,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(5)+up+1,locplusup,pt)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup,pt)
    else                                ! at boundary w/ normal [0,1]
           call putpcgc(inp(1,3,3),ptindx(4)+up-1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(1,2,3),ptindx(4)+up,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(4)+up+1,locplusup,pt)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup,pt)
    end if
    call putpcgc(inp(1,2,2),ptindx(1)+up,locplusup,pt)
  end if

  if( normal(2) == 0.d0 )then            ! at boundary normal to x, centered diffs in y 
        if( normal(1) == -1.d0 )then     ! at boundary w/ normal [-1,0]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup,pt)
           call putpcgc( inp(2,3,3)+inp(2,1,2),ptindx(3)+up,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup,pt)
        else                                 ! at boundary w/ normal [1,0]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup,pt)
           call putpcgc( inp(2,3,3)+inp(2,3,2),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup,pt)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup,pt)
    end if
    call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
  end if

  ! at corners where only one-side diffs. apply
  if( normal(1) > 0.d0 .and. normal(2) /= 0.d0 )then
    if( normal(2) > 0.d0 )then      ! corner w/ normal [ 1/sqrt(2), 1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(2,3,2)+inp(1,2,3),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup,pt)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup,pt)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup,pt)
    else                                 ! corner w/ normal [ 1/sqrt(2), -1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(1,2,1)+inp(2,3,2),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup,pt)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup,pt)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup,pt)
    end if
  end if

  if( normal(1) < 0.d0 .and. normal(2) /= 0.d0 )then
    if( normal(2) > 0.d0 )then       ! corner w/ normal [ -1/sqrt(2), 1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(1,2,3)+inp(2,1,2),ptindx(3)+up,locplusup,pt)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup,pt)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup,pt)
    else                                  ! corner w/ normal [ -1/sqrt(2), -1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(2,1,2)+inp(1,2,1),ptindx(3)+up,locplusup,pt)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup,pt)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup,pt)
    end if
  end if

  return

end subroutine fillsprsebndy

!***********************************************************************

subroutine getlatboundinfo( ew, ns, up, ewn, nsn, upn,  &
                           thck, loc_array,             &
                           fwdorbwd, normal, loc_latbc)

  ! Calculate map plane normal vector at 45 deg. increments
  ! for regions of floating ice
  implicit none

  integer, intent(in) :: ew, ns, up
  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(ewn,nsn), intent(in) :: loc_array
  real(dp), dimension(3,3), intent(in) :: thck

  real(dp), dimension(2), intent(out) :: fwdorbwd, normal
  integer, dimension(6), intent(out) :: loc_latbc

  real(dp), dimension(3,3) :: mask, maskcorners
  real(dp), dimension(3,3) :: thckmask
  real(dp), dimension(3) :: testvect
  real(dp) :: phi, deg2rad

  deg2rad = 3.141592654d0 / 180.d0
  loc_latbc = 0; phi = 0
  mask(:,1) = (/ 0.d0, 180.d0, 0.d0 /)
  mask(:,2) = (/ 270.d0, 0.d0, 90.d0 /)
  mask(:,3) = (/ 0.d0, 360.d0, 0.d0 /)
  maskcorners(:,1) = (/ 225.d0, 0.d0, 135.d0 /)
  maskcorners(:,2) = (/ 0.d0, 0.d0, 0.d0 /)
  maskcorners(:,3) = (/ 315.d0, 0.d0, 45.d0 /)

  ! specify new value of 'loc' vector such that fwd/bwd diffs. are set up correctly in sparse matrix
  ! when function 'fillsprsebndy' is called. Also, specify appropriate values for the vectors 'normal'
  ! and 'fwdorbwd', which specify the orientation of the boundary normal and the direction of forward or
  ! backward differencing to be done in the lateral boundary condition functions 'normhorizmainbc_lat'
  ! and 'crosshorizmainbc_lat'

  ! following is algorithm for calculating boundary normal at 45 deg. increments, based on arbitray
  ! boundary shape (based on initial suggestions by Anne LeBrocq)
  where( thck /= 0.d0 )
        thckmask = 0.d0
  elsewhere( thck == 0.d0 )
        thckmask = 1.d0
  endwhere

  testvect = sum( thckmask * mask, 1 )

    ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 
    ! 90 deg. at 3 O'clock, etc.
    if( sum( sum( thckmask, 1 ) ) == 1.d0 )then
        phi = sum( sum( thckmask * maskcorners, 1 ) )
    else
        if( any( testvect == 360.d0 ) )then
            if( sum( testvect ) == 450.d0 )then
                phi = 45.d0
            elseif( sum( testvect ) == 630.d0 )then
                phi = 315.d0
            else
                phi = 0.d0
            end if
        elseif( all( testvect /= 360 ) )then
            phi = sum( testvect ) / sum( testvect/testvect, testvect /= 0.d0 )
        end if
    end if

    ! define normal vectors and change definition of loc_array based on this angle
    if( phi == 0.d0 )then
         loc_latbc(1) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         loc_latbc(2) = loc_array(ew+1,ns); loc_latbc(3) = loc_array(ew-1,ns)
         normal = (/ 0.d0, 1.d0 /); fwdorbwd = (/ -1.d0, -1.d0 /)
    elseif( phi == 45.d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(6) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         normal = (/ 1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0) /); fwdorbwd = (/ -1.d0, -1.d0 /)
    elseif( phi == 90.d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(4) = loc_array(ew,ns+1); loc_latbc(5) = loc_array(ew,ns-1)
         normal = (/ 1.d0, 0.d0 /); fwdorbwd = (/ -1.d0, -1.d0 /)
    elseif( phi == 135.d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(6) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns)
         normal = (/ 1.d0/sqrt(2.d0), -1.d0/sqrt(2.d0) /); fwdorbwd = (/ -1.d0, 1.d0 /)
    elseif( phi == 180.d0 )then
         loc_latbc(1) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns)
         loc_latbc(2) = loc_array(ew+1,ns); loc_latbc(3) = loc_array(ew-1,ns)
         normal = (/ 0.d0, -1.d0 /); fwdorbwd = (/ 1.d0, 1.d0 /)
    elseif( phi == 225.d0 )then
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(6) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns);
         normal = (/ -1.d0/sqrt(2.d0), -1.d0/sqrt(2.d0) /); fwdorbwd = (/ 1.d0, 1.d0 /)
    elseif( phi == 270.d0 )then
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(4) = loc_array(ew,ns+1); loc_latbc(5) = loc_array(ew,ns-1)
         normal = (/ -1.d0, 0.d0 /); fwdorbwd = (/ 1.d0, 1.d0 /)
    else
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(6) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         normal = (/ -1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0) /); fwdorbwd = (/ 1.d0, -1.d0 /)
    end if

  return

end subroutine getlatboundinfo

!***********************************************************************

function indshift( which, ew, ns, up, ewn, nsn, upn, loc_array, thck )

  ! Subroutine to rearrange indices slightly at sfc,bed, and lateral boundaries,
  ! so that values one index inside of the domain are used for, e.g. eff. visc.

  ! Function output is a vector containing necessary index shifts for portions of 'othervel' and 'efvs' 
  ! extracted near domain boundaries. NOTE that this contains duplication of some of the code in the 
  ! subroutine "getlatboundinfo", and the two could be combined at some point.
  implicit none

  integer, intent(in) :: which
  integer, intent(in) :: ew, ns, up, ewn, nsn, upn
  integer, dimension(ewn,nsn), intent(in) :: loc_array
  real(dp), dimension(3,3), intent(in) :: thck

  integer, dimension(3) :: indshift
  integer :: upshift = 0, ewshift = 0, nsshift = 0

  real(dp), dimension(3,3) :: mask, maskcorners
  real(dp), dimension(3,3) :: thckmask
  real(dp), dimension(3) :: testvect
  real(dp) :: phi, deg2rad

  deg2rad = 3.141592654d0 / 180.d0
  mask(:,1) = (/ 0.d0, 180.d0, 0.d0 /)
  mask(:,2) = (/ 270.d0, 0.d0, 90.d0 /)
  mask(:,3) = (/ 0.d0, 360.d0, 0.d0 /)
  maskcorners(:,1) = (/ 225.d0, 0.d0, 135.d0 /)
  maskcorners(:,2) = (/ 0.d0, 0.d0, 0.d0 /)
  maskcorners(:,3) = (/ 315.d0, 0.d0, 45.d0 /)

  if( up == 1 )then   !! first treat bed/sfc, which aren't complicated
      upshift = 1
  elseif( up == upn )then
      upshift = -1
  else
      upshift = 0
  end if

  select case(which)

      case(0)   !! internal to lateral boundaries; no shift to ew,ns indices

          ewshift = 0; nsshift = 0;

      case(1)   !! at lateral boundaries; shift to ew,ns may be non-zero

          where( thck /= 0.d0 )
            thckmask = 0.d0
          elsewhere( thck == 0.d0 )
            thckmask = 1.d0
          endwhere

          testvect = sum( thckmask * mask, 1 )

        ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 90 deg. at 3 O'clock, etc.
        if( sum( sum( thckmask, 1 ) ) == 1.d0 )then
            phi = sum( sum( thckmask * maskcorners, 1 ) )
        else
            if( any( testvect == 360.d0 ) )then
                if( sum( testvect ) == 450.d0 )then
                    phi = 45.d0
                elseif( sum( testvect ) == 630.d0 )then
                    phi = 315.d0
                else
                    phi = 0.d0
                end if
            elseif( all( testvect /= 360 ) )then
                phi = sum( testvect ) / sum( testvect/testvect, testvect /= 0.d0 )
            end if
        end if

        ! define shift to indices based on this angle 
        if( phi == 0.d0 )then
            nsshift = -1; ewshift = 0
        elseif( phi == 45.d0 )then
            nsshift = -1; ewshift = -1
        elseif( phi == 90.d0 )then
            nsshift = 0; ewshift = -1
        elseif( phi == 135.d0 )then
            nsshift = 1; ewshift = -1
        elseif( phi == 180.d0 )then
            nsshift = 1; ewshift = 0
        elseif( phi == 225.d0 )then
            nsshift = 1; ewshift = 1
        elseif( phi == 270.d0 )then
            nsshift = 0; ewshift = 1
        elseif( phi == 315.d0 )then
            nsshift = -1; ewshift = 1
        end if

  end select

  indshift = (/ upshift, ewshift, nsshift /)

  return

end function indshift

!***********************************************************************

subroutine calcbetasquared (whichbabc,               &
                            dew,         dns,        &
                            ewn,         nsn,        &
                            lsrf,        topg,       &
                            thck,                    &
                            thisvel,     othervel,   &
                            minTauf, beta,           &
                            betasquared, betafile)

  ! subroutine to calculate map of betasquared sliding parameter, based on 
  ! user input ("whichbabc" flag, from config file as "which_ho_babc").
  implicit none

  integer, intent(in) :: whichbabc
  integer, intent(in) :: ewn, nsn

  real(dp), intent(in) :: dew, dns
  real(dp), intent(in), dimension(:,:) :: lsrf, topg, thck
  real(dp), intent(in), dimension(:,:) :: thisvel, othervel, minTauf, beta

  real(dp), intent(out), dimension(ewn-1,nsn-1) :: betasquared

  character (len=30), intent(in), optional :: betafile
  real(dp) :: smallnum = 1.0d-2
  real(dp), dimension(ewn) :: grounded
  real(dp) :: alpha, dx, thck_gl, betalow, betahigh, roughness
  integer :: ew, ns

  ! Note that the dimensional scale (tau0 / vel0 / scyr) is used here for making the basal traction coeff.
  ! betasquared dimensional, within the subroutine, and then non-dimensional again before being sent back out
  ! for use in the code. This scale is the same as scale2d_beta defined in libglimmer/glimmer_scales.F90.

!TODO - Remove scaling here?

  select case(whichbabc)

    case(0)     ! constant value; useful for debugging and test cases

      betasquared = 10.d0

    case(1)     ! simple pattern; also useful for debugging and test cases
                ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      betasquared = 1.d4

!TODO - Should these 5's be hardwired?  Is 10.d1 correct?  (Change to 100.d0?)
      do ns=5, nsn-5; do ew=1, ewn-1; 
        betasquared(ew,ns) = 10.d1 
      end do; end do


    case(2)     ! take input value for till yield stress and force betasquared to be implemented such
                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).

      !!! NOTE: Eventually, this option will provide the till yield stress as calculate from the basal processes
      !!! submodel. Currently, to enable sliding over plastic till, simple specify the value of "betasquared" as 
      !!! if it were the till yield stress (in units of Pascals).
!      betasquared = minTauf*tau0 / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )

      betasquared = ( beta * ( tau0 / vel0 / scyr  ) ) &
                    / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )

    case(3)     ! circular ice shelf: set B^2 ~ 0 except for at center, where B^2 >> 0 to enforce u,v=0 there

      betasquared = 1.d-5
      betasquared( (ewn-1)/2:(ewn-1)/2+1, (nsn-1)/2:(nsn-1)/2+1 ) = 1.d10

    case(4)    ! frozen (u=v=0) ice-bed interface

      betasquared = 1.d10

    case(5)    ! use value passed in externally from CISM (NOTE not dimensional when passed in) 

!TODO - Careful with scaling here.
      ! scale CISM input value to dimensional units of (Pa yrs 1/m)
      betasquared = beta * ( tau0 / vel0 / scyr )

      ! this is a check for NaNs, which indicate, and are replaced by no slip
      where ( betasquared /= betasquared )
        betasquared = 1.d10
      end where

    ! NOTE: cases (6) and (7) are handled external to this subroutine

  end select

  ! convert whatever the specified value is to dimensional units of (Pa s m^-1 ) 
  ! and then non-dimensionalize using PP dyn core specific scaling.
  betasquared = betasquared / ( tau0 / vel0 / scyr )

end subroutine calcbetasquared

!***********************************************************************

subroutine plasticbediteration( )  

    !!! NOTE: This subroutine is under development, in support of  Newton-type iteration on a plastic-till
    !!!       basal boundary condition. It will be supported in a future release of CISM.

end subroutine plasticbediteration

!***********************************************************************

!TODO - Might be cleaner to just inline the vertical loop wherever this function is called.

function vertintg(upn, sigma, in)

  implicit none

  integer, intent(in) :: upn
  real(dp), dimension(:), intent(in) :: sigma
  real(dp), dimension(:), intent(in) :: in
  real(dp) :: vertintg

  integer :: up

  vertintg = 0.d0

  do up = upn-1, 1, -1
    vertintg = vertintg + sum(in(up:up+1)) * dups(up)
  end do

  vertintg = vertintg / 2.d0

  return

end function vertintg

!***********************************************************************

subroutine geom2derscros(dew,  dns,   &
                         ipvr, stagthck, opvrewns)

  ! geometric (2nd) cross-deriv. for generic input variable 'ipvr', output as 'opvr'       

  implicit none

  real(dp), intent(in) :: dew, dns
  real(dp), intent(out), dimension(:,:) :: opvrewns
  real(dp), intent(in), dimension(:,:) :: ipvr, stagthck

!TODO - Should these be loops over locally owned velocity points? I.e. (ilo-1:ihi, jlo-1:jhi).
 
  ! consider replacing by a loop over ewn, nsn?
  where (stagthck /= 0.d0)
    opvrewns = (eoshift(eoshift(ipvr,1,0.d0,2),1,0.d0,1) + ipvr   &
               - eoshift(ipvr,1,0.d0,1) - eoshift(ipvr,1,0.d0,2)) / (dew*dns)
  elsewhere
    opvrewns = 0.d0
  end where

  return

end subroutine geom2derscros


!***********************************************************************

subroutine geom2ders(ewn,    nsn,  &
                     dew,    dns,  &
                     ipvr,   stagthck,  &
                     opvrew, opvrns)

  ! geometric 1st deriv. for generic input variable 'ipvr', 
  ! output as 'opvr' (includes 'upwinding' for boundary values)

  implicit none

  integer, intent(in) :: ewn, nsn
  real(dp), intent(in) :: dew, dns
  real(dp), intent(out), dimension(:,:) :: opvrew, opvrns
  real(dp), intent(in), dimension(:,:) :: ipvr, stagthck

  integer :: ew, ns
  real(dp) :: dewsq4, dnssq4

  integer :: pt(2)

  dewsq4 = 4.d0 * dew * dew
  dnssq4 = 4.d0 * dns * dns

!TODO - I think these loops should be over locally owned velocity points: (ilo-1:ihi, jlo-1:jhi).
!       Provided nhalo >= 2, we should have enough points to compute a centered difference.

  do ns = 2, nsn-2
  do ew = 2, ewn-2
    if (stagthck(ew,ns) > 0.d0) then
      opvrew(ew,ns) = centerew(ew,ns,ipvr,dewsq4)
      opvrns(ew,ns) = centerns(ew,ns,ipvr,dnssq4)
    else
      opvrew(ew,ns) = 0.d0
      opvrns(ew,ns) = 0.d0
    end if
  end do
  end do

  ! *** 2nd order boundaries using upwinding

!TODO - If nhalo = 2, then I'm not clear on why upwinding is needed.
!       Where are these values used in the computation?
!       I don't think they should be used for any interior halo cells.
!       Are they needed at the global boundaries?  If so, then need to use the correct indices for global boundaries.
!       Would be easier if we could set global halos in a way that gives reasonable 2nd derivs
!        without a special case.

  do ew = 1, ewn-1, ewn-2

    pt = whichway(ew)

    do ns = 2, nsn-2
      if (stagthck(ew,ns) > 0.d0) then
        opvrew(ew,ns) = boundyew(ns,pt,ipvr,dewsq4)
        opvrns(ew,ns) = centerns(ew,ns,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.d0
        opvrns(ew,ns) = 0.d0
      end if
    end do

  end do

!TODO - If nhalo = 2, then I'm not clear on why upwinding is needed.
!       Where are these values used in the computation?

  do ns = 1, nsn-1, nsn-2

    pt = whichway(ns)

    do ew = 2, ewn-2
      if (stagthck(ew,ns) > 0.d0) then
        opvrew(ew,ns) = centerew(ew,ns,ipvr,dewsq4)
        opvrns(ew,ns) = boundyns(ew,pt,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.d0
        opvrns(ew,ns) = 0.d0
      end if
    end do

  end do

!TODO - If nhalo = 2, then I'm not clear on why upwinding is needed.
!       Where are these values used in the computation?

  do ns = 1, nsn-1, nsn-2
    do ew = 1, ewn-1, ewn-2
      if (stagthck(ew,ns) > 0.d0) then
        pt = whichway(ew)
        opvrew(ew,ns) = boundyew(ns,pt,ipvr,dewsq4)
        pt = whichway(ns)
        opvrns(ew,ns) = boundyns(ew,pt,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.d0
        opvrns(ew,ns) = 0.d0
      end if
    end do
  end do

end subroutine geom2ders

!***********************************************************************

  function centerew(ew, ns, ipvr, dewsq4)
 
    implicit none

    integer, intent(in) :: ew, ns 
    real(dp), intent(in) :: ipvr(:,:)
    real(dp), intent(in) :: dewsq4
    real(dp) :: centerew
 
    centerew = (sum(ipvr(ew+2,ns:ns+1)) + sum(ipvr(ew-1,ns:ns+1)) - &
                sum(ipvr(ew+1,ns:ns+1)) - sum(ipvr(ew,ns:ns+1))) / dewsq4
 
    return
    
  end function centerew 
 
!***********************************************************************

  function centerns(ew, ns, ipvr, dnssq4)
 
    implicit none
 
    integer, intent(in) :: ew, ns 
    real(dp), intent(in) :: ipvr(:,:)
    real(dp), intent(in) :: dnssq4
    real(dp) :: centerns
 
    centerns = (sum(ipvr(ew:ew+1,ns+2)) + sum(ipvr(ew:ew+1,ns-1)) - &
                sum(ipvr(ew:ew+1,ns+1)) - sum(ipvr(ew:ew+1,ns))) / dnssq4
 
    return
    
  end function centerns 
 
!***********************************************************************

  function boundyew(ns,pt,ipvr,dewsq4)
 
    implicit none

    integer, intent(in) :: ns
    integer, intent(in) :: pt(2)
    real(dp), intent(in) :: ipvr(:,:)
    real(dp), intent(in) :: dewsq4
    real(dp) :: boundyew
 
    boundyew = pt(1) * (3.d0 * sum(ipvr(pt(2),ns:ns+1)) - 7.d0 * sum(ipvr(pt(2)+pt(1),ns:ns+1)) + &
               5.d0 * sum(ipvr(pt(2)+2*pt(1),ns:ns+1)) - sum(ipvr(pt(2)+3*pt(1),ns:ns+1))) / dewsq4
 
    return
 
  end function boundyew
 
!***********************************************************************

  function boundyns(ew,pt,ipvr,dnssq4)
 
    implicit none
 
    integer, intent(in) :: ew
    integer, intent(in) :: pt(2)
    real(dp), intent(in) :: ipvr(:,:)
    real(dp), intent(in) :: dnssq4
    real(dp) :: boundyns
 
    boundyns = pt(1) * (3.d0 * sum(ipvr(ew:ew+1,pt(2))) - 7.d0 * sum(ipvr(ew:ew+1,pt(2)+pt(1))) + &
               5.d0 * sum(ipvr(ew:ew+1,pt(2)+2*pt(1))) - sum(ipvr(ew:ew+1,pt(2)+3*pt(1)))) / dnssq4
 
    return
 
  end function boundyns
 
!***********************************************************************

  function whichway(i)
 
    implicit none
 
    integer, intent(in) :: i
    integer :: whichway(2) 
 
    if (i == 1) then 
      whichway = (/1,1/)
    else
      whichway = (/-1,i+1/)
    end if
 
    return
 
  end function whichway
 

!***********************************************************************

    function hsum(inp) 
 
      implicit none
 
      real(dp), dimension(:,:,:), intent(in) :: inp
      real(dp), dimension(size(inp,dim=1)) :: hsum
 
      hsum = sum(sum(inp(:,:,:),dim=3),dim=2)
 
      return 
 
    end function hsum

!***********************************************************************

subroutine putpcgc(value,col,row,pt) 
 
  implicit none
 
  integer, intent(in) :: row, col
  integer, intent(in), optional :: pt
  real(dp), intent(in) :: value 

   !*sfp* For now, ignoring the possibility of using JFNK w/ Trilinos ...
   if( nonlinear == HO_NONLIN_PICARD )then

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
        ! Option to load entry into Triad sparse matrix format
        if (value /= 0.d0) then
          pcgval(ct) = value
          pcgcol(ct) = col
          pcgrow(ct) = row
          ct = ct + 1
        end if
#ifdef TRILINOS
    else
        ! Option to load entry directly into Trilinos sparse matrix 
        if (value /= 0.d0) then
           !AGS: If we find that sparsity changes inside a time step,
           !     consider adding entry even for value==0.
           call putintotrilinosmatrix(row, col, value) 

           !JEFF: Verify matrix matches for globalIDs case
           ! call verify_trilinos_rowcolval(row, col, value)
        end if
#endif
    end if
 
 
   !*sfp* if using JFNK, store the main block diagonal coeffs and off diag coeffs 
   elseif ( nonlinear == HO_NONLIN_JFNK )then

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then      ! if using Triad format to store matrix entries

      if ( .not. storeoffdiag ) then ! block diag coeffs in normal storage space
          ! load entry into Triad sparse matrix format
          if (value /= 0.d0) then
            pcgval(ct) = value
            pcgcol(ct) = col
            pcgrow(ct) = row
            ct = ct + 1
          end if
      else if ( storeoffdiag ) then ! off-block diag coeffs in other storage
          ! load entry into Triad sparse matrix format
          if( pt == 1 )then ! store uv coeffs 
              if (value /= 0.d0) then
                pcgvaluv(ct2) = value
                pcgcoluv(ct2) = col
                pcgrowuv(ct2) = row
                ct2 = ct2 + 1
              end if
          else if( pt == 2 )then ! store vu coeffs
              if (value /= 0.d0) then
                pcgvalvu(ct2) = value
                pcgcolvu(ct2) = col
                pcgrowvu(ct2) = row
                ct2 = ct2 + 1
              end if
          end if
      end if

#ifdef TRILINOS
    else    ! if storing matrix entires directly in Trilinos sparse format

      ! *sfp* NOTE: that this option does not allow for storing offidag terms at present
      ! (I assume that Andy can grab these when we know they are correct)
      if (.not. storeoffdiag) then
        if (value /= 0.d0) then
           !AGS: If we find that sparsity changes inside a time step,
           !     consider adding entry even for value==0.
         call putintotrilinosmatrix(row, col, value) 
        end if
      end if
#endif

    end if  ! end of "if using Triad or Trilinos storage format" construct
 
   end if   ! end of "if using Picard or JFNK for nonlinear solve" construct
   
  return
 
end subroutine putpcgc 

!***********************************************************************

  subroutine distributed_create_partition(ewn, nsn, upstride, indxmask, mySize, myIndices, myX, myY, myZ)
  	! distributed_create_partition builds myIndices ID vector for Trilinos using (ns,ew) coordinates in indxmask
  	! upstride is the total number of vertical layers including any ghosts
  	! indxmask is ice mask with non-zero values for cells with ice.
  	! mySize is number of elements in myIndices
  	! myIndices is integer vector in which IDs are def
	  use parallel

 	  implicit none

	  integer, intent(in) :: ewn, nsn, upstride
	  integer, intent(in), dimension(:,:) :: indxmask
	  integer, intent(in) :: mySize
	  integer, intent(out), dimension(:) :: myIndices
	  real(dp),  intent(out), dimension(:) :: myX, myY, myZ

	  integer :: ew, ns, pointno
	  integer :: glblID, upindx, slnindx

!TODO - Loop over locally owned velocity points?

      ! Step through indxmask, but exclude halo
          do ns = 1+staggered_shalo,size(indxmask,2)-staggered_nhalo
             do ew = 1+staggered_whalo,size(indxmask,1)-staggered_ehalo
	        if ( indxmask(ew,ns) /= 0 ) then
	          pointno = indxmask(ew,ns)  ! Note that pointno starts at value 1.  If we step through correctly then consecutive values
	          ! write(*,*) "pointno = ", pointno
	          ! first layer ID is set from parallel_globalID, rest by incrementing through layers
	          glblID = parallel_globalID(ns, ew, upstride)
	          ! write(*,*) "global ID (ew, ns) = (", ew, ",", ns, ") ", glblID
	          upindx = 0
	          do slnindx = (pointno - 1) * upstride + 1, pointno * upstride
	              ! slnindx is offset into myIndices for current ice cell's layers. upindx is offset from current globalID.
  	              myIndices(slnindx) = glblID + upindx
                      ! Return coordinates for nodes. Assumes structured with dx=1,dy=1,dz=1.0e6
                      myX(slnindx) = (ewlb+ew) * 1.0
                      myY(slnindx) = (nslb+ns) * 1.0
                      myZ(slnindx) = upindx * 1.0e-6
  	              upindx = upindx + 1
  	              ! write(*,*) "myIndices offset = ", slnindx
	          end do
	        endif
	      end do
	  end do

	  return
  end subroutine

!***********************************************************************

  function distributed_globalID_to_localindex(globalID)
  	! distributed_globalID_to_localindex converts a globalID to its position in the solution vector. 
        ! It is a utility function that is not currently used, but retained for future debugging capability.
  	! The function searches loc2_array(:,:,1) for the globalID closest to the 
        ! given globalID, then uses this difference and loc2_array(:,:,2) for the same ew,ns coordinates
        ! to calculate (and return) the corresponding index.
        ! Result is checked using myIndices.
  	! loc2_array is assumed to be a module-level variable set by the routine getlocationarray.
  	! myIndices is assumed to be a module-level variable which holds the local processor's ID partition list.
  	! This function will work for both globalIDs and regular partitions.
  	! In the latter case it is redundant, because the ID will be at the same index, so it is just an identity function.
  	! Original implementation using myIndices, and then fast inverse, by JEFF 11/2010 and 11/2011
        ! Current loc2_array-based implementation by PW 12/2011
	  use parallel

 	  implicit none

	  integer, intent(in) :: globalID

	  integer :: distributed_globalID_to_localindex

#ifdef globalIDs
      !JEFF integer :: GlobalIDsGet ! C++ function with return value
          integer :: ew, ns
          integer :: minew, minns
          integer :: curdiff, mindiff
          integer :: lindex

      ! loc2_array-based search
          minew = 1
          minns = 1
          mindiff = globalID
!          do ns = 1+staggered_shalo,size(loc2_array,2)-staggered_nhalo
!            do ew = 1+staggered_whalo,size(loc2_array,1)-staggered_ehalo
          ! loc2_array(:,:,1) defined for all ew,ns, 
          ! while loc2_array(:,:,2) == 0 for halos and ice-free loactions
          do ns = 1,size(loc2_array,2)
            do ew = 1,size(loc2_array,1)
              curdiff = globalID-loc2_array(ew,ns,1)
              if ((curdiff >= 0) .and. (curdiff < mindiff)) then
                mindiff = globalID-loc2_array(ew,ns,1)
                minew = ew
                minns = ns
              endif
            enddo
          enddo
          lindex = loc2_array(minew,minns,2) + mindiff

          if ( myIndices(lindex) == globalID ) then
            distributed_globalID_to_localindex = lindex
            return
          else
            write(*,*) "Error in distributed_globalID_to_localindex()."
            write(*,*) "GlobalID to match = ", globalID
            write(*,*) "GlobalID found = ", myIndices(lindex), "(lindex = ",lindex,")"
            stop
          endif

      ! linear search from beginning of myIndices.
      ! Inefficient.  There could be some ordering of myIndices that would enable us to us a binary search.  Not certain at this time.
      !JEFF    do lindex = 1, size(myIndices)
      !JEFF	     if ( myIndices(lindex) == globalID ) then
      !JEFF	     	distributed_globalID_to_localindex = lindex
      !JEFF	     	return
      !JEFF	     endif
      !JEFF	  end do

#else
      distributed_globalID_to_localindex = globalID
      return
#endif

  end function distributed_globalID_to_localindex

!***********************************************************************

  subroutine verify_trilinos_rowcolval(row, col, value)
     ! Translates back globalID row and col values to their original grid values and outputs the set
     ! For verification of the matrix passed to Trilinos.
     ! JEFF November 2010
     integer, intent(in) :: row, col
     real(dp), intent(in) :: value
     integer :: locrow, loccol

#ifdef globalIDs
     locrow = distributed_globalID_to_localindex(row)
     loccol = distributed_globalID_to_localindex(col)
#else
     locrow = row
     loccol = col
#endif

     write (*,*) "Row = ", locrow, " Col = ", loccol, " Value = ", value
  end subroutine verify_trilinos_rowcolval

!***********************************************************************

function scalebasalbc( coeffblock, bcflag, lateralboundry, beta, efvs )

  ! *sfp* This function is used to scale the matrix coeffs and rhs vector coeff
  ! of the basal boundary condition when using JFNK for the nonlinear iteration
  ! (iteration on viscosity). 
  implicit none

  integer, dimension(2), intent(in) :: bcflag         
  logical :: lateralboundry
  real(dp), dimension(:,:,:), intent(in) :: coeffblock 
  real(dp), dimension(:,:,:), intent(in) :: efvs       
  real(dp), intent(in) :: beta

  real(dp) :: scale, scalebasalbc 

    if( nonlinear == 1 )then
        if( bcflag(1) == 1 )then

           ! use the dominant terms in the coeff associated with the velocity under consideration
           !scale = beta / ( sum( efvs(2,:,:) ) / 4.d0 ) * (len0 / thk0)

           ! Use the magnitude of the coeff associated with the vert stress gradients. 
           ! NOTE that relevant coeffs are stored in diff parts of block depending 
           ! on type of boudnary     
           if( lateralboundry )then
               scale = abs( coeffblock(3,3,3) );  
           else
               scale = abs( coeffblock(3,2,2) );     
           end if                

           if( scale <= 0.d0 )then
            scale = 1.d0
           end if

        else
            scale = 1.d0
        end if

    else
        scale = 1.d0
    end if

    scalebasalbc = scale

  return

end function scalebasalbc   

!***********************************************************************

subroutine init_resid_type(resid_object, model, uindx, umask, &
     d2thckdewdns, d2usrfdewdns, pcgsize, gx_flag, matrixA, matrixC, L2norm, ewn, nsn)
  
  use iso_c_binding 
  use glide_types, only : glide_global_type, pass_through
  use glimmer_sparse_type, only : sparse_matrix_type
  
  implicit none
  
  type(glide_global_type)  ,intent(in) :: model
  type(sparse_matrix_type) ,intent(in) :: matrixA, matrixC
  
  integer :: i, j
  integer                   ,intent(in) :: ewn, nsn
  integer, dimension(2)     ,intent(in) :: pcgsize
  integer                   ,intent(in) :: gx_flag(2*pcgsize(1)) ! 0 :reg cell
  integer                   ,intent(in) :: uindx(ewn-1,nsn-1), umask(ewn-1,nsn-1)
  real(dp)          ,intent(in) :: L2norm
  real(dp)          ,intent(in) :: d2thckdewdns(ewn-1,nsn-1), d2usrfdewdns(ewn-1,nsn-1)
  
  type(pass_through)     ,intent(out) :: resid_object
  
  allocate(resid_object%ui(ewn-1,nsn-1) )
  allocate(resid_object%um(ewn-1,nsn-1) ) 
  allocate(resid_object%d2thckcross(ewn-1,nsn-1) )
  allocate(resid_object%d2usrfcross(ewn-1,nsn-1) ) 
  allocate(resid_object%gxf( 2*pcgsize(1) ) )
  
  resid_object%model%general%ewn = model%general%ewn
  resid_object%model%general%nsn = model%general%nsn
  resid_object%model%general%upn = model%general%upn
  resid_object%model%numerics%dew = model%numerics%dew
  resid_object%model%numerics%dns = model%numerics%dns
  resid_object%model%numerics%sigma => model%numerics%sigma(:)
  resid_object%model%numerics%stagsigma => model%numerics%stagsigma(:)
  resid_object%model%geometry%thck => model%geometry%thck(:,:)
  resid_object%model%geometry%lsrf => model%geometry%lsrf(:,:)
  resid_object%model%geometry%topg => model%geometry%topg(:,:)
  resid_object%model%geomderv%dthckdew => model%geomderv%dthckdew(:,:)
  resid_object%model%geomderv%dthckdns => model%geomderv%dthckdns(:,:)
  resid_object%model%geomderv%dusrfdew => model%geomderv%dusrfdew(:,:)
  resid_object%model%geomderv%dusrfdns => model%geomderv%dusrfdns(:,:)
  resid_object%model%geomderv%dlsrfdew => model%geomderv%dlsrfdew(:,:)
  resid_object%model%geomderv%dlsrfdns => model%geomderv%dlsrfdns(:,:)
  resid_object%model%geomderv%d2thckdew2 => model%geomderv%d2thckdew2(:,:) 
  resid_object%model%geomderv%d2thckdns2 => model%geomderv%d2thckdns2(:,:) 
  resid_object%model%geomderv%d2usrfdew2 => model%geomderv%d2usrfdew2(:,:) 
  resid_object%model%geomderv%d2usrfdns2 => model%geomderv%d2usrfdns2(:,:) 
  resid_object%model%geomderv%stagthck => model%geomderv%stagthck(:,:)
  resid_object%model%temper%flwa => model%temper%flwa(:,:,:)
  resid_object%model%basalproc%minTauf => model%basalproc%minTauf(:,:)
  resid_object%model%velocity%btraction => model%velocity%btraction(:,:,:)
  resid_object%model%options%which_ho_babc = model%options%which_ho_babc
  resid_object%model%options%which_ho_efvs = model%options%which_ho_efvs
  resid_object%model%options%which_ho_sparse = model%options%which_ho_sparse
  resid_object%model%velocity%beta => model%velocity%beta(:,:)
  do i = 1, ewn-1 
   do j = 1, nsn-1 
    resid_object%ui(i,j)  = uindx(i,j)
    resid_object%um(i,j)  = umask(i,j)
    resid_object%d2thckcross(i,j) = d2thckdewdns(i,j) 
    resid_object%d2usrfcross(i,j) = d2usrfdewdns(i,j) 
   end do
  end do
  resid_object%pcgsize = pcgsize
  do i = 1, 2*pcgsize(1)
   resid_object%gxf(i) = gx_flag(i)
  end do
  resid_object%L2norm = L2norm
  resid_object%matrixA = matrixA
  resid_object%matrixC = matrixC
  resid_object%model%stress%efvs => model%stress%efvs(:,:,:)
  resid_object%model%velocity%uvel => model%velocity%uvel(:,:,:)
  resid_object%model%velocity%vvel => model%velocity%vvel(:,:,:)

end subroutine init_resid_type


!***********************************************************************************************
!OSBC: Below here are older subroutines that have been replaced by newer or slightly altered ones
! to allow for the new one-sided bc implementatation at the sfc and bed.
!***********************************************************************************************

function normhorizmainbc(dew,       dns,        &
                         dusrfdew,  dusrfdns,   &
                         dsigmadew, dsigmadns,  &
                         which,     bcflag,     &
                         dup,                   &
                         oneorfour, fourorone)

    ! Determines higher-order surface and basal boundary conditions for LHS of equation.
    ! Gives 3x3x3 coeff. array for either u or v component of velocity, depending on the 
    ! value of the flag 'which'. Example of function call:
    !
    !  g = normhorizmainbc(dusrfew(ew,ns),dusrfnx(ew,ns),dsigmadew(up),dsigmadns(up),which,up,bcflag)   
    !
    ! ... where g is a 3x3x3 array.
    !
    ! 'bcflag' is a 1 x 2 vector to indicate (1) which b.c. is being solved for (surface or bed) and 
    ! (2), if solving for the bed b.c., which type of b.c. to use. For example, bcflag = [ 0, 0 ] 
    ! denotes free sfc bc; bcflag = [ 1, 0 ] denotes basal bc w/ u=v=0, etc. (see also subroutine
    ! "bodyset"). "fourorone" and "oneorfour" are given by vectors: fourorone = [ 4 1 ]; oneorfour = [ 1 4 ].
    ! A single value is chosen from each vector and applied to the calculation of coefficients below.
    ! The "correct" value needed to satisfy the expression is chosen based on the "which" flag, which
    ! takes on a value of 1 for calculations in the x direction and a value of 2 for calculations in 
    ! the y direction. 

    implicit none

    real(dp), intent(in) :: dew, dns
    real(dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real(dp), intent(in), dimension(2) :: oneorfour, fourorone
    real(dp), dimension(3,3,3) :: normhorizmainbc
    real(dp), dimension(3,3,3) :: g
    real(dp) :: c

    integer, intent(in) :: which
    integer, intent(in), dimension(2) :: bcflag

    c = 0.d0
    g(:,:,:) = 0.d0

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this handles the case for specified stress at the bed as well, as we 
    ! simply pass in a different value for the normal vector (slope) components (still
    ! called "dusrfdns", "dusrfdew" here, but args passed in are different).
    if( bcflag(1) == 1 )then

           ! first, coeff. that go with du/dsigma, and thus are associated 
           ! with u(1,2,2) and u(3,2,2) ...
           c = ( fourorone(which) * dusrfdew * dsigmadew   &
               + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup)
           g(3,2,2) = -c
           g(1,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
           c = fourorone(which) * dusrfdew / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

           c = oneorfour(which) * dusrfdns / (2*dns)
           g(2,2,3) = c
           g(2,2,1) = -c

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! note that this requires that rhs(up) be set to 0 as well ...
    else if( bcflag(1) == 0 )then

           g(:,:,:) = 0.d0
           g(2,2,2) = 1.d0;

    end if

    normhorizmainbc = g

    return

end function normhorizmainbc

!***********************************************************************

function croshorizmainbc(dew,       dns,       &
                         dusrfdew,  dusrfdns,  &
                         dsigmadew, dsigmadns, &
                         which,     bcflag,    &
                         dup,       local_othervel,  &
                         efvs,                       &
                         oneortwo,  twoorone,        &
                         g_cros, velbc, plastic_coeff )

    ! As described for "normhorizmainbc" above. The vectors "twoorone" and 
    ! "oneortwo" are given by: twoorone = [ 2 1 ]; oneortwo = [ 1 2 ];

    implicit none

    integer, intent(in) :: which
    integer, intent(in), dimension(:) :: bcflag

    real(dp), intent(in) :: dew, dns
    real(dp), intent(in), dimension(:) :: oneortwo, twoorone
    real(dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real(dp), intent(in), dimension(:,:,:) :: local_othervel
    real(dp), intent(in), dimension(:,:,:) :: efvs
    real(dp), intent(in), optional :: velbc, plastic_coeff
    real(dp), intent(out),dimension(:,:,:) :: g_cros


    real(dp), dimension(3,3,3) :: g, croshorizmainbc
    real(dp) :: c
    integer :: nz

    c = 0.d0
    g(:,:,:) = 0.d0
    g_cros = g
    nz = 0

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this handles the case for specified stress at the bed as well, as we 
    ! simply pass in a different value for the normal vector (slope) components (still
    ! called "dusrfdns", "dusrfdew" here, but args passed in are different).
    if( bcflag(1) == 1 )then

           ! first, coeff. that go with du/dsigma, and thus are associated
           ! with u(1,2,2) and u(3,2,2) ...
           c = ( - twoorone(which) * dusrfdew * dsigmadns   &
                 - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup)
           g(3,2,2) = -c
           g(1,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
           c = - oneortwo(which) * dusrfdns / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

           c = - twoorone(which) * dusrfdew / (2*dns)
           g(2,2,3) = c
           g(2,2,1) = -c

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! This forces the multiplication by 'local_otherval' in the main program 
    ! to result in a value of 1, thus leaving the boundary vel. unchanged
    ! ... conditional makes sure there is no div by zero if the bc value IS also zero
    else if( bcflag(1) == 0 )then

        g(:,:,:) = 0.d0

        where( local_othervel /= 0.d0 )
            g = 1
        elsewhere
            g = 0.d0
        endwhere

        nz = sum( g )
        g(:,:,:) = 0.d0

        where( local_othervel /= 0.d0 )
            g = ( velbc / nz ) / local_othervel
        elsewhere
            g = 0.d0
        endwhere

     end if

     ! NOTE: here we define 'g_cros' FIRST, because we want the value w/o the plastic
     ! bed coeff. included (needed for estimate of basal traction in plastic bed iteration)
     g_cros = g

     if( bcflag(2) == 2 )then        ! add on coeff. associated w/ plastic bed iteration

         g(2,2,2) = g(2,2,2) + plastic_coeff / ( sum( efvs(2,:,:) ) / 4.d0 ) * (len0 / thk0)

     end if

    croshorizmainbc = g

    return

end function croshorizmainbc

!***********************************************************************************************
!OSBC:ABOVE here are older subroutines that have been replaced by newer or slightly altered ones
!***********************************************************************************************


end module glam_strs2

!***********************************************************************
