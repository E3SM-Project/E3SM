!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glam_strs2.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! 3d velocity calculation based on Blatter/Pattyn, 1st-order equations, by Tony Payne (Univ.
! of Bristol) and Steve Price (Univ. of Bristol / Los Alamos Nat. Lab.). Boundary conditions
! available include periodic (lateral), free surface, zero slip at bed, specified basal 
! traction at bed, and specified basal yield stress at bed (all three of which are implemented
! through various verions of the specified traction b.c.)
! include macros for glide mask definitions
#include "glide_mask.inc"
#include "config.inc"

!TODO - Get rid of the globalIDs option.
!       Make it the default for Trilinos, else not used.

!GlobalIDs are for distributed TRILINOS variable IDs
#ifdef TRILINOS
#define globalIDs
#endif

!TODO:  In this module there are chunks of code that are used more than once, for Picard as well as JFNK.
!       It would be better to combine these chunks of code into subroutines that can be called
!        from multiple places in the code--or even better, to remove the extra chunks of code
!        if they are no longer needed.
! KJE looked into creating a generic initialization solver routine but most of init is passing 
! variables, so its not worth it IMHO

!***********************************************************************
module glam_strs2
!***********************************************************************

use iso_c_binding
use glimmer_paramets, only : dp
use glimmer_physcon,  only : gn, rhoi, rhoo, grav, pi, scyr

!TODO:     Remove scaling parameters tau0, vis0, evs0, etc. from this module.
!          Note: if thk0 = 1, then tau0 = rhoi*grav
          
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

  !SCALING - This corresponds to an effective min strain rate of 1.0d-20 s^(-1).
  real(dp), parameter :: effstrminsq = (1.0d-20 * tim0)**2
  real(dp) :: homotopy = 0.d0

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

  real(dp), dimension(:,:,:,:), allocatable :: ghostbvel

  ! variables for use in sparse matrix calculation
  real(dp), dimension(:), allocatable :: pcgval, rhsd, rhsx
  integer, dimension(:), allocatable :: pcgcol, pcgrow
  integer, dimension(2) :: pcgsize
  integer :: ct_nonzero  ! number of nonzero matrix entries

!*SFP* NOTE: these redefined here so that they are "in scope" and can avoid being passed as args
  integer :: whatsparse ! needed for putpgcg()
  integer :: nonlinear  ! flag for indicating type of nonlinar iteration (Picard vs. JFNK)

  logical, save :: inisoln = .false.      ! true only if a converged solution (velocity fields) exists

  real(dp) :: linearSolveTime = 0.d0
  real(dp) :: totalLinearSolveTime = 0.d0 ! total linear solve time

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

! WJS: The following routine doesn't compile on gnu; commenting it out for now
! subroutine dumpvels(name, uvel, vvel)
!     !JEFF routine to track the uvel and vvel calculations in Picard Iteration for debugging
!     !3/28/11
!     use parallel
!     implicit none

!     character(*) :: name
!     real(dp), dimension(:,:,:), intent(inout) :: uvel, vvel  ! horiz vel components: u(z), v(z)

!     if (distributed_execution()) then
!        if (this_rank == 0) then
!            write(*,*) name, "Proc 0 uvel & vvel (1,7:8,16:17)", uvel(1,7:8,16:17), vvel(1,7:8,16:17)
!        else
!            write(*,*) name, "Proc 1 uvel & vvel (1,7:8,0:1)", uvel(1,7:8,0:1), vvel(1,7:8,0:1)
!        endif
!     else
!        write(*,*) name, "Parallel uvel & vvel (1,5:6,15:16)", uvel(1,5:6,15:16), vvel(1,5:6,15:16)
!     endif 
! end subroutine dumpvels


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

    allocate(ghostbvel(2,3,ewn-1,nsn-1))        !! for saving the fictious basal vels at the bed !!

    ghostbvel(:,:,:,:) = 0.d0

    flwafact = 0.d0

     ! define constants used in various FD calculations associated with the 
     ! subroutine 'findcoefst'   
     call calccoeffsinit(upn, dew, dns)

    dups = (/ (sigma(up+1) - sigma(up), up=1,upn-1), 0.d0 /)

end subroutine glam_velo_init


!***********************************************************************

! This is the driver subroutine, called from subroutine glissade_velo_driver in
! module glissade_velo.F90. 

subroutine glam_velo_solver(ewn,      nsn,    upn,  &
                            dew,      dns,          &
                            sigma,    stagsigma,    &
                            thck,     usrf,         &
                            lsrf,     topg,         &
                            dthckdew, dthckdns,     &
                            dusrfdew, dusrfdns,     &
                            dlsrfdew, dlsrfdns,     &
                            stagthck, flwa,         &
                            bwat,     mintauf,      &
                            btraction,              &
                            umask,                  &
                            whichbabc,              &
                            whichefvs,              &
                            whichresid,             &
                            whichnonlinear,         &
                            whichsparse,            &
                            beta,                   &
                            uvel,     vvel,         &
                            uflx,     vflx,         &
                            efvs )

  use parallel
!!  use glimmer_horiz_bcs, only: horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, horiz_bcs_unstag_scalar
  use glimmer_paramets, only: GLC_DEBUG

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
  real(dp), dimension(:,:),   intent(in)  :: bwat                   ! thickness of basal water layer
  real(dp), dimension(:,:),   intent(in)  :: mintauf                ! till yield stress
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

  real(dp), dimension(:,:,:), intent(inout) :: uvel, vvel  ! horiz vel components: u(z), v(z)
  real(dp), dimension(:,:),   intent(out) :: uflx, vflx  ! horiz fluxs: u_bar*H, v_bar*H
  real(dp), dimension(:,:,:), intent(out) :: efvs        ! effective viscosity

  integer :: ew, ns, up     ! counters for horiz and vert do loops

  real(dp), parameter :: minres = 1.0d-4    ! assume vel fields converged below this resid 
  real(dp), parameter :: NL_tol = 1.0d-6    ! to have same criterion than with JFNK
  real(dp), save, dimension(2) :: resid     ! vector for storing u resid and v resid 

  integer, parameter :: cmax = 100                  ! max no. of iterations
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
!TODO - Can we get rid of whatsparse and use only whichsparse?
  whatsparse = whichsparse

  ! assign value for nonlinear iteration flag
  nonlinear = whichnonlinear

!TODO - Note: d2usrfdew2 and d2usrfdns2 are needed at all locally owned velocity points.
!       I am not sure where and why the upwind 2nd derivatives are computed.
!TODO MJH These 2nd derivatives are already calculated in subroutine geometry_derivs(model) in glide_thck.  
!These calls could either be deleted and just use those previous calculations, or possibly use that module here.  
!First it needs to be determined that they are making the same (or not) calculation!

  ! calc geometric 2nd deriv. for generic input variable 'ipvr', returns 'opvr'
  call geom2ders(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdew2, d2usrfdns2)
  call geom2ders(ewn, nsn, dew, dns, thck, stagthck, d2thckdew2, d2thckdns2)

  ! calc geometric 2nd cross-deriv. for generic input variable 'ipvr', returns 'opvr'
  call geom2derscros(ewn, nsn, dew, dns, thck, stagthck, d2thckdewdns)
  call geom2derscros(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdewdns)

  allocate(uindx(ewn-1,nsn-1))

  ! If a point from the 2d array 'mask' is associated with a non-zero ice thickness
  ! assign it a unique number. If not assign a zero.             
  uindx = indxvelostr(ewn, nsn, upn, umask,pcgsize(1))

!!!!!!!!!! Boundary conditions HACKS section !!!!!!!!!!!!!

!TODO - Remove this commented-out code if no longer needed.

!! A hack of the boundary condition mask needed for the Ross Ice Shelf exp.
!! The quick check of whether or not this is the Ross experiment is to look
!! at the domain size.
! if( ewn == 151 .and. nsn == 115 )then
!    call not_parallel(__FILE__, __LINE__)
!    do ns=1,nsn-1; do ew=1,ewn-1
!        if( umask(ew,ns) == 21 .or. umask(ew,ns) == 5 )then
!            umask(ew,ns) = 73
!        endif
!    end do; end do
! end if

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

     ! Now send this partition to Trilinos initialization routines
     call inittrilinos(20, mySize, myIndices, myX, myY, myZ, comm) 

     ! Set if need full solution vector returned or just owned portion

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
  L2norm = 1.d20

  ! intialize outer loop test vars
  outer_it_criterion = 1.d0
  outer_it_target = 0.d0

  if (main_task) then
     ! print some info to the screen to update on iteration progress
     print *, ' '
     print *, 'Running Payne/Price higher-order dynamics solver'
     print *, ' '
     if( whichresid == HO_RESID_L2NORM ) then
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
  outer_it_criterion = 1.d0
  outer_it_target = 0.d0

  do while ( outer_it_criterion >= outer_it_target .and. counter < cmax)    ! use L2 norm for resid calculation
 call t_startf("PICARD_in_iter")

  ! choose outer loop stopping criterion
  if( counter > 1 )then
    if( whichresid == HO_RESID_L2NORM )then
      outer_it_criterion = L2norm
      outer_it_target = NL_target
    else
      outer_it_criterion = maxval(resid)
      outer_it_target = minres
    end if
  else
    outer_it_criterion = 1.d10
    outer_it_target = 1.d-12
  end if

  ! WJS: commenting out the following block, because it leads to lots of extra files,
  ! which is undesirable even when GLC_DEBUG=.true.
  ! if (GLC_DEBUG) then
  !   !JEFF Debugging Output to see what differences in final vvel and tvel.
  !   write(loopnum,'(i3.3)') counter
  !   write(Looptime, '(i3.3)') overallloop
  !   loopnum = trim(loopnum)  ! Trying to get rid of spaces in name.
  !   Looptime = trim(Looptime)
  !   call distributed_print("uvela_ov"//Looptime//"_pic"//loopnum//"_tsk", uvel)

  !   call distributed_print("vvela_ov"//Looptime//"_pic"//loopnum//"_tsk", vvel)

  !   ! call dumpvels("Before findefvsstr", uvel, vvel)

  !   ! call distributed_print("preefvs_ov"//Looptime//"_pic"//loopnum//"_tsk", efvs)
  ! end if

 call t_startf("PICARD_findefvsstr")
    ! calc effective viscosity using previously calc vel. field
    call findefvsstr(ewn,  nsn,  upn,      &
                     stagsigma,  counter,  &
                     whichefvs,  efvs,     &
                     uvel,       vvel,     &
                     flwa,       thck,     &
                     dusrfdew,   dthckdew, &
                     dusrfdns,   dthckdns, &
                     umask)
 call t_stopf("PICARD_findefvsstr")

 call t_startf("PICARD_findcoefstr1")
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
                     mintauf,     flwa,           &
                     beta,        btraction,      &
                     bwat,        0 )
 call t_stopf("PICARD_findcoefstr1")

 call t_startf("PICARD_solver_pre1")
    ! put vels and coeffs from 3d arrays into sparse vector format
    call solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, vvel )
 call t_stopf("PICARD_solver_pre1")

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

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 call t_startf("PICARD_findcoefstr2")
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
                     mintauf,     flwa,           &
                     beta,        btraction,      &
                     bwat,        0 )
 call t_stopf("PICARD_findcoefstr2")

 call t_startf("PICARD_solver_pre2")
    ! put vels and coeffs from 3d arrays into sparse vector format
    call solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, uvel )
 call t_stopf("PICARD_solver_pre2")

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

 call t_startf("PICARD_findcoefstr3")
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
                     mintauf,     flwa,           &
                     beta,        btraction,      &
                     bwat,        1 )

    call findcoefstr(ewn,  nsn,   upn,            &
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
                     mintauf,     flwa,           &
                     beta,        btraction,      &
                     bwat,        1 )

 call t_stopf("PICARD_findcoefstr3")

    ! apply unstable manifold correction to converged velocities

 call t_startf("PICARD_mindcrsh")

    call mindcrshstr(1,whichresid,uvel,counter,resid(1))

    vvel = tvel
    call mindcrshstr(2,whichresid,vvel,counter,resid(2))

 call t_stopf("PICARD_mindcrsh")

!HALO - I'm pretty sure these updates *are* needed.
!       
 call t_startf("PICARD_halo_upds")
    ! coordinate halos for updated uvel and vvel
    call staggered_parallel_halo(uvel)
!    call horiz_bcs_stag_vector_ew(uvel)
    call staggered_parallel_halo(vvel)
!    call horiz_bcs_stag_vector_ns(vvel)
 call t_stopf("PICARD_halo_upds")

    !call dumpvels("After mindcrsh", uvel, vvel)

    if (this_rank == 0) then

        !TODO - Does this comment still apply, or is parallel_single defunct?

        ! Can't use main_task flag because main_task is true for all processors in case of parallel_single
        ! output the iteration status: iteration number, max residual, and location of max residual
        ! (send output to the screen or to the log file, per whichever line is commented out) 

        if( whichresid == HO_RESID_L2NORM ) then
            print '(i4,3g20.6)', counter, L2norm, NL_target    ! Output when using L2norm for convergence
            !print '(a,i4,3g20.6)', "sup-norm uvel, vvel=", counter, resid(1), resid(2), minres
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

!TODO - I don't think uflx and vflx are needed; they are not used by the remapping subroutine.

  do ns = 1+staggered_lhalo, size(umask,2)-staggered_uhalo
      do ew = 1+staggered_lhalo, size(umask,1)-staggered_uhalo
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

!TODO - Do we need halo updates for btraction and efvs?
!       I think we don't need an update for efvs, because it is already computed in a layer of halo cells.
!       And I think we don't need an update for btraction, because it is computed in bodyset for all
!        locally owned velocity points.

  call parallel_halo(efvs)
!  call horiz_bcs_unstag_scalar(efvs)

  call staggered_parallel_halo(btraction)

  !TODO - Pretty sure we don't need these updates; uflx and vflx are not used elsewhere.
  call staggered_parallel_halo(uflx)
!  call horiz_bcs_stag_vector_ew(uflx)

  call staggered_parallel_halo(vflx)
!  call horiz_bcs_stag_vector_ns(vflx)

  if (GLC_DEBUG) then
     !JEFF Debugging Output to see what differences in final vvel and tvel.
     ! write(CurrTimeLoopStr, '(i3.3)') CurrTimeLoop
     ! call distributed_print("uvel_post_ov"//CurrTimeLoopStr//"_tsk", uvel)
     !
     ! call distributed_print("vvel_post_ov"//CurrTimeLoopStr//"_tsk", vvel)
  end if

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
  deallocate(uk_1, vk_1, F, g_flag)

  !JEFF Debugging output
  overallloop = overallloop + 1
 call t_stopf("PICARD_post")

  return

end subroutine glam_velo_solver

!***********************************************************************

subroutine JFNK_velo_solver  (model,umask)

  use parallel
!!  use glimmer_horiz_bcs, only: horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns,  horiz_bcs_unstag_scalar
  use glimmer_paramets, only: GLC_DEBUG

  use iso_c_binding 
  use glide_types, only : glide_global_type

  implicit none

  type(glide_global_type) ,target, intent(inout) :: model

  !TODO - Can we make the mask intent in?

  integer, dimension(:,:),   intent(inout)  :: umask  !*SFP* replaces the prev., internally calc. mask
                                                      ! ... 'inout' status allows for a minor alteration
                                                      ! to cism defined mask, which don't necessarily 
                                                      ! associate all/any boundaries as a unique mask value.

  type(glide_global_type) ,pointer :: fptr=>NULL()
  type(c_ptr)                 :: c_ptr_to_object

  integer(c_int) :: xk_size
  real(dp), dimension(:), allocatable :: xk_1
  integer ,dimension(:) ,allocatable :: gx_flag

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
  real(dp), dimension(:,:)   ,pointer :: mintauf
  real(dp), dimension(:,:,:) ,pointer :: btraction            ! consistent basal traction array
  
!TODO - Anything to update here?
  !*SFP* This is the betasquared field from CISM (externally specified), and should eventually
  ! take the place of the subroutine 'calcbetasquared' below (for now, using this value instead
  ! will simply be included as another option within that subroutine) 
  real(dp), dimension(:,:)  ,pointer :: beta 

  real(dp), dimension(:,:)  ,pointer :: bwat 

  integer :: whichbabc
  integer :: whichefvs
  integer :: whichresid
  integer :: whichsparse
  integer :: whichnonlinear

!TODO - Should the following be passed out explicitly?
! intent(out)
  real(dp), dimension(:,:,:) ,pointer :: uvel, vvel
  real(dp), dimension(:,:)   ,pointer :: uflx, vflx
  real(dp), dimension(:,:,:) ,pointer :: efvs

  integer :: ew, ns, up, nele
  real(dp), parameter :: NL_tol = 1.0d-6

! currently needed to assess whether basal traction is updated after each nonlinear iteration
!  integer :: k 
!TODO: "k" is not needed in order to calculate basal traction; note that new subroutine calls
! at lines 1175 below pass in a dummy value for this variable. In the long run, we can likely remove
! this argument altogether - it was originally passed in to aid in stabilization
! of the ice shelf boundary conditions but may no longer be needed (grep for the variable "cc" within
! the subroutine "bodyset" to see where it is currently used)

  character(len=100) :: message

!*SFP* needed to incorporate generic wrapper to solver
  type(sparse_matrix_type) :: matrixA, matrixC, matrixtp, matrixAuv, matrixAvu
  real(dp) :: L2norm

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
  mintauf => model%basalproc%mintauf(:,:)
  btraction => model%velocity%btraction(:,:,:)
  whichbabc = model%options%which_ho_babc
  whichefvs = model%options%which_ho_efvs
  whichresid = model%options%which_ho_resid
  whichsparse = model%options%which_ho_sparse
  whichnonlinear = model%options%which_ho_nonlinear
  beta => model%velocity%beta(:,:)
  bwat => model%temper%bwat(:,:)
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
!TODO MJH: can we put these derivative calculations in the diagnostic solve part where the other derivatives are calculated?

  ! *SFP* geometric 1st deriv. for generic input variable 'ipvr',
  !      output as 'opvr' (includes 'upwinding' for boundary values)
  call geom2ders(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdew2, d2usrfdns2)
  call geom2ders(ewn, nsn, dew, dns, thck, stagthck, d2thckdew2, d2thckdns2)

  ! *SFP* geometric (2nd) cross-deriv. for generic input variable 'ipvr', output as 'opvr'
  call geom2derscros(ewn, nsn, dew, dns, thck, stagthck, d2thckdewdns)
  call geom2derscros(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdewdns)

!TODO - Do these derivatives have to go in the model derived type and the residual object?
  model%geomderv%d2thckdew2 = d2thckdew2
  model%geomderv%d2thckdns2 = d2thckdns2
  model%geomderv%d2usrfdew2 = d2usrfdew2
  model%geomderv%d2usrfdns2 = d2usrfdns2

  ! *SFP* make a 2d array identifying if the associated point has zero thickness,
  !      has non-zero thickness and is interior, or has non-zero thickness
  !      and is along a boundary

  !*SFP* This subroutine has been altered from its original form (was a function, still included
  ! below w/ subroutine but commented out) to allow for a tweak to the CISM calculated mask (adds
  ! in an unique number for ANY arbitrary boundary, be it land, water, or simply at the edge of
  ! the calculation domain). 

  allocate(uindx(ewn-1,nsn-1))

  ! *SFP* if a point from the 2d array 'mask' is associated with non-zero ice thickness,
  !      either a boundary or interior point, give it a unique number. If not, give it a zero			 
  uindx = indxvelostr(ewn, nsn, upn, umask, pcgsize(1))

  L2norm = 1.0d20
 
  ! *SFP* an initial guess at the size of the sparse matrix
  pcgsize(2) = pcgsize(1) * 20

  ! Structure to become NOX implementation for JFNK solve
  xk_size=2*pcgsize(1)

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

#ifdef TRILINOS
  if (whatsparse == STANDALONE_TRILINOS_SOLVER) then
     if (main_task) write(*,*) "Using GlobalIDs..."
	 ! JEFF: Define myIndices in terms of globalIDs
     allocate(myIndices(pcgsize(1)))  ! myIndices is an integer vector with a unique ID for each layer for ice grid points
     allocate(myX(pcgsize(1))) ! Coordinates of nodes, used by ML preconditioner
     allocate(myY(pcgsize(1))) 
     allocate(myZ(pcgsize(1))) 
     call distributed_create_partition(ewn, nsn, (upn + 2) , uindx, pcgsize(1), myIndices, myX, myY, myZ)  ! Uses uindx mask to determine ice grid points.
     mySize = pcgsize(1)  ! Set variable for inittrilinos

     if (GLC_DEBUG) then
        write(*,*) "GlobalIDs myIndices..."
        write(*,*) "pcgsize = ", pcgsize(1)
        write(*,*) "myIndices = ", myIndices
        !call parallel_stop(__FILE__, __LINE__)
     end if

     call inittrilinos(25, mySize, myIndices, myX, myY, myZ, comm)   !re: Why 25 not 20 for PIC? needed the mem space

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

  allocate( xk_1(2*pcgsize(1)), gx_flag(2*pcgsize(1)) )

  ! *SFP* allocate space matrix variables
  allocate (pcgrow(pcgsize(2)),pcgcol(pcgsize(2)), rhsd(pcgsize(1)), rhsx(2*pcgsize(1)), &
            pcgval(pcgsize(2)))
  allocate(matrixA%row(pcgsize(2)), matrixA%col(pcgsize(2)), &
            matrixA%val(pcgsize(2))) 
  allocate(matrixC%row(pcgsize(2)), matrixC%col(pcgsize(2)), &
            matrixC%val(pcgsize(2)))
  allocate(matrixtp%row(pcgsize(2)), matrixtp%col(pcgsize(2)), &
            matrixtp%val(pcgsize(2)))

  allocate(model%solver_data%ui(ewn-1,nsn-1) )
  allocate(model%solver_data%um(ewn-1,nsn-1) ) 
  allocate(model%solver_data%d2thckcross(ewn-1,nsn-1) )
  allocate(model%solver_data%d2usrfcross(ewn-1,nsn-1) ) 
  allocate(model%solver_data%gxf( 2*pcgsize(1) ) )
  
  call assign_resid(model, uindx, umask, d2thckdewdns, d2usrfdewdns, &
                       pcgsize, gx_flag, matrixA, matrixC, L2norm, ewn, nsn)

  fptr => model
  c_ptr_to_object = c_loc(fptr)

  call ghost_preprocess_jfnk( ewn, nsn, upn, uindx, ughost, vghost, &
                         xk_1, uvel, vvel, gx_flag, pcgsize(1)) ! jfl_20100430

if (main_task) then
  print *, ' '
  print *, 'Running Payne/Price higher-order dynamics with JFNK solver' 
end if

 call t_stopf("JFNK_pre")

#ifdef TRILINOS 

!==============================================================================
! Newton loop Using Trilinos NOX. Solves F(x) = 0 for x where x = [v, u] and
!                                                       F = [Fv(u,v), Fu(u,v)] 
!==============================================================================
 
 call t_startf("JFNK_noxinit")
  call noxinit(xk_size, xk_1, comm, c_ptr_to_object)
 call t_stopf("JFNK_noxinit")

 call t_startf("JFNK_noxsolve")
  call noxsolve(xk_size, xk_1, c_ptr_to_object)
 call t_stopf("JFNK_noxsolve")

 call t_startf("JFNK_noxfinish")
  call noxfinish()
 call t_stopf("JFNK_noxfinish")

!TODO remove since not needed?
! k = 0

#else

!TODO Is the slapsolve code still used? 

!==============================================================================
! SLAP JFNK loop: calculate F(u^k-1,v^k-1)
!==============================================================================

 call t_startf("JFNK_SLAP")
 call slapsolve(xk_1, xk_size, c_ptr_to_object, NL_tol, pcgsize)
 call t_stopf("JFNK_SLAP")

! k = 1

#endif  

 call t_startf("JFNK_post")

! need to update these values from fptr%uvel,vvel,stagthck etc
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

! KJE this is now outside the loop of both JFNK methods (and has been for while) 
! appears to be redundant, but leaving commented for a while in case an unknown issues pops up
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
                     mintauf,     flwa,           &
                     beta,        btraction,      &
                     bwat,        1 )

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
                     mintauf,     flwa,           &
                     beta,        btraction,      &
                     bwat,        1 )

  inisoln = .true.

  if (GLC_DEBUG) then
     print*,"Solution vector norm after JFNK = " ,sqrt(DOT_PRODUCT(xk_1,xk_1))
  end if

!TODO - The remaining code in this subroutine is cut and pasted from above.
!       Can we encapsulate this repeated code in a subroutine?

!TODO - I don't think uflx and vflux are needed.

!LOOP - Locally owned velocity points
  do ns = 1+staggered_lhalo, size(umask,2)-staggered_uhalo
      do ew = 1+staggered_lhalo, size(umask,1)-staggered_uhalo
      ! *SFP* calc. fluxes from converged vel. fields (for input to thickness evolution subroutine)
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

  ! *SFP* de-allocation of sparse matrix solution variables 
  deallocate(uindx)
  deallocate(pcgval,pcgrow,pcgcol,rhsd, rhsx)
  deallocate(matrixA%row, matrixA%col, matrixA%val)
  deallocate(matrixC%row, matrixC%col, matrixC%val)
  deallocate(matrixtp%row, matrixtp%col, matrixtp%val)
  deallocate(gx_flag )
  deallocate(model%solver_data%ui)
  deallocate(model%solver_data%um)
  deallocate(model%solver_data%d2thckcross)
  deallocate(model%solver_data%d2usrfcross)
  deallocate(model%solver_data%gxf)
  
!TODO - Not sure whether these are needed.  Where does JFNK do its parallel halo updates for uvel, vvel?

 !PW following are needed for glam_velo_fordsiapstr - putting here until can be convinced
 !   that they are not needed (or that they should be delayed until later)
  call staggered_parallel_halo(uvel)
!  call horiz_bcs_stag_vector_ew(uvel)
  call staggered_parallel_halo(vvel)
!  call horiz_bcs_stag_vector_ns(vvel)

!TODO - Not sure we need these two updates
!       I think we do not need an update for efvs, because it is already computed in a layer of halo cells.
!       And I think we don't need an update for btraction, because it is computed in bodyset for all
!        locally owned velocity points.

  call parallel_halo(efvs)
!  call horiz_bcs_unstag_scalar(efvs)

  call staggered_parallel_halo(btraction)

!TODO - Probably do not need these two updates
  call staggered_parallel_halo(uflx)
!  call horiz_bcs_stag_vector_ew(uflx)
  call staggered_parallel_halo(vflx)
!  call horiz_bcs_stag_vector_ns(vflx)

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

  do ns = 1+staggered_lhalo, size(mask,2)-staggered_uhalo
     do ew = 1+staggered_lhalo, size(mask,1)-staggered_uhalo
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

  ! calculate the effective viscosity    

  use parallel
  use glimmer_paramets, only: GLC_DEBUG
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

     ! effstrminsq = (1.0d-20 * tim0)**2

     if (GLC_DEBUG) then

        !  if (main_task) then
        !    print *, 'nsn=', nsn
        !    print *, 'ewn=', ewn
        !    print *, 'uvel shape =', shape(uvel)
        !    print *, 'vvel shape =', shape(vvel)
        !    print *, 'thck shape =', shape(thck)
        !    print *, 'efvs shape =', shape(efvs)
        !    print *, 'flwafact shape =', shape(flwafact)
        !  endif

     end if

!TODO - If we are not supporting glam_strs2 together with the old Glimmer temperature routines,
!       then we can assume that temp and flwa live on the staggered vertical grid.

     if (size(flwa,1)==upn-1) then   ! temperature and flwa live on staggered vertical grid

        !Note: To avoid parallel halo calls for efvs within glam_strs2, we need to compute efvs in one layer of halo cells
        !      surrounding the locally owned velocity cells.

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

  case(HO_EFVS_CONSTANT)       ! set the eff visc to a constant value

   do ns = 2,nsn-1
     do ew = 2,ewn-1
        if (thck(ew,ns) > 0.d0) then
           ! Steve recommends 10^6 to 10^7 Pa yr
           efvs(1:upn-1,ew,ns) = 1.d7  * scyr/tim0 / tau0    ! tau0 = rhoi*grav*thk0   
        else        
           efvs(:,ew,ns) = effstrminsq ! if the point is associated w/ no ice, set to min value
        endif
     enddo
   enddo

  case(HO_EFVS_FLOWFACT)    ! set the eff visc to a value based on the flow rate factor 

!   *SFP* changed default setting for linear viscosity so that the value of the rate
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

  case(HO_EFVS_NONLINEAR)      ! calculate eff. visc. using eff. strain rate

!TODO - This code may not work correctly if nhalo = 1.  
!       In that case we would need a halo update of efvs to make sure we have the correct value
!        in all neighbors of locally owned velocity cells.
 
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

  end select

! JEFF Halo does NOT verify, because used a staggered array to hold unstaggered data.  
! The unstaggered data fits because remove two rows and columns of data.
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

  integer :: k

!WHL - Rewriting to get code to run on Mac with array bounds checking
!!  vertideriv(1:upn-1) = dupm * (varb(2:upn) - varb(1:upn-1)) / thck

  do k = 1, upn-1 
     vertideriv(k) = dupm(k) * (varb(k+1) - varb(k)) / thck
  enddo

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

!! WHL - Testing whether this function will work for single-processor parallel runs
!!       with solvers other than trilinos

function getlocationarray(ewn, nsn, upn, mask, indxmask, return_global_IDs)
!function getlocationarray(ewn, nsn, upn, mask, indxmask)

  use parallel

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:), intent(in) :: mask
  integer, dimension(:,:), intent(in) :: indxmask
  logical, intent(in), optional :: return_global_IDs

  integer, dimension(ewn,nsn,2) :: getlocationarray

  logical :: return_globalIDs  ! set to return_global_IDs, if present 

  integer :: ew, ns
  integer, dimension(ewn,nsn) :: temparray
  integer :: cumsum

  if (present(return_global_IDs)) then
     if (return_global_IDs) then
        return_globalIDs = .true.
     else
        return_globalIDs = .false.
     endif
  else
     return_globalIDs = .true.
  endif

!TODO - Make this if which_ho_sparse = 4 instead (or ifdef Trilinos?)
#ifdef globalIDs
  ! Returns in (:,:,1) the global ID bases for each grid point, including 
  ! halos and those without ice.
  ! Since the code checks elsewhere whether ice occurs at a given grid point, 
  ! this information is not encoded here. For the local indices (see below)
  ! the mask information is used since ice-free grid points are not indexed
  ! locally

!WHL - debug
!   print*, 'In getlocationarray, ifdef globalIDs' 
!   print*, 'return_globalIDs =', return_globalIDs

!LOOP TODO - Not sure if these loops are correct.
!       Is the input mask on the scalar (ice) grid? 
!SFP: Need to check indices here - getlocationarray should exist on the velocity grid, not the thickness (scalar) grid

!WHL - added this conditional

  if (return_globalIDs) then

     do ns = 1,nsn
        do ew = 1,ewn
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

     ! Step through indxmask, but exclude halo

     do ns = 1+staggered_lhalo, size(indxmask,2)-staggered_uhalo
        do ew = 1+staggered_lhalo, size(indxmask,1)-staggered_uhalo
           if ( indxmask(ew,ns) /= 0 ) then
              getlocationarray(ew,ns,2) = (indxmask(ew,ns) - 1) * (upn+2) + 1
           endif
        end do
     end do

!TODO - Clean this up, so we always use this procedure when solving without Trilinos.

  else  ! use the procedure below under #else

     ! initialize to zero
     cumsum = 0
     temparray = 0
     getlocationarray = 0

     do ns=1+staggered_lhalo, size(mask,2)-staggered_uhalo
        do ew=1+staggered_lhalo, size(mask,1)-staggered_uhalo
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

     getlocationarray(:,:,1) = ( getlocationarray(:,:,1) + 1 ) - temparray(:,:)
     getlocationarray(:,:,2) = getlocationarray(:,:,1)

  endif   ! return_globalIDs

#else

  ! initialize to zero
  cumsum = 0
  temparray = 0
  getlocationarray = 0

  do ns=1+staggered_lhalo, size(mask,2)-staggered_uhalo
    do ew=1+staggered_lhalo, size(mask,1)-staggered_uhalo
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

  getlocationarray(:,:,1) = ( getlocationarray(:,:,1) + 1 ) - temparray(:,:)
  getlocationarray(:,:,2) = getlocationarray(:,:,1)

#endif

  return

end function getlocationarray

!***********************************************************************

!TODO - Remove this function?  I don't think it is called anywhere.

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

  pcgsize(2) = ct_nonzero - 1

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

!LOOP TODO - Are loop bounds OK? Since this is for the serial SLAP solver, I think so.

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

  pcgsize(2) = ct_nonzero - 1

  matrix%order = pcgsize(1)
  matrix%nonzeros = pcgsize(2)
  matrix%symmetric = .false.

  matrix%row = pcgrow
  matrix%col = pcgcol
  matrix%val = pcgval

  ! Initial estimate for vel. field; take from 3d array and put into
  ! the format of a solution vector.

  do ns = 1+staggered_lhalo, size(uindx,2)-staggered_uhalo
   do ew = 1+staggered_lhalo, size(uindx,1)-staggered_uhalo
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

  do ns = 1+staggered_lhalo, size(uindx,2)-staggered_uhalo
      do ew = 1+staggered_lhalo, size(uindx,1)-staggered_uhalo
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

   do ns = 1+staggered_lhalo, size(uindx,2)-staggered_uhalo
       do ew = 1+staggered_lhalo, size(uindx,1)-staggered_uhalo
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

subroutine resvect_postprocess_jfnk( ewn, nsn, upn, uindx, pcg1, answrapped, ansunwrappedv, &
                                    ansunwrappedu, ansunwrappedmag )   
! Unwrap the jfnk residual vector from the solution vector and place into a 3d array.
   use parallel

   implicit none

   integer :: pcg1
   integer, intent(in) :: ewn, nsn, upn
   integer, dimension(:,:), intent(in) :: uindx
   real(dp), dimension(:), intent(in) :: answrapped
   real(dp), dimension(upn,ewn-1,nsn-1), intent(out), optional :: ansunwrappedv, ansunwrappedu, ansunwrappedmag

   integer, dimension(2) :: loc
   integer :: ew, ns

   do ns = 1+staggered_lhalo, size(uindx,2)-staggered_uhalo
       do ew = 1+staggered_lhalo, size(uindx,1)-staggered_uhalo
           if (uindx(ew,ns) /= 0) then
             loc = getlocrange(upn, uindx(ew,ns))
             ansunwrappedv(:,ew,ns) = answrapped(loc(1):loc(2))
             ansunwrappedu(:,ew,ns) = answrapped(pcg1+loc(1):pcg1+loc(2))
           else
             ansunwrappedv(:,ew,ns) = 0.d0
             ansunwrappedu(:,ew,ns) = 0.d0
           end if
       end do
   end do

   ansunwrappedmag = dsqrt( ansunwrappedu**2.d0 + ansunwrappedv**2.d0 )

end subroutine resvect_postprocess_jfnk

!***********************************************************************

subroutine form_matrix( matrix ) ! for JFNK solver

  ! Puts sparse matrix variables in SLAP triad format into "matrix" 
  ! derived type. Similar to solver_preprocess but does not form answer vector

  implicit none

!  integer, intent(in) :: ewn, nsn, upn
  type(sparse_matrix_type), intent(inout) :: matrix

  pcgsize(2) = ct_nonzero - 1

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
      expo      = 2.d0

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
  type(glide_global_type) ,pointer        :: fptr=>NULL()
  type(c_ptr) ,intent(inout)         :: c_ptr_to_object

  integer :: nu1, nu2, whichsparse
  integer :: iter
  type(sparse_matrix_type) :: matrixA, matrixC
  real(dp), dimension(xk_size) :: wk1
  real(dp), dimension(xk_size) :: wk2
  real(dp), allocatable, dimension(:) :: answer, vectp
  real(dp) :: err

  call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr= model

  matrixA = fptr%solver_data%matrixA
  matrixC = fptr%solver_data%matrixC
  whichsparse = fptr%options%which_ho_sparse
  pcgsize = fptr%solver_data%pcgsize
           
  nu1 = pcgsize(1)
  nu2 = 2*pcgsize(1)
  allocate ( answer(nu1) )
  allocate ( vectp(nu1) )
  wk1  = wk1_nox

! ID as a test
!  wk2_nox  = wk1

! precondition v component 
       
      answer = 0.d0 ! initial guess
      vectp(:) = wk1(1:nu1) ! rhs for precond v
 call t_startf("nox_precond_v")
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
 call t_stopf("nox_precond_v")
      wk2(1:nu1) = answer(:)

! precondition u component 
       
      answer = 0.d0 ! initial guess
      vectp(:) = wk1(nu1+1:nu2) ! rhs for precond u
 call t_startf("nox_precond_u")
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
 call t_stopf("nox_precond_u")
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
  if (esm_factor > 1.0d-10) then
    homotopy = 10.0**( esm_factor - 9.0 )
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
  use glide_types ,only : glide_global_type
  use parallel
!!  use glimmer_horiz_bcs, only: horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns

  implicit none

   integer(c_int) ,intent(in) ,value  :: xk_size
! ispert is 0 for base calculations, 1 for perturbed calculations
   integer(c_int) ,intent(in) ,value  :: ispert 
   real(c_double)  ,intent(in)        :: xtp(xk_size)
   real(c_double)  ,intent(out)       :: F(xk_size)
   type(glide_global_type) ,pointer        :: fptr=>NULL()
   type(c_ptr) ,intent(inout)         :: c_ptr_to_object

  integer :: ewn, nsn, upn, counter, whichbabc, whichefvs, i
  integer  ,dimension(2)   :: pcgsize
  integer  ,dimension(:) ,allocatable :: gxf ! 0 :reg cell
  integer  ,dimension(:,:) ,allocatable :: ui, um
  real(dp) :: dew, dns
  real(dp), dimension(:)  ,pointer :: sigma, stagsigma
  real(dp), dimension(:,:) ,pointer :: thck, dusrfdew, dthckdew, dusrfdns, dthckdns, &
                                         dlsrfdew, dlsrfdns, stagthck, lsrf, topg, mintauf, beta, bwat
  real(dp), dimension(:,:) ,pointer ::  d2usrfdew2, d2thckdew2, d2usrfdns2, d2thckdns2
  real(dp), dimension(:,:,:) ,pointer :: efvs, btraction
  real(dp), dimension(:,:,:) ,pointer :: uvel, vvel, flwa
!  real(dp), dimension(:,:,:) ,pointer :: ures, vres, magres  !! used for output of residual fields  
  type(sparse_matrix_type) :: matrixA, matrixC
  real(dp), dimension(:) ,allocatable :: vectx
  real(dp), dimension(:) ,allocatable :: vectp

  real(dp) :: L2square
!  real(dp), intent(inout):: L2norm

!  real(dp) :: Ft(xk_size)   !! used for output of residual fields (ures,vres,magres) 
                            ! storage for "F" vector when using F to output residual fields for plotting (because it 
                            ! during the process of calc. the resid. and unwrapping it and we don't want to alter the
                            ! actual F vector)
  real(dp) :: L2norm

 call t_startf("Calc_F")
  call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr= model
           
  ewn = fptr%general%ewn
  nsn = fptr%general%nsn
  upn = fptr%general%upn
  whichbabc = fptr%options%which_ho_babc
  whichefvs = fptr%options%which_ho_efvs
  dew = fptr%numerics%dew
  dns = fptr%numerics%dns
  sigma => fptr%numerics%sigma(:)
  stagsigma => fptr%numerics%stagsigma(:)
  thck => fptr%geometry%thck(:,:)
  lsrf => fptr%geometry%lsrf(:,:)
  topg => fptr%geometry%topg (:,:)
  stagthck => fptr%geomderv%stagthck(:,:)
  dthckdew => fptr%geomderv%dthckdew(:,:)
  dthckdns => fptr%geomderv%dthckdns(:,:)
  dusrfdew => fptr%geomderv%dusrfdew(:,:)
  dusrfdns => fptr%geomderv%dusrfdns(:,:)
  dlsrfdew => fptr%geomderv%dlsrfdew(:,:)
  dlsrfdns => fptr%geomderv%dlsrfdns(:,:)
  d2thckdew2 => fptr%geomderv%d2thckdew2(:,:)
  d2thckdns2 => fptr%geomderv%d2thckdns2(:,:)
  d2usrfdew2 => fptr%geomderv%d2usrfdew2(:,:)
  d2usrfdns2 => fptr%geomderv%d2usrfdns2(:,:)
  mintauf => fptr%basalproc%mintauf(:,:)
  beta => fptr%velocity%beta(:,:)
  bwat => fptr%temper%bwat(:,:)
!intent (inout) terms
  btraction => fptr%velocity%btraction(:,:,:)
  flwa => fptr%temper%flwa(:,:,:)
  efvs => fptr%stress%efvs(:,:,:)
  uvel => fptr%velocity%uvel(:,:,:)
  vvel => fptr%velocity%vvel(:,:,:)
!  ures => fptr%velocity%ures(:,:,:)         !! used for output of residual fields
!  vres => fptr%velocity%vres(:,:,:)         !! used for output of residual fields
!  magres => fptr%velocity%magres(:,:,:)     !! used for output of residual fields
  L2norm = fptr%solver_data%L2norm

  allocate( ui(ewn-1,nsn-1), um(ewn-1,nsn-1) )
  ui= fptr%solver_data%ui
  um = fptr%solver_data%um

  pcgsize = fptr%solver_data%pcgsize
  allocate( gxf(2*pcgsize(1)) )

  gxf = fptr%solver_data%gxf
! temporary to test JFNK -  need to take out
  counter = 1

  d2usrfdewdns = fptr%solver_data%d2usrfcross
  d2thckdewdns = fptr%solver_data%d2thckcross

  matrixA = fptr%solver_data%matrixA
  matrixC = fptr%solver_data%matrixC
  allocate( vectp( pcgsize(1)) )
  allocate( vectx(2*pcgsize(1)) )

  call solver_postprocess_jfnk( ewn, nsn, upn, ui, &
                                xtp, vvel, uvel, ghostbvel, pcgsize(1) )

    ! coordinate halos for updated uvel and vvel
 call t_startf("Calc_F_uvhalo_upd")
    call staggered_parallel_halo(uvel)
!    call horiz_bcs_stag_vector_ew(uvel)
    call staggered_parallel_halo(vvel)
!    call horiz_bcs_stag_vector_ns(vvel)
 call t_stopf("Calc_F_uvhalo_upd")

 call t_startf("Calc_F_findefvsstr")
    call findefvsstr(ewn,  nsn,  upn,       &
                     stagsigma,  counter,  &
                     whichefvs,  efvs,     &
                     uvel,       vvel,     &
                     flwa,       thck,     &
                     dusrfdew,   dthckdew, &
                     dusrfdns,   dthckdns, &
                     um)
 call t_stopf("Calc_F_findefvsstr")

!==============================================================================
! jfl 20100412: residual for v comp: Fv= A(utp,vtp)vtp - b(utp,vtp)  
!==============================================================================

    ! *SFP* calculation of coeff. for stress balance calc. 
 call t_startf("Calc_F_findcoefstr1")
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
                     mintauf,     flwa,           &
                     beta,        btraction,      &
                     bwat,        0 )

 call t_stopf("Calc_F_findcoefstr1")

    rhsx(1:pcgsize(1)) = rhsd ! Fv

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
 call t_startf("Calc_F_form_matrix1")
      call form_matrix ( matrixA ) ! to get A(utp,vtp)
 call t_stopf("Calc_F_form_matrix1")
#ifdef TRILINOS
    else
      if (ispert == 0) then
 call t_startf("Calc_F_savetrilinos1")
        call savetrilinosmatrix(0); 
 call t_stopf("Calc_F_savetrilinos1")
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

 call t_startf("Calc_F_findcoefstr2")

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
                     mintauf,     flwa,           &
                     beta,        btraction,      &
                     bwat,        0 )

 call t_stopf("Calc_F_findcoefstr2")

    rhsx(pcgsize(1)+1:2*pcgsize(1)) = rhsd ! Fv

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
 call t_startf("Calc_F_form_matrix2")
      call form_matrix ( matrixC ) ! to get C(utp,vtp)
 call t_stopf("Calc_F_form_matrix2")
#ifdef TRILINOS
    else
      if (ispert == 0) then
 call t_startf("Calc_F_savetrilinos2")
        call savetrilinosmatrix(1); 
 call t_stopf("Calc_F_savetrilinos2")
      endif
#endif
    end if
    
    vectp(1:pcgsize(1)) = xtp(pcgsize(1)+1:2*pcgsize(1))

 call t_startf("Calc_F_res_vect")
    call res_vect(matrixC, vectp, rhsd, pcgsize(1), gxf, L2square, whatsparse)
 call t_stopf("Calc_F_res_vect")
    L2norm = sqrt(L2norm + L2square)

    F(pcgsize(1)+1:2*pcgsize(1)) = vectp(1:pcgsize(1)) 

!TODO: Older code that doesn't seem to be needed anymore? Note that "res_vect_jfnk" sits inside of "res_vect.F90"
! and should NOT be removed. It is still useful, as per below where it can be used during debug/perf. testing to
! output the 3d residual fields.
!
!   vectx = xtp 
!   call res_vect_jfnk(matrixA, matrixC, vectx, rhsx, pcgsize(1), 2*pcgsize(1), gxf, L2square, whatsparse)
!   L2norm = L2square
!   F = vectx 

    call solver_postprocess_jfnk( ewn, nsn, upn, ui, xtp, vvel, uvel, ghostbvel, pcgsize(1) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This section used and active only if / for output of residual fields !! 
!  Ft = F  !! need a temp variable to pass in here because "res_vect_jfnk" alters the value of "F"
!  call res_vect_jfnk(matrixA, matrixC, Ft, rhsx, pcgsize(1), 2*pcgsize(1), gxf, L2square, whatsparse)
!  call resvect_postprocess_jfnk( ewn, nsn, upn, ui, pcgsize(1), Ft, vres, ures, magres )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fptr%velocity%btraction => btraction(:,:,:)
  fptr%temper%flwa => flwa(:,:,:)
  fptr%stress%efvs => efvs(:,:,:)
  fptr%velocity%uvel => uvel(:,:,:)
  fptr%velocity%vvel => vvel(:,:,:)

!  fptr%velocity%ures => ures(:,:,:)     !! used for output of residual fields 
!  fptr%velocity%vres => vres(:,:,:)     !! used for output of residual fields
!  fptr%velocity%magres => magres(:,:,:) !! used for output of residual fields

  fptr%solver_data%L2norm = L2norm
  fptr%solver_data%matrixA = matrixA
  fptr%solver_data%matrixC = matrixC
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

  do ns = 1+staggered_lhalo, size(uindx,2)-staggered_uhalo
   do ew = 1+staggered_lhalo, size(uindx,1)-staggered_uhalo
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

   do ns = 1+staggered_lhalo, size(uindx,2)-staggered_uhalo
    do ew = 1+staggered_lhalo, size(uindx,1)-staggered_uhalo
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

  do ns = 1+staggered_lhalo, size(uindx,2)-staggered_uhalo
      do ew = 1+staggered_lhalo, size(uindx,1)-staggered_uhalo
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

   do ns = 1+staggered_lhalo, size(uindx,2)-staggered_uhalo
       do ew = 1+staggered_lhalo, size(uindx,1)-staggered_uhalo
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
  use glimmer_paramets, only: GLC_DEBUG

  implicit none

  real(dp), intent(inout), dimension(:,:,:) :: vel
  integer, intent(in) :: counter, pt, whichresid

  real(dp), intent(out) :: resid

!TODO - critlimit is never used
!TODO - SCALING - Does 'small' need a velocity scale factor?
  real(dp), parameter :: ssthres = 5.d0 * pi / 6.d0, &
                         critlimit = 10.d0 / (scyr * vel0), &
                         small = 1.0d-16

  real(dp) :: temp_vel

  integer, dimension(2), save :: new = 1, old = 2
  !JEFF integer :: locat(3)
  integer ew, ns, nr

  integer, dimension(size(vel,1),size(vel,2),size(vel,3)) :: vel_ne_0
  real(dp) :: sum_vel_ne_0

!WHL - debug (to print out intermediate terms in equations)
!!  real(dp) :: alpha, theta

! Note: usav and corr initialized to zero upon allocation; following probably
! not necessary, but occurs only once (per nonlinear solve)
  if (counter == 1) then
    usav(:,:,:,pt) = 0.d0
    corr(:,:,:,old(pt),pt) = 0.d0
  end if

  ! RESIDUAL CALCULATION

  !TODO - Remove hardwired numbers for whichresid (see glide_types.F90)

  select case (whichresid)
  ! options for residual calculation method, as specified in configuration file
  ! (see additional notes in "higher-order options" section of documentation)
  ! case(0): use max of abs( vel_old - vel ) / vel )
  ! case(1): use max of abs( vel_old - vel ) / vel ) but ignore basal vels
  ! case(2): use mean of abs( vel_old - vel ) / vel )
  ! case(3): use max of abs( vel_old - vel ) / vel ) (in addition to L2 norm calculated externally)

   case(HO_RESID_MAXU)

    ! resid = maxval( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel /= 0.d0)
    resid = 0.d0

    do ns = 1 + staggered_lhalo, size(vel,3) - staggered_uhalo
      do ew = 1 + staggered_lhalo, size(vel,2) - staggered_uhalo
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

   case(HO_RESID_MAXU_NO_UBAS)
    ! nr = size( vel, dim=1 ) ! number of grid points in vertical ...
    ! resid = maxval( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ), MASK = vel /= 0.d0)
    resid = 0.d0

    do ns = 1 + staggered_lhalo, size(vel,3) - staggered_uhalo
      do ew = 1 + staggered_lhalo, size(vel,2) - staggered_uhalo
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

   case(HO_RESID_MEANU)
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

   case(HO_RESID_L2NORM)

!! SFP - the L2norm option is handled entirely external to this subroutine. That is, if the L2norm option
!! for the residul is specified (it is currently the default), the residual is calculated as the L2norm of 
!! the system residul, r = Ax - b (rather than defining the residual according to the velocity update, as
!! is done in all the parts of this subroutine). If the L2norm option is active, the value of "residual" 
!! passed out of this subroutine is NOT used for determining when to halt iterations on the velocity solution.
!! The original code that was here for this option has been removed.

  end select

  if (GLC_DEBUG) then
     ! Additional debugging line, useful when trying to determine if convergence is being consistently
     ! help up by the residual at one or a few particular locations in the domain.
     ! print '("* ",i3,g20.6,3i6,g20.6)', counter, resid, locat, vel(locat(1),locat(2),locat(3))*vel0
  end if

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

    do ns = 1 + staggered_lhalo, size(vel,3) - staggered_uhalo
      do ew = 1 + staggered_lhalo, size(vel,2) - staggered_uhalo
        do nr = 1, size(vel, 1)
          temp_vel = vel(nr,ew,ns)

          if (acos((corr(nr,ew,ns,new(pt),pt) * corr(nr,ew,ns,old(pt),pt)) / &
               (abs(corr(nr,ew,ns,new(pt),pt)) * abs(corr(nr,ew,ns,old(pt),pt)) + small)) > &
              ssthres .and. corr(nr,ew,ns,new(pt),pt) - corr(nr,ew,ns,old(pt),pt) /= 0.d0) then

               ! theta and alpha are intermediate terms that might be useful to print out
!!             theta = acos((corr(nr,ew,ns,new(pt),pt) * corr(nr,ew,ns,old(pt),pt)) / &
!!                    (abs(corr(nr,ew,ns,new(pt),pt)) * abs(corr(nr,ew,ns,old(pt),pt)) + small))

!!             alpha = abs(corr(nr,ew,ns,old(pt),pt)) / &
!!                     abs(corr(nr,ew,ns,new(pt),pt) - corr(nr,ew,ns,old(pt),pt))

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

! TODO take out the older one?

  !*SFP* Old version
  ! if (new(pt) == 1) then; old(pt) = 1; new(pt) = 2; else; old(pt) = 1; new(pt) = 2; end if

  !*SFP* correction from Carl Gladdish
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
    
   if (theta  < (1.d0/8.d0)*pi) then
        mindcrshstr2 = usav(:,:,:,pt) + cvg_accel * corr(:,:,:,new(pt),pt)
!        print *, theta/pi, 'increased correction'
   else if(theta < (19.d0/20.d0)*pi) then
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

   case(HO_RESID_MAXU)
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
!    sig_rel_diff = sqrt( sum((rel_diff - mean_rel_diff) ** 2.d0 )/ sum(vel_ne_0) )
!    print *, 'mean', mean_rel_diff, 'sig', sig_rel_diff

    !write(*,*) 'locat', locat
    !call write_xls('resid1.txt',abs((usav(1,:,:,pt) - mindcrshstr2(1,:,:)) / (mindcrshstr2(1,:,:) + 1e-20)))

   case(HO_RESID_MAXU_NO_UBAS)
    !**cvg*** should replace vel by mindcrshstr2 in the following lines, I belive
    nr = size( vel, dim=1 ) ! number of grid points in vertical ...
    resid = maxval( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
                        MASK = vel /= 0.d0)
    locat = maxloc( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
            MASK = vel /= 0.d0)

   case(HO_RESID_MEANU)
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
                       mintauf,     flwa,           &
                       beta,        btraction,      &
                       bwat,        assembly )

  ! Main subroutine for determining coefficients that go into the LHS matrix A 
  ! in the expression Au = b. Calls numerous other subroutines, including boundary
  ! condition subroutines, which determine "b".

  use parallel
!!  use glimmer_horiz_bcs, only: ghost_shift

  implicit none

  integer, intent(in) :: ewn, nsn, upn, assembly
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

  real(dp), dimension(:,:), intent(in) :: mintauf
  real(dp), dimension(:,:), intent(in) :: beta
  real(dp), dimension(:,:), intent(in) :: bwat 
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

  logical :: comp_bound

  ct_nonzero = 1        ! index to count the number of non-zero entries in the sparse matrix

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
                        mintauf, beta,          &
                        betasquared, mask,      &
                        bwat )
  ! intent(out) betasquared

  ! Note loc2_array is defined only for non-halo ice grid points.
  ! JEFFLOC returns an array with starting indices into solution vector for each ice grid point.
 
  allocate(loc2_array(ewn,nsn,2))

!WHL - Using a different procedure depending on whether or not we are using trilinos.
!      This is needed to avoid an error when using the SLAP solver in a
!        single-processor parallel run.
!TODO: Find a more elegant solution?
               
  loc2_array = getlocationarray(ewn, nsn, upn, mask, uindx)

  if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
     loc2_array = getlocationarray(ewn, nsn, upn, mask, uindx, &
                                   return_global_IDs = .false.)
  else
     loc2_array = getlocationarray(ewn, nsn, upn, mask, uindx)
  endif

!WHL - debug
!  print*, ' '
!  print*, 'loc2_array(1)'
!  do ns = nsn, 1, -1
!     write(6,'(34i6)') loc2_array(:,ns,1)
!  enddo

!  print*, ' '
!  print*, 'loc2_array(2)'
!  do ns = nsn, 1, -1
!     write(6,'(34i6)') loc2_array(:,ns,2)
!  enddo

  !  !!!!!!!!! useful for debugging !!!!!!!!!!!!!!
  !    print *, 'loc2_array = '
  !    print *, loc2_array
  !    pause
  
  ! Note: With nhalo = 2, efvs has been computed in a layer of halo cells, 
  !       so we have its value in all neighbors of locally owned velocity points.

  do ns = 1+staggered_lhalo, size(mask,2)-staggered_uhalo
    do ew = 1+staggered_lhalo, size(mask,1)-staggered_uhalo

      !Theoretically, this should just be .false. to remove it from the if statements and let the ghost cells
      !take over. However, with only one process, this give an exception error when calc_F calls savetrilinosmatrix(0).
      !Therefore, it will currently revert back to the old BC's when using only one task for now. I am working to
      !debug and fix this case, but for now, it does no harm for the original BC's.

!      comp_bound = ( nslb < 1          .and. ns <              staggered_lhalo+1+ghost_shift ) .or. &
!                   ( ewlb < 1          .and. ew <              staggered_lhalo+1+ghost_shift ) .or. &
!                   ( nsub > global_nsn .and. ns > size(mask,2)-staggered_uhalo  -ghost_shift ) .or. &
!                   ( ewub > global_ewn .and. ew > size(mask,1)-staggered_uhalo  -ghost_shift )

      comp_bound = .false.

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
      if ( GLIDE_HAS_ICE(mask(ew,ns)) .and. .not. &
           comp_bound .and. .not. &
           GLIDE_IS_MARGIN(mask(ew,ns)) .and. .not. &
           GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) .and. .not. &
           GLIDE_IS_CALVING(mask(ew,ns) ) .and. .not. &
           GLIDE_IS_THIN(mask(ew,ns) ) ) then
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

          !HALO - Note that ew and ns below are locally owned velocity points.
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
        enddo  ! upn
        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !TODO - Not sure COMP_DOMAIN_BND condition is needed
      elseif ( GLIDE_IS_CALVING( mask(ew,ns) ) .and. .not. &
               comp_bound .and. .not. &
               GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) .and. .not. &
               GLIDE_IS_THIN(mask(ew,ns) ) ) then
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
                       abar=flwabar)
        enddo
        lateralboundry = .false.
        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !TODO - Here we deal with cells on the computational domain boundary.
        !       Currently the velocity is always set to a specified value on this boundary.
        !       With open (non-Dirichlet) BCs, we might want to solve for these velocities,
        !        using the code above to compute the matrix elements.
      elseif ( GLIDE_HAS_ICE(mask(ew,ns)) .and. ( GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) .or. &
               comp_bound ) .or. GLIDE_IS_LAND_MARGIN(mask(ew,ns)) .or. &
               GLIDE_IS_THIN(mask(ew,ns)) ) then
        !    print*, ' '
        !    print*, 'At a NON-SHELF boundary ... ew, ns = ', ew, ns
        !    print*, 'LAND_MARGIN =', GLIDE_IS_LAND_MARGIN(mask(ew,ns))
        !    print*, 'MASK(ew,ns) =', mask(ew,ns)
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
         !call valueset( 0.d0 )                ! vel at margin set to 0 
        enddo
      endif
    enddo      ! ew 
  enddo        ! ns

  deallocate(loc2_array)

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
                   abar)

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

  logical :: fons, foew     ! true when geom. requires using 1st-order one sided diffs. at floating ice boundary
                            ! (default is 2nd-order, which requires larger stencil)

  loc2plusup = loc2(1,:) + up

  if( lateralboundry )then

  ! *********************************************************************************************
  ! lateral boundary conditions 

  ! if at sfc or bed, source due to seawater pressure is 0 and bc normal vector
  ! should contain sfc/bed slope components, e.g. (-ds/dx, -ds/dy, 1) or (db/dx, db/dy, -1)
     source = 0.d0

     call getlatboundinfo( ew,  ns,  up,                                 &
                           ewn, nsn, upn,                                &
                           stagthck(ew-2:ew+2, ns-2:ns+2),               &
                           loc2_array(:,:,1), fwdorbwd, normal,          & 
                           loc_latbc, foew, fons)

     if( up == 1 .or. up == upn )then

        if( up == 1 )then                ! specify necessary variables and flags for free sfc
           bcflag = (/1,0/)
           loc2plusup = loc2(1,:) + up - 1   ! reverse the sparse matrix / rhs vector row index by 1 ...
           slopex = -dusrfdew(ew,ns); slopey = -dusrfdns(ew,ns); nz = 1.d0
        else                             ! specify necessary variables and flags for basal bc
   
           if( whichbabc == HO_BABC_NO_SLIP )then
                bcflag = (/0,0/)             ! flag for u=v=0 at bed; doesn't work well so commented out here...
                                             ! better to specify very large value for betasquared below
           elseif( whichbabc == HO_BABC_CONSTANT     .or. whichbabc == HO_BABC_SIMPLE         .or.  &
                   whichbabc == HO_BABC_YIELD_PICARD .or. whichbabc == HO_BABC_BETA_BWAT .or.  &
                   whichbabc == HO_BABC_LARGE_BETA   .or. whichbabc == HO_BABC_EXTERNAL_BETA) then
                bcflag = (/1,1/)              ! flag for specififed stress at bed: Tau_zx = betasquared * u_bed,
                                              ! where betasquared is MacAyeal-type traction parameter
           end if   

           loc2plusup = loc2(1,:) + up + 1   ! advance the sparse matrix / rhs row vector index by 1 ...
           slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns); nz = -1.d0

        end if

!TODO: conduct realistic test cases with and w/o this hack
!        !! Hack to avoid bad sfc and basal bc normal vectors !!
!        slopex = 0.d0; slopey = 0.d0

        ! get coeffs. associated with horiz. normal stresses lateral boundary
        g = normhorizmainbc_lat(dew,           dns,             &
                                slopex,        slopey,          &
                                dsigmadew(up), dsigmadns(up),   &
                                pt,            2,               &
                                dup(up),       local_efvs,      &
                                oneorfour,     fourorone,       &
                                onesideddiff,                   &
                                normal,        fwdorbwd,        &
                                foew,          fons    )

        ! add on coeffs. associated with vertical shear stresses
        g(:,3,3) = g(:,3,3) &
                 + vertimainbc( stagthck(ew,ns), bcflag, dup(up),              &
                                local_efvs,      betasquared,   g_vert,  nz ) 

        !! scale basal bc coeffs when using JFNK solver 
        scalebabc = scalebasalbc( g, bcflag, lateralboundry, betasquared, local_efvs )
        g = g / scalebabc ! put the coeff. for the b.c. equation in the same place as the prev. equation
        ! (w.r.t. cols), on a new row ...
        call fillsprsebndy( g, loc2plusup(1), loc_latbc, up, normal, pt )


        ! get coeffs. for horiz shear stress terms, multiply by other vel and put into RHS vector

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
                                                   normal,        fwdorbwd,      &
                                                   foew,          fons )         & 
                                                 * local_othervel ) /scalebabc

    end if     ! up = 1 or up = upn (IF at lateral boundary and IF at surface or bed)

    ! If in main body and at ice/ocean boundary, calculate depth-averaged stress
    ! due to sea water, bc normal vector components should be boundary normal 
    loc2plusup = loc2(1,:) + up

    ! for this bc, the normal vector components are not the sfc/bed slopes but are taken
    ! from a normal to the shelf front in map view (x,y plane); slopex,slopey are simply renamed here 
    slopex = normal(1)
    slopey = normal(2)

    ! There are two options here for the source term associated with the boundary condition for 
    ! floating ice:
    !
    !  (1) use the 1d solution that involves the rate factor (not accurate for 
    !      3d domains, but can be more robust and stable) 
    !  (2) use the more general solution that involves the eff. visc. and normal
    !      vector orientation at lateral boundary
    !
    ! Only one of these options should be active at a time (i.e. comment the other lines out)
    ! The default setting is (2), which is the more general case that should also work for 1d problems. 

    ! In some cases, the two options can be used together to improve performance, e.g. for the Ross
    ! ice shelf experiment, a number of early iterations could use the more simple bc (option 1) and then
    ! when the solution has converged a bit, we switch to the more realistic implementation (option 2).
    ! This has the advantage of "conditioning" the eff. visc. in the source term a bit before turning 
    ! the source term dependence on the eff. visc. "on". 

    ! NOTE that the newer sfc, basal, and lateral bc subroutines keep the eff. visc. terms with the LHS
    ! matrix coeffs. In this case, they do not have any affect on the source term for floating ice bcs
    ! and the considerations in the above paragraph do not apply (w.r.t. adversely affecting the source term).  

!    ! --------------------------------------------------------------------------------------
!    ! (1) source term (strain rate at shelf/ocean boundary) from Weertman's analytical solution 
!    ! This is primarily of use for debugging purposes, e.g. when a 1d test case is run. Also useful
!    ! if one wants to turn "off" the eff. visc. dependence in the matrix coeffs. that go with this
!    ! boundary condition, since this form of it has no eff. visc. terms.
!    ! --------------------------------------------------------------------------------------
!    ! See eq. 2, Pattyn+, 2006, JGR v.111; eq. 8, Vieli&Payne, 2005, JGR v.110). Note that this 
!    ! contains the 1d assumption that ice is not spreading lateraly !(assumes dv/dy = 0 for u along flow)
!
!    source = abar * vis0 * ( 1.d0/4.d0 * rhoi * grav * stagthck(ew,ns)*thk0 * ( 1.d0 - rhoi/rhoo))**3.d0
!
!    ! multiply by 4 so that case where v=0, du/dy = 0, LHS gives: du/dx = du/dx|_shelf 
!    ! (i.e. LHS = 4*du/dx, requires 4*du/dx_shelf)
!    source = source * 4.d0
!
!    ! split source based on the boundary normal orientation and non-dimensinoalize
!    ! Note that it is not really appropriate to apply option (1) to 2d flow, since terms other than du/dx in 
!    ! eff. strain rate are ignored. For 2d flow, should use option (2) below. 
!     source = source * normal(pt)
!     source = source * tim0 ! make source term non-dim
!    ! --------------------------------------------------------------------------------------

    ! --------------------------------------------------------------------------------------
    ! (2) source term (strain rate at shelf/ocean boundary) from MacAyeal depth-ave solution. 
    ! --------------------------------------------------------------------------------------

    source = (rhoi*grav*stagthck(ew,ns)*thk0) / tau0 / 2.d0 * ( 1.d0 - rhoi / rhoo )

    source = source * normal(pt) ! partition according to normal vector at lateral boundary
                                 ! NOTE that source term is already non-dim here 
    ! --------------------------------------------------------------------------------------

    ! get matrix coefficients that go with horiz normal stresses at a floating ice boundary
    g = normhorizmainbc_lat(dew,           dns,        &
                            slopex,        slopey,     &
                            dsigmadew(up), dsigmadns(up),  &
                            pt,            1,          &
                            dup(up),       local_efvs, &
                            oneorfour,     fourorone,  &
                            onesideddiff,              &
                            normal,        fwdorbwd,   &
                            foew,          fons    )

    ! NOTE that for lateral floating ice boundary, we assume u_sfc ~ u_bed and stress free bc
    ! at both upper and lower sfc boundaries, so that there are no coeffs. for vert. shear stresses 

    ! put the coeff. for the b.c. equation in the same place as the prev. equation
    ! (w.r.t. cols), on a new row ...

!TODO: is above comment correct or is this now just a normal scatter of coeffs. into the matrix?
    call fillsprsebndy( g, loc2plusup(1), loc_latbc, up, normal, pt )


    ! get matrix coefficients that go with the horiz shear stresses at a floating ice
    ! boundary, multiply by their respective "other" velocity and put into RHS vector
    
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
                                               normal,        fwdorbwd,       &
                                               foew,          fons )          & 
                                              * local_othervel ) + source

  else   ! NOT at a lateral boundary 

! *********************************************************************************************
! normal discretization for points inside of lateral boundary and inside main body of ice sheet

    ! This if construct skips the normal discretization for the RHS and LHS for the sfc and basal indices
    ! because these vertical levels are handled by different subroutines.
    if( up /= upn .and. up /= 1 )then

       g = normhorizmain(pt,up,local_efvs)                      ! normal stress grad coeffs

       g(:,2,2) = g(:,2,2) + vertimain(hsum(local_efvs),up)     ! add vert stress grad coeffs

    ! NOTE that version of 'fillspremain' for one-sided bcs needs additional argument to specify a
    ! column shift of coeffs. of rows in LHS matrix. That is the "0" past last here (no shift for internal bcs)
       call fillsprsemain(g,loc2plusup(1),loc2(:,1),up,pt,0)

    ! NOTE that in the following expression, the "-" sign on the crosshoriz terms,
    ! which results from moving them from the LHS over to the RHS, is explicit and
    ! hast NOT been moved inside of "croshorizmin" (as is the case for the analogous
    ! boundary condition routines).
       rhsd(loc2plusup(2)) = thisdusrfdx(ew,ns) - &             ! shear stress grad coeffs into RHS vector
                             sum(croshorizmain(pt,up,local_efvs) * local_othervel)
    end if
 
   ! The follow two if constructs set the ghost cell storage to have ones on the martrix diag and zeros 
   ! on the rhs, enforcing a zero vel bc for the ghost cells. Eventually, the capacity allowing for ghost
   ! cells can probably be removed but keeping here for now for backward compatibility.
    if( up == upn  )then
       loc2plusup = loc2(1,:) + upn + 1    ! basal ghost cells
       call valueset(0.d0, loc2plusup)
    endif
    if( up == 1  )then
       loc2plusup = loc2(1,:)              ! sfc ghost cells
       call valueset(0.d0, loc2plusup)
    endif
 
  end if

! *********************************************************************************************
! higher-order sfc and bed boundary conditions in main body of ice sheet (NOT at lat. boundry)

  if(  ( up == upn  .or. up == 1 ) .and. .not. lateralboundry) then

        if( up == 1 )then                ! specify necessary variables and flags for free sfc
           bcflag = (/1,0/)
           loc2plusup = loc2(1,:) + up - 1   ! reverse the sparse matrix / rhs vector row index by 1 ...
           slopex = -dusrfdew(ew,ns); slopey = -dusrfdns(ew,ns); nz = 1.d0
        else                             ! specify necessary variables and flags for basal bc

           if( whichbabc == HO_BABC_NO_SLIP )then
                bcflag = (/0,0/)             ! flag for u=v=0 at bed; doesn't work well so commented out here...
                                             ! better to specify very large value for betasquared below

           elseif( whichbabc == HO_BABC_CONSTANT     .or. whichbabc == HO_BABC_SIMPLE         .or.  &
                   whichbabc == HO_BABC_YIELD_PICARD .or. whichbabc == HO_BABC_BETA_BWAT .or.  &
                   whichbabc == HO_BABC_LARGE_BETA   .or. whichbabc == HO_BABC_EXTERNAL_BETA) then
                bcflag = (/1,1/)              ! flag for specififed stress at bed: Tau_zx = betasquared * u_bed,
                                              ! where betasquared is MacAyeal-type traction parameter
           end if
           
           loc2plusup = loc2(1,:) + up + 1   ! advance the sparse matrix / rhs row vector index by 1 ...
           slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns); nz = -1.d0
        
        end if

        ! get matrix coefficients that go with normal stresses at sfc or basal boundary 
        g = normhorizmainbcos(dew,         dns,            &
                         slopex,        slopey,         &
                         dsigmadew(up), dsigmadns(up),  &
                         pt,            bcflag,         &
                         dup(up),       local_efvs,     &
                         oneorfour,     fourorone)

        g_norm = g              ! save these coeffs, as needed for basal traction calculation


        ! get matrix coefficients that go with vertical stresses at sfc or basal boundary 
        g(:,2,2) = g(:,2,2)   &
                 + vertimainbcos( stagthck(ew,ns),bcflag,dup(up),local_efvs, &
                                    betasquared, g_vert, nz )

        !! scale basal bc coeffs when using JFNK solver 
        scalebabc = scalebasalbc( g, bcflag, lateralboundry, betasquared, local_efvs )
        g = g / scalebabc
   
        loc2plusup = loc2(1,:) + up  ! Need to reset this index since we want the bc on the actual row
                                     ! coinciding with the boundary at up=1
   
        ! Replace ghost cells w/ one-sided diffs at sfc/basal indices. This section shifts the LHS matrix coeffs for the sfc
        ! and basal bcs back on to the main diagonal, as opposed to staggered off the diag, which was necessary for the ghost 
        ! cell implementation.
        if( up == 1 .or. up == upn )then
           if( up == 1 )then
             call fillsprsemain(g,loc2plusup(1),loc2(:,1),up,pt,1)
           else if( up == upn )then
             call fillsprsemain(g,loc2plusup(1),loc2(:,1),up,pt,-1)
           end if
        end if

        ! calc shear stress coeffs., multiply by other vel and move to RHS vector
        rhsd(loc2plusup(2)) = sum( croshorizmainbcos(dew,          dns,            &
                                   slopex,        slopey,         &
                                   dsigmadew(up), dsigmadns(up),  &
                                   pt,            bcflag,         &
                                   dup(up),       local_othervel, &
                                   local_efvs,                    &
                                   oneortwo, twoorone, g_cros )   &
                                   * local_othervel ) / scalebabc 

        ! The following calculates the basal traction AFTER an updated solution is obtain by passing the new
        ! values of uvel, vvel back to the matrix assembly routines, and thus obtaining updated values of the 
        ! relevant coefficients. The if construct allows the assembly routines to be called for only the vert
        ! layers that are needed to cacluate the basal traction (as opposed to all vert levels 1:upn).
        if( assembly == 1 )then

           g_vel_lhs = local_thisvel
           g_vel_rhs = local_othervel

!HALO - Since ew and ns are locally owned velocity points, we will have btraction at all such points.
           btraction(pt,ew,ns) = sum( (g_norm+g_vert)*g_vel_lhs*thk0/len0 ) &
                                  - sum( g_cros*g_vel_rhs*thk0/len0 )
        end if

  end if   ! (up = 1 or up = upn) and lateralboundry = F

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

function vertimainbc(thck, bcflag, dup, efvs, betasquared, g_vert, nz )

! altered form of 'vertimain' that calculates coefficients for higher-order
! b.c. that go with the 'normhorizmain' term: -(X/H)^2 * dsigma/dzhat * du/dsigma 

    implicit none

    real(dp), intent(in) :: dup, thck, betasquared 
    real(dp), intent(in) :: nz                      ! sfc normal vect comp in z-dir
    real(dp), intent(in), dimension(2,2,2) :: efvs
    real(dp), intent(out), dimension(3,3,3) :: g_vert
    integer, intent(in), dimension(2) :: bcflag

    real(dp) :: c
    real(dp), dimension(3) :: vertimainbc

    c = 0.d0
    g_vert = 0.d0

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    if( bcflag(1) == 1 )then

           c = nz / thck / (2.d0*dup) * (len0**2 / thk0**2)   ! value of coefficient

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

function vertimainbcos(thck, bcflag, dup, efvs, betasquared, g_vert, nz )

! altered form of 'vertimain' that calculates coefficients for higher-order
! b.c. that go with the 'normhorizmain' term: -(X/H)^2 * dsigma/dzhat * du/dsigma

    implicit none

    real (dp), intent(in) :: dup, thck, betasquared
    real (dp), intent(in) :: nz                      ! sfc normal vect comp in z-dir
    real (dp), intent(in), dimension(2,2,2) :: efvs
    real (dp), intent(out), dimension(3,3,3) :: g_vert
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

           c = nz / thck / (2.d0*dup) * (len0**2 / thk0**2) * efvsbar_sfc ! value of coefficient

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

           c = nz / thck / (2.d0*dup) * (len0**2 / thk0**2) * efvsbar_bed ! value of coefficient

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
!               + oneorfour(which) * dusrfdns * dsigmadns )/(2.d0*dup)
           c = ( fourorone(which) * dusrfdew * dsigmadew   &
               + oneorfour(which) * dusrfdns * dsigmadns )/(2.d0*dup) * ( sum( efvs(1,:,:) ) / 4.d0 )

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
           c = fourorone(which) * dusrfdew / (2.d0*dew) * ( sum( efvs(2,:,:) ) / 4.d0 )
           g(3,3,2) = c
           g(3,1,2) = -c

!           c = oneorfour(which) * dusrfdns / (2*dns)
           c = oneorfour(which) * dusrfdns / (2.d0*dns) * ( sum( efvs(2,:,:) ) / 4.d0 )
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
                         g_cros, velbc )

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
    real (kind = dp), intent(in), optional :: velbc
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
!                 - oneortwo(which) * dusrfdns * dsigmadew )/(2.d0*dup)
           c = ( - twoorone(which) * dusrfdew * dsigmadns   &
                 - oneortwo(which) * dusrfdns * dsigmadew )/(2.d0*dup) * ( sum( efvs(1,:,:) ) / 4.d0 )

           g(1,2,2) = 3.d0*c
           g(2,2,2) = -4.d0*c
           g(3,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
!           c = - oneortwo(which) * dusrfdns / (2*dew)
           c = - oneortwo(which) * dusrfdns / (2.d0*dew) * ( sum( efvs(1,:,:) ) / 4.d0 )
           g(1,3,2) = c
           g(1,1,2) = -c

!           c = - twoorone(which) * dusrfdew / (2*dns)
           c = - twoorone(which) * dusrfdew / (2.d0*dns) * ( sum( efvs(1,:,:) ) / 4.d0 )
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
                 - oneortwo(which) * dusrfdns * dsigmadew )/(2.d0*dup) * ( sum( efvs(2,:,:) ) / 4.d0 )

           g(1,2,2) = -1.d0*c
           g(2,2,2) = 4.d0*c
           g(3,2,2) = -3.d0*c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
!           c = - oneortwo(which) * dusrfdns / (2*dew)
           c = - oneortwo(which) * dusrfdns / (2.d0*dew) * ( sum( efvs(2,:,:) ) / 4.d0 )
           g(3,3,2) = c
           g(3,1,2) = -c


!           c = - twoorone(which) * dusrfdew / (2*dns)
           c = - twoorone(which) * dusrfdew / (2.d0*dns) * ( sum( efvs(2,:,:) ) / 4.d0 )
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
            g = 1.d0
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
                             normal,    fwdorbwd,  &
                             foew,      fons )

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

    logical, intent(in) :: fons, foew   ! true when geom. requires 1st-order one sided diffs for shelf bcs

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
    !efvsbar = 1.0d0; 

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

           if( foew )then
               c =  -1.d0 * fwdorbwd(1) * fourorone(which) * dusrfdew / dew * efvsbar
           else
               c = fourorone(which) * fwdorbwd(1) * onesideddiff(1) * dusrfdew / (2.d0*dew) * efvsbar
           endif
           g(2,2-int(fwdorbwd(1)),2) = c

           if( foew )then
               c = fwdorbwd(1)*fourorone(which) * dusrfdew / dew * efvsbar
           else
               c = fourorone(which) * fwdorbwd(1) * onesideddiff(2) * dusrfdew / (2.d0*dew) * efvsbar
           endif
           g(2,2,2) = c

           if( foew )then
               c = 0.d0
           else
               c = fourorone(which) * fwdorbwd(1) * onesideddiff(3) * dusrfdew / (2.d0*dew) * efvsbar
           endif
           g(2,2+int(fwdorbwd(1)),2) = c

    end if

    if( normal(2) == 0.d0 ) then   ! centered in y ... 
                                       ! (NOTE that y coeff. are stored in g(1,:,:) )

           c = oneorfour(which) * dusrfdns / (2*dns) * efvsbar
           g(1,2,3) = c
           g(1,2,1) = -c

    elseif( normal(2) /= 0.d0) then ! forward/backward in y ...

           if( fons )then
               c =  -1.d0 * fwdorbwd(2) * oneorfour(which) * dusrfdns / dns * efvsbar
           else
               c = oneorfour(which) * fwdorbwd(2) * onesideddiff(1) * dusrfdns / (2.d0*dns) * efvsbar
           endif
           g(1,2,2-int(fwdorbwd(2))) = c

           if( fons )then
               c = fwdorbwd(2)*oneorfour(which) * dusrfdns / dns * efvsbar
           else
               c = oneorfour(which) * fwdorbwd(2) * onesideddiff(2) * dusrfdns / (2.d0*dns) * efvsbar
           endif
           g(1,2,2) = c

           if( fons )then
               c = 0.d0
           else
               c = oneorfour(which) * fwdorbwd(2) * onesideddiff(3) * dusrfdns / (2.d0*dns) * efvsbar
           endif
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
                              normal,    fwdorbwd,  &
                              foew,      fons )

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

    logical, intent(in) :: fons, foew   ! true when geom. requires 1st-order one sided diffs for shelf bcs

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
              - oneortwo(which) * dusrfdns * dsigmadew )/(2.d0*dup) * efvsbar
    gvert(3) = -c * whichbc(what)
    gvert(1) = c * whichbc(what)

    if( normal(1) == 0.d0 )then        ! centered in x ...

           c = -oneortwo(which) * dusrfdns / (2.d0*dew) * efvsbar
           g(2,3,2) = c
           g(2,1,2) = -c

    elseif( normal(1) /= 0.d0 )then    ! forward/backward in x ...
                                           ! (NOTE that x coeff. are stored in g(2,:,:) )

           if( foew )then
               c =  oneortwo(which) * fwdorbwd(1) * dusrfdns / dew * efvsbar
           else
               c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(1) * dusrfdns / (2.d0*dew) * efvsbar
           endif
           g(2,2-int(fwdorbwd(1)),2) = c

           if( foew )then
               c = -oneortwo(which) * fwdorbwd(1) * dusrfdns / dew * efvsbar
           else
               c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(2) * dusrfdns / (2.d0*dew) * efvsbar
           endif
           g(2,2,2) = c

           if( foew )then
               c = 0.d0
           else
               c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(3) * dusrfdns / (2.d0*dew) * efvsbar
           endif
           g(2,2+int(fwdorbwd(1)),2) = c

    end if

    if( normal(2) == 0.d0 )then    ! centered in y ...
                                       ! (NOTE that y coeff. are stored in g(1,:,:) )

           c = -twoorone(which) * dusrfdew / (2.d0*dns) * efvsbar
           g(1,2,3) = c
           g(1,2,1) = -c

    elseif( normal(2) /= 0.d0 )then ! forward/backward in y ...

           if( fons )then
               c =  twoorone(which) * fwdorbwd(2) * dusrfdew / dns * efvsbar
           else
               c = -twoorone(which) * fwdorbwd(2) * onesideddiff(1) * dusrfdew / (2.d0*dns) * efvsbar
           endif
           g(1,2,2-int(fwdorbwd(2))) = c

           if( fons )then
               c = -twoorone(which) * fwdorbwd(2) * dusrfdew / dns * efvsbar
           else
               c = -twoorone(which) * fwdorbwd(2) * onesideddiff(2) * dusrfdew / (2.d0*dns) * efvsbar
           endif
           g(1,2,2) = c

           if( fons )then
               c = 0.d0
           else
               c = -twoorone(which) * fwdorbwd(2) * onesideddiff(3) * dusrfdew / (2.d0*dns) * efvsbar
           endif
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

subroutine getlatboundinfo( ew, ns, up, ewn, nsn, upn,    &
                           thckin, loc_array,             &
                           fwdorbwd, normal, loc_latbc,   &
                           foew, fons)

  ! Calculate map plane normal vector at 45 deg. increments
  ! for regions of floating ice
  implicit none

  integer, intent(in) :: ew, ns, up
  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(ewn,nsn), intent(in) :: loc_array

  real(dp), dimension(5,5), intent(in) :: thckin

  real(dp), dimension(2), intent(out) :: fwdorbwd, normal
  integer, dimension(6), intent(out) :: loc_latbc

  logical, intent(out) :: fons, foew

  real(dp), dimension(3,3) :: mask, maskcorners

  integer, dimension(5,5) :: thckinmask

  real(dp), dimension(3,3) :: thckmask, thck
  real(dp), dimension(3) :: testvect
  real(dp) :: phi, deg2rad

  thck(:,:) = thckin(2:4,2:4)
  thckinmask = 0

!  deg2rad = 3.141592654d0 / 180.d0
  deg2rad = pi / 180.d0
  loc_latbc = 0; phi = 0.d0
  mask(:,1) = (/ 0.d0, 180.d0, 0.d0 /)
  mask(:,2) = (/ 270.d0, 0.d0, 90.d0 /)
  mask(:,3) = (/ 0.d0, 360.d0, 0.d0 /)
  maskcorners(:,1) = (/ 225.d0, 0.d0, 135.d0 /)
  maskcorners(:,2) = (/ 0.d0, 0.d0, 0.d0 /)
  maskcorners(:,3) = (/ 315.d0, 0.d0, 45.d0 /)

  !! first section below contains logic to ID where 1st-order one-sided diffs are needed
  where( thckin /= 0.d0 )
        thckinmask = 1
  endwhere
  !! check if 1st-order one sided diffs. are needed in n/s direction
  if( (thckinmask(3,3)+thckinmask(3,4)+thckinmask(3,5)) < 3 .and. (thckinmask(3,1)+thckinmask(3,2)) < 2 )then
        !print *, '1st-order one-sided diffs. in N/S direction at ew,ns = ', ew, ns
        fons = .true. 
  elseif( (thckinmask(3,1)+thckinmask(3,2)+thckinmask(3,3)) < 3 .and. (thckinmask(3,4)+thckinmask(3,5)) < 2 )then
        !print *, '1st-order one-sided diffs. in N/S direction at ew,ns = ', ew, ns
        fons = .true. 
  else 
        fons = .false.
  endif
  !! check if 1st-order one sided diffs. are needed in n/s direction
  if( (thckinmask(3,3)+thckinmask(4,3)+thckinmask(5,3)) < 3 .and. (thckinmask(1,3)+thckinmask(2,3)) < 2 )then
        !print *, '1st-order one-sided diffs. in E/W direction at ew,ns = ', ew, ns
        foew = .true. 
  elseif( (thckinmask(1,3)+thckinmask(2,3)+thckinmask(3,3)) < 3 .and. (thckinmask(4,3)+thckinmask(5,3)) < 2 )then
        !print *, '1st-order one-sided diffs. in E/W direction at ew,ns = ', ew, ns
        foew = .true. 
  else 
        foew = .false.
  endif

  ! specify new value of 'loc' vector such that fwd/bwd diffs. are set up correctly in sparse matrix
  ! when function 'fillsprsebndy' is called. Also, specify appropriate values for the vectors 'normal'
  ! and 'fwdorbwd', which specify the orientation of the boundary normal and the direction of forward or
  ! backward differencing to be done in the lateral boundary condition functions 'normhorizmainbc_lat'
  ! and 'croshorizmainbc_lat'

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

!TODO: This function does not use loc_array.  Remove from argument list?

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

!  deg2rad = 3.141592654d0 / 180.d0
  deg2rad = pi / 180.d0
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

  !TODO - Remove hardwiring of case numbers?
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
                            mintauf, beta,           &
                            betasquared, mask,       &
                            bwat,                    &
                            betafile)

  ! subroutine to calculate map of betasquared sliding parameter, based on 
  ! user input ("whichbabc" flag, from config file as "which_ho_babc").

  implicit none

  integer, intent(in) :: whichbabc
  integer, intent(in) :: ewn, nsn

  real(dp), intent(in) :: dew, dns
  real(dp), intent(in), dimension(:,:) :: lsrf, topg, thck
  real(dp), intent(in), dimension(:,:) :: thisvel, othervel, mintauf, beta, bwat

  integer, intent(in), dimension(:,:) :: mask 

  real(dp), intent(out), dimension(ewn-1,nsn-1) :: betasquared

  character (len=30), intent(in), optional :: betafile
  real(dp) :: smallnum = 1.0d-2
  real(dp), dimension(ewn) :: grounded
  real(dp) :: alpha, dx, thck_gl, betalow, betahigh, roughness 
  integer :: ew, ns

  ! SFP added for making beta a function of basal water flux 
  real(dp), dimension(:,:), allocatable :: unstagbetasquared
  real(dp) :: C, m

  ! Note that the dimensional scale (tau0 / vel0 / scyr ) is used here for making the basal traction coeff.
  ! betasquared dimensional, within the subroutine (mainly for debugging purposes), and then non-dimensional 
  ! again before being sent back out for use in the code. This scale is the same as "scale_beta" defined in 
  ! libglimmer/glimmer_scales.F90. See additional notes where that scale is defined. In general, below, it is
  ! assumed that values for betasquared being accessed from inside the code are already dimensionless and 
  ! any hardwired values have units of Pa yr/m. 

!TODO - Remove scaling here?

  select case(whichbabc)

    case(HO_BABC_CONSTANT)  ! constant value; useful for debugging and test cases

      betasquared(:,:) = 10.d0       ! Pa yr/m

    case(HO_BABC_SIMPLE)    ! simple pattern; also useful for debugging and test cases
                            ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      betasquared(:,:) = 1.d4        ! Pa yr/m

!TODO - Should these 5's be hardwired?  This will give strange results in parallel.
!       Could fix by applying small value of betasquared on global domain and scattering to local.  
!TODO - Is 10.d1 correct?  (Change to 100.d0?)
      do ns=5, nsn-5
      do ew=1, ewn-1
        betasquared(ew,ns) = 10.d1      ! Pa yr/m
      end do
      end do

    case(HO_BABC_YIELD_PICARD)  ! take input value for till yield stress and force betasquared to be implemented such
                                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).

      !!! NOTE: Eventually, this option will provide the till yield stress as calculate from the basal processes
      !!! submodel. Currently, to enable sliding over plastic till, simple specify the value of "betasquared" as 
      !!! if it were the till yield stress (in units of Pascals).

       !TODO - Will we ever use mintauf?  It is passed throughout the code but never used
!      betasquared = mintauf*tau0 / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )

      betasquared(:,:) = ( beta(:,:) * ( tau0 / vel0 / scyr ) ) &     ! Pa yr/m
                         / dsqrt( (thisvel(:,:)*vel0*scyr)**2 + (othervel(:,:)*vel0*scyr)**2 + (smallnum)**2 )

    case(HO_BABC_BETA_BWAT)  ! set value of beta as proportional to value of bwat                                         

      C = 10.d0
      m = 1.d0

      allocate(unstagbetasquared(ewn,nsn))

      unstagbetasquared(:,:) = 200.d0       

      where ( bwat > 0.d0 .and. unstagbetasquared > 200.d0 )
          unstagbetasquared = C / ( bwat**m )
      endwhere

      ! average betas from unstag grid onto stag grid
      betasquared = 0.5d0 * ( unstagbetasquared(1:ewn-1,:) + unstagbetasquared(2:ewn,:) )
      betasquared = 0.5d0 * ( unstagbetasquared(:,1:nsn-1) + unstagbetasquared(:,2:nsn) )
   
      deallocate(unstagbetasquared) 

    case(HO_BABC_LARGE_BETA)      ! frozen (u=v=0) ice-bed interface

      betasquared(:,:) = 1.d10       ! Pa yr/m

    case(HO_BABC_EXTERNAL_BETA)   ! use value passed in externally from CISM (NOTE not dimensional when passed in) 

      ! scale CISM input value to dimensional units of (Pa yr/m)

      betasquared(:,:) = beta(:,:) * ( tau0 / vel0 / scyr )

      ! this is a check for NaNs, which indicate, and are replaced by no slip
      !TODO: Not sure I follow the logic of this ... keep/omit? Added by the UMT crew at some point

      do ns=1, nsn-1
      do ew=1, ewn-1 
        if( betasquared(ew,ns) /= betasquared(ew,ns) )then
          betasquared(ew,ns) = 1.d10     ! Pa yr/m
        endif 
      end do
      end do

      ! check for areas where ice is floating or grounded and make sure beta in these regions is 0  

      !TODO: Ideally, these mask values should not be hardwired, but keeping it this way for now until
      ! we decide which mask values to keep/remove

      do ns=1, nsn-1
      do ew=1, ewn-1 
        !if( ( mask(ew,ns) >= 21 .and. mask(ew,ns) <= 23 ) .or. ( mask(ew,ns) >= 41 .and. mask(ew,ns) <= 57 ) &
        !! less agressive than apply beta = 0 at g.l., which will make some test cases fail (e.g. circ. shelf)
        !! because of lack of fully grounded area.
        if( ( mask(ew,ns) >= 41 .and. mask(ew,ns) <= 43 ) &     
             .or. mask(ew,ns) == 9 .or. mask(ew,ns) == 11 )then
           betasquared(ew,ns) = 0.d0
        endif
      end do
      end do

    ! NOTE: cases (HO_BABC_NO_SLIP) and (HO_BABC_YIELD_NEWTON) are handled external to this subroutine

  end select

  ! convert the dimensional value of betasquared to non-dimensional units by dividing by scale factor.
  betasquared(:,:) = betasquared(:,:) / ( tau0 / vel0 / scyr )    !! scale in parentheses is: Pa * sec/m * yr/sec = Pa yr/m

end subroutine calcbetasquared

!***********************************************************************

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

subroutine geom2derscros(ewn,  nsn,   &
                         dew,  dns,   &
                         ipvr, stagthck, opvrewns)

  ! geometric (2nd) cross-deriv. for generic input variable 'ipvr', output as 'opvr'       

  implicit none

  integer, intent(in) :: ewn, nsn
  real(dp), intent(in) :: dew, dns
  real(dp), intent(out), dimension(:,:) :: opvrewns
  real(dp), intent(in), dimension(:,:) :: ipvr, stagthck

  integer :: ew, ns
  real(dp) :: dewdns

  dewdns = dew*dns
 
! TODO: Check this over and if ok remove old code !!
!  *SFP* OLD method; replaced (below) w/ loops and logic for compatibility w/ gnu compilers
!  where (stagthck /= 0.d0)
!    opvrewns = (eoshift(eoshift(ipvr,1,0.d0,2),1,0.d0,1) + ipvr   &
!               - eoshift(ipvr,1,0.d0,1) - eoshift(ipvr,1,0.d0,2)) / (dewdns)
!  elsewhere
!    opvrewns = 0.d0
!  end where

!  *SFP* NEW method

  opvrewns = ( ipvr(2:ewn,2:nsn) - ipvr(2:ewn,1:nsn-1) - ipvr(1:ewn-1,2:nsn) + ipvr(1:ewn-1,1:nsn-1) ) / dewdns

  do ns = 1, nsn-1
      do ew = 1, ewn-1
        if (stagthck(ew,ns) == 0.d0) then
           opvrewns(ew,ns) = 0.d0
        end if
      end do
  end do

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


!LOOP TODO - Please confirm that these are the right boundaries
! Note: Provided nhalo >= 2, we should have enough points to compute a centered difference.
!       Not sure what happens if nhalo = 1
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

   !*SFP*for now, ignoring the possibility of using JFNK w/ Trilinos ...
   if( nonlinear == HO_NONLIN_PICARD )then

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
        ! Option to load entry into Triad sparse matrix format
        if (value /= 0.d0) then
          pcgval(ct_nonzero) = value
          pcgcol(ct_nonzero) = col
          pcgrow(ct_nonzero) = row
          ct_nonzero = ct_nonzero + 1
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
 
 
   !*SFP* if using JFNK, store the main block diagonal coeffs and off diag coeffs 
   elseif ( nonlinear == HO_NONLIN_JFNK )then

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then      ! if using Triad format to store matrix entries

          ! load entry into Triad sparse matrix format
          if (value /= 0.d0) then
            pcgval(ct_nonzero) = value
            pcgcol(ct_nonzero) = col
            pcgrow(ct_nonzero) = row
            ct_nonzero = ct_nonzero + 1
          end if

#ifdef TRILINOS
    else    ! if storing matrix entires directly in Trilinos sparse format

        if (value /= 0.d0) then
           !AGS: If we find that sparsity changes inside a time step,
           !     consider adding entry even for value==0.
           call putintotrilinosmatrix(row, col, value) 
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

      ! Step through indxmask, but exclude halo

!    SFP: debug line below
!    print *, 'mySize = ', mySize

      do ns = 1+staggered_lhalo, size(indxmask,2)-staggered_uhalo
         do ew = 1+staggered_lhalo, size(indxmask,1)-staggered_uhalo
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

  end subroutine distributed_create_partition

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

!LOOP TODO: Please confirm that these are the correct loop bounds.
         ! loc2_array-based search
          minew = 1
          minns = 1
          mindiff = globalID
!          do ns = 1+staggered_lhalo,size(loc2_array,2)-staggered_uhalo
!            do ew = 1+staggered_lhalo,size(loc2_array,1)-staggered_uhalo
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

  ! *SFP* This function is used to scale the matrix coeffs and rhs vector coeff
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

subroutine assign_resid(model, uindx, umask, &
     d2thckdewdns, d2usrfdewdns, pcgsize, gx_flag, matrixA, matrixC, L2norm, ewn, nsn)
  
  
  use iso_c_binding 
  use glide_types, only : glide_global_type
  use glimmer_sparse_type, only : sparse_matrix_type
  
  implicit none
  
  type(glide_global_type)  ,intent(inout) :: model
  type(sparse_matrix_type) ,intent(in) :: matrixA, matrixC
  
  integer :: i, j
  integer                   ,intent(in) :: ewn, nsn
  integer, dimension(2)     ,intent(in) :: pcgsize
  integer                   ,intent(in) :: gx_flag(2*pcgsize(1)) ! 0 :reg cell
  integer                   ,intent(in) :: uindx(ewn-1,nsn-1), umask(ewn-1,nsn-1)
  real(dp)          ,intent(in) :: L2norm
  real(dp)          ,intent(in) :: d2thckdewdns(ewn-1,nsn-1), d2usrfdewdns(ewn-1,nsn-1)
  
!LOOP TODO: Would it be sufficient to loop only over locally owned velocity points?      
!LOOP TODO - Switch i and j to reduce strides?

  do i = 1, ewn-1 
   do j = 1, nsn-1 
    model%solver_data%ui(i,j)  = uindx(i,j)
    model%solver_data%um(i,j)  = umask(i,j)
    model%solver_data%d2thckcross(i,j) = d2thckdewdns(i,j) 
    model%solver_data%d2usrfcross(i,j) = d2usrfdewdns(i,j) 
   end do
  end do

  model%solver_data%pcgsize = pcgsize
  do i = 1, 2*pcgsize(1)
   model%solver_data%gxf(i) = gx_flag(i)
  end do
  model%solver_data%L2norm = L2norm
  model%solver_data%matrixA = matrixA
  model%solver_data%matrixC = matrixC

end subroutine assign_resid

!-------------------------------------------------------------------

!  uvec is either u^k-1 or v^k-1 on input and Av-b or Cu-d on output

subroutine res_vect ( matrix, uvec, bvec, nu, g_flag, L2square, whatsparse)

use parallel

use glimmer_paramets, only : dp
use glimmer_sparse_type
use glimmer_sparse
use glide_mask
use profile

implicit none

integer :: i, j, nu, nele, whatsparse ! nu: size of uvec and bvec
integer, dimension(nu), intent(in) :: g_flag ! 0 :reg cell
                                             ! 1 :top ghost, 2 :base ghost

type(sparse_matrix_type),  intent(in) :: matrix

real(dp), dimension(nu), intent(in) :: bvec
real(dp), dimension(nu), intent(inout) :: uvec
real(dp), dimension(nu) :: Au_b_wig
real(dp), intent(out) :: L2square
! 
real(dp) :: scale_ghosts = 0.0d0

! calculate residual vector of the u OR v component

      Au_b_wig = 0d0 ! regular+ghost cells

call t_startf("res_vect_matvec")
      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then

        do nele = 1, matrix%nonzeros 

           i = matrix%row(nele)
           j = matrix%col(nele)
           Au_b_wig(i) = Au_b_wig(i) + matrix%val(nele) * uvec(j)

        enddo

#ifdef TRILINOS
      else 
        call matvecwithtrilinos(uvec, Au_b_wig);
#endif
      endif 
call t_stopf("res_vect_matvec")

      do i = 1, nu
         Au_b_wig(i) = Au_b_wig(i) - bvec(i)
      enddo

      uvec = Au_b_wig

! AGS: Residual norm includes scaling to decrease importance of ghost values
! By calling it a redefinition of an inner product, it is kosher.
      L2square = 0.d0
      do i = 1, nu
         if (g_flag(i) == 0) then
            L2square = L2square + Au_b_wig(i) * Au_b_wig(i)
         else
            L2square = L2square + scale_ghosts * Au_b_wig(i) * Au_b_wig(i)
         endif
      end do

      !JEFF Sum L2square across nodes
call t_startf("res_vect_reduce")
      L2square = parallel_reduce_sum(L2square)
call t_stopf("res_vect_reduce")

      return

end subroutine res_vect

!-------------------------------------------------------------------

subroutine res_vect_jfnk ( matrixA, matrixC, uvec, bvec, nu1, nu2, g_flag, L2square, whatsparse)

! similar to res_vect, but state vector uvec and rhs vector bvec are now both velocities 
! A and C matrices are separate, but eventually could be combined

use glimmer_paramets, only : dp
use glimmer_sparse_type
use glimmer_sparse
use glide_mask

implicit none

integer :: i, j, nu1, nu2, nele, whatsparse ! nu2: size of uvec and bvec, size of u, v within

type(sparse_matrix_type),  intent(in) :: matrixA, matrixC

integer, dimension(nu2) :: g_flag  ! 0=reg cell, 1: top ghost, 2, base ghost
real(dp), dimension(nu2), intent(in) :: bvec
real(dp), dimension(nu2), intent(inout) :: uvec
real(dp), dimension(nu1) :: Au_b_wig, Cv_d_wig
real(dp), intent(out) :: L2square
! 
real(dp) :: scale_ghosts = 0.0d0

! calculate residual vector of the u and v component

      Au_b_wig = 0d0 ! regular+ghost cells
      Cv_d_wig = 0d0 ! regular+ghost cells

      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then

        do nele = 1, matrixA%nonzeros

           i = matrixA%row(nele)
           j = matrixA%col(nele)
           Au_b_wig(i) = Au_b_wig(i) + matrixA%val(nele) * uvec(j)

        enddo

        do nele = 1, matrixC%nonzeros

           i = matrixC%row(nele)
           j = matrixC%col(nele)
           Cv_d_wig(i) = Cv_d_wig(i) + matrixC%val(nele) * uvec(nu1+j)

        enddo

#ifdef TRILINOS
      else

        call matvecwithtrilinos(uvec(1:nu1), Au_b_wig);
        call matvecwithtrilinos(uvec(nu1+1:nu2), Cv_d_wig);
#endif
      endif

      do i = 1, nu1

         Au_b_wig(i) = Au_b_wig(i) - bvec(i)
         Cv_d_wig(i) = Cv_d_wig(i) - bvec(nu1+i)

      enddo

! to do: combine A and C

      do i = 1, nu1

         uvec(i)    = Au_b_wig(i)
         uvec(nu1+i) = Cv_d_wig(i)

      enddo

! AGS: Residual norm includes scaling to decrease importance of ghost values
! By calling it a redefinition of an inner product, it is kosher.
!      L2square = 0.0
!      do i = 1, nu1
!         if (g_flag(i) == 0) then
!            L2square = L2square + Au_b_wig(i) * Au_b_wig(i)
!         else
!            L2square = L2square + scale_ghosts * Au_b_wig(i) * Au_b_wig(i)
!         endif
!      end do
!
!      do i = 1, nu1
!         if (g_flag(nu1+i) == 0) then
!            L2square = L2square + Cv_d_wig(i) * Cv_d_wig(i)
!         else
!            L2square = L2square + scale_ghosts * Cv_d_wig(i) * Cv_d_wig(i)
!         endif
!      end do
! when the combined version is used, convergence wrong
!TODO (KJE) what is the comment above. What is wrong?

      do i = 1, nu2
         if (g_flag(i) == 0) then
            L2square = L2square + uvec(i) * uvec(i)
         else
            L2square = L2square + scale_ghosts * uvec(i) * uvec(i)
         endif
      end do


      return

end subroutine res_vect_jfnk

!-------------------------------------------------------------------

!TODO - Is this subroutine still needed?

subroutine slapsolve(xk_1, xk_size, c_ptr_to_object, NL_tol, pcgsize)

  use iso_c_binding  
  use glimmer_paramets, only : dp
  use glide_types ,only : glide_global_type
  use parallel

  implicit none

  real(dp), dimension(:), intent(out) :: xk_1
  integer(c_int) ,intent(in) ,value  :: xk_size
  type(c_ptr) ,intent(inout)         :: c_ptr_to_object
  real(dp) ,intent(in) :: NL_tol
  integer, dimension(2) :: pcgsize

  type(glide_global_type) ,pointer        :: fptr=>NULL()

  real(dp), dimension(:), allocatable :: xk_1_plus, vectx
  real(dp), dimension(:), allocatable :: dx, F, F_plus
  real(dp), dimension(:), allocatable :: wk1, wk2, rhs
  real(dp), dimension(:,:), allocatable :: vv, wk
  real(dp) :: L2norm_wig, tol, gamma_l, epsilon, NL_target
  integer :: tot_its, itenb, maxiteGMRES, iout, icode
  integer, parameter :: img = 20, img1 = img+1, kmax = 500
  integer :: k

  type(sparse_matrix_type) :: matrixA, matrixC
  real(dp) :: L2norm

  allocate( vectx(2*pcgsize(1)), xk_1_plus(2*pcgsize(1)) )
  allocate( F(2*pcgsize(1)), F_plus(2*pcgsize(1)), dx(2*pcgsize(1)) )
  allocate( wk1(2*pcgsize(1)), wk2(2*pcgsize(1)), rhs(2*pcgsize(1)) )
  allocate( vv(2*pcgsize(1),img1), wk(2*pcgsize(1),img) )

! Iteration loop

  do k = 1, kmax

    call calc_F (xk_1, F, xk_size, c_ptr_to_object, 0)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    L2norm = fptr%solver_data%L2norm
    matrixA = fptr%solver_data%matrixA
    matrixC = fptr%solver_data%matrixC

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

 10 CONTINUE
! icode = 0 means that fgmres has finished and sol contains the app. solution
      
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

!------------------------------------------------------------------------
! Update solution vectors (x^k = x^k-1 + dx) and 3D fields 
!------------------------------------------------------------------------
      xk_1 = xk_1 + dx(1:2*pcgsize(1))

 end do   ! k = 1, kmax 

  deallocate(dx, vectx, xk_1_plus)
  deallocate(F, F_plus, rhs)
  deallocate(wk1, wk2)
  deallocate(vv, wk)

end subroutine slapsolve

!-----------------------------------------------------------------------

!TODO - Is this subroutine still needed?  It's called from slapsolve above.

  subroutine fgmres (n,im,rhs,sol,i,vv,w,wk1, wk2, &
                     eps,maxits,iout,icode,its) 

! JFL to be removed

!-----------------------------------------------------------------------
! jfl Dec 1st 2006. We modified the routine so that it is double precison.
! Here are the modifications:
! 1) implicit real (a-h,o-z) becomes implicit real*8 (a-h,o-z) 
! 2) real bocomes real*8
! 3) subroutine scopy.f has been changed for dcopy.f
! 4) subroutine saxpy.f has been changed for daxpy.f
! 5) function sdot.f has been changed for ddot.f
! 6) 1e-08 becomes 1d-08
!
! Be careful with the dcopy, daxpy and ddot code...there is a slight 
! difference with the single precision versions (scopy, saxpy and sdot).
! In the single precision versions, the array are declared sightly differently.
! It is written for single precision:
!
! modified 12/3/93, array(1) declarations changed to array(*)
!-----------------------------------------------------------------------

      implicit double precision (a-h,o-z) !jfl modification
      integer n, im, maxits, iout, icode
      double precision rhs(*), sol(*), vv(n,im+1),w(n,im)
      double precision wk1(n), wk2(n), eps
!-----------------------------------------------------------------------
! flexible GMRES routine. This is a version of GMRES which allows a 
! a variable preconditioner. Implemented with a reverse communication 
! protocole for flexibility -
! DISTRIBUTED VERSION (USES DISTDOT FOR DDOT) 
! explicit (exact) residual norms for restarts  
! written by Y. Saad, modified by A. Malevsky, version February 1, 1995
!-----------------------------------------------------------------------
! This Is A Reverse Communication Implementation. 
!------------------------------------------------- 
! USAGE: (see also comments for icode below). FGMRES
! should be put in a loop and the loop should be active for as
! long as icode is not equal to 0. On return fgmres will
!    1) either be requesting the new preconditioned vector applied
!       to wk1 in case icode==1 (result should be put in wk2) 
!    2) or be requesting the product of A applied to the vector wk1
!       in case icode==2 (result should be put in wk2) 
!    3) or be terminated in case icode == 0. 
! on entry always set icode = 0. So icode should be set back to zero
! upon convergence.
!-----------------------------------------------------------------------
! Here is a typical way of running fgmres: 
!
!      icode = 0
! 1    continue
!      call fgmres (n,im,rhs,sol,i,vv,w,wk1, wk2,eps,maxits,iout,icode)
!
!      if (icode == 1) then
!         call  precon(n, wk1, wk2)    <--- user's variable preconditioning
!         goto 1
!      else if (icode >= 2) then
!         call  matvec (n,wk1, wk2)    <--- user's matrix vector product. 
!         goto 1
!      else 
!         ----- done ---- 
!         .........
!-----------------------------------------------------------------------
! list of parameters 
!------------------- 
!
! n     == integer. the dimension of the problem
! im    == size of Krylov subspace:  should not exceed 50 in this
!          version (can be reset in code. looking at comment below)
! rhs   == vector of length n containing the right hand side
! sol   == initial guess on input, approximate solution on output
! vv    == work space of size n x (im+1)
! w     == work space of length n x im 
! wk1,
! wk2,  == two work vectors of length n each used for the reverse
!          communication protocole. When on return (icode \= 1)
!          the user should call fgmres again with wk2 = precon * wk1
!          and icode untouched. When icode==1 then it means that
!          convergence has taken place.
!          
! eps   == tolerance for stopping criterion. process is stopped
!          as soon as ( ||.|| is the euclidean norm):
!          || current residual||/||initial residual|| <= eps
!
! maxits== maximum number of iterations allowed
!
! iout  == output unit number number for printing intermediate results
!          if (iout <= 0) no statistics are printed.
! 
! icode = integer. indicator for the reverse communication protocole.
!         ON ENTRY : icode should be set to icode = 0.
!         ON RETURN: 
!       * icode == 1 value means that fgmres has not finished
!         and that it is requesting a preconditioned vector before
!         continuing. The user must compute M**(-1) wk1, where M is
!         the preconditioing  matrix (may vary at each call) and wk1 is
!         the vector as provided by fgmres upun return, and put the 
!         result in wk2. Then fgmres must be called again without
!         changing any other argument. 
!       * icode == 2 value means that fgmres has not finished
!         and that it is requesting a matrix vector product before
!         continuing. The user must compute  A * wk1, where A is the
!         coefficient  matrix and wk1 is the vector provided by 
!         upon return. The result of the operation is to be put in
!         the vector wk2. Then fgmres must be called again without
!         changing any other argument. 
!       * icode == 0 means that fgmres has finished and sol contains 
!         the approximate solution.
!         comment: typically fgmres must be implemented in a loop
!         with fgmres being called as long icode is returned with 
!         a value \= 0. 
!-----------------------------------------------------------------------
!     local variables -- !jfl modif
      double precision hh(201,200),c(200),s(200),rs(201),t,ro,ddot,sqrt 
!
!-------------------------------------------------------------
!     arnoldi size should not exceed 50 in this version..
!     to reset modify sizes of hh, c, s, rs       
!-------------------------------------------------------------

      save
      data epsmac/1.d-16/

      !WHL - added integer declarations
      integer :: i, its, i1, ii, j, jj, k, k1, n1
!     
!     computed goto 
!     
      goto (100,200,300,11) icode +1
 100  continue
      n1 = n + 1
      its = 0
!-------------------------------------------------------------
!     **  outer loop starts here..
!--------------compute initial residual vector --------------
! 10   continue
      call dcopy (n, sol, 1, wk1, 1) !jfl modification
      icode = 3
      return
 11   continue
      do j=1,n
         vv(j,1) = rhs(j) - wk2(j) 
      enddo
 20   ro = ddot(n, vv, 1, vv,1) !jfl modification
      ro = sqrt(ro)
      if (ro == 0.0d0) goto 999 
      t = 1.0d0/ ro 
      do j=1, n
         vv(j,1) = vv(j,1)*t 
      enddo
      if (its == 0) eps1=eps
      if (its == 0) r0 = ro
      if (iout > 0) write(*, 199) its, ro!&
!           print *,'chau',its, ro !write(iout, 199) its, ro
!     
!     initialize 1-st term  of rhs of hessenberg system..
!     
      rs(1) = ro
      i = 0
 4    i=i+1
      its = its + 1
      i1 = i + 1
      do k=1, n
         wk1(k) = vv(k,i) 
      enddo
!     
!     return
!     
      icode = 1

      return
 200  continue
      do k=1, n
         w(k,i) = wk2(k) 
      enddo
!     
!     call matvec operation
!     
      icode = 2
      call dcopy(n, wk2, 1, wk1, 1) !jfl modification
!
!     return
!     
      return
 300  continue
!     
!     first call to ope corresponds to intialization goto back to 11.
!     
!      if (icode == 3) goto 11
      call  dcopy (n, wk2, 1, vv(1,i1), 1) !jfl modification
!     
!     modified gram - schmidt...
!     
      do j=1, i
         t = ddot(n, vv(1,j), 1, vv(1,i1), 1) !jfl modification
         hh(j,i) = t
         call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1) !jfl modification
      enddo
      t = sqrt(ddot(n, vv(1,i1), 1, vv(1,i1), 1)) !jfl modification
      hh(i1,i) = t
      if (t == 0.0d0) goto 58
      t = 1.0d0 / t
      do k=1,n
         vv(k,i1) = vv(k,i1)*t
      enddo
!     
!     done with modified gram schimd and arnoldi step. 
!     now  update factorization of hh
!     
 58   if (i == 1) goto 121
!     
!     perfrom previous transformations  on i-th column of h
!     
      do k=2,i
         k1 = k-1
         t = hh(k1,i)
         hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
         hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
      enddo
 121  gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
      if (gam == 0.0d0) gam = epsmac
!-----------#determine next plane rotation  #-------------------
      c(i) = hh(i,i)/gam
      s(i) = hh(i1,i)/gam
      rs(i1) = -s(i)*rs(i)
      rs(i) =  c(i)*rs(i)
!     
!     determine res. norm. and test for convergence-
!     
      hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
      ro = abs(rs(i1))
      if (iout > 0) &
           write(*, 199) its, ro
      if (i < im .and. (ro > eps1))  goto 4
!     
!     now compute solution. first solve upper triangular system.
!     
      rs(i) = rs(i)/hh(i,i)
      do ii=2,i
         k=i-ii+1
         k1 = k+1
         t=rs(k)
         do j=k1,i
            t = t-hh(k,j)*rs(j)
         enddo
         rs(k) = t/hh(k,k)
      enddo
!     
!     done with back substitution..
!     now form linear combination to get solution
!     
      do j=1, i
         t = rs(j)
         call daxpy(n, t, w(1,j), 1, sol,1) !jfl modification
      enddo
!     
!     test for return 
!     
      if (ro <= eps1 .or. its >= maxits) goto 999
!     
!     else compute residual vector and continue..
!     
!       goto 10

     do j=1,i
        jj = i1-j+1
        rs(jj-1) = -s(jj-1)*rs(jj)
        rs(jj) = c(jj-1)*rs(jj)
     enddo
     do j=1,i1
        t = rs(j)
        if (j == 1)  t = t-1.0d0
        call daxpy (n, t, vv(1,j), 1,  vv, 1)
     enddo
!     
!     restart outer loop.
!     
     goto 20
 999  icode = 0

 199  format('   -- fmgres its =', i4, ' res. norm =', d26.16)
!     
      return 

   end subroutine fgmres
!-----------------------------------------------------------------------

!***********************************************************************************************
!BELOW here are deprecated boundary condition subroutines that have been replaced by newer 
! ones (using one sided differences) or slightly altered ones.
!***********************************************************************************************

!***********************************************************************************************
!NOTE: This subroutine has been deprecated because it is has been replaced by
! 'normhorizmainbcos', where the "os" stands for one-sided difference.
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
!NOTE: This subroutine has been deprecated because it is has been replaced by
! 'croshorizmainbcos', where the "os" stands for one-sided difference.
function croshorizmainbc(dew,       dns,       &
                         dusrfdew,  dusrfdns,  &
                         dsigmadew, dsigmadns, &
                         which,     bcflag,    &
                         dup,       local_othervel,  &
                         efvs,                       &
                         oneortwo,  twoorone,        &
                         g_cros, velbc )

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
    real(dp), intent(in), optional :: velbc
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

    croshorizmainbc = g

    return

end function croshorizmainbc

!***********************************************************************************************
!ABOVE here are deprecated boundary condition subroutines that have been replaced by newer 
! ones (using one sided differences) or slightly altered ones.
!***********************************************************************************************
 

end module glam_strs2

!!!***********************************************************************
