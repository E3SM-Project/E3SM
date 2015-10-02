module fv_prints
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: fv_prints --- print maxima and minima of dycore varibles
!
! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
      use perf_mod
      use cam_logfile,  only: iulog
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC     fv_out
!
! !DESCRIPTION:
!
!   This module provides basic utilities to evaluate the dynamics state
!
! !REVISION HISTORY:
!   00.08.01   Lin     Creation
!   01.01.05   Boville Modifications
!   01.03.26   Sawyer  Added ProTex documentation
!   03.04.17   Sawyer  Bug fix: pls=pls/2*plon instead of 2*plat (Boville)
!   05.07.06   Sawyer  Simplified interface with grid
!   06.02.21   Sawyer  Converted to XY decomposition
!   06.07.01   Sawyer  Transitioned tracers q3 to T_TRACERS
!   06.09.10   Sawyer  Isolated magic numbers with F90 parameters
!   08.07.03   Worley  Introduced repro_sum logic
!   12.10.29   Santos  repro_sum_mod is now shr_reprosum_mod
!
!EOP
!-------------------------------------------------------------------------

private
  real(r8), parameter ::  D0_0                    =   0.0_r8
  real(r8), parameter ::  D0_01                   =   0.01_r8
  real(r8), parameter ::  D1_0                    =   1.0_r8
  real(r8), parameter ::  D2_0                    =   2.0_r8
  real(r8), parameter ::  D864_0                  = 864.0_r8
  real(r8), parameter ::  G_EARTH                 = 9.80616_r8
  real(r8), parameter ::  SECS_PER_1000_DAYS      = 86400000.0_r8

CONTAINS

!-------------------------------------------------------------------------
!BOP
! !IROUTINE: fv_out --- Write out maxima and minima of dynamics state
!
! !INTERFACE: 
  subroutine  fv_out( grid,   pk,    pt, ptop,       ps,                  &
                      tracer, delp,  pe, surf_state, phys_state,          &
                      ncdate, ncsec, full_phys  )

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use dynamics_vars,  only : T_FVDYCORE_GRID
    use ppgrid,         only: begchunk, endchunk, pcols, pver
    use phys_grid,      only: get_ncols_p
    use physics_types,  only: physics_state
    use camsrfexch,     only: cam_out_t
    use constituents,   only: cnst_name
#if defined( SPMD )
    use parutilitiesmodule, only : sumop, parcollective
    use mpishorthand, only: mpicom
#endif
    use shr_reprosum_mod, only : shr_reprosum_calc, shr_reprosum_tolExceeded
                              
    use phys_gmean,    only : gmean

    implicit none

! !INPUT PARAMETERS:
    type (T_FVDYCORE_GRID), intent(in) :: grid

    integer ncdate                      ! Date
    integer ncsec                       ! Time

    real(r8) :: ptop                       ! Pressure at top
! Surface pressure
    real(r8) :: ps(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
! Pe**kappa
    real(r8) :: pk(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km+1)
! Potential temperature
    real(r8) :: pt(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
! Layer thickness (pint(k+1) - pint(k))
    real(r8) :: delp(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
! Tracers
    real(r8), intent(inout) ::   &
        tracer(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km,grid%ntotq)
! Edge pressure
    real(r8) ::  pe(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy)

    type(cam_out_t),     intent(in), dimension(begchunk:endchunk) :: surf_state

    type(physics_state), intent(in), dimension(begchunk:endchunk) :: phys_state
    logical full_phys                   ! Full physics on?

!
! !DESCRIPTION:
!
!   Determine maxima and minima of dynamics state and write them out
!
! !REVISION HISTORY:
!   00.08.01   Lin     Creation
!   01.01.05   Boville Modifications
!   01.03.26   Sawyer  Added ProTex documentation
!   01.06.27   Mirin   Converted to 2D yz decomposition
!   01.12.18   Mirin   Calculate average height (htsum) metric
!   02.02.13   Eaton   Pass precc and precl via cam_out_t type
!   05.07.06   Sawyer  Simplified interface with grid
!   06.02.21   Sawyer  Converted to XY decomposition
!   06.07.01   Sawyer  Transitioned tracers q3 to T_TRACERS
!   08.07.03   Worley  Introduced repro_sum and gmean logic
!   12.10.2=   Santos  repro_sum is now shr_reprosum_mod
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    integer i, j, k, ic, nj, lchnk, nck, ncol
    real(r8), dimension(begchunk:endchunk)    :: pmax, tmax, umax, vmax, wmax
    real(r8), dimension(begchunk:endchunk)    :: pmin, tmin, umin, vmin, wmin
    real(r8), dimension(pcols,begchunk:endchunk,1) :: precc ! convective precip rate
    real(r8), dimension(pcols,begchunk:endchunk,1) :: precl ! large-scale precip rate
    real(r8), dimension(begchunk:endchunk)    :: preccmax, preclmax
    real(r8), dimension(begchunk:endchunk)    :: preccmin, preclmin
    real(r8) :: fac, precmax, precmin
    real(r8) :: pcon(1), pls(1)
    real(r8) :: p1, p2, dtmp, apcon, htsum(1)
    real(r8), pointer :: qtmp(:,:,:)
    real(r8) :: htg(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
    real(r8) :: rel_diff(2)
    
    integer :: im, jm, km, ifirstxy, ilastxy, jfirstxy, jlastxy
    integer :: itot, jtot, ltot

    integer :: ntotq                     ! No. of total tracers
    integer :: iam

    integer n, nhmsf

    logical  :: write_warning, exceeded

! statement function for hour minutes seconds of day
    nhmsf(n)  = n/3600*10000 + mod(n,3600 )/ 60*100 + mod(n, 60)

! Initialize variables from grid (for convenience)

    im      = grid%im
    jm      = grid%jm
    km      = grid%km
    ifirstxy= grid%ifirstxy
    ilastxy = grid%ilastxy
    jfirstxy= grid%jfirstxy
    jlastxy = grid%jlastxy
    ntotq   = grid%ntotq

    itot    = (ilastxy-ifirstxy) + 1
    jtot    = (jlastxy-jfirstxy) + 1
    ltot    = itot*jtot

    iam     = grid%iam

    if (iam == 0) then
       write(iulog,*) ' '
       write(iulog,*) nhmsf(ncsec), ncdate
    endif

!
! Check total air and dry air mass.

    call dryairm( grid, .true.,  ps,    tracer,  delp,     &
                  pe,   .true.)

!$omp parallel do private(lchnk, ncol)
    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       pmax(lchnk) = maxval(phys_state(lchnk)%ps(1:ncol))
       pmin(lchnk) = minval(phys_state(lchnk)%ps(1:ncol))
       tmax(lchnk) = maxval(phys_state(lchnk)%t(1:ncol,1:pver))
       tmin(lchnk) = minval(phys_state(lchnk)%t(1:ncol,1:pver))
       umax(lchnk) = maxval(phys_state(lchnk)%u(1:ncol,1:pver))
       umin(lchnk) = minval(phys_state(lchnk)%u(1:ncol,1:pver))
       vmax(lchnk) = maxval(phys_state(lchnk)%v(1:ncol,1:pver))
       vmin(lchnk) = minval(phys_state(lchnk)%v(1:ncol,1:pver))
       wmax(lchnk) = maxval(phys_state(lchnk)%omega(1:ncol,1:pver))
       wmin(lchnk) = minval(phys_state(lchnk)%omega(1:ncol,1:pver))
    end do

#if defined( SPMD )
    nck = endchunk - begchunk + 1
    call pmaxmin2('PS',         pmin, pmax, nck, D0_01, mpicom)
    call pmaxmin2('U ',         umin, umax, nck, D1_0, mpicom)
    call pmaxmin2('V ',         vmin, vmax, nck, D1_0, mpicom)
    call pmaxmin2('T ',         tmin, tmax, nck, D1_0, mpicom)
    call pmaxmin2('W (mb/day)', wmin, wmax, nck, D864_0, mpicom)
#endif

#if 0
!
! This code is currently inactive:  the maxima and minima were not
! being used
!
    nj = (jlastxy - jfirstxy + 1) * (ilastxy - ifirstxy + 1)
    do ic=1,ntotq
       qtmp => tracer(:,:,:,ic)
       call pmaxmin(cnst_name(ic), qtmp, p1, p2, nj, km, D1_0, grid%commxy)
!
! Do something with p1 and p2?
!
    end do
#endif

!
! Calculate the vertically integrated heights 
!
    htg(:,:) = D0_0
    apcon = D1_0/G_EARTH

!$omp parallel do private(i, j, k)
    do j=jfirstxy,jlastxy
      do k=1,km
        do i=ifirstxy,ilastxy
          htg(i,j) = htg(i,j) + apcon * pt(i,j,k) * (pk(i,j,k+1)-pk(i,j,k))
        enddo
      enddo
    enddo

!$omp parallel do private(i, j, k)
    do j=jfirstxy,jlastxy
       do i=ifirstxy,ilastxy
          htg(i,j) = htg(i,j)*grid%cosp(j)
       enddo
    enddo

    call t_startf("fv_out_reprosum")
    call shr_reprosum_calc(htg, htsum, ltot, ltot, 1, gbl_count=im*jm, &
                   commid=grid%commxy, rel_diff=rel_diff)
    call t_stopf("fv_out_reprosum")

    ! check that "fast" reproducible sum is accurate enough.
    ! NOTE: not recomputing if difference too large. This
    !  value is output only, so does not feed back into the
    !  simulation
    write_warning = .false.
    if (iam == 0) write_warning = .true.
    exceeded = shr_reprosum_tolExceeded('fv_out', 1, write_warning, &
                                      iulog, rel_diff)
    
    if (iam == 0) then
      htsum(1) = htsum(1) / (D2_0*im)
      write(iulog,*) 'Average Height (geopotential units) = ', htsum(1)
    endif

    if ( .not. full_phys ) return

! Global means:

    fac = SECS_PER_1000_DAYS                     ! convert to mm/day

!$omp parallel do private(lchnk, ncol)
    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       precc(:ncol,lchnk,1) = surf_state(lchnk)%precc(:ncol)
       precl(:ncol,lchnk,1) = surf_state(lchnk)%precl(:ncol)
       preccmax(lchnk) = maxval(precc(1:ncol,lchnk,1))
       preccmin(lchnk) = minval(precc(1:ncol,lchnk,1))
       preclmax(lchnk) = maxval(precl(1:ncol,lchnk,1))
       preclmin(lchnk) = minval(precl(1:ncol,lchnk,1))
    end do

#if defined( SPMD )
    nck = endchunk - begchunk + 1
    call pmaxmin2('PRECC', preccmin, preccmax, nck, fac, mpicom)
    call pmaxmin2('PRECL', preclmin, preclmax, nck, fac, mpicom)
#endif

    call gmean(precc,pcon,1)
    call gmean(precl,pls,1)

    if (iam == 0) then
       pcon(1) = pcon(1) * fac
       pls(1)  = pls(1)  * fac
       write(iulog,*) 'Total precp=',pcon(1)+pls(1), &
                      ' CON=', pcon(1),' LS=',pls(1)
       write(iulog,*) ' '
    endif

!EOC
  end subroutine fv_out
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pmaxmin --- Find and print the maxima and minima of a field
!
! !INTERFACE: 
  subroutine pmaxmin( qname, a, pmin, pmax, im, jm, fac, commun )

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid==0)
    use parutilitiesmodule, only : gid, maxop, parcollective
#else
#define CPP_PRT_PREFIX
#endif
    implicit none

! !INPUT PARAMETERS:
    character*(*)  qname             ! Name of field
    integer  im                      ! Total longitudes
    integer  jm                      ! Total latitudes
    integer commun                   ! Communicator
    real(r8) a(im,jm)                ! 2D field
    real(r8) fac                     ! multiplication factor

! !OUTPUT PARAMETERS:
    real(r8) pmax                    ! Field maximum
    real(r8) pmin                    ! Field minimum

! !DESCRIPTION:
!
!   Parallelized utility routine for computing/printing global 
!   max/min from input lists of max/min's (usually for each latitude).  
! 
! !REVISION HISTORY:
!   00.03.01   Lin     Creation
!   00.05.01   Mirin   Coalesce variables to minimize collective ops
!   01.08.05   Sawyer  Modified to use parcollective
!   01.03.26   Sawyer  Added ProTex documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    integer  i, j
    real(r8) qmin(jm), qmax(jm)
    real(r8) pm(2)

!$omp  parallel do default(shared) private(i,j, pmax, pmin)

    do j=1,jm
       pmax = a(1,j)
       pmin = a(1,j)
       do i=2,im
          pmax = max(pmax, a(i,j))
          pmin = min(pmin, a(i,j))
       enddo
       qmax(j) = pmax
       qmin(j) = pmin
    enddo
!
! Now find max/min of qmax/qmin
!
    pmax = qmax(1)
    pmin = qmin(1)
    do j=2,jm
       pmax = max(pmax, qmax(j))
       pmin = min(pmin, qmin(j))
    enddo

#if defined( SPMD )
    pm(1) = pmax
    pm(2) = -pmin
    call parcollective( commun, maxop, 2, pm )
    pmax = pm(1)
    pmin = -pm(2)
#endif

    CPP_PRT_PREFIX write(iulog,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

    return
!EOC
  end subroutine pmaxmin
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pmaxmin2 --- Find and print the maxima and minima of 1-D array
!
! !INTERFACE: 
  subroutine pmaxmin2( qname, qmin, qmax, nj, fac, commun )

! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid==0)
    use parutilitiesmodule, only : gid, maxop, parcollective
#else
#define CPP_PRT_PREFIX
#endif
    implicit none

! !INPUT PARAMETERS:
    character*(*)  qname
    integer nj
    integer commun
    real(r8), intent(in), dimension(nj) :: qmax, qmin      ! Fields
    real(r8) fac                     ! multiplication factor

! !DESCRIPTION:
!
!   Parallelized utility routine for computing/printing global max/min from 
!   input lists of max/min's (usually for each latitude). The primary purpose 
!   is to allow for the original array and the input max/min arrays to be 
!   distributed across nodes.
! 
! !REVISION HISTORY:
!   00.10.01   Lin     Creation from pmaxmin
!   01.03.26   Sawyer  Added ProTex documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    real(r8) pm(2)
    real(r8) pmin, pmax

    pmax = maxval(qmax)
    pmin = minval(qmin)

#if defined( SPMD )
    pm(1) = pmax
    pm(2) = -pmin
    call parcollective( commun, maxop, 2, pm )
    pmax = pm(1)
    pmin = -pm(2)
#endif

    CPP_PRT_PREFIX write(iulog,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

    return
!EOC
  end subroutine pmaxmin2
!-----------------------------------------------------------------------

end module fv_prints

