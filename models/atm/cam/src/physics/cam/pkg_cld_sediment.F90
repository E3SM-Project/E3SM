#undef OLDLIQSED
module pkg_cld_sediment

!---------------------------------------------------------------------------------
! Purpose:
!
! Contains routines to compute tendencies from sedimentation of cloud liquid and 
! ice particles
!
! Author: Byron Boville  Sept 19, 2002 from code by P. J. Rasch
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use spmd_utils,    only: masterproc
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit, latvap, latice, rair, rhoh2o
  use cldwat,        only: icritc
  use pkg_cldoptics, only: reitab, reltab
  use abortutils,    only: endrun
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: cld_sediment_readnl, cld_sediment_vel, cld_sediment_tend


  real (r8), parameter :: vland  = 1.5_r8            ! liquid fall velocity over land  (cm/s)
  real (r8), parameter :: vocean = 2.8_r8            ! liquid fall velocity over ocean (cm/s)
  real (r8), parameter :: mxsedfac   = 0.99_r8       ! maximum sedimentation flux factor

  logical,   parameter :: stokes = .true.         ! use Stokes velocity instead of McFarquhar and Heymsfield

! parameter for modified McFarquhar and Heymsfield
  real (r8), parameter :: vice_small = 1._r8         ! ice fall velocity for small concentration (cm/s)

! parameters for Stokes velocity
  real (r8), parameter :: eta =  1.7e-5_r8           ! viscosity of air (kg m / s)
  real (r8), parameter :: r40 =  40._r8              !  40 micron radius
  real (r8), parameter :: r400= 400._r8              ! 400 micron radius
  real (r8), parameter :: v400= 1.00_r8              ! fall velocity of 400 micron sphere (m/s)
  real (r8)            :: v40 ! = (2._r8/9._r8) * rhoh2o * gravit/eta * r40**2 * 1.e-12_r8  
                                                     ! Stokes fall velocity of 40 micron sphere (m/s)
  real (r8)            :: vslope !  = (v400 - v40)/(r400 -r40) ! linear slope for large particles m/s/micron

  ! namelist variables
  real(r8) :: cldsed_ice_stokes_fac = huge(1._r8)    ! factor applied to the ice fall velocity computed from 
                                                     ! stokes terminal velocity

!===============================================================================
contains
!===============================================================================

subroutine cld_sediment_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cld_sediment_readnl'

   namelist /cldsed_nl/ cldsed_ice_stokes_fac
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cldsed_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cldsed_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      write(iulog,*) subname//': cldsed_ice_stokes_fac = ', cldsed_ice_stokes_fac

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(cldsed_ice_stokes_fac, 1, mpir8, 0, mpicom)
#endif

end subroutine cld_sediment_readnl

!===============================================================================

  subroutine cld_sediment_vel (ncol,                               &
       icefrac , landfrac, ocnfrac , pmid    , pdel    , t       , &
       cloud   , cldliq  , cldice  , pvliq   , pvice   , landm, snowh)

!----------------------------------------------------------------------

! Compute gravitational sedimentation velocities for cloud liquid water
! and ice, based on Lawrence and Crutzen (1998).

! LIQUID

! The fall velocities assume that droplets have a gamma distribution
! with effective radii for land and ocean as assessed by Han et al.;
! see Lawrence and Crutzen (1998) for a derivation.

! ICE

! The fall velocities are based on data from McFarquhar and Heymsfield
! or on Stokes terminal velocity for spheres and the effective radius.

! NEED TO BE CAREFUL - VELOCITIES SHOULD BE AT THE *LOWER* INTERFACE
! (THAT IS, FOR K+1), FLUXES ARE ALSO AT THE LOWER INTERFACE (K+1), 
! BUT MIXING RATIOS ARE AT THE MIDPOINTS (K)...

! NOTE THAT PVEL IS ON PVERP (INTERFACES), WHEREAS VFALL IS FOR THE CELL
! AVERAGES (I.E., MIDPOINTS); ASSUME THE FALL VELOCITY APPLICABLE TO THE 
! LOWER INTERFACE (K+1) IS THE SAME AS THAT APPLICABLE FOR THE CELL (V(K))

!-----------------------------------------------------------------------
!     MATCH-MPIC version 2.0, Author: mgl, March 1998
! adapted by P. J. Rasch
!            B. A. Boville, September 19, 2002
!            P. J. Rasch    May 22, 2003 (added stokes flow calc for liquid
!                                         drops based on effect radii)
!-----------------------------------------------------------------------


! Arguments
    integer, intent(in) :: ncol                     ! number of colums to process

    real(r8), intent(in)  :: icefrac (pcols)        ! sea ice fraction (fraction)
    real(r8), intent(in)  :: landfrac(pcols)        ! land fraction (fraction)
    real(r8), intent(in)  :: ocnfrac (pcols)        ! ocean fraction (fraction)
    real(r8), intent(in)  :: pmid  (pcols,pver)     ! pressure of midpoint levels (Pa)
    real(r8), intent(in)  :: pdel  (pcols,pver)     ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: cloud (pcols,pver)     ! cloud fraction (fraction)
    real(r8), intent(in)  :: t     (pcols,pver)     ! temperature (K)
    real(r8), intent(in)  :: cldliq(pcols,pver)     ! cloud water, liquid (kg/kg)
    real(r8), intent(in)  :: cldice(pcols,pver)     ! cloud water, ice    (kg/kg)
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)

    real(r8), intent(out) :: pvliq (pcols,pverp)    ! vertical velocity of cloud liquid drops (Pa/s)
    real(r8), intent(out) :: pvice (pcols,pverp)    ! vertical velocity of cloud ice particles (Pa/s)
    real(r8), intent(in) :: landm(pcols)            ! land fraction ramped over water
! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

! Local variables
    real (r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real (r8) :: vfall                              ! settling velocity of cloud particles (m/s)
    real (r8) :: icice                              ! in cloud ice water content (kg/kg)
    real (r8) :: iciwc                              ! in cloud ice water content in g/m3
    real (r8) :: icefac
    real (r8) :: logiwc

    real (r8) :: rei(pcols,pver)                    ! effective radius of ice particles (microns)
    real (r8) :: rel(pcols,pver)                    ! effective radius of liq particles (microns)
    real(r8)  pvliq2 (pcols,pverp)    ! vertical velocity of cloud liquid drops (Pa/s)

    integer i,k

    real (r8) :: lbound, ac, bc, cc

!-----------------------------------------------------------------------
!------- initialize linear ramp variables for fall velocity ------------
!-----------------------------------------------------------------------

  v40 = (2._r8/9._r8) * rhoh2o * gravit/eta * r40**2 * 1.e-12_r8  
  vslope = (v400 - v40)/(r400 -r40)

!-----------------------------------------------------------------------
!--------------------- liquid fall velocity ----------------------------
!-----------------------------------------------------------------------

    ! compute air density
    rho(:ncol,:) = pmid(:ncol,:) / (rair * t(:ncol,:))

    pvliq(:ncol,:) = 0._r8

    ! get effective radius of liquid drop
    call reltab(ncol, t, landfrac, landm, icefrac, rel, snowh)

    do k = 1,pver
       do i = 1,ncol
          if (cloud(i,k) > 0._r8 .and. cldliq(i,k) > 0._r8) then

#ifdef OLDLIQSED
! oldway
             ! merge the liquid fall velocities for land and ocean (cm/s)
             ! SHOULD ALSO ACCOUNT FOR ICEFRAC
             vfall = vland*landfrac(i) + vocean*(1._r8-landfrac(i))
!!$          vfall = vland*landfrac(i) + vocean*ocnfrac(i) + vseaice*icefrac(i)

             ! convert the fall speed to pressure units, but do not apply the traditional
             ! negative convention for pvel.
             pvliq(i,k+1) = vfall     &
                  * 0.01_r8                 & ! cm to meters
                  * rho(i,k)*gravit        ! meters/sec to pascals/sec
#else

! newway
             if (rel(i,k) < 40._r8 ) then
                vfall = 2._r8/9._r8 * rhoh2o * gravit * rel(i,k)**2 / eta  * 1.e-12_r8  ! micons^2 -> m^2
             else
                vfall = v40 + vslope * (rel(i,k)-r40)      ! linear above 40 microns
             end if
             ! convert the fall speed to pressure units
             ! but do not apply the traditional
             ! negative convention for pvel.
!             pvliq2(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
             pvliq(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
#endif
          end if
       end do
    end do

    !-----------------------------------------------------------------------
    !--------------------- ice fall velocity -------------------------------
    !-----------------------------------------------------------------------

    pvice(:ncol,:) = 0._r8

    if (stokes) then

       !-----------------------------------------------------------------------
       !--------------------- stokes terminal velocity < 40 microns -----------
       !-----------------------------------------------------------------------

       ! get effective radius
       call reitab(ncol, t, rei)

       do k = 1,pver
          do i = 1,ncol
             if (cloud(i,k) > 0._r8 .and. cldice(i,k) > 0._r8) then
                if (rei(i,k) < 40._r8 ) then
                   vfall = 2._r8/9._r8 * rhoh2o * gravit * rei(i,k)**2 / eta  * 1.e-12_r8  ! micons^2 -> m^2
                   vfall = vfall * cldsed_ice_stokes_fac
                else
                   vfall = v40 + vslope * (rei(i,k)-r40)      ! linear above 40 microns
                end if

                ! convert the fall speed to pressure units, but do not apply the traditional
                ! negative convention for pvel.
                pvice(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
             end if
          end do
       end do

    else

       !-----------------------------------------------------------------------
       !--------------------- McFarquhar and Heymsfield > icritc --------------
       !-----------------------------------------------------------------------

       ! lower bound for iciwc

       cc = 128.64_r8 
       bc = 53.242_r8
       ac = 5.4795_r8
       lbound = (-bc + sqrt(bc*bc-4*ac*cc))/(2*ac)
       lbound = 10._r8**lbound

       do k = 1,pver
          do i = 1,ncol
             if (cloud(i,k) > 0._r8 .and. cldice(i,k) > 0._r8) then

                ! compute the in-cloud ice concentration (kg/kg)
                icice = cldice(i,k) / cloud(i,k)

                ! compute the ice water content in g/m3
                iciwc = icice * rho(i,k) * 1.e3_r8

                ! set the fall velocity (cm/s) to depend on the ice water content in g/m3,
                if (iciwc > lbound) then ! need this because of log10
                   logiwc = log10(iciwc)
                   !          Median - 
                   vfall = 128.64_r8 + 53.242_r8*logiwc + 5.4795_r8*logiwc**2
                   !          Average - 
                !!$             vfall = 122.63 + 44.111*logiwc + 4.2144*logiwc**2
                else
                   vfall = 0._r8
                end if

                ! set ice velocity to 1 cm/s if ice mixing ratio < icritc, ramp to value
                ! calculated above at 2*icritc
                if (icice <= icritc) then
                   vfall = vice_small
                else if(icice < 2*icritc) then
                   icefac = (icice-icritc)/icritc
                   vfall = vice_small * (1._r8-icefac) + vfall * icefac
                end if

                ! bound the terminal velocity of ice particles at high concentration
                vfall = min(100.0_r8, vfall)

                ! convert the fall speed to pressure units, but do not apply the traditional
                ! negative convention for pvel.
                pvice(i,k+1) = vfall     &
                   * 0.01_r8                 & ! cm to meters
                   * rho(i,k)*gravit        ! meters/sec to pascals/sec
             end if
          end do
       end do

    end if

 end subroutine cld_sediment_vel


!===============================================================================
  subroutine cld_sediment_tend (ncol, dtime  ,               &
       pint   , pmid   , pdel   , t      ,                   &
       cloud  , cldliq , cldice , pvliq  , pvice  ,          &
       liqtend, icetend, wvtend , htend  , sfliq  , sfice   )

!----------------------------------------------------------------------
!     Apply Cloud Particle Gravitational Sedimentation to Condensate
!----------------------------------------------------------------------


! Arguments
    integer,  intent(in)  :: ncol                      ! number of colums to process

    real(r8), intent(in)  :: dtime                     ! time step
    real(r8), intent(in)  :: pint  (pcols,pverp)       ! interfaces pressure (Pa)
    real(r8), intent(in)  :: pmid  (pcols,pver)        ! midpoint pressures (Pa)
    real(r8), intent(in)  :: pdel  (pcols,pver)        ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: cloud (pcols,pver)        ! cloud fraction (fraction)
    real(r8), intent(in)  :: t     (pcols,pver)        ! temperature (K)
    real(r8), intent(in)  :: cldliq(pcols,pver)        ! cloud liquid water (kg/kg)
    real(r8), intent(in)  :: cldice(pcols,pver)        ! cloud ice water    (kg/kg)
    real(r8), intent(in)  :: pvliq (pcols,pverp)       ! vertical velocity of liquid drops  (Pa/s)
    real(r8), intent(in)  :: pvice (pcols,pverp)       ! vertical velocity of ice particles (Pa/s)
! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

    real(r8), intent(out) :: liqtend(pcols,pver)       ! liquid condensate tend
    real(r8), intent(out) :: icetend(pcols,pver)       ! ice condensate tend
    real(r8), intent(out) :: wvtend (pcols,pver)       ! water vapor tend
    real(r8), intent(out) :: htend  (pcols,pver)       ! heating rate
    real(r8), intent(out) :: sfliq  (pcols)            ! surface flux of liquid (rain, kg/m/s)
    real(r8), intent(out) :: sfice  (pcols)            ! surface flux of ice    (snow, kg/m/s)

! Local variables
    real(r8) :: fxliq(pcols,pverp)                     ! fluxes at the interfaces, liquid (positive = down)
    real(r8) :: fxice(pcols,pverp)                     ! fluxes at the interfaces, ice    (positive = down)
    real(r8) :: cldab(pcols)                           ! cloud in layer above
    real(r8) :: evapliq                                ! evaporation of cloud liquid into environment
    real(r8) :: evapice                                ! evaporation of cloud ice into environment
    real(r8) :: cldovrl                                ! cloud overlap factor

    integer :: i,k
!----------------------------------------------------------------------

! initialize variables
    fxliq  (:ncol,:) = 0._r8 ! flux at interfaces (liquid)
    fxice  (:ncol,:) = 0._r8 ! flux at interfaces (ice)
    liqtend(:ncol,:) = 0._r8 ! condensate tend (liquid)
    icetend(:ncol,:) = 0._r8 ! condensate tend (ice)
    wvtend(:ncol,:)  = 0._r8 ! environmental moistening
    htend(:ncol,:)   = 0._r8 ! evaporative cooling
    sfliq(:ncol)     = 0._r8 ! condensate sedimentation flux out bot of column (liquid)
    sfice(:ncol)     = 0._r8 ! condensate sedimentation flux out bot of column (ice)

! fluxes at interior points
    call getflx(ncol, pint, cldliq, pvliq, dtime, fxliq)
    call getflx(ncol, pint, cldice, pvice, dtime, fxice)

! calculate fluxes at boundaries
    do i = 1,ncol
       fxliq(i,1) = 0._r8
       fxice(i,1) = 0._r8
! surface flux by upstream scheme
       fxliq(i,pverp) = cldliq(i,pver) * pvliq(i,pverp) * dtime
       fxice(i,pverp) = cldice(i,pver) * pvice(i,pverp) * dtime
    end do

! filter out any negative fluxes from the getflx routine
! (typical fluxes are of order > 1e-3 when clouds are present)
    do k = 2,pver
       fxliq(:ncol,k) = max(0._r8, fxliq(:ncol,k))
       fxice(:ncol,k) = max(0._r8, fxice(:ncol,k))
    end do

! Limit the flux out of the bottom of each cell to the water content in each phase.
! Apply mxsedfac to prevent generating very small negative cloud water/ice
! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
! ***Should we include the flux in the top, to allow for thin surface layers?
! ***Requires simple treatment of cloud overlap, already included below.
    do k = 1,pver
       do i = 1,ncol
          fxliq(i,k+1) = min( fxliq(i,k+1), mxsedfac * cldliq(i,k) * pdel(i,k) )
          fxice(i,k+1) = min( fxice(i,k+1), mxsedfac * cldice(i,k) * pdel(i,k) )
!!$        fxliq(i,k+1) = min( fxliq(i,k+1), cldliq(i,k) * pdel(i,k) + fxliq(i,k))
!!$        fxice(i,k+1) = min( fxice(i,k+1), cldice(i,k) * pdel(i,k) + fxice(i,k))
!!$        fxliq(i,k+1) = min( fxliq(i,k+1), cloud(i,k) * cldliq(i,k) * pdel(i,k) )
!!$        fxice(i,k+1) = min( fxice(i,k+1), cloud(i,k) * cldice(i,k) * pdel(i,k) )
       end do
    end do

! Now calculate the tendencies assuming that condensate evaporates when
! it falls into environment, and does not when it falls into cloud.
! All flux out of the layer comes from the cloudy part.
! Assume maximum overlap for stratiform clouds
!  if cloud above < cloud,  all water falls into cloud below
!  if cloud above > cloud,  water split between cloud  and environment
    do k = 1,pver
       cldab(:ncol) = 0._r8
       do i = 1,ncol
! cloud overlap cloud factor
          cldovrl  = min( cloud(i,k) / (cldab(i)+.0001_r8), 1._r8 )
          cldab(i) = cloud(i,k)
! evaporation into environment cause moistening and cooling
          evapliq = fxliq(i,k) * (1._r8-cldovrl) / (dtime * pdel(i,k))  ! into env (kg/kg/s)
          evapice = fxice(i,k) * (1._r8-cldovrl) / (dtime * pdel(i,k))  ! into env (kg/kg/s)
          wvtend(i,k) = evapliq + evapice                          ! evaporation into environment (kg/kg/s)
          htend (i,k) = -latvap*evapliq -(latvap+latice)*evapice   ! evaporation (W/kg)
! net flux into cloud changes cloud liquid/ice (all flux is out of cloud)
          liqtend(i,k)  = (fxliq(i,k)*cldovrl - fxliq(i,k+1)) / (dtime * pdel(i,k))
          icetend(i,k)  = (fxice(i,k)*cldovrl - fxice(i,k+1)) / (dtime * pdel(i,k))
       end do
    end do

! convert flux out the bottom to mass units Pa -> kg/m2/s
    sfliq(:ncol) = fxliq(:ncol,pverp) / (dtime*gravit)
    sfice(:ncol) = fxice(:ncol,pverp) / (dtime*gravit)

    return
  end subroutine cld_sediment_tend

!===============================================================================
  subroutine getflx(ncol, xw, phi, vel, deltat, flux)

!.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
!....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
!....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
!.........phi1......phi2.......phi3.....phi4.......phi5.......



    integer, intent(in) :: ncol                      ! number of colums to process

    integer i
    integer k

    real (r8), intent(in) :: vel(pcols,pverp)
    real (r8) flux(pcols,pverp)
    real (r8) xw(pcols,pverp)
    real (r8) psi(pcols,pverp)
    real (r8), intent(in) :: phi(pcols,pverp-1)
    real (r8) fdot(pcols,pverp)
    real (r8) xx(pcols)
    real (r8) fxdot(pcols)
    real (r8) fxdd(pcols)

    real (r8) psistar(pcols)
    real (r8) deltat

    real (r8) xxk(pcols,pver)

    do i = 1,ncol
!        integral of phi
       psi(i,1) = 0._r8
!        fluxes at boundaries
       flux(i,1) = 0._r8
       flux(i,pverp) = 0._r8
    end do

!     integral function
    do k = 2,pverp
       do i = 1,ncol
          psi(i,k) = phi(i,k-1)*(xw(i,k)-xw(i,k-1)) + psi(i,k-1)
       end do
    end do


!     calculate the derivatives for the interpolating polynomial
    call cfdotmc_pro (ncol, xw, psi, fdot)

!  NEW WAY
!     calculate fluxes at interior pts
    do k = 2,pver
       do i = 1,ncol
          xxk(i,k) = xw(i,k)-vel(i,k)*deltat
       end do
    end do
    do k = 2,pver
       call cfint2(ncol, xw, psi, fdot, xxk(1,k), fxdot, fxdd, psistar)
       do i = 1,ncol
          flux(i,k) = (psi(i,k)-psistar(i))
       end do
    end do


    return
  end subroutine getflx



!##############################################################################

  subroutine cfint2 (ncol, x, f, fdot, xin, fxdot, fxdd, psistar)



! input
    integer ncol                      ! number of colums to process

    real (r8) x(pcols, pverp)
    real (r8) f(pcols, pverp)
    real (r8) fdot(pcols, pverp)
    real (r8) xin(pcols)

! output
    real (r8) fxdot(pcols)
    real (r8) fxdd(pcols)
    real (r8) psistar(pcols)

    integer i
    integer k
    integer intz(pcols)
    real (r8) dx
    real (r8) s
    real (r8) c2
    real (r8) c3
    real (r8) xx
    real (r8) xinf
    real (r8) psi1, psi2, psi3, psim
    real (r8) cfint
    real (r8) cfnew
    real (r8) xins(pcols)

!     the minmod function 
    real (r8) a, b, c
    real (r8) minmod
    real (r8) medan
    logical found_error

    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do i = 1,ncol
       xins(i) = medan(x(i,1), xin(i), x(i,pverp))
       intz(i) = 0
    end do

! first find the interval 
    do k =  1,pverp-1
       do i = 1,ncol
          if ((xins(i)-x(i,k))*(x(i,k+1)-xins(i)).ge.0) then
             intz(i) = k
          endif
       end do
    end do

    found_error=.false.
    do i = 1,ncol
       if (intz(i).eq.0._r8) found_error=.true.
    end do
    if(found_error) then
       do i = 1,ncol
          if (intz(i).eq.0._r8) then
             write(iulog,*) ' interval was not found for col i ', i
             call endrun('CFINT2')
          endif
       end do
    endif

! now interpolate
    do i = 1,ncol
       k = intz(i)
       dx = (x(i,k+1)-x(i,k))
       s = (f(i,k+1)-f(i,k))/dx
       c2 = (3*s-2*fdot(i,k)-fdot(i,k+1))/dx
       c3 = (fdot(i,k)+fdot(i,k+1)-2*s)/dx**2
       xx = (xins(i)-x(i,k))
       fxdot(i) =  (3*c3*xx + 2*c2)*xx + fdot(i,k)
       fxdd(i) = 6*c3*xx + 2*c2
       cfint = ((c3*xx + c2)*xx + fdot(i,k))*xx + f(i,k)

!        limit the interpolant
       psi1 = f(i,k)+(f(i,k+1)-f(i,k))*xx/dx
       if (k.eq.1) then
          psi2 = f(i,1)
       else
          psi2 = f(i,k) + (f(i,k)-f(i,k-1))*xx/(x(i,k)-x(i,k-1))
       endif
       if (k+1.eq.pverp) then
          psi3 = f(i,pverp)
       else
          psi3 = f(i,k+1) - (f(i,k+2)-f(i,k+1))*(dx-xx)/(x(i,k+2)-x(i,k+1))
       endif
       psim = medan(psi1, psi2, psi3)
       cfnew = medan(cfint, psi1, psim)
       if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36_r8)  .gt..03_r8) then
!     CHANGE THIS BACK LATER!!!
!     $        .gt..1) then


!     UNCOMMENT THIS LATER!!!
!            write(iulog,*) ' cfint2 limiting important ', cfint, cfnew


       endif
       psistar(i) = cfnew
    end do

    return
  end subroutine cfint2



!##############################################################################

  subroutine cfdotmc_pro (ncol, x, f, fdot)

!     prototype version; eventually replace with final SPITFIRE scheme

!     calculate the derivative for the interpolating polynomial
!     multi column version



! input
    integer ncol                      ! number of colums to process

    real (r8) x(pcols, pverp)
    real (r8) f(pcols, pverp)
! output
    real (r8) fdot(pcols, pverp)          ! derivative at nodes

! assumed variable distribution
!     x1.......x2.......x3.......x4.......x5.......x6     1,pverp points
!     f1.......f2.......f3.......f4.......f5.......f6     1,pverp points
!     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
!     .........d2.......d3.......d4.......d5.........     2,pver points
!     .........s2.......s3.......s4.......s5.........     2,pver points
!     .............dh2......dh3......dh4.............     2,pver-1 points
!     .............eh2......eh3......eh4.............     2,pver-1 points
!     ..................e3.......e4..................     3,pver-1 points
!     .................ppl3......ppl4................     3,pver-1 points
!     .................ppr3......ppr4................     3,pver-1 points
!     .................t3........t4..................     3,pver-1 points
!     ................fdot3.....fdot4................     3,pver-1 points


! work variables


    integer i
    integer k

    real (r8) a                    ! work var
    real (r8) b                    ! work var
    real (r8) c                    ! work var
    real (r8) s(pcols,pverp)             ! first divided differences at nodes
    real (r8) sh(pcols,pverp)            ! first divided differences between nodes
    real (r8) d(pcols,pverp)             ! second divided differences at nodes
    real (r8) dh(pcols,pverp)            ! second divided differences between nodes
    real (r8) e(pcols,pverp)             ! third divided differences at nodes
    real (r8) eh(pcols,pverp)            ! third divided differences between nodes
    real (r8) pp                   ! p prime
    real (r8) ppl(pcols,pverp)           ! p prime on left
    real (r8) ppr(pcols,pverp)           ! p prime on right
    real (r8) qpl
    real (r8) qpr
    real (r8) ttt
    real (r8) t
    real (r8) tmin
    real (r8) tmax
    real (r8) delxh(pcols,pverp)


!     the minmod function 
    real (r8) minmod
    real (r8) medan
    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do k = 1,pver


!        first divided differences between nodes
       do i = 1, ncol
          delxh(i,k) = (x(i,k+1)-x(i,k))
          sh(i,k) = (f(i,k+1)-f(i,k))/delxh(i,k)
       end do

!        first and second divided differences at nodes
       if (k.ge.2) then
          do i = 1,ncol
             d(i,k) = (sh(i,k)-sh(i,k-1))/(x(i,k+1)-x(i,k-1))
             s(i,k) = minmod(sh(i,k),sh(i,k-1))
          end do
       endif
    end do

!     second and third divided diffs between nodes
    do k = 2,pver-1
       do i = 1, ncol
          eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
          dh(i,k) = minmod(d(i,k),d(i,k+1))
       end do
    end do

!     treat the boundaries
    do i = 1,ncol
       e(i,2) = eh(i,2)
       e(i,pver) = eh(i,pver-1)
!        outside level
       fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1)  &
            - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
       fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
       fdot(i,pverp) = sh(i,pver) + d(i,pver)*delxh(i,pver)  &
            + eh(i,pver-1)*delxh(i,pver)*(x(i,pverp)-x(i,pver-1))
       fdot(i,pverp) = minmod(fdot(i,pverp),3*sh(i,pver))
!        one in from boundary
       fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1) - eh(i,2)*delxh(i,1)*delxh(i,2)
       fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
       fdot(i,pver) = sh(i,pver) - d(i,pver)*delxh(i,pver)   &
            - eh(i,pver-1)*delxh(i,pver)*delxh(i,pver-1)
       fdot(i,pver) = minmod(fdot(i,pver),3*s(i,pver))
    end do


    do k = 3,pver-1
       do i = 1,ncol
          e(i,k) = minmod(eh(i,k),eh(i,k-1))
       end do
    end do



    do k = 3,pver-1

       do i = 1,ncol

!           p prime at k-0.5
          ppl(i,k)=sh(i,k-1) + dh(i,k-1)*delxh(i,k-1)  
!           p prime at k+0.5
          ppr(i,k)=sh(i,k)   - dh(i,k)  *delxh(i,k)

          t = minmod(ppl(i,k),ppr(i,k))

!           derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
          pp = sh(i,k-1) + d(i,k)*delxh(i,k-1) 

!           quartic estimate of fdot
          fdot(i,k) = pp                            &
               - delxh(i,k-1)*delxh(i,k)            &
               *(  eh(i,k-1)*(x(i,k+2)-x(i,k  ))    &
               + eh(i,k  )*(x(i,k  )-x(i,k-2))      &
               )/(x(i,k+2)-x(i,k-2))

!           now limit it
          qpl = sh(i,k-1)       &
               + delxh(i,k-1)*minmod(d(i,k-1)+e(i,k-1)*(x(i,k)-x(i,k-2)), &
               d(i,k)  -e(i,k)*delxh(i,k))
          qpr = sh(i,k)         &
               + delxh(i,k  )*minmod(d(i,k)  +e(i,k)*delxh(i,k-1),        &
               d(i,k+1)+e(i,k+1)*(x(i,k)-x(i,k+2)))

          fdot(i,k) = medan(fdot(i,k), qpl, qpr)

          ttt = minmod(qpl, qpr)
          tmin = min(0._r8,3*s(i,k),1.5_r8*t,ttt)
          tmax = max(0._r8,3*s(i,k),1.5_r8*t,ttt)

          fdot(i,k) = fdot(i,k) + minmod(tmin-fdot(i,k), tmax-fdot(i,k))

       end do

    end do

    return
  end subroutine cfdotmc_pro
end module pkg_cld_sediment
