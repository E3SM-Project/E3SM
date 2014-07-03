!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_bwater.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!TODO - Much of this module consists of Jesse's water-routing code,
!        which is serial only.
!       Add a subroutine that is supported for parallel code?

module glide_bwater

   use glimmer_global, only: rk, sp
   use glide_types

   implicit none

contains

  subroutine calcbwat(model, which, bmlt, bwat, bwatflx, thck, topg, btem, floater, wphi)

    use parallel
    use glimmer_global, only : dp 
    use glimmer_paramets, only : thk0
    use glide_grid_operators, only: stagvarb

    implicit none

    type(glide_global_type),intent(inout) :: model
    integer, intent(in) :: which
    real(dp), dimension(:,:), intent(inout) :: bwat, wphi, bwatflx
    real(dp), dimension(:,:), intent(in) :: bmlt, thck, topg, btem
    logical, dimension(:,:), intent(in) :: floater

    real(dp), dimension(2), parameter :: &
         blim = (/ 0.00001 / thk0, 0.001 / thk0 /)

    real(dp) :: dwphidew, dwphidns, dwphi, pmpt, bave

    integer :: t_wat,ns,ew

    real(dp),  dimension(model%general%ewn,model%general%nsn) :: N, flux, lakes
    real(dp) :: c_effective_pressure,c_flux_to_depth,p_flux_to_depth,q_flux_to_depth

    real(dp), parameter :: const_bwat = 10.d0   ! constant value for basal water depth (m)

! TODO: move these declarations into a parameters derived type?
    c_effective_pressure = 0.0d0       ! For now estimated with c/w
    c_flux_to_depth = 1./(1.8d-3*12.0d0)  ! 
    p_flux_to_depth = 2.0d0            ! exponent on the depth
    q_flux_to_depth = 1.0d0            ! exponent on the potential gradient

    select case (which)

    ! which = BWATER_NONE Nothing, basal water depth = 0.
    ! which = BWATER_LOCAL Completely local, bwat_new = c1 * melt_rate + c2 * bwat_old
    ! which = BWATER_FLUX Flux based calculation
    ! which = BWATER_BASAL_PROC, till water content in the basal processes module

    case(BWATER_LOCAL)

       ! model%tempwk%c(1) =  model%tempwk%dt_wat
       !              c(2) =  1.0d0 - 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim
       !              c(3) =  1.0d0 + 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim

       do t_wat = 1, model%tempwk%nwat

          !LOOP - For glissade, loop should be over locally owned cells (ilo:ihi,jo:jhi).

          do ns = 1,model%general%nsn
             do ew = 1,model%general%ewn

                if (model%numerics%thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   bwat(ew,ns) = (model%tempwk%c(1) * bmlt(ew,ns) + model%tempwk%c(2) * bwat(ew,ns)) / &
                        model%tempwk%c(3)
                   if (bwat(ew,ns) < blim(1)) then
                      bwat(ew,ns) = 0.0d0
                   end if
                else
                   bwat(ew,ns) = 0.0d0
                end if

             end do
          end do
       end do

       model%tempwk%smth = 0.
       do ns = 2,model%general%nsn-1
          do ew = 2,model%general%ewn-1
             call smooth_bwat(ew-1,ew,ew+1,ns-1,ns,ns+1)
          end do
       end do

       ! apply periodic BC
       if (model%options%periodic_ew) then
          do ns = 2,model%general%nsn-1
             call smooth_bwat(model%general%ewn-1,1,2,ns-1,ns,ns+1)
             call smooth_bwat(model%general%ewn-1,model%general%ewn,2,ns-1,ns,ns+1)
          end do
       end if

       bwat(1:model%general%ewn,1:model%general%nsn) = &
                       model%tempwk%smth(1:model%general%ewn,1:model%general%nsn)

    ! Case added by Jesse Johnson 11/15/08
    ! Steady state routing of basal water using flux calculation

    case(BWATER_FLUX)

       !TODO - Test this option

       call not_parallel(__FILE__,__LINE__)

       call effective_pressure(bwat,c_effective_pressure,N)
       call pressure_wphi(thck,topg,N,wphi,model%numerics%thklim,floater)
       call route_basal_water(wphi,bmlt,model%numerics%dew,model%numerics%dns,bwatflx,lakes)
       call flux_to_depth(bwatflx,wphi,c_flux_to_depth,p_flux_to_depth,q_flux_to_depth,model%numerics%dew,model%numerics%dns,bwat)

    case(BWATER_CONST)

       ! Use a constant thickness of water, to force Tpmp.
       bwat(:,:) = const_bwat / thk0

!!    case(BWATER_BASAL_PROC)   ! not currently supported

       ! Normalized basal water 

!!     bwat = model%basalproc%Hwater / thk0

    case default   ! includes BWATER_NONE

       bwat(:,:) = 0.0d0

    end select

    ! now also calculate basal water in velocity (staggered) coord system
    call stagvarb(model%temper%bwat, &
                  model%temper%stagbwat ,&
                  model%general%ewn, &
                  model%general%nsn)

  contains

    ! Internal subroutine for smoothing
    subroutine smooth_bwat(ewm,ew,ewp,nsm,ns,nsp)
      ! smoothing basal water distrib
      implicit none
      integer, intent(in) :: ewm,ew,ewp,nsm,ns,nsp

      if (bwat(ew,ns) > blim(2)) then
         model%tempwk%smth(ew,ns) = bwat(ew,ns) + model%paramets%bwat_smooth * &
              (bwat(ewm,ns) + bwat(ewp,ns) + bwat(ew,nsm) + bwat(ew,nsp) - 4.0d0 * bwat(ew,ns))
      else 
         model%tempwk%smth(ew,ns) = bwat(ew,ns)
      end if   
    end subroutine smooth_bwat

  end subroutine calcbwat
  
  subroutine find_dt_wat(dttem,estimate,dt_wat,nwat)
    
    implicit none
    
    real(dp), intent(out) :: dt_wat
    integer, intent(out) :: nwat
    real(dp), intent(in) :: dttem, estimate
    
    nwat = int(dttem/estimate) + 1
    dt_wat = dttem / nwat

  end subroutine find_dt_wat

  ! Note: This routing is supported in serial code only.

  subroutine route_basal_water(wphi,melt,dx,dy,flux,lakes)
    !*FD Routes water from melt field to its destination, recording flux
    !*FD of water along the route. Water flow direction is determined according
    !*FD to the gradient of a wphi elevation field. For the algorithm to 
    !*FD function properly depressions in the wphi surface must be filled.
    !*FD this results in the lakes field, which is the difference between the
    !*FD filled surface and the original wphi.
    !*FD The method used is by Quinn et. al. (1991).
    !*FD
    !*FD 12/9/05 Jesse Johnson based on code from the glimmer_routing file
    !*FD by Ian Rutt.

    implicit none

    real(rk),dimension(:,:),intent(in)  :: wphi    !*FD Input potential surface
    real(rk),dimension(:,:),intent(in)  :: melt    !*FD Input melting field
    real(rk),               intent(in)  :: dx      !*FD Input $x$ grid-spacing
    real(rk),               intent(in)  :: dy      !*FD Input $y$ grid-spacing
    real(rk),dimension(:,:),intent(out) :: flux    !*FD Output flux field
    real(rk),dimension(:,:),intent(out) :: lakes   !*FD Output lakes field

    ! Internal variables --------------------------------------

    integer :: nx,ny,k,nn,cx,cy,px,py,x,y
    integer, dimension(:,:),allocatable :: mask    !*FD Masked points
    integer, dimension(:,:),allocatable :: sorted
    real(rk),dimension(:,:),allocatable :: flats,potcopy
    real(rk),dimension(-1:1,-1:1) :: slopes
    real(rk),dimension(-1:1,-1:1) :: dists
    logical :: flag

    ! Set up grid dimensions ----------------------------------

    nx=size(wphi,1) ; ny=size(wphi,2)
    nn=nx*ny

    ! Change these distances for slope determination

    dists(-1,:)=(/sqrt(dx**2+dy**2),dy,sqrt(dx**2+dy**2)/)
    dists(0,:)=(/dx,0d0,dx/)
    dists(1,:)=dists(-1,:)

    ! Allocate internal arrays and copy data ------------------

    allocate(sorted(nn,2),flats(nx,ny),potcopy(nx,ny),mask(nx,ny))
    potcopy=wphi
    mask=1

    ! Fill holes in data, and sort heights --------------------

    call fillholes(potcopy,flats,mask)
    call heights_sort(potcopy,sorted)

    lakes=potcopy-wphi

    ! Initialise flux with melt, which will then be --------
    ! redistributed. Multiply by area, so volumes are found.---

    flux=melt * dx * dy

    ! Begin loop over points, highest first -------------------

    do k=nn,1,-1
    
      ! Get location of current point -------------------------

      x=sorted(k,1)
      y=sorted(k,2)

      ! Only propagate down slope positive values 
      if (melt(x,y) > 0) then

        ! Reset flags and slope arrays --------------------------

        flag=.true.
        slopes=0.0

        ! Loop over adjacent points, and calculate slopes -------

        do cx=-1,1,1
          do cy=-1,1,1
            ! If this is the centre point, ignore
            if (cx==0.and.cy==0) continue
            ! Otherwise do slope calculation 
            px=x+cx ; py=y+cy
            if (px > 0 .and. px<=nx .and. py > 0 .and. py <= ny) then
                ! Only allow flow to points that are melted or freezing.
                ! Testing relax this condition (Hell, Frank does).
                !if (potcopy(px,py)<potcopy(x,y) .and. melt(px,py)/=0.0) then
                if (potcopy(px,py)<potcopy(x,y)) then
                  slopes(cx,cy)=(potcopy(x,y)-potcopy(px,py))/dists(cx,cy)
                endif
            endif
          enddo
        enddo

        ! If there are places for the water to drain to, --------
        ! distribute it accordingly -----------------------------

        if (sum(slopes)/=0.0) then
          slopes=slopes/sum(slopes)
          do cx=-1,1
            do cy=-1,1
              px=x+cx ;py=y+cy
              if (slopes(cx,cy)/=0.0) then
                flux(px,py)=flux(px,py)+flux(x,y)*slopes(cx,cy)
              endif
            enddo
          enddo

        ! Note that sources are not zeroed in this case.---------

        endif

      ! End test for positive melt rate.-------------------------
      endif

      ! End of main loop ----------------------------------------

    enddo

    ! Tidy up -------------------------------------------------

    deallocate(sorted,flats)

  end subroutine route_basal_water

!==============================================================

  subroutine flux_to_depth(flux,wphi,c,p,q,dew,dns,bwat)

  !*FD Assuming that the flow is steady state, this function simply solves
  !*FD              flux = depth * velocity
  !*FD for the depth, assuming that the velocity is a function of depth,
  !*FD and pressure potential. This amounts to assuming a Weertman film,
  !*FD or Manning flow, both of which take the form of a constant times water
  !*FD depth to a power, times pressure wphi to a power.

    use glam_grid_operators, only: df_field_2d      ! Find grad_wphi
    use glimmer_physcon, only : scyr                ! Seconds per year

    real(dp),dimension(:,:),intent(in) :: flux      ! Basal water flux
    real(dp),dimension(:,:),intent(in) :: wphi      ! Pressure wphi
    real(dp)               ,intent(in) ::  c        ! Constant of proportionality
    real(dp)               ,intent(in) ::  p        ! Exponent of the water depth
    real(dp)               ,intent(in) ::  q        ! Exponent of the pressure pot.
    real(dp)               ,intent(in) ::  dew      ! Grid spacing, ew direction
    real(dp)               ,intent(in) ::  dns      ! Grid spacing, ns direction
    real(dp),dimension(:,:),intent(out)::  bwat     ! Water Depth

    ! Internal variables 
    real(rk),dimension(:,:),allocatable :: grad_wphi, dwphidx, dwphidy

    integer nx,ny,nn

    ! Set up grid dimensions ----------------------------------
    nx=size(flux,1) ; ny=size(flux,2)
    nn=nx*ny

    ! Allocate internal arrays and copy data ------------------
    allocate(dwphidx(nx,ny),dwphidy(nx,ny),grad_wphi(nx,ny))

    ! Compute the gradient of the potential field.
    call df_field_2d(wphi,dew,dns,dwphidx,dwphidy)

    grad_wphi = sqrt(dwphidx**2 + dwphidy**2)

    where (grad_wphi /= 0.d0) 
        bwat = ( flux / (c * scyr *  dns * grad_wphi ** q) ) ** (1./(p+1.))
    elsewhere
        bwat = 0.d0
    endwhere
       

  end subroutine flux_to_depth

!==============================================================

  subroutine effective_pressure(bwat,c,N)
    real(dp),dimension(:,:),intent(in) ::  bwat! Water depth
    real(dp)               ,intent(in) ::  c   ! Constant of proportionality
    real(dp),dimension(:,:),intent(out) :: N   ! Effective pressure

    where (bwat > 0.d0)
        N = c / bwat
    elsewhere
        N = 0.d0
    endwhere
  end subroutine effective_pressure

!==============================================================

  subroutine pressure_wphi(thck,topg,N,wphi,thicklim,floater)
  !*FD Compute the pressure wphi at the base of the ice sheet according to
  !*FD ice overburden plus bed height minus effective pressure.
  !*FD
  !*FD whpi/(rhow*g) = topg + bwat * rhoi / rhow * thick - N / (rhow * g)

    use glimmer_physcon, only : rhoi,rhow,grav
    implicit none
    real(dp),dimension(:,:),intent(in) :: thck      ! Thickness
    real(dp),dimension(:,:),intent(in) :: topg      ! Bed elevation
    real(dp),dimension(:,:),intent(in) :: N         ! Effective pressure
    logical,dimension(:,:),intent(in)  :: floater   ! Mask of floating ice
    real(dp),intent(in)                :: thicklim  ! Minimal ice thickness
    real(dp),dimension(:,:),intent(out) :: wphi     ! Pressure wphi


    where (thck > thicklim .and. .not. floater)
      wphi = thck + rhow/rhoi * topg - N / (rhow * grav)
    elsewhere
      wphi = max(topg *rhow/rhoi,0.0d0)
    end where

  end subroutine pressure_wphi
  
!==============================================================
! Internal subroutines
!==============================================================

  subroutine fillholes(phi,flats,mask)

    implicit none

    real(rk),dimension(:,:),intent(inout) :: phi
    real(rk),dimension(:,:),intent(inout) :: flats
    integer, dimension(:,:),intent(in)    :: mask

    ! Internal variables --------------------------------------

    real(rk),allocatable,dimension(:,:) :: old_phi
    integer, allocatable,dimension(:,:) :: pool

    real(rk) :: pvs(9), max_val
    real(rk), parameter :: null = 1e+20
    integer :: flag,nx,ny,i,j

    ! ---------------------------------------------------------

    nx=size(phi,1) ; ny=size(phi,2)

    allocate(pool(nx,ny),old_phi(nx,ny))

    flag = 1

    ! ---------------------------------------------------------

    do while (flag == 1)

       flag = 0

       old_phi = phi

       do i=2,nx-1
          do j=2,ny-1

             flats(i,j) = 0

             if (mask(i,j) == 1) then

                if (any(old_phi(i-1:i+1,j-1:j+1) < old_phi(i,j))) then
                   pool(i,j) = 0
                else
                   pool(i,j) = 1
                end if

                if (pool(i,j) == 1) then

                   flag = 1

                   pvs = (/ old_phi(i-1:i+1,j-1), old_phi(i-1:i+1,j+1), old_phi(i-1:i+1,j) /)

                   where (pvs == old_phi(i,j))
                      pvs = null
                   end where

                   max_val = minval(pvs)

                   if (max_val /= null) then
                      phi(i,j) = max_val
                   else
                      flag = 0
                      flats(i,j) = 1
                   end if

                end if

             end if
          end do
       end do

    end do

    deallocate(pool,old_phi)

  end subroutine fillholes

!==============================================================

  subroutine heights_sort(wphi,sorted)

    real(rk),dimension(:,:) :: wphi
    integer,dimension(:,:) :: sorted

    integer :: nx,ny,nn,i,j,k
    real(rk),dimension(:),allocatable :: vect
    integer,dimension(:),allocatable :: ind

    nx=size(wphi,1) ; ny=size(wphi,2)
    nn=size(sorted,1)

    allocate(vect(nn),ind(nn)) 

    if (nn/=nx*ny.or.size(sorted,2) /= 2) then
      print*,'Wrong dimensions'
      stop
    endif

    k=1

    do i=1,nx
      do j=1,ny
        vect(k)=wphi(i,j)
        k=k+1
      enddo
    enddo

    call indexx(vect,ind)

    do k=1,nn
      sorted(k,1)=floor(real(ind(k)-1)/real(ny))+1
      sorted(k,2)=mod(ind(k)-1,ny)+1
    enddo

    do k=1,nn
      vect(k)=wphi(sorted(k,1),sorted(k,2))
    enddo
    
  end subroutine heights_sort

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ! The following two subroutines perform an index-sort of an array. 
  ! They are a GPL-licenced replacement for the Numerical Recipes routine indexx. 
  ! They are not derived from any NR code, but are based on a quicksort routine by
  ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
  ! in C, and issued under the GNU General Public License. The conversion to 
  ! Fortran 90, and modification to do an index sort was done by Ian Rutt.
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine indexx(array,index)

    use glimmer_log

    !*FD Performs an index sort of \texttt{array} and returns the result in
    !*FD \texttt{index}. The order of elements in \texttt{array} is unchanged.
    !*FD
    !*FD This is a GPL-licenced replacement for the Numerical Recipes routine indexx. 
    !*FD It is not derived from any NR code, but are based on a quicksort routine by
    !*FD Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !*FD in C, and issued under the GNU General Public License. The conversion to 
    !*FD Fortran 90, and modification to do an index sort was done by Ian Rutt.

    real(rk),dimension(:) :: array !*FD Array to be indexed.
    integer, dimension(:) :: index !*FD Index of elements of \texttt{array}.
    integer :: i

    if (size(array) /= size(index)) then
      call write_log('ERROR: INDEXX size mismatch.',GM_FATAL,__FILE__,__LINE__)
    endif

    do i=1,size(index)
       index(i)=i
    enddo

    call q_sort_index(array,index,1,size(array))

  end subroutine indexx

!==============================================================

  recursive subroutine q_sort_index(numbers,index,left,right)

    !*FD This is the recursive subroutine actually used by \texttt{indexx}. 
    !*FD
    !*FD This is a GPL-licenced replacement for the Numerical Recipes routine indexx. 
    !*FD It is not derived from any NR code, but are based on a quicksort routine by
    !*FD Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !*FD in C, and issued under the GNU General Public License. The conversion to 
    !*FD Fortran 90, and modification to do an index sort was done by Ian Rutt.

    implicit none

    real(rk),dimension(:) :: numbers !*FD Numbers being sorted
    integer, dimension(:) :: index   !*FD Returned index
    integer :: left, right           !*FD Limit of sort region

    integer :: ll,rr
    integer :: pv_int,l_hold, r_hold,pivpos
    real(rk) :: pivot

    ll=left
    rr=right

    l_hold = ll
    r_hold = rr
    pivot = numbers(index(ll))
    pivpos=index(ll)

    do
       if (.not.(ll < rr)) exit

       do 
          if  (.not.((numbers(index(rr)) >= pivot) .and. (ll < rr))) exit
          rr=rr-1
       enddo

       if (ll /= rr) then
          index(ll) = index(rr)
          ll=ll+1
       endif

       do
          if (.not.((numbers(index(ll)) <= pivot) .and. (ll < rr))) exit
          ll=ll+1
       enddo

       if (ll /= rr) then
          index(rr) = index(ll)
          rr=rr-1
       endif
    enddo

    index(ll) = pivpos
    pv_int = ll
    ll = l_hold
    rr = r_hold
    if (ll < pv_int)  call q_sort_index(numbers, index,ll, pv_int-1)
    if (rr > pv_int)  call q_sort_index(numbers, index,pv_int+1, rr)

  end subroutine q_sort_index

end module glide_bwater
