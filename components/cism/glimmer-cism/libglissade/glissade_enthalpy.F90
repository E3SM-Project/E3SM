!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   glissade_enthalpy.F90 - part of the Community Ice Sheet Model (CISM)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This module computes temperature diffusion, strain heating, and local
!  melting and refreezing in each ice column using enthalpy as the
!  primary state variable.

#include "glide_mask.inc"

!TODO - The glissade_enthalpy module needs to be tested.
!       Also, we may want to reorganize a bit to reduce the amount of duplicate code.
module glissade_enthalpy

    use glimmer_global, only : dp
    use glide_types
    use glimmer_log 

    implicit none

    private
    public :: glissade_init_enthalpy, glissade_enthalpy_findvtri, &
              enth2temp, temp2enth, glissade_enthalpy_calcbmlt

    !WHL - debug                                                                                                                          
    integer :: itest, jtest, rtest

!--------------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------------

  subroutine glissade_init_enthalpy (model)

    ! initialization for the enthalpy scheme
    use glimmer_physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow, trpt
    use glimmer_paramets, only : thk0, tim0

    type(glide_global_type),intent(inout) :: model    ! ice model parameters

    !BDM - dups(1) and dups(2) are done in glissade_init_temp

    ! write out cons(1) which is coefficient before diffusion terms 
    ! thk0 is just a scaling parameter, not actual thickness at ew, ns
    model%tempwk%cons(1) = tim0 * model%numerics%dttem/ (2.0d0 * thk0**2)

  end subroutine glissade_init_enthalpy

!--------------------------------------------------------------------------------

  subroutine glissade_enthalpy_findvtri (model, ew,   ns,          &
                                         subd,  diag, supd, rhsd,  &
                                         float, alpha_enth)

    ! solve for tridiagonal entries of sparse matrix
    use glimmer_physcon,  only : rhoi, shci, lhci, rhow, grav, coni
    use glimmer_paramets, only : thk0, tim0

    !WHL - debug
    use parallel, only: this_rank

    ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
    ! whereas temperature/enthalpy is indexed from 0 to upn.
    ! The first row of the matrix is the equation for enthalpy(0,ew,ns),
    ! the last row is the equation for enthalpy(upn,ew,ns), and so on.

    !I/O variables
    type(glide_global_type), intent(inout) :: model
    integer, intent(in) :: ew, ns
    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd
    logical, intent(in) :: float
    real(dp), dimension(:), intent(out) :: alpha_enth  ! half-node diffusivity (m^2/s) for enthalpy
	                                               ! located halfway between temperature points
    ! local variables
    real(dp) :: dsigbot  ! bottom layer thicknes in sigma coords.
    real(dp) :: alphai ! cold ice diffusivity
    real(dp) :: alpha0 ! temperate ice diffusivity
    real(dp) :: fact ! coefficient in tridiag
    integer  :: up, k
    real(dp), dimension(0:model%general%upn) :: enthalpy   ! specific enthalpy (J/m^3)
    real(dp), dimension(0:model%general%upn) :: pmptemp    ! pmptemp all nodes (deg C)
!!    real(dp), dimension(1:model%general%upn) :: alpha_half ! half-node diffusivity (m^2/s)
    real(dp), dimension(0:model%general%upn) :: enth_T     ! specific enthalpy (J/m^3)
    real(dp) :: denth    ! enthalpy difference between adjacent layers
    real(dp) :: denth_T  ! difference in temperature component of enthalpy between adjacent layers
    real(dp) :: alpha_fact  ! factor for averaging diffusivity, 0 <= fact <= 1

    logical, parameter :: &
         alpha_harmonic_avg  = .false.  ! if true, take harmonic average of alpha in adjacent layers
                                        ! if false, take arithmetic average

    !WHL - debug
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ! define diffusivities alpha_i and alpha_0
    alphai = coni / rhoi / shci
    alpha0 = alphai / 100.0d0
	
    !WHL - Moved temp2enth call to temperature driver
    ! Convert model%temper%temp and model%temper%waterfrac to enthalpy
    ! using temp and waterfrac from last timestep.  For interior and boundary nodes.
    ! BDM enthalpy will be size 0:upn
!    call temp2enth(enthalpy(0:model%general%upn),                       &
!                   model%temper%temp(0:model%general%upn,ew,ns),        &
!                   model%temper%waterfrac(1:model%general%upn-1,ew,ns), &
!                   model%geometry%thck(ew,ns),                          &
!                   model%numerics%stagsigma(1:model%general%upn-1))

    !WHL - Copy enthalpy to model derived type so that it's available on output

    ! Set local 1D enthalpy variable
    enthalpy(:) = model%temper%enthalpy(:,ew,ns)

    ! find pmptemp for this column (interior nodes and boundary)
    pmptemp(0) = 0.0d0
    call glissade_calcpmpt(pmptemp(1:model%general%upn-1), model%geometry%thck(ew,ns), &
                           model%numerics%stagsigma(1:model%general%upn-1))
							
    call glissade_calcpmpt_bed(pmptemp(model%general%upn), model%geometry%thck(ew,ns))

    !WHL - debug                                                                                                                                          
    if (ew==itest .and. ns==jtest) then
       print*, ' '
       print*, 'Starting enthalpy calc, i, j =', ew, ns
       print*, 'k, temp, wfrac, enthalpy/(rhoi*ci), pmpt:'
       k = 0
       print*, k, model%temper%temp(k,ew,ns), 0.d0, enthalpy(k)/(rhoi*shci)
       do k = 1, model%general%upn-1
          print*, k, model%temper%temp(k,ew,ns), model%temper%waterfrac(k,ew,ns), &
               enthalpy(k)/(rhoi*shci), pmptemp(k)
       enddo
       k = model%general%upn
       print*, k, model%temper%temp(k,ew,ns), 0.d0, enthalpy(k)/(rhoi*shci), pmptemp(k)
    endif

    !WHL - Commenting out the following and replacing it with a new way of computing alpha.
    !      The commented-out code can result in sudden large changes in alpha that
    !       lead to oscillations in the thickness, temperature and velocity fields.
    !      These oscillations have a period of ~1 yr or more, spatial scale of
    !       many grid cells, and amplitude of ~10 m in thickness, 1 deg in temperature,
    !       and 2 m/s in velocity.

    ! create a column vector of size (0:upn) of diffusivity based on 
    ! previous timestep's temp.  Boundary nodes need a value so half-node
    ! diffusivity can be calculated at interior nodes (1:upn-1)

!    do up = 0,model%general%upn
!       if (model%temper%temp(up,ew,ns) < pmptemp(up)) then
!          alpha(up) = alphai
!       else
!          alpha(up) = alpha0
!       endif
!    end do
    
    ! Find half-node diffusivity using harmonic average between nodes.
    ! The vector will be size (1:upn) - the first value is the half-node
    ! between nodes 0 and 1, the last value is the half-node between
    ! nodes upn-1 and upn. 
    
!    do up = 1,model%general%upn
!       alpha_enth(up) = 2.d0 / ((1.d0/alpha(up-1)) + (1.d0/alpha(up)))
!    end do
       
    !--------------------------------------------------------------------
    !WHL - Trying a different approach to the diffusivity at layer interfaces.
    ! Let d(enth)/dz = the gradient of enthalpy
    ! Can write 
    !    d(enth)/dz = d(enth_T)/dz + d(enth_w)/dz,
    ! where
    !    enth_T = (1-phi_w) * rhoi*ci*T
    !    enth_w =    phi_w  * rhow*(L + ci*Tpmp)
    !
    ! Now let f = d(enth_T)/z / d(enth)/dz
    !   (f -> 0 if f is computed to be negative)
    ! For cold ice, f = 1 and alpha = alphai
    ! For temperate ice, f ~ 0 and alpha = alpha0
    ! At the interface between cold and temperate ice,
    !  f ~ 0 if the temperate ice has large phi_w, but
    !  f ~ 1 if the temperate ice has close to zero phi_w.
    ! Two ways to average:
    ! (1) arithmetic average:  alpha = f*alphai + (1-f)*alpha0
    ! (2) harmonic average:    alpha = 1 / (f/alphai + (1-f)/alpha0).
    ! Both methods have the same asymptotic values at f = 0 or 1,
    !  but the arithmetic average gives greater diffusivity for
    !  intermediate values.
    !
    ! Still to be determined which is more accurate.
    ! The harmonic average allows large temperature gradients between the 
    !  bottom layer and the next layer up; the arithmetic average gives
    !  smoother gradients.
    !--------------------------------------------------------------------
    !
    ! At each temperature point, compute the temperature part of the enthalpy.
    ! enth_T = enth for cold ice, enth_T < enth for temperate ice

    do up = 0, model%general%upn
       enth_T(up) = (1.d0 - model%temper%waterfrac(up,ew,ns)) * rhoi*shci*model%temper%temp(up,ew,ns)
    enddo

!WHL - debug
    if (ew==itest .and. ns==jtest) then
       print*, ' '
       print*, 'k, denth_T/(rhoi*shci), denth/(rhoi*shci), alpha_fact, alpha_enth(up):'
    endif

    ! Compute factors relating the temperature gradient to the total enthalpy gradient.
    ! Use these factors to average the diffusivity between adjacent temperature points.
    do up = 1,model%general%upn
       denth   = enthalpy(up) - enthalpy(up-1)
       denth_T = enth_T(up) - enth_T(up-1)   ! = denth in cold ice, < denth in temperate ice
       if (abs(denth) > 1.d-20 * rhow*lhci) then
          alpha_fact = max(0.d0, denth_T/denth)
          alpha_fact = min(1.d0, alpha_fact)
       else
          alpha_fact = 0.d0
       endif

       if (alpha_harmonic_avg) then  ! take a harmonic average
                                     ! This gives slower cooling of temperate layers and allows
                                     !  large temperature gradients between cold and temperate layers
          alpha_enth(up) = 1.d0 / ((alpha_fact/alphai) + (1.d0-alpha_fact)/alpha0)
       else   ! take an arithmetic average
              ! This gives faster cooling of temperate layers and smaller gradients 
          alpha_enth(up) = alpha_fact*alphai + (1.d0-alpha_fact)*alpha0
       endif

!WHL - debug
       if (ew==itest .and. ns==jtest) then
          print*, up, denth_T/(rhoi*shci), denth/(rhoi*shci), alpha_fact, alpha_enth(up)
       endif

    end do

    !WHL - debug
    if (ew==itest .and. ns==jtest) then
!       print*, ' '
!       print*, 'alphai, alpha0 =', alphai, alpha0
!       print*, ' '
!       print*, 'k, alpha_enth(up), alpha(up-1), alpha(up), i, j =', ew, ns
!       do up = 1, model%general%upn
!          print*, up, alpha_enth(up), alpha(up-1), alpha(up)
!       enddo
    endif
	
    ! Compute subdiagonal, diagonal, and superdiagonal matrix elements
    
    ! upper boundary: set to surface air temperature
    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0
    rhsd(1) = dmin1(0.0d0,dble(model%climate%artm(ew,ns))) * rhoi * shci
  
    ! RJH - Multiplied fact by a factor of 2 to become EB coefficients
    ! cons(1) is (dt * tim0) / (2 * thk0^2)
    ! fact is (dt * tim0) / (2 * H^2 * thk0^2)
    fact = 2 * model%tempwk%cons(1) / model%geometry%thck(ew,ns)**2

    ! RJH - Altered rhsd to become fully implicit (backward Euler). 
    !       This included deleting subd and supd terms and dropping enthalpy(1:model%general%upn-1) coefficient

    ! ice interior. layers 1:upn-1  (matrix elements 2:upn)
    subd(2:model%general%upn) = -fact * alpha_enth(1:model%general%upn-1)        &
                                * model%tempwk%dups(1:model%general%upn-1,1)
    supd(2:model%general%upn) = -fact * alpha_enth(2:model%general%upn)          &
                                * model%tempwk%dups(1:model%general%upn-1,2)
    diag(2:model%general%upn) = 1.0d0 - subd(2:model%general%upn)                &
                                - supd(2:model%general%upn)
    rhsd(2:model%general%upn) = enthalpy(1:model%general%upn-1)                 &
                              + model%temper%dissip(1:model%general%upn-1,ew,ns) * rhoi * shci

    ! BDM I'm assuming that model%temper%dissip has units of phi/rhoi/shci.
    ! For an enthalpy calc, we want just phi, so model%temper%dissip * rhoi * shci
	
    ! basal boundary:
    ! for grounded ice, a heat flux is applied
    ! for floating ice, the basal temperature is held constant

    !NOTE: This lower BC is different from the one in standard glide_temp.
    !      If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
    !       That is, net heat flux at lower boundary must equal zero.
    !      If T(upn) >= Tpmp, then set T(upn) = Tpmp

    if (float) then
       supd(model%general%upn+1) = 0.0d0 
       subd(model%general%upn+1) = 0.0d0
       diag(model%general%upn+1) = 1.0d0
       rhsd(model%general%upn+1) = enthalpy(model%general%upn)
    
    else    ! grounded ice

       !WHL - debug
       if (ew==itest .and. ns==jtest) then
          up = model%general%upn-1
          print*, 'temp(upn-1), pmptemp(upn-1):', model%temper%temp(up,ew,ns), pmptemp(up)
          up = model%general%upn
          print*, 'temp(upn), pmptemp(upn):', model%temper%temp(up,ew,ns), pmptemp(up)
       endif

    !Positive-Thickness Basal Temperate Boundary Layer

    !WHL - Not sure whether this condition is right.  
    !      It implies that the enthalpy at the bed (upn) = enthalpy in layer (upn-1). 
       if (abs(model%temper%temp(model%general%upn-1,ew,ns) -       &
            pmptemp(model%general%upn-1)) < 0.001d0) then   
       
          subd(model%general%upn+1) = -1.0d0
          supd(model%general%upn+1) =  0.0d0 
          diag(model%general%upn+1) = 1.0d0 
          rhsd(model%general%upn+1) = 0.0d0

          !WHL - debug
          if (ew==itest .and. ns==jtest) then
             print*, 'basal BC: branch 1 (finite-thck BL)'
          endif

          !Zero-Thickness Basal Temperate Boundary Layer
       elseif (abs(model%temper%temp(model%general%upn,ew,ns) -       &
            pmptemp(model%general%upn)) < 0.001d0) then  ! melting
          
          ! hold basal temperature at pressure melting point
          supd(model%general%upn+1) = 0.0d0
          subd(model%general%upn+1) = 0.0d0
          diag(model%general%upn+1) = 1.0d0
          rhsd(model%general%upn+1) = pmptemp(model%general%upn) * rhoi * shci
          
          !WHL - debug
          if (ew==itest .and. ns==jtest) then
             print*, 'basal BC: branch 2 (zero-thck BL)'
          endif
          
       else  
          
          !WHL - debug
          if (ew==itest .and. ns==jtest) then
             print*, 'basal BC: branch 3 (cold ice)'
          endif
          
          ! frozen at bed
          ! maintain balance of heat sources and sinks
          ! (conductive flux, geothermal flux, and basal friction)
          ! Note: Heat fluxes are positive down, so slterm <= 0 and bheatflx <= 0.
          
          ! Note: The heat source due to basal sliding (bfricflx) is computed in subroutine calcbfric.
          ! Also note that bheatflx is generally <= 0, since defined as positive down.
          
          ! calculate dsigma for the bottom layer between the basal boundary and the temp. point above
          dsigbot = (1.0d0 - model%numerics%stagsigma(model%general%upn-1))                                                                  
          
          ! =====Backward Euler flux basal boundary condition=====
          ! MJH: If Crank-Nicolson is desired for the b.c., it is necessary to
          ! ensure that the i.c. temperature for the boundary satisfies the
          ! b.c. - otherwise oscillations will occur because the C-N b.c. only
          ! specifies the basal flux averaged over two consecutive time steps.
          subd(model%general%upn+1) = -1.0d0
          supd(model%general%upn+1) =  0.0d0 
          diag(model%general%upn+1) = 1.0d0 
          rhsd(model%general%upn+1) = (model%temper%bfricflx(ew,ns) -          &
                                       model%temper%bheatflx(ew,ns)) *          & 
                                       dsigbot * model%geometry%thck(ew,ns) *   &
                                       thk0 * rhoi * shci / coni
          ! BDM temp approach should work out to be dT/dsigma, so enthalpy approach
          ! should just need dT/dsigma * rhoi * shci for correct units

          ! =====Basal boundary using heat equation with specified flux====
          ! MJH: These coefficients are based on those used in the old temperature code 
          ! (eqns. 3.60-3.62 in the documentation).
          ! The implementation assumes the basal fluxes are the same at both time steps (lagged).
          ! The flux b.c. above was determined to be preferable, but this is left
          ! as an alternative.  It gives similar, but slightly different results.
          ! Because this formulation uses C-N time averaging, it results
          ! in a slight oscillation.
          !subd(model%general%upn+1) = -fact / dsigbot**2                                                                                     
          !supd(model%general%upn+1) =  0.0d0                                                                                                 
          !diag(model%general%upn+1) = 1.0d0 + fact / dsigbot**2       
          !model%tempwk%inittemp(model%general%upn,ew,ns) =    &   
          !       model%temper%temp(model%general%upn-1,ew,ns) * fact / dsigbot**2  &
          !       + model%temper%temp(model%general%upn,  ew,ns)  &
          !       * (1.0d0 - fact/dsigbot**2)   &   
          !       - fact *2.0d0 * & 
          !       model%geometry%thck(ew,ns) * thk0 / coni / dsigbot *  &
          !       (model%temper%bheatflx(ew,ns) & ! geothermal (H/k)*G
          !       - model%temper%bfricflx(ew,ns) )  ! sliding (H/k)*taub*ub.
          !rhsd(model%general%upn+1) = model%tempwk%inittemp(model%general%upn,ew,ns)
          
       endif   ! melting or frozen
       
    end if     ! floating or grounded

  end subroutine glissade_enthalpy_findvtri

!--------------------------------------------------------------------------------

  subroutine enth2temp (enthalpy, temp, waterfrac, thck, stagsigma)

    ! BDM convert from specific enthalpy to ice temperature and water content
    ! takes a vertical column of size enthalpy(dup-1,1,1) and converts to
    ! temp(1:dup-1,1,1) and waterfrac(1:dup-1,1,1)
    ! enthalpy(1:dup-1,1,1)
    ! temp(1:dup-1,1,1)
    ! waterfrac(1:dup-1,1,1)
    ! thck(1,1,1)
    ! stagsigma(1:dup-1,1,1)
	
    use glimmer_physcon, only : rhoi, shci, lhci, rhow

    ! I/O variables
    real(dp), dimension(0:), intent(inout)               :: enthalpy !enthalpy is (0:upn)
    real(dp), intent(in)                                 :: thck
    real(dp), dimension(:), intent(in)                   :: stagsigma !stagsigma is (1:upn-1)
    real(dp), dimension(0:size(enthalpy)-1), intent(out) :: temp !temp is (0:upn)
    real(dp), dimension(1:size(enthalpy)-2), intent(out) :: waterfrac !waterfrac is (1:upn-1)

    ! local variables
    real(dp), dimension(0:size(enthalpy)-1)              :: pmptemp ! (0:upn)
    real(dp), dimension(0:size(enthalpy)-1)              :: pmpenthalpy ! (0:upn)
    integer                                              :: upn !used for convenience
    integer                                              :: up
	
    upn = size(enthalpy)-1
	
    ! Find pmpenthalpy(0:upn)
    pmptemp(0) = 0.0d0
    call glissade_calcpmpt(pmptemp(1:upn-1), thck, stagsigma(1:upn-1))
    call glissade_calcpmpt_bed(pmptemp(upn), thck)
	
    pmpenthalpy = pmptemp * rhoi * shci
	
    !solve for temp and waterfrac
    if(enthalpy(0) >= pmpenthalpy(0)) then ! temperate ice
       temp(0) = pmptemp(0)                ! temperate ice
       !WHL - Resetting enthalpy so that it's consistent with the new temperature
       !      This is consistent with energy conservation because the top surface
       !       is infinitesimally thin.
       enthalpy(0) = pmpenthalpy(0)
    else
       temp(0) = enthalpy(0) / (rhoi*shci) ! temp is cold
    endif
	
    do up = 1, upn-1
       if(enthalpy(up) >= pmpenthalpy(up)) then ! temperate ice
          temp(up) = pmptemp(up)                ! temp = pmptemp
          waterfrac(up) = (enthalpy(up)-pmpenthalpy(up)) /                 &
                          ((rhow-rhoi) * shci * pmptemp(up) + rhow * lhci)
       else ! cold ice

          !WHL - debug
          if (waterfrac(up) > 0.d0) then
             print*, 'Zeroing out waterfrac: k, waterfrac =', up, waterfrac(up)
          endif

          temp(up) = enthalpy(up) / (rhoi*shci) ! temp is cold
          waterfrac(up) = 0.0d0                 ! waterfrac = 0
       endif
    end do
	
    if(enthalpy(upn) >= pmpenthalpy(upn)) then  ! temperate ice
       temp(upn) = pmptemp(upn)                 ! temp = pmptemp
    else
       temp(upn) = enthalpy(upn) / (rhoi*shci)  ! temp is cold
       !WHL - Resetting enthalpy so that it's consistent with the new temperature
       !      This is consistent with energy conservation because the top surface
       !       is infinitesimally thin.
       enthalpy(upn) = pmpenthalpy(upn)
    endif

  end subroutine enth2temp

!--------------------------------------------------------------------------------

  subroutine temp2enth (enthalpy, temp, waterfrac, thck, stagsigma)

    ! BDM convert from temperature and water content and converts to specific enthalpy
    ! takes a vertical column of size temp(0:dup) and converts to enthalpy(0:dup,1,1)
    ! waterfrac is only size(1:dup-1), so will assume no waterfrac at boundaries

    use glimmer_physcon, only : rhoi, shci, lhci, rhow

    ! I/O variables
    real(dp), dimension(0:), intent(out)        :: enthalpy !enthalpy is (0:upn)
    real(dp), intent(in)                        :: thck
    real(dp), dimension(:), intent(in)          :: stagsigma !stagsigma is (1:upn-1)
    real(dp), dimension(0:), intent(in)         :: temp !temp is (0:upn)
    real(dp), dimension(:), intent(in)          :: waterfrac !waterfrac is (1:upn-1)

    ! local variables
    real(dp), dimension(0:size(temp)-1)         :: pmptemp !(0:upn)
    real(dp), dimension(0:size(temp)-1)         :: pmpenthalpy !(0:upn)
    integer                                     :: up
    integer                                     :: upn !used for convenience
		
    upn = size(temp)-1
	
    ! Find pmpenthalpy(0:dup,1,1)
    pmptemp(0) = 0.0d0
    call glissade_calcpmpt(pmptemp(1:upn-1), thck, stagsigma(1:upn-1))
    call glissade_calcpmpt_bed(pmptemp(upn), thck)
    
    !WHL - This variable is not used below
    pmpenthalpy = rhoi * shci * pmptemp
	
    ! solve for enthalpy
    ! assume waterfrac = 0 at upper and lower ice surfaces
    enthalpy(0) = temp(0) * rhoi * shci
    do up = 1, upn-1
       enthalpy(up) = ((1 - waterfrac(up)) * rhoi * shci * temp(up))          &
                      + waterfrac(up) * rhow * ((shci * pmptemp(up)) + lhci)
    end do
    enthalpy(upn) = temp(upn) * rhoi * shci
	
  end subroutine temp2enth

!--------------------------------------------------------------------------------

  subroutine glissade_enthalpy_calcbmlt(model,                  &
                                        temp,        waterfrac, &
                                        stagsigma,   thck,      &
                                        bmlt_ground, floater)

    ! Compute the amount of basal melting.
    ! The basal melting computed here is applied to the ice thickness
    !  by glissade_transport_driver, conserving mass and energy.
    !
    ! This is done with the enthalpy formulation by taking any internal
    ! water content above 1% and draining it to the bed.
    ! 
    ! Note: Since this module is deprecated, the calculation has not been updated
    !       to include bmlt_float.

    use glimmer_physcon, only: shci, rhoi, lhci
    use glimmer_paramets, only : thk0, tim0

    type(glide_global_type) :: model

    real(dp), dimension(0:,:,:), intent(inout) :: temp
    real(dp), dimension(1:,:,:), intent(inout) :: waterfrac
    real(dp), dimension(0:),     intent(in) :: stagsigma
    real(dp), dimension(:,:),    intent(in) :: thck
    real(dp), dimension(:,:),    intent(out):: bmlt_ground  ! scaled melt rate (m/s * tim0/thk0)
                                                            ! > 0 for melting, < 0 for freeze-on
    logical,  dimension(:,:),    intent(in) :: floater

    real(dp), dimension(size(stagsigma))    :: pmptemp   ! pressure melting point temperature
    real(dp) :: bflx    ! heat flux available for basal melting (W/m^2)
    real(dp) :: hmlt    ! depth of internal melting (m)
    real(dp) :: internal_melt_rate   ! rate of internal melt sent to the bed (m/s)
    integer :: up, ew, ns

    bmlt_ground(:,:) = 0.0d0

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1

          if (thck(ew,ns) > model%numerics%thklim_temp .and. .not. floater(ew,ns)) then

             ! Basal friction term is computed above in subroutine glissade_calcbfric

             ! Compute basal melting
             ! Note: bmlt > 0 for melting, < 0 for freeze-on
             !       bfricflx >= 0 by definition
             !       bheatflx is positive down, so usually bheatflx < 0 (with negative values contributing to melt)
             !       lcondflx is positive down, so lcondflx < 0 for heat is flowing from the bed toward the surface

             !TODO - This equation allows for freeze-on (bmlt < 0) if the conductive term 
             !       (lcondflx, positive down) is carrying enough heat away from the boundary.  
             !       But freeze-on requires a local water supply, bwat > 0.
             !       What should we do if bwat = 0?

             bflx = model%temper%bfricflx(ew,ns) + model%temper%lcondflx(ew,ns) - model%temper%bheatflx(ew,ns)
             bmlt_ground(ew,ns) = bflx * model%tempwk%f(2)   ! f(2) = tim0 / (thk0 * lhci * rhoi)

             ! Add internal melting associated with waterfrac > waterfrac_max (1%)
             ! Note: glissade_calcpmpt does not compute pmpt at the top surface or the bed.

             call glissade_calcpmpt(pmptemp(:), thck(ew,ns),   &
                                    stagsigma(:) )

             !WHL - Any correction for rhoi/rhow here?
             do up = 1, model%general%upn-1
                if (waterfrac(up,ew,ns) > 0.01d0) then
                   hmlt = (waterfrac(up,ew,ns) - 0.01d0) * (model%geometry%thck(ew,ns) * thk0)  &
                        * (model%numerics%sigma(up+1) - model%numerics%sigma(up))     ! m
                   internal_melt_rate = hmlt / (model%numerics%dttem * tim0)          ! m/s
                   bmlt_ground(ew,ns) = bmlt_ground(ew,ns) + internal_melt_rate * tim0/thk0
                   waterfrac(up,ew,ns) = 0.01d0
                endif
             enddo

             ! Reset basal temp to pmptemp, if necessary
             !WHL - Is this necessary for enthalpy code?

             up = model%general%upn
             call glissade_calcpmpt_bed(pmptemp(up), thck(ew,ns))
             temp(up,ew,ns) = min (temp(up,ew,ns), pmptemp(up))

             ! If freeze-on was computed above (bmlt < 0) and Tbed = Tpmp but no basal water is present, then set T(upn) < Tpmp.
             ! Note: In subroutine findvtri, we solve for Tbed (instead of holding it at Tpmp) when Tbed < 0.001.
             !       With an offset here of 0.01, we will solve for T_bed at the next timestep.
             ! Note: Energy is not exactly conserved here.

             up = model%general%upn  ! basal level
             if (bmlt_ground(ew,ns) < 0.d0 .and. model%temper%bwat(ew,ns)==0.d0 .and. temp(up,ew,ns) >= pmptemp(up)) then
                temp(up,ew,ns) = pmptemp(up) - 0.01d0
             endif

          endif   ! thk > thklim_temp

       enddo
    enddo

  end subroutine glissade_enthalpy_calcbmlt

!-----------------------------------------------------------------------------------

  !TODO - Remove or inline these subroutines?  They are copies of subroutines in glissade_temp.F90.
  subroutine glissade_calcpmpt(pmptemp, thck, stagsigma)

    ! Compute the pressure melting point temperature in the column
    ! (but not at the surface or bed).
    ! Note: pmptemp and stagsigma should have dimensions (1:upn-1).

    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    real(dp), dimension(:), intent(out) :: pmptemp  ! pressure melting point temperature (deg C)
    real(dp), intent(in) :: thck                    ! ice thickness
    real(dp), intent(in), dimension(:) :: stagsigma ! staggered vertical coordinate
                                                    ! (defined at layer midpoints)

    pmptemp(:) = - grav * rhoi * pmlt * thk0 * thck * stagsigma(:)

  end subroutine glissade_calcpmpt

!-----------------------------------------------------------------------

  subroutine glissade_calcpmpt_bed(pmptemp_bed, thck)

    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    real(dp), intent(out) :: pmptemp_bed ! pressure melting point temp at bed (deg C)
    real(dp), intent(in) :: thck         ! ice thickness

    pmptemp_bed = - grav * rhoi * pmlt * thk0 * thck 

  end subroutine glissade_calcpmpt_bed

!-------------------------------------------------------------------

end module glissade_enthalpy

!--------------------------------------------------------------------------------
