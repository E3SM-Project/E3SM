!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   glissade_enthalpy.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)
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
!
! This module computes temperature diffusion, strain heating, and local
!  melting and refreezing in each ice column using enthalpy as the
!  primary state variable.

#include "glide_mask.inc"

module glissade_enthalpy

    use glimmer_global, only : dp
    use glide_types
    use glimmer_log 

    implicit none

    private
    public :: glissade_init_enthalpy, glissade_enthalpy_findvtri, &
              enth2temp, temp2enth, glissade_enthalpy_calcbmlt

!--------------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------------

  subroutine glissade_init_enthalpy (model)

    ! initialization for the enthalpy scheme
    use glimmer_physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow, trpt
    use glimmer_paramets, only : thk0, tim0

    type(glide_global_type),intent(inout) :: model    ! ice model parameters
    integer :: up

    !BDM - dups(1) and dups(2) are done in glissade_init_temp

    ! write out cons(1) which is coefficient before diffusion terms 
    ! thk0 is just a scaling parameter, not actual thickness at ew, ns
    ! Usually thk0 = 1?
    model%tempwk%cons(1) = tim0 * model%numerics%dttem/ (2.0d0 * thk0**2)

  end subroutine glissade_init_enthalpy

!--------------------------------------------------------------------------------

  subroutine glissade_enthalpy_findvtri (model, ew,   ns,          &
                                         subd,  diag, supd, rhsd,  &
                                         float)

    ! solve for tridiagonal entries of sparse matrix
    use glimmer_physcon,  only : rhoi, shci, lhci, rhow, grav, coni
    use glimmer_paramets, only : thk0, tim0

    ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
    ! whereas temperature/enthalpy is indexed from 0 to upn.
    ! The first row of the matrix is the equation for enthalpy(0,ew,ns),
    ! the last row is the equation for enthalpy(upn,ew,ns), and so on.

    !I/O variables
    type(glide_global_type), intent(inout) :: model
    integer, intent(in) :: ew, ns
    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd
    logical, intent(in) :: float
	
    ! local variables
    real(dp) :: dsigbot  ! bottom layter thicknes in sigma coords.
    real(dp) :: alphai ! cold ice diffusivity
    real(dp) :: alpha0 ! temperate ice diffusivity
    real(dp) :: fact ! coefficient in tridiag
    integer  :: up
    real(dp), dimension(0:model%general%upn) :: alpha !diffusivity vector
    real(dp), dimension(0:model%general%upn) :: enthalpy !specific enthalpy
    real(dp), dimension(0:model%general%upn) :: pmptemp !pmptemp all nodes
    real(dp), dimension(1:model%general%upn) :: alpha_half !half-node diffusivity
	
    ! define diffusivities alpha_i and alpha_0
    alphai = coni / rhoi / shci
    alpha0 = alphai / 100.0d0
	
    ! convert model%temper%temp and model%temper%waterfrac to enthalpy
    ! using temp from last timestep.  For interior and boundary nodes.
    ! BDM enthalpy will be size 0:upn
    call temp2enth(enthalpy(0:model%general%upn),                       &
                   model%temper%temp(0:model%general%upn,ew,ns),        &
                   model%temper%waterfrac(1:model%general%upn-1,ew,ns), &
                   model%geometry%thck(ew,ns),                          &
                   model%numerics%stagsigma(1:model%general%upn-1))	

    ! find pmptemp for this column (interior nodes and boundary)
    pmptemp(0) = 0.0d0
    call glissade_calcpmpt(pmptemp(1:model%general%upn-1), model%geometry%thck(ew,ns), &
                           model%numerics%stagsigma(1:model%general%upn-1))
							
    call glissade_calcpmpt_bed(pmptemp(model%general%upn), model%geometry%thck(ew,ns))

    ! create a column vector of size (0:upn) of diffusivity based on 
    ! previous timestep's temp.  Boundary nodes need a value so half-node
    ! diffusivity can be calculated at interior nodes (1:upn-1)
    do up = 0,model%general%upn
       if (model%temper%temp(up,ew,ns) < pmptemp(up)) then
          alpha(up) = alphai
       else
          alpha(up) = alpha0
       endif
    end do
	
    ! Find half-node diffusivity using harmonic average between nodes.
    ! The vector will be size (1:upn) - the first value is the half-node
    ! between nodes 0 and 1, the last value is the half-node between
    ! nodes upn-1 and upn. 
    do up = 1,model%general%upn	
       alpha_half(up) = 2 / ((1/alpha(up-1)) + (1/alpha(up)))
    end do
	
    ! Compute subdiagonal, diagonal, and superdiagonal matrix elements
    
    ! upper boundary: set to surface air temperature
    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0
    rhsd(1) = dmin1(0.0d0,dble(model%climate%artm(ew,ns))) * rhoi * shci
  
    ! cons(1) is (dt * tim0) / (2 * thk0^2)
    ! fact is (dt * tim0) / (2 * H^2 * thk0^2)
    fact = model%tempwk%cons(1) / model%geometry%thck(ew,ns)**2

    ! ice interior. layers 1:upn-1  (matrix elements 2:upn)
    subd(2:model%general%upn) = -fact * alpha_half(1:model%general%upn-1)        &
                                * model%tempwk%dups(1:model%general%upn-1,1)
    supd(2:model%general%upn) = -fact * alpha_half(2:model%general%upn)          &
                                * model%tempwk%dups(1:model%general%upn-1,2)
    diag(2:model%general%upn) = 1.0d0 - subd(2:model%general%upn)                &
                                - supd(2:model%general%upn)
    rhsd(2:model%general%upn) =                                                  &
           enthalpy(1:model%general%upn-1) * (2.0d0 - diag(2:model%general%upn)) &
         - enthalpy(0:model%general%upn-2) * subd(2:model%general%upn)           &
         - enthalpy(2:model%general%upn)   * supd(2:model%general%upn)           & 
         + model%tempwk%dissip(1:model%general%upn-1,ew,ns) * rhoi * shci
    ! BDM I'm assuming that model%tempwk%dissip has units of phi/rhoi/shci.
    ! For an enthalpy calc, we want just phi, so model%tempwk%dissip * rhoi * shci
	
    ! basal boundary:
    ! for grounded ice, a heat flux is applied
    ! for floating ice, the basal temperature is held constant

    !WHL - This lower BC is different from the one in standard glide_temp.
    !      If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
    !       That is, net heat flux at lower boundary must equal zero.
    !      If T(upn) >= Tpmp, then set T(upn) = Tpmp

    if (float) then
       supd(model%general%upn+1) = 0.0d0 
       subd(model%general%upn+1) = 0.0d0
       diag(model%general%upn+1) = 1.0d0
       rhsd(model%general%upn+1) = enthalpy(model%general%upn)
    
    else    ! grounded ice

    !Positive-Thickness Basal Temperate Boundary Layer
    if (abs(model%temper%temp(model%general%upn-1,ew,ns) -       &
        pmptemp(model%general%upn-1)) < 0.001d0) then   
       
       subd(model%general%upn+1) = -1.0d0
       supd(model%general%upn+1) =  0.0d0 
       diag(model%general%upn+1) = 1.0d0 
       rhsd(model%general%upn+1) = 0.0d0

    !Zero-Thickness Basal Temperate Boundary Layer
    elseif (abs(model%temper%temp(model%general%upn,ew,ns) -       &
        pmptemp(model%general%upn)) < 0.001d0) then  ! melting

       ! hold basal temperature at pressure melting point
       supd(model%general%upn+1) = 0.0d0
       subd(model%general%upn+1) = 0.0d0
       diag(model%general%upn+1) = 1.0d0
       rhsd(model%general%upn+1) = pmptemp(model%general%upn) * rhoi * shci

    else  
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
    real(dp), dimension(0:), intent(in)                  :: enthalpy !enthalpy is (0:upn)
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
    if(enthalpy(0) >= pmpenthalpy(0)) then !temperate ice
       temp(0) = pmptemp(0) !temperate ice
    else	
       temp(0) = enthalpy(0) / rhoi / shci !temp is cold
    endif		
	
    do up = 1, upn-1
       if(enthalpy(up) >= pmpenthalpy(up)) then !temperate ice
          temp(up) = pmptemp(up) !temp = pmptemp
          waterfrac(up) = (enthalpy(up)-pmpenthalpy(up)) /                 &
                          ((rhow-rhoi) * shci * pmptemp(up) + rhow * lhci)
       else !cold ice
          temp(up) = enthalpy(up) / rhoi / shci !temp is cold
          waterfrac(up) = 0.0d0 !waterfrac = 0
       endif
    end do
	
    if(enthalpy(upn) >= pmpenthalpy(upn)) then !temperate ice
       temp(upn) = pmptemp(upn) !temp = pmptemp
    else	
       temp(upn) = enthalpy(upn) / rhoi / shci !temp is cold
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
	
    pmpenthalpy = rhoi * shci * pmptemp
	
    !solve for enthalpy
    enthalpy(0) = temp(0) * rhoi * shci
    do up = 1, upn-1
       enthalpy(up) = ((1 - waterfrac(up)) * rhoi * shci * temp(up))          &
                      + waterfrac(up) * rhow * ((shci * pmptemp(up)) + lhci)
    end do
    enthalpy(upn) = temp(upn) * rhoi * shci
	
  end subroutine temp2enth

!--------------------------------------------------------------------------------

  subroutine glissade_enthalpy_calcbmlt(model,                &
                                        temp,                 &
                                        waterfrac, stagsigma, &
                                        thck,      stagthck,  &
                                        bmlt,      floater)

    ! Compute the amount of basal melting.
    ! The basal melting computed here is applied to the ice thickness
    !  by glissade_transport_driver, conserving mass and energy.
    !
    ! This is done with the enthalpy formulation by taking any internal
    ! water content above 1% and draining it to the bed

    use glimmer_physcon, only: shci, rhoi, lhci
    use glimmer_paramets, only : thk0, tim0

    type(glide_global_type) :: model

    real(dp), dimension(0:,:,:), intent(inout) :: temp
    real(dp), dimension(1:,:,:), intent(inout) :: waterfrac
    real(dp), dimension(0:),     intent(in) :: stagsigma
    real(dp), dimension(:,:),    intent(in) :: thck,  stagthck  
    real(dp), dimension(:,:),    intent(out):: bmlt    ! scaled melt rate (m/s * tim0/thk0)
                                                       ! > 0 for melting, < 0 for freeze-on
    logical,  dimension(:,:),    intent(in) :: floater

    real(dp), dimension(size(stagsigma))    :: pmptemp   ! pressure melting point temperature
    real(dp) :: bflx    ! heat flux available for basal melting (W/m^2)
    real(dp) :: hmlt    ! scaled depth of internal melting (m/thk0)
    integer :: up, ew, ns

    bmlt(:,:) = 0.0d0

    !LOOP TODO - This loop should be over locally owned cells? (ilo:ihi,jlo:jhi)

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1

          if (thck(ew,ns) > model%numerics%thklim .and. .not. floater(ew,ns)) then

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
             bmlt(ew,ns) = bflx * model%tempwk%f(2)   ! f(2) = tim0 / (thk0 * lhci * rhoi)

             ! Add internal melting associated with waterfrac > waterfrac_max (1%)
             ! Note: glissade_calcpmpt does not compute pmpt at the top surface or the bed.

             call glissade_calcpmpt(pmptemp(:), thck(ew,ns),   &
                                    stagsigma(:) )

             do up = 1, model%general%upn-1
                if (waterfrac(up,ew,ns) > 0.01) then
                   hmlt = (waterfrac(up,ew,ns) - 0.01) * (model%numerics%sigma(up+1) - &
                          model%numerics%sigma(up))
                   bmlt(ew,ns) = bmlt(ew,ns) + (hmlt * tim0 / model%numerics%dttem)
                   waterfrac(up,ew,ns) = 0.01
                endif
             enddo

             ! Reset basal temp to pmptemp, if necessary

             up = model%general%upn
             call glissade_calcpmpt_bed(pmptemp(up), thck(ew,ns))
             temp(up,ew,ns) = min (temp(up,ew,ns), pmptemp(up))

          endif   ! thk > thklim

       enddo
    enddo

  end subroutine glissade_enthalpy_calcbmlt

!-----------------------------------------------------------------------------------

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

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp(:) = fact * thck * stagsigma(:)

  end subroutine glissade_calcpmpt

!-----------------------------------------------------------------------

  subroutine glissade_calcpmpt_bed(pmptemp_bed, thck)

    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    real(dp), intent(out) :: pmptemp_bed ! pressure melting point temp at bed (deg C)
    real(dp), intent(in) :: thck         ! ice thickness

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp_bed = fact * thck 

  end subroutine glissade_calcpmpt_bed

  !-------------------------------------------------------------------

end module glissade_enthalpy

!--------------------------------------------------------------------------------
