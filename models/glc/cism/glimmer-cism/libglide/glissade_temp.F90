!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_temp.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

! This module computes temperature diffusion and strain heating
!  in a local column without doing horizontal or vertical advection.
! Temperature advection is done separately, e.g. using the incremental 
!  remapping transport scheme.
! It is assumed here that temperature values are staggered in the
!  vertical compared to the velocity.  That is, the temperature lives
!  at the midpoint of each layer instead of at layer interfaces.
! As in the unstaggered case, the temperature is also defined at the 
!  upper and lower surfaces with appropriate boundary conditions.   

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade_temp

    use glimmer_global, only : dp 
    use glide_types
    use glimmer_log

    implicit none

    private
    public :: glissade_init_temp, glissade_temp_driver, glissade_calcflwa

contains

!****************************************************    

  subroutine glissade_init_temp (model)

    ! initialization subroutine for the case that temperature lives on the
    ! vertically staggered grid (i.e., at layer centers)

!TODO - Remove scaling parameters from this module
!       Note: if thk0 = 1, then tau0 = rhoi*grav

    use glimmer_physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow, trpt
    use glimmer_paramets, only : tim0, thk0, len0, vis0, vel0, tau0
    use glide_bwater, only: find_dt_wat
    use parallel, only: lhalo, uhalo
    use glissade_enthalpy, only: glissade_init_enthalpy

    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.

    integer, parameter :: p1 = gn + 1  
    integer up, ns, ew
    real(dp) :: estimate

    !TODO - Should these allocations be done in glide_allocarr?
    !TODO -  Make sure the arrays allocated here are deallocated at the end of the run.
    !        Might want to move allocation/deallocation to subroutines in glide_types.

    ! Note vertical dimensions here.  Dissipation is computed for each of (upn-1) layers.
    ! Temperature is defined at midpoint of each layer, plus upper and lower surfaces.
    allocate(model%tempwk%dups(model%general%upn+1,2))
    allocate(model%tempwk%inittemp(model%general%upn+1,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%dissip  (model%general%upn-1,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%compheat(model%general%upn-1,model%general%ewn,model%general%nsn))
    model%tempwk%compheat = 0.0d0

    allocate(model%tempwk%c1(model%general%upn-1))   ! upn-1 for staggered grid

    allocate(model%tempwk%dupa(model%general%upn))
    allocate(model%tempwk%dupb(model%general%upn))
    allocate(model%tempwk%dupc(model%general%upn))

    allocate(model%tempwk%smth(model%general%ewn,model%general%nsn))  !TODO - Is this used for glissade?
    allocate(model%tempwk%wphi(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatu(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatv(model%general%ewn,model%general%nsn))

    model%tempwk%dups = 0.0d0

!WHL - Note that the 'dups' grid coefficients are not the same as for unstaggered temperatures.

    up = 1
    model%tempwk%dups(up,1) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                    (model%numerics%stagsigma(up) - model%numerics%sigma(up)) )
    do up = 2, model%general%upn-1
       model%tempwk%dups(up,1) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                       (model%numerics%stagsigma(up) - model%numerics%stagsigma(up-1)) )
    enddo

    do up = 1, model%general%upn-2
       model%tempwk%dups(up,2) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                       (model%numerics%stagsigma(up+1) - model%numerics%stagsigma(up)) )
    end do
    up = model%general%upn-1
    model%tempwk%dups(up,2) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                    (model%numerics%sigma(up+1) - model%numerics%stagsigma(up)) )

    model%tempwk%zbed = 1.0d0 / thk0
    model%tempwk%dupn = model%numerics%sigma(model%general%upn) - model%numerics%sigma(model%general%upn-1)

!WHL - We need only two of these (cons(1) and cons(5)) for the staggered case with no advection.
!!    model%tempwk%cons = (/ 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2), &
!!         model%numerics%dttem / 2.0d0, &
!!         VERT_DIFF*2.0d0 * tim0 * model%numerics%dttem / (thk0 * rhoi * shci), &
!!         VERT_ADV*tim0 * acc0 * model%numerics%dttem / coni, &
!!         0.0d0, &   !WHL - no vertical advection
!!         ( tau0 * vel0 / len0 ) / ( rhoi * shci ) * ( model%numerics%dttem * tim0 ) /)  
         !*sfp* added last term to vector above for use in HO & SSA dissip. cacl

!TODO - Give these constants (tempwk%cons, tempwk%f, tempwk%c) better names.

!WHL - The factor of 2 in the numerator in the original code can be traced to a missing factor of 0.5
!      in the denominator of the dups coefficients.  On the vertically staggered grid, there is no
!      factor of 0.5 in the dups coefficients, so there is no factor of 2 here.
!WHL - The factor of 2 in the denominator is a Crank-Nicolson averaging factor.

!!    model%tempwk%cons(1) = 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)
    model%tempwk%cons(1) = tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)

    model%tempwk%cons(5) = (tau0 * vel0 / len0 ) / (rhoi * shci) * (model%numerics%dttem * tim0)

    !Note: stagsigma here instead of sigma
    model%tempwk%c1(1:model%general%upn-1) =   &
                             (model%numerics%stagsigma(1:model%general%upn-1)   &
                             * rhoi * grav * thk0**2 / len0)**p1 * 2.0d0 * vis0              &
                             * model%numerics%dttem * tim0 / (16.0d0 * rhoi * shci)

    model%tempwk%f = (/ tim0 * coni / (thk0**2 * lhci * rhoi), &
                        tim0 / (thk0 * lhci * rhoi), &
                        tim0 * thk0 * rhoi * shci /  (thk0 * tim0 * model%numerics%dttem * lhci * rhoi), &
                        tim0 * thk0**2 * vel0 * grav * rhoi / (4.0d0 * thk0 * len0 * rhoi * lhci), &
                        tim0 * vel0 * tau0 / (4.0d0 * thk0 * rhoi * lhci) /)      
                        !*sfp* added the last term in the vect above for HO and SSA dissip. calc. 

    select case(model%options%whichbwat)

       case(BWATER_LOCAL)
          model%paramets%hydtim = tim0 / (model%paramets%hydtim * scyr)
          estimate = 0.2d0 / model%paramets%hydtim
          
          model%tempwk%c = (/ model%tempwk%dt_wat,   &
                              1.0d0 - 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, &
                              1.0d0 + 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, &
                              0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /) 

       case(BWATER_FLUX) ! steady-state routing using flux calculation

         !TODO - Test this option for one-processor runs.
         !       It has not been parallelized.

          model%tempwk%watvel = model%paramets%hydtim * tim0 / (scyr * len0)
          estimate = (0.2d0 * model%tempwk%watvel) / min(model%numerics%dew,model%numerics%dns)
          call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat)

          !print *, model%numerics%dttem*tim0/scyr, model%tempwk%dt_wat*tim0/scyr, model%tempwk%nwat

          model%tempwk%c = (/ rhow * grav, rhoi * grav, 2.0d0 * model%numerics%dew, 2.0d0 * model%numerics%dns, &
	         0.25d0 * model%tempwk%dt_wat / model%numerics%dew, 0.25d0 * model%tempwk%dt_wat / model%numerics%dns, &
               0.5d0 * model%tempwk%dt_wat / model%numerics%dew, 0.5d0 * model%tempwk%dt_wat / model%numerics%dns /)
          
    end select
 
      !==== Initialize ice temperature.============
      !This block of code is similar to that in glide_init_temp

    ! Five possibilities:
    ! (1) Set ice temperature to 0 C everywhere in column (TEMP_INIT_ZERO)
    ! (2) Set ice temperature to surface air temperature everywhere in column (TEMP_INIT_ARTM)
    ! (3) Set up a linear temperature profile, with T = artm at the surface and T <= Tpmp
    !     at the bed (TEMP_INIT_LINEAR). 
    !     A parameter (pmpt_offset) controls how far below Tpmp the initial bed temp is set.
    ! (4) Read ice temperature from an initial input file.
    ! (5) Read ice temperature from a restart file.
    !
    ! If restarting, we always do (5).
    ! If not restarting and the temperature field is present in the input file, we do (4).
    ! If (4) or (5), then the temperature field should already have been read from a file,
    !  and the rest of this subroutine will do nothing.
    ! Otherwise, the initial temperature is controlled by model%options%temp_init,
    !  which can be read from the config file.
    !
    !TODO - For reading from restart or input file, make sure that halo values are correct.

    if (model%options%is_restart == RESTART_TRUE) then

       ! Temperature has already been initialized from a restart file. 
       ! (Temperature is always a restart variable.)

       call write_log('Initializing ice temperature from the restart file')

    elseif ( minval(model%temper%temp(1:model%general%upn, &
                    1+lhalo:model%general%ewn-lhalo, 1+uhalo:model%general%nsn-uhalo)) > &
                    (-1.0d0 * trpt) ) then    ! trpt = 273.15 K
                                              ! Default initial temps in glide_types are -999

       ! Temperature has already been initialized from an input file.
       ! (We know this because the default initial temps of -999 have been overwritten.)

       call write_log('Initializing ice temperature from an input file')

    else   ! not reading temperature from restart or input file
           ! initialize it here basee on model%options%temp_init

       ! First set T = 0 C everywhere

       model%temper%temp(:,:,:) = 0.0d0                                              
                                                    
       if (model%options%temp_init == TEMP_INIT_ZERO) then

          call write_log('Initializing ice temperature to 0 deg C')

          ! No call is needed to glissade_init_temp_column because the
          ! ice temperature has been set to zero above

       elseif (model%options%temp_init == TEMP_INIT_ARTM) then

          ! Initialize ice column temperature to min(artm, 0 C).

          !Note: Old glide sets temp = artm everywhere without regard to whether ice exists in a column.
          !TODO - Verify that this makes no difference for model results.

          call write_log('Initializing ice temperature to the surface air temperature')

          !TODO - Locally owned cells only?
          do ns = 1, model%general%nsn
             do ew = 1, model%general%ewn

                call glissade_init_temp_column(model%options%temp_init,         &
                                               model%numerics%stagsigma(:),     &
                                               dble(model%climate%artm(ew,ns)), &
                                               model%geometry%thck(ew,ns),      &
                                               model%temper%temp(:,ew,ns) )
             end do
          end do

       elseif (model%options%temp_init == TEMP_INIT_LINEAR) then

          ! Initialize ice column temperature with a linear profile:
          ! T = artm at the surface, and T <= Tpmp at the bed.
          ! Loop over physical cells where artm is defined (not temperature halo cells)

          call write_log('Initializing ice temperature to a linear profile in each column')

          !TODO - Locally owned cells only?
          do ns = 1, model%general%nsn
             do ew = 1, model%general%ewn

                call glissade_init_temp_column(model%options%temp_init,         &
                                               model%numerics%stagsigma(:),     &
                                               dble(model%climate%artm(ew,ns)), &
                                               model%geometry%thck(ew,ns),      &
                                               model%temper%temp(:,ew,ns) )

             end do
          end do

       endif ! model%options%temp_init

    endif    ! restart file, input file, or other options

!WHL - Removed glissade_calcflwa call here; now computed in glissade_diagnostic_variable_solve.
!TODO - If only locally owned values are filled here, then make sure there is a halo update in glissade_initialise. 

    !BDM - make call to glissade_enthalpy_init if using enthalpy approach
    if (model%options%whichtemp == TEMP_ENTHALPY) then
       call glissade_init_enthalpy(model)
    end if

  end subroutine glissade_init_temp

!****************************************************    

  subroutine glissade_init_temp_column(temp_init,                 &
                                       stagsigma,   artm,         &
                                       thck,        temp)

  ! Initialize temperatures in a column based on the value of temp_init.
  ! Threee possibilities:
  ! (1) Set ice temperature in column to 0 C (TEMP_INIT_ZERO)
  ! (2) Set ice temperature in column to surface air temperature (TEMP_INIT_ARTM)
  ! (3) Set up a linear temperature profile, with T = artm at the surface and T <= Tpmp
  !     at the bed (TEMP_INIT_LINEAR). 
  !     A local parameter (pmpt_offset) controls how far below Tpmp the initial bed temp is set.
  !
  ! This subroutine is functionally equivalent to glide_init_temp_column.
  ! The only difference is that temperature is staggered in the vertical
  !  (i.e., located at layer midpoints as well as the top and bottom surfaces).

  ! In/out arguments
 
  integer, intent(in) :: temp_init          ! option for temperature initialization

  real(dp), dimension(:), intent(in)     :: stagsigma  ! staggered vertical coordinate
                                                       ! includes layer midpoints, but not top and bottom surfaces
  real(dp), intent(in)                   :: artm   ! surface air temperature (deg C)
                                                   ! Note: artm should be passed in as double precision
  real(dp), intent(in)                   :: thck   ! ice thickness
  real(dp), dimension(0:), intent(inout) :: temp   ! ice column temperature (deg C)
                                                   ! Note first index of zero
                                                   
  ! Local variables and parameters

  real(dp) :: pmptb                              ! pressure melting point temp at the bed
  real(dp), dimension(size(stagsigma)) :: pmpt   ! pressure melting point temp thru the column
  integer :: upn                                 ! number of vertical levels (deduced from temp array)

!TODO - Define elsewhere? (glimmer_paramets?)
  real(dp), parameter :: pmpt_offset = 2.d0  ! offset of initial Tbed from pressure melting point temperature (deg C)
                                             ! Note: pmtp_offset is positive for T < Tpmp

  upn = size(temp) - 1     ! temperature array has dimension (0:model%general%upn)

  ! Set the temperature in the column

  select case(temp_init)

  case(TEMP_INIT_ZERO)     ! set T = 0 C

     temp(:) = 0.d0

  case(TEMP_INIT_ARTM)     ! initialize ice-covered areas to the min of artm and 0 C
                           ! set ice-free areas to T = 0 C

     if (thck > 0.0d0) then
        temp(:) = dmin1(0.0d0, artm)  !TODO - dmin1 --> min?
     else
        temp(:) = 0.d0
     endif

  case(TEMP_INIT_LINEAR)

     ! Tsfc = artm, Tbed = Tpmp = pmpt_offset, linear profile in between 

     temp(0) = artm

     call glissade_calcpmpt_bed (pmptb, thck)
     temp(upn) = pmptb - pmpt_offset

     temp(1:upn-1) = temp(0) + (temp(upn) - temp(0))*stagsigma(:)
                               
     ! Make sure T <= Tpmp - pmpt_offset in column interior
     ! TODO: Change condition to T <= Tpmp?

     call glissade_calcpmpt(pmpt(:), thck, stagsigma(:))
     temp(1:upn-1) = min(temp(1:upn-1), pmpt(1:upn-1) - pmpt_offset)

  end select

  end subroutine glissade_init_temp_column

!****************************************************    

  subroutine glissade_temp_driver(model, whichtemp)

    ! Calculates the ice temperature 

    use glimmer_utils,  only : tridiag
    use glimmer_paramets, only : thk0, tim0
    use glimmer_physcon, only: shci, coni, rhoi
    use glide_mask
    use glide_bwater
    use glissade_enthalpy

    !TODO - Use glam_grid_operators instead?
    use glide_grid_operators, only: stagvarb

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_global_type),intent(inout) :: model       ! Ice model parameters.
    integer,                intent(in)    :: whichtemp    ! Flag to choose method.

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    integer :: ew, ns, up, upn
    character(len=100) :: message

    ! These arrays have the same size as the vertical dimension of temperature:
    !  upn+1 on the staggered grid.
    real(dp), dimension(size(model%temper%temp,1)) :: subd, diag, supd, rhsd

    ! These have the same dimensions as staggered temperature
    real(dp),dimension(0:model%general%upn) :: Tstagsigma, prevtemp_stag, enthalpy

    ! for energy conservation check
    real(dp) :: einit, efinal, delta_e, dTtop, dTbot

    upn = model%general%upn

    select case(whichtemp)

    case(TEMP_SURFACE_AIR_TEMP)  ! Set column to surface air temperature ------------------

       ! JEFF - OK for distributed since using air temperature at grid point to initialize.

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(model%climate%artm(ew,ns)))
          end do
       end do

    case(TEMP_PROGNOSTIC) ! Local column calculation (with advection done elsewhere)

            ! No horizontal or vertical advection; vertical diffusion and strain heating only.
            ! Temperatures are vertically staggered relative to velocities.  
            ! That is, the temperature is defined at the midpoint of each layer 
            ! (and at the top and bottom surfaces).

       !TODO - Change to stagwbndsigma
       ! Set Tstagsigma (= stagsigma except that it has values at the top and bottom surfaces).

       Tstagsigma(0) = 0.d0
       Tstagsigma(1:model%general%upn-1) = model%numerics%stagsigma(1:model%general%upn-1)
       Tstagsigma(model%general%upn) = 1.d0

       model%tempwk%inittemp = 0.0d0

       ! Calculate interior heat dissipation -------------------------------------

       call glissade_finddisp( model,                   &
                               model%geometry%thck,     &
                               model%options%which_disp,&
                               model%stress%efvs,       &
                               model%geomderv%stagthck, &
                               model%geomderv%dusrfdew, &
                               model%geomderv%dusrfdns, &
                               model%temper%flwa)

       ! Calculate heating from basal friction -----------------------------------

       !TODO - dusrfdew/dns are not needed as inputs
       call glissade_calcbfric( model,                        &
                                model%geometry%thck,          &
                                model%velocity%btraction,     &
                                model%velocity%ubas,          &
                                model%velocity%vbas,          &
                                GLIDE_IS_FLOAT(model%geometry%thkmask) )

       ! Note: No iteration is needed here since we are doing a local tridiagonal solve without advection.

       !LOOP TODO: Change to 1+lhalo, nsn-uhalo?
       !           If so, don't forget to do halo updates later for temp and flwa.

       do ns = 2,model%general%nsn-1
       do ew = 2,model%general%ewn-1
          if(model%geometry%thck(ew,ns) > model%numerics%thklim) then

             ! compute initial internal energy in column (for energy conservation check)
             einit = 0.0d0
             do up = 1, upn-1
                einit = einit + model%temper%temp(up,ew,ns) *  &
                               (model%numerics%sigma(up+1) -   &
                                model%numerics%sigma(up) )
             enddo
             einit = einit * rhoi * shci * model%geometry%thck(ew,ns)*thk0

             ! compute matrix elements

             call glissade_findvtri( model, ew,   ns,             &
                                     subd,  diag, supd, rhsd,     &            
                                     GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))
                
             prevtemp_stag(:) = model%temper%temp(:,ew,ns)

             ! solve the tridiagonal system

             ! Note: Temperature is indexed from 0 to upn, with indices 1 to upn-1 colocated
             !  with stagsigma values of the same index.
             ! However, the matrix elements are indexed 1 to upn+1, with the first row
             !  corresponding to the surface temperature, temp(0,:,:).

             call tridiag(subd(1:model%general%upn+1), &
                          diag(1:model%general%upn+1), &
                          supd(1:model%general%upn+1), &
                          model%temper%temp(0:model%general%upn,ew,ns), &
                          rhsd(1:model%general%upn+1))

             ! Check that the net input of energy to the column is equal to the difference
             !  between the initial and final internal energy.
             !TODO - Make this check optional and/or move it to a subroutine.

             ! compute the final internal energy

             efinal = 0.0d0
             do up = 1, upn-1
                efinal = efinal + model%temper%temp(up,ew,ns) *  &
                                 (model%numerics%sigma(up+1) - model%numerics%sigma(up))
             enddo
             efinal = efinal * rhoi*shci * model%geometry%thck(ew,ns)*thk0

             ! compute net heat flux to the column

             ! conductive flux = (k/H * dT/dsigma) at upper and lower surfaces; positive down

             dTtop = 0.5d0 * ( model%temper%temp(1,ew,ns) - model%temper%temp(0,ew,ns) &
                             +     prevtemp_stag(1)       -     prevtemp_stag(0) )
             model%temper%ucondflx(ew,ns) = (-coni / (model%geometry%thck(ew,ns)*thk0) )         &
                                           * dTtop / (Tstagsigma(1) - Tstagsigma(0))


             dTbot = 0.5d0 * ( model%temper%temp(upn,ew,ns) - model%temper%temp(upn-1,ew,ns) &
                            +     prevtemp_stag(upn)       -   prevtemp_stag(upn-1) )
             model%temper%lcondflx(ew,ns) = (-coni / (model%geometry%thck(ew,ns)*thk0) )         &
                                           * dTbot / (Tstagsigma(upn) - Tstagsigma(upn-1))

             ! total dissipation in column

             model%temper%dissipcol(ew,ns) = 0.0d0
             do up = 1, upn-1
                model%temper%dissipcol(ew,ns) = model%temper%dissipcol(ew,ns) + &
                                              model%tempwk%dissip(up,ew,ns)  &
                                           * (model%numerics%sigma(up+1) - model%numerics%sigma(up))  
             enddo 
             model%temper%dissipcol(ew,ns) = model%temper%dissipcol(ew, ns)     &
                                     * thk0*model%geometry%thck(ew,ns)*rhoi*shci / (tim0*model%numerics%dttem)  

             ! Verify that the net input of energy into the column is equal to the change in
             ! internal energy.  

             delta_e = (model%temper%ucondflx(ew,ns) - model%temper%lcondflx(ew,ns)  &
                      + model%temper%dissipcol(ew,ns)) * tim0*model%numerics%dttem

             if ( abs((efinal-einit-delta_e)/(tim0*model%numerics%dttem)) > 1.0d-8 ) then
                write(message,*) 'WARNING: Energy conservation error, ew, ns =', ew, ns
                call write_log(message)
! Can uncomment the following for diagnostics
!                write(50,*) 'Interior fluxes:'
!                write(50,*) 'ftop (pos up)=', -model%temper%ucondflx(ew,ns) 
!                write(50,*) 'fbot (pos up)=', -model%temper%lcondflx(ew,ns)
!                write(50,*) 'fdissip =',       model%temper%dissipcol(ew,ns)
!                write(50,*) 'Net flux =', delta_e/(tim0*model%numerics%dttem)
!                write(50,*) ' '
!                write(50,*) 'delta_e =', delta_e
!                write(50,*) 'einit =',  einit
!                write(50,*) 'efinal =', efinal
!                write(50,*) 'einit + delta_e =', einit + delta_e
!                write(50,*) ' '
!                write(50,*) 'Energy imbalance =', efinal - einit - delta_e
!                write(50,*) ' '
!                write(50,*) 'Basal fluxes:'
!                write(50,*) 'ffric =', model%temper%bfricflx(ew,ns)
!                write(50,*) 'fgeo =', -model%temper%bheatflx(ew,ns)
!                write(50,*) 'flux for bottom melting =', model%temper%bfricflx(ew,ns)   &
!                                                       - model%temper%bheatflx(ew,ns)   &
!                                                       + model%temper%lcondflx(ew,ns)
             endif

!WHL - No call here to corrpmpt.  Temperatures above pmpt are set to pmpt 
!      in glissade_calcbmlt (conserving energy).

          endif  ! thck > thklim
       end do    ! ew
       end do    ! ns

       ! set temperature of thin ice to the air temperature and set ice-free nodes to zero

      !TODO - Loop over locally owned cells only?
       do ns = 1, model%general%nsn
          do ew = 1, model%general%ewn

             if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
                model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
             else if (model%geometry%thkmask(ew,ns) < 0) then
                model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
             !else if (model%geometry%thkmask(ew,ns) < -1) then
             !   model%temper%temp(:,ew,ns) = 0.0d0
             end if

!TODO - Maybe it should be done in the following way, so that the temperature profile for thin ice
!       is consistent with the temp_init option, with T = 0 for ice-free cells.

             ! NOTE: Calling this subroutine will maintain a sensible temperature profile
             !        for thin ice, but in general does *not* conserve energy.
             !       To conserve energy, we need either thklim = 0, or some additional
             !        energy accounting and correction.
 
             !TODO - Why not simply 0 < thck < thklim?
!             if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
!                call glissade_init_temp_column(model%options%temp_init,         &
!                                               model%numerics%stagsigma(:),     &
!                                               dble(model%climate%artm(ew,ns)), &
!                                               model%geometry%thck(ew,ns),      &
!                                               model%temper%temp(:,ew,ns) )
!             else if (model%geometry%thkmask(ew,ns) < 0) then
!                model%temper%temp(:,ew,ns) = 0.d0
!             end if

          end do
       end do

       ! Calculate basal melt rate
       ! Temperature above the pressure melting point are reset to Tpmp,
       !  with excess heat contributing to melting.

!TODO - Some of these arguments are not needed:
!       stagthck, dusrfdew, dusrfdns, ubas, vbas

       call glissade_calcbmlt( model,                     &
                               model%temper%temp,         &
                               Tstagsigma,                &
                               model%geometry%thck,       &
                               model%geomderv%stagthck,   &
                               model%geomderv%dusrfdew,   &
                               model%geomderv%dusrfdns,   &
                               model%velocity%ubas,       &
                               model%velocity%vbas,       &
                               model%temper%bmlt,         &
                               GLIDE_IS_FLOAT(model%geometry%thkmask))

       ! Calculate basal water depth ------------------------------------------------

       !TODO - Is it necessary to pass 'model'?

       call calcbwat( model,                     &
                      model%options%whichbwat,   &
                      model%temper%bmlt,         &
                      model%temper%bwat,         &
                      model%temper%bwatflx,      &
                      model%geometry%thck,       &
                      model%geometry%topg,       &
                      model%temper%temp(model%general%upn,:,:), &
                      GLIDE_IS_FLOAT(model%geometry%thkmask),   &
                      model%tempwk%wphi)

!WHL - Here I restored some Glimmer calls that were removed earlier.
!      We need stagbpmp for one of the basal traction cases.

!TODO - Think about whether Glissade will support the same basal traction cases as Glide.
!TODO - Use a staggered difference routine from glam_grid_operators?

       ! Transform basal temperature and pressure melting point onto velocity grid

       call stagvarb(model%temper%temp(model%general%upn, 1:model%general%ewn, 1:model%general%nsn), &
                     model%temper%stagbtemp ,&
                     model%general%  ewn, &
                     model%general%  nsn)
       
       call glissade_calcbpmp(model, &
                              model%geometry%thck,  &
                              model%temper%bpmp)

       call stagvarb(model%temper%bpmp, &
                     model%temper%stagbpmp ,&
                     model%general%  ewn, &
                     model%general%  nsn)

       !WHL - Removed glissade_calcflwa call here; moved it to glissade_diagnostic_variable_solve.

   case(TEMP_STEADY)! do nothing

   case(TEMP_ENTHALPY)! BDM Local column calculation (with advection done elsewhere)

      ! No horizontal or vertical advection; vertical diffusion and strain heating only.
      ! Enthalpy is vertically staggered relative to velocities.  
      ! That is, enthalpy is defined at the midpoint of each layer 
      ! (and at the top and bottom surfaces).

      ! BDM Enthalpy Gradient Method is used here to solve for temp. and water content.

      ! BDM If I'm using TEMP_ENTHALPY, should I make a call to glissade_init_enthalpy here?

      !TODO - Change to stagwbndsigma
      ! Set Tstagsigma (= stagsigma except that it has values at the top and bottom surfaces).

      Tstagsigma(0) = 0.d0
      Tstagsigma(1:model%general%upn-1) = model%numerics%stagsigma(1:model%general%upn-1)
      Tstagsigma(model%general%upn) = 1.d0

      ! Calculate interior heat dissipation -------------------------------------

      call glissade_finddisp(  model,                   &
                               model%geometry%thck,     &
                               model%options%which_disp,&
                               model%stress%efvs,       &
                               model%geomderv%stagthck, &
                               model%geomderv%dusrfdew, &
                               model%geomderv%dusrfdns, &
                               model%temper%flwa)

      ! Calculate heating from basal friction -----------------------------------

      call glissade_calcbfric( model,                        &
                               model%geometry%thck,          &
                               model%velocity%btraction,     &
                               model%velocity%ubas,          &
                               model%velocity%vbas,          &
                               GLIDE_IS_FLOAT(model%geometry%thkmask) )

      ! Note: No iteration is needed here since we are doing a local tridiagonal solve without advection.

      !LOOP TODO: Change to 1+lhalo, nsn-uhalo?
      !           If so, don't forget to do halo updates later for temp and flwa.

      do ns = 2,model%general%nsn-1
         do ew = 2,model%general%ewn-1
            if(model%geometry%thck(ew,ns) > model%numerics%thklim) then

               ! BDM compute matrix elements using Enthalpy Gradient Method

               call glissade_enthalpy_findvtri(model, ew, ns, subd, diag, supd, rhsd, &
                                               GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))
					
               ! BDM leave as prevtemp because it's only used for dT/dsigma at top and bottom boundaries,
               ! we don't want this as enthalpy
               prevtemp_stag(:) = model%temper%temp(:,ew,ns)

               ! solve the tridiagonal system
               ! Note: Temperature is indexed from 0 to upn, with indices 1 to upn-1 colocated
               ! with stagsigma values of the same index.
               ! However, the matrix elements are indexed 1 to upn+1, with the first row
               ! corresponding to the surface temperature, temp(0,:,:).

               call tridiag(subd(1:model%general%upn+1),   &
                            diag(1:model%general%upn+1),   &
                            supd(1:model%general%upn+1),   &
                            enthalpy(0:model%general%upn), &
                            rhsd(1:model%general%upn+1))
							  
               ! BDM convert back to temperature and water content
               call enth2temp(enthalpy(0:model%general%upn),                       &
                              model%temper%temp(0:model%general%upn,ew,ns),        &
                              model%temper%waterfrac(1:model%general%upn-1,ew,ns), &
                              model%geometry%thck(ew,ns),                          &
                              model%numerics%stagsigma(1:model%general%upn-1))	

            endif  ! thck > thklim
         end do    ! ew
      end do    ! ns

      ! set temperature of thin ice to the air temperature and set ice-free nodes to zero

      !TODO - Loop over locally owned cells only?
      do ns = 1, model%general%nsn
         do ew = 1, model%general%ewn

            if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
               model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
            else if (model%geometry%thkmask(ew,ns) < 0) then
               model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
            !else if (model%geometry%thkmask(ew,ns) < -1) then
            !   model%temper%temp(:,ew,ns) = 0.0d0
            end if

!TODO - Maybe it should be done in the following way, so that the temperature profile for thin ice
!       is consistent with the temp_init option, with T = 0 for ice-free cells.

! NOTE: Calling this subroutine will maintain a sensible temperature profile
!        for thin ice, but in general does *not* conserve energy.
!       To conserve energy, we need either thklim = 0, or some additional
!        energy accounting and correction.
 
!TODO - Why not simply 0 < thck < thklim?
!             if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
!                call glissade_init_temp_column(model%options%temp_init,         &
!                                               model%numerics%stagsigma(:),     &
!                                               dble(model%climate%artm(ew,ns)), &
!                                               model%geometry%thck(ew,ns),      &
!                                               model%temper%temp(:,ew,ns) )
!             else if (model%geometry%thkmask(ew,ns) < 0) then
!                model%temper%temp(:,ew,ns) = 0.d0
!             end if

         end do !ew
      end do !ns

      ! BDM since there will be no temps above PMP, need new subroutine to calculate basal melt
      ! and basal water depth
      call glissade_enthalpy_calcbmlt(model,                     &
                                      model%temper%temp,         &
                                      model%temper%waterfrac,    &
                                      Tstagsigma,                &
                                      model%geometry%thck,       &
                                      model%geomderv%stagthck,   &
                                      model%temper%bmlt,         &
                                      GLIDE_IS_FLOAT(model%geometry%thkmask))

      call calcbwat( model,                                    &
                     model%options%whichbwat,                  &
                     model%temper%bmlt,                        &
                     model%temper%bwat,                        &
                     model%temper%bwatflx,                     &
                     model%geometry%thck,                      &
                     model%geometry%topg,                      &
                     model%temper%temp(model%general%upn,:,:), &
                     GLIDE_IS_FLOAT(model%geometry%thkmask),   &
                     model%tempwk%wphi)

!WHL - Here I restored some Glimmer calls that were removed earlier.
!      We need stagbpmp for one of the basal traction cases.

!TODO - Think about whether Glissade will support the same basal traction cases as Glide.
!TODO - Use a staggered difference routine from glam_grid_operators?

       ! Transform basal temperature and pressure melting point onto velocity grid

      call stagvarb(model%temper%temp(model%general%upn, 1:model%general%ewn, 1:model%general%nsn), &
                    model%temper%stagbtemp ,                                                        &
                    model%general%  ewn,                                                            &
                    model%general%  nsn)
       
      call glissade_calcbpmp(model,                &
                             model%geometry%thck,  &
                             model%temper%bpmp)

      call stagvarb(model%temper%bpmp,      &
                    model%temper%stagbpmp,  &
                    model%general%ewn,      &
                    model%general%nsn)

      !WHL - Removed glissade_calcflwa call here; moved it to glissade_diagnostic_variable_solve.

    end select

  end subroutine glissade_temp_driver

  !-------------------------------------------------------------------------

  subroutine glissade_findvtri (model, ew,   ns,          &
                                subd,  diag, supd, rhsd,  &
                                float)

    ! compute matrix elements for the tridiagonal solve

    use glimmer_paramets, only : thk0
    use glimmer_physcon,  only : rhoi, grav, coni

    ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
    !             whereas temperature is indexed from 0 to upn.
    !            The first row of the matrix is the equation for temp(0,ew,ns),
    !             the second row is the equation for temp(1,ew,ns), and so on.

    type(glide_global_type), intent(inout) :: model
    integer, intent(in) :: ew, ns
    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd
    logical, intent(in) :: float

    ! local variables

    real(dp) :: pmptempb  ! pressure melting temp at bed
    real(dp) :: fact
    real(dp) :: dsigbot  ! bottom layter thicknes in sigma coords.

    ! set surface temperature

    model%temper%temp(0,ew,ns) = dble(model%climate%artm(ew,ns))

    ! Compute subdiagonal, diagonal, and superdiagonal matrix elements

    ! upper boundary: set to surface air temperature

    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0
    rhsd(1) = model%temper%temp(0,ew,ns)

    ! ice interior. layers 1:upn-1  (matrix elements 2:upn)

    ! model%tempwk%cons(1) = 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)

    fact = model%tempwk%cons(1) / model%geometry%thck(ew,ns)**2
    subd(2:model%general%upn) = -fact * model%tempwk%dups(1:model%general%upn-1,1)
    supd(2:model%general%upn) = -fact * model%tempwk%dups(1:model%general%upn-1,2)
    diag(2:model%general%upn) = 1.0d0 - subd(2:model%general%upn)     &
                                      - supd(2:model%general%upn)

    model%tempwk%inittemp(1:model%general%upn-1,ew,ns) =   &
           model%temper%temp(1:model%general%upn-1,ew,ns) * (2.0d0 - diag(2:model%general%upn)) &
         - model%temper%temp(0:model%general%upn-2,ew,ns) * subd(2:model%general%upn) &
         - model%temper%temp(2:model%general%upn,  ew,ns) * supd(2:model%general%upn) & 
         + model%tempwk%dissip(1:model%general%upn-1,ew,ns)
    
    rhsd(2:model%general%upn) = model%tempwk%inittemp(1:model%general%upn-1,ew,ns)

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

       model%tempwk%inittemp(model%general%upn,ew,ns) = model%temper%temp(model%general%upn,ew,ns) 
       rhsd(model%general%upn+1) = model%temper%temp(model%general%upn,ew,ns)

    else    ! grounded ice

!TODO - This call (and those below) could be inlined.
       call glissade_calcpmpt_bed(pmptempb, model%geometry%thck(ew,ns))

       if (abs(model%temper%temp(model%general%upn,ew,ns) - pmptempb) < 0.001d0) then  ! melting

          ! hold basal temperature at pressure melting point

          supd(model%general%upn+1) = 0.0d0
          subd(model%general%upn+1) = 0.0d0
          diag(model%general%upn+1) = 1.0d0

          model%tempwk%inittemp(model%general%upn,ew,ns) = pmptempb
          rhsd(model%general%upn+1) = pmptempb

       else   ! frozen at bed
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

          model%tempwk%inittemp(model%general%upn,ew,ns) =    &
             (model%temper%bfricflx(ew,ns)  - model%temper%bheatflx(ew,ns)) &
             * dsigbot * model%geometry%thck(ew,ns) * thk0 / coni
          rhsd(model%general%upn+1) = model%tempwk%inittemp(model%general%upn,ew,ns)
          
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

  end subroutine glissade_findvtri

  !-----------------------------------------------------------------------

  subroutine glissade_calcbfric (model,                 &
                                 thck,     btraction,   &
                                 ubas,     vbas,        &
                                 float)

    ! compute frictional heat source due to sliding at the bed

    use glimmer_physcon,  only: rhoi, grav
    use glimmer_paramets, only: thk0, vel0, vel_scale

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(in) :: thck
    real(dp), dimension(:,:), intent(in) :: ubas, vbas
    real(dp), dimension(:,:,:), intent(in) :: btraction
    logical, dimension(:,:), intent(in) :: float

    real(dp) :: slterm       ! sliding friction
 
    integer :: ewp, nsp, ew, ns
    integer :: slide_count   ! number of neighbor cells with nonzero sliding

       ! compute heat source due to basal friction
       ! Note: slterm and bfricflx are defined to be >= 0

       !LOOP TODO - This loop should be over locally owned cells? (ilo:ihi,jlo:jhi)
 
       do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1

          slterm = 0.d0
          slide_count = 0


             !WHL - copied Steve Price's formulation from calcbmlt
             ! btraction is computed in glam_strs2.F90

             !TODO - Make sure we have btraction for all locally owned velocity cells.

             if (thck(ew,ns) > model%numerics%thklim .and. .not. float(ew,ns)) then
                do nsp = ns-1,ns
                do ewp = ew-1,ew

!SCALING - WHL: Multiplied ubas by vel0/vel_scale so we get the same result in these two cases:
!           (1) Old Glimmer with scaling:         vel0 = vel_scale = 500/scyr, and ubas is non-dimensional
!           (2) New Glimmer-CISM without scaling: vel0 = 1, vel_scale = 500/scyr, and ubas is in m/s.

!!!                   if (abs(model%velocity%ubas(ewp,nsp)) > 1.0d-6 .or.   &
!!!                       abs(model%velocity%vbas(ewp,nsp)) > 1.0d-6) then
                   if ( abs(model%velocity%ubas(ewp,nsp))*(vel0/vel_scale) > 1.0d-6 .or.   &
                        abs(model%velocity%vbas(ewp,nsp))*(vel0/vel_scale) > 1.0d-6 ) then
                      slide_count = slide_count + 1
                      slterm = slterm + model%velocity%btraction(1,ewp,nsp) * &
                                        model%velocity%uvel(model%general%upn,ewp,nsp) &
                                      + model%velocity%btraction(2,ewp,nsp) * &
                                        model%velocity%vvel(model%general%upn,ewp,nsp) 
                   end if
                end do
                end do

             endif  ! thk > thklim, not floating

          ! include sliding contrib only if temperature node is surrounded by sliding velo nodes
          !TODO - This may result in non-conservation of energy.
          !       Why not include all nonzero terms? 

          if (slide_count >= 4) then
             slterm = 0.25d0 * slterm
          else
             slterm = 0.0d0
          end if

          model%temper%bfricflx(ew,ns) = slterm

       enddo    ! ns
       enddo    ! ew

  end subroutine glissade_calcbfric

  !-----------------------------------------------------------------------------------

!TODO - Some of these arguments are not needed:
!       stagthck, dusrfdew, dusrfdns, ubas, vbas

  subroutine glissade_calcbmlt( model,                   &
                                temp,     stagsigma,     &
                                thck,     stagthck,      &
                                dusrfdew, dusrfdns,      &
                                ubas,     vbas,          &
                                bmlt,     floater)

    ! Compute the amount of basal melting.
    ! The basal melting computed here is applied to the ice thickness
    !  by glissade_transport_driver, conserving mass and energy.
    !
    ! Any internal temperatures above the pressure melting point are reset to the
    !  pmp temperature, with excess energy applied toward basal melting.
    !  Hopefully this is rare.
    ! TODO: Moving all internal melting to the basal surface is not very realistic 
    !       and should be revisited.

    use glimmer_physcon, only: shci, rhoi, lhci

    type(glide_global_type) :: model

    real(dp), dimension(0:,:,:), intent(inout) :: temp
    real(dp), dimension(0:),     intent(in) :: stagsigma
    real(dp), dimension(:,:),    intent(in) :: thck,  stagthck, dusrfdew, dusrfdns, ubas, vbas  
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

            ! Add internal melting associated with temp > pmptemp
            ! Note: glissade_calcpmpt does not compute pmpt at the top surface or the bed.

             call glissade_calcpmpt(pmptemp(:), thck(ew,ns),   &
                                    stagsigma(:) )

             do up = 1, model%general%upn-1
                 if (temp(up,ew,ns) > pmptemp(up)) then
                    hmlt = (shci * thck(ew,ns) * (temp(up,ew,ns) - pmptemp(up))) / (rhoi * lhci) 
                    !BDM adding in what I think should be correct hmlt and bmlt for temp. based
                    !hmlt = (rhoi * shci * (model%numerics%sigma(up+1) - model%numerics%sigma(up))&
                    !       * (temp(up,ew,ns) - pmptemp(up))) / (rhow * lhci * thk0)
                    bmlt(ew,ns) = bmlt(ew,ns) + hmlt / model%numerics%dttem 
                    !bmlt(ew,ns) = bmlt(ew,ns) + hmlt * tim0 / (model%numerics%dttem)
                    temp(up,ew,ns) = pmptemp(up)
                 endif
             enddo

             ! Reset basal temp to pmptemp, if necessary

             up = model%general%upn
             call glissade_calcpmpt_bed(pmptemp(up), thck(ew,ns))
             temp(up,ew,ns) = min (temp(up,ew,ns), pmptemp(up))

          endif   ! thk > thklim

       enddo
    enddo

  end subroutine glissade_calcbmlt

!-------------------------------------------------------------------
 
!TODO - The HO version does not need dusrfdew, dusrfdns as inputs.
!       Since the HO version always computes temp on a staggered vertical grid, can remove some code.

  subroutine glissade_finddisp (model,     &
                                thck,      &
                                whichdisp, &
                                efvs,      &
                                stagthck,  &
                                dusrfdew,dusrfdns,  &
                                flwa)

    ! Compute the dissipation source term associated with strain heating.
    ! Note that the dissipation is computed in the same way on either a staggered or an
    !  unstaggered vertical grid.  
    ! Note also that dissip and flwa must have the same vertical dimension 
    !  (1:upn on an unstaggered vertical grid, or 1:upn-1 on a staggered vertical grid).
    
    use glimmer_physcon, only : gn

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(in) :: thck, stagthck, dusrfdew, dusrfdns
    real(dp), dimension(:,:,:), intent(in) :: flwa, efvs
    integer, intent(in) :: whichdisp

    integer, parameter :: p1 = gn + 1  
    integer :: ew, ns
    integer :: iew, ins    !*sfp* for HO and SSA dissip. calc.

    real(dp) :: c2 

    !*sfp* The next 2 declarations needed for HO and SSA dissip. calc. ... only needed 
    ! for internal work, so not clear if it is necessary to declare/allocate them elsewhere 
    real(dp) :: c4                         
    real(dp), dimension(model%general%upn) :: c5     
    
    select case( whichdisp ) 

!TODO - Do we need to support SIA dissipation for glissade code?
!       I am leaving this for now, because the standard dome config files do not specify whichdisp,
!        and therefore it uses the default SIA value (whichdisp = 0).
!       I think it would be better to remove the whichdisp option, and simply use the option
!        appropriate to the dycore.

    case( SIA_DISP )

    !*sfp* 0-order SIA case only 
    ! two methods of doing this. 
    ! 1. find dissipation at u-pts and then average
    ! 2. find dissipation at H-pts by averaging quantities from u-pts
    ! 2. works best for eismint divide (symmetry) but 1 likely to be better for full expts

    model%tempwk%dissip(:,:,:) = 0.0d0

    !LOOP TODO: Locally owned cells only
    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim) then
             
             c2 = (0.25*sum(stagthck(ew-1:ew,ns-1:ns)) * dsqrt((0.25*sum(dusrfdew(ew-1:ew,ns-1:ns)))**2 &
                  + (0.25*sum(dusrfdns(ew-1:ew,ns-1:ns)))**2))**p1
             
             model%tempwk%dissip(:,ew,ns) = c2 * model%tempwk%c1(:) * ( &
                  flwa(:,ew-1,ns-1) + flwa(:,ew-1,ns+1) + flwa(:,ew+1,ns+1) + flwa(:,ew+1,ns-1) + &
                  2*(flwa(:,ew-1,ns)+flwa(:,ew+1,ns)+flwa(:,ew,ns-1)+flwa(:,ew,ns+1)) + &
                  4*flwa(:,ew,ns)) 

          end if
       end do
    end do

    !the compensatory heating(compheat) is initialized to zero and allows
    !for modifying the calculated dissip.  This is needed for exact verification tests,
    !model%tempwk%dissip = model%tempwk%dissip + model%tempwk%compheat 

    case( FIRSTORDER_DISP )

    !*sfp* 1st-order, NON-depth integrated SIA case only (Pattyn, Payne-Price models) 
    ! NOTE: this needs tau and efvs (3d arrays), which are the eff. stress and the eff. visc. calculated
    ! from and/or consistent with the HO model. For simplicity, tau can be calculated from: tau = 2*efvs*eps_eff,
    ! where eps_eff is the eff. strain rate. Further, eps_eff can be calculated from the efvs according to a 
    ! re-arrangement of: efvs = 1/2 * ( 1 / A(T) )^(1/n) * eps_eff^((1-n)/n), in which case only the efvs and rate
    ! factor arrays need to be passed in for this calculation.

    model%tempwk%dissip(:,:,:) = 0.0d0

    if (size(model%tempwk%dissip,1) /= model%general%upn-1) then  ! staggered vertical grid
        !TODO - Write an error message and exit gracefully
    endif

    !LOOP TODO: Locally owned cells only

    do ns = 1, model%general%nsn
       do ew = 1, model%general%ewn
          if (thck(ew,ns) > model%numerics%thklim) then

             c5(:) = 0.0d0

             if ( sum( efvs(:,ew,ns) ) /= 0.0d0) then

                ! Use space in c5 vector to store dissip terms that apply at layer midpoints 
                ! (i.e. on staggered vertical grid).  No vertical averaging is needed, since
                ! temp and dissip are colocated with eff stress and eff viscosity.

                 c5(1:model%general%upn-1) = c5(1:model%general%upn-1)                      &
                                            +  model%stress%tau%scalar(:,ew,ns)**2 /  &
                                                efvs(1:model%general%upn-1,ew,ns)
             endif

             !Note: model%tempwk%cons(5) = (tau0*vel0/len0) / (rhoi*shci) * (model%numerics%dttem*tim0)

             model%tempwk%dissip(:,ew,ns) = c5(:) * model%tempwk%cons(5)

          endif
       enddo
    enddo

!   case( SSA_DISP )     !!! Waiting for an SSA solver !!!
!    !*sfp* 1st-order, depth-integrated case only (SSA model) 
!    ! NOTE: this needs taus and efvss (2d arrays), which are depth-integrated and averaged 
!    ! effective stress and effective viscosity fields calculated from and/or consistent
!    ! with the SSA model.
!
!    model%tempwk%dissip = 0.0d0
!    do ns = 2, model%general%nsn-1
!       do ew = 2, model%general%ewn-1
!          if (thck(ew,ns) > model%numerics%thklim) then
!            c4 = 0.0d0
!            do ins = ns-1,ns; do iew = ew-1,ew; 
!                if (efvss(iew,ins)  /=  0.0d0) then                     
!                    c4 = c4 + taus(iew,ins)**2 / efvss(iew,ins)
!                end if; 
!            end do; end do
!            model%tempwk%dissip(:,ew,ns) = c4 * model%tempwk%cons(5)
!          end if
!       end do
!    end do

    end select

  end subroutine glissade_finddisp

  !-----------------------------------------------------------------------------------
 
!TODO - Inline these subroutines above?

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

  subroutine glissade_calcbpmp(model,thck,bpmp)

    ! Calculate the pressure melting point at the base of the ice sheet

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(in)  :: thck
    real(dp), dimension(:,:), intent(out) :: bpmp

    integer :: ew,ns

    bpmp = 0.d0

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1
          !TODO - Inline this code?
          call glissade_calcpmpt_bed(bpmp(ew,ns),thck(ew,ns))
       end do
    end do

  end subroutine glissade_calcbpmp

  !-------------------------------------------------------------------

  subroutine glissade_calcpmpt_bed(pmptemp_bed, thck)

    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    real(dp), intent(out) :: pmptemp_bed ! pressure melting point temp at bed (deg C)
    real(dp), intent(in) :: thck         ! ice thickness

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp_bed = fact * thck 

  end subroutine glissade_calcpmpt_bed

  !-------------------------------------------------------------------

  subroutine glissade_calcflwa(stagsigma, thklim, flwa, temp, thck, flow_factor, &
                               default_flwa_arg, flag, waterfrac)

    !*FD Calculates Glen's $A$ over the three-dimensional domain,
    !*FD using one of three possible methods.
    !*FD
    !*FD The primary method is to use this equation from \emph{Paterson and Budd} [1982]:
    !*FD \[
    !*FD A(T^{*})=a \exp \left(\frac{-Q}{RT^{*}}\right)
    !*FD \]
    !*FD This is equation 9 in {\em Payne and Dongelmans}. $a$ is a constant of proportionality,
    !*FD $Q$ is the activation energy for for ice creep, and $R$ is the universal gas constant.
    !*FD The pressure-corrected temperature, $T^{*}$ is given by:
    !*FD \[
    !*FD T^{*}=T-T_{\mathrm{pmp}}+T_0
    !*FD \] 
    !*FD \[
    !*FD T_{\mathrm{pmp}}=T_0-\sigma \rho g H \Phi
    !*FD \]
    !*FD $T$ is the ice temperature, $T_{\mathrm{pmp}}$ is the pressure melting point 
    !*FD temperature, $T_0$ is the triple point of water, $\rho$ is the ice density, and 
    !*FD $\Phi$ is the (constant) rate of change of melting point temperature with pressure.

    use glimmer_physcon
    use glimmer_paramets, only : thk0, vis0

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

!   Note: The flwa, temp, and stagsigma arrays should have the same vertical dimension
!         (1:upn-1 on the staggered vertical grid).

    real(dp),dimension(:),      intent(in)    :: stagsigma ! vertical coordinate at layer midpoints
    real(dp),                   intent(in)    :: thklim    ! thickness threshold
    real(dp),dimension(:,:,:),  intent(in)    :: temp      ! 3D temperature field
    real(dp),dimension(:,:),    intent(in)    :: thck      ! ice thickness
    real(dp)                                  :: flow_factor ! fudge factor in Arrhenius relationship
    real(dp),                   intent(in)    :: default_flwa_arg ! Glen's A to use in isothermal case 
                                                                  ! Units: Pa^{-n} yr^{-1} 
    integer,                    intent(in)    :: flag      !*FD Flag to select the method
                                                           !*FD of calculation
    real(dp),dimension(:,:,:),  intent(out)   :: flwa      !*FD The calculated values of $A$
    real(dp),dimension(:,:,:),  intent(in), optional :: waterfrac!internal water content fraction, 0 to 1

    !*FD \begin{description}
    !*FD \item[0] {\em Paterson and Budd} relationship.
    !*FD \item[1] {\em Paterson and Budd} relationship, with temperature set to -5$^{\circ}$C.
    !*FD \item[2] Set to prescribed constant value.
    !*FD \end{description}

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp) :: default_flwa
    integer :: ew, ns, up, ewn, nsn, uflwa
    real(dp) :: tempcor

    real(dp), parameter :: fact = grav * rhoi * pmlt * thk0
    real(dp), parameter :: const_temp = -5.0d0

    real(dp),dimension(4), parameter ::  &
       arrfact = (/ arrmlh / vis0,      &   ! Value of A when T* is above -263K
                    arrmll / vis0,      &   ! Value of A when T* is below -263K
                   -actenh / gascon,    &   ! Value of -Q/R when T* is above -263K
                   -actenl / gascon/)       ! Value of -Q/R when T* is below -263K
    
    !------------------------------------------------------------------------------------ 
   
    uflwa=size(flwa,1) ; ewn=size(flwa,2) ; nsn=size(flwa,3)

    ! Check that the temperature array has the desired vertical dimension

    if (size(temp,1) /= size(flwa,1)) then
! debug
!       print*, 'upn =', upn
!       print*, 'size(temp,1) =', size(temp,1)
!       print*, 'size(flwa,1) =', size(flwa,1)
       call write_log('glissade_calcflwa: temp and flwa must have the same vertical dimensions', &
                       GM_FATAL)
    endif

    ! Scale the default rate factor (default value has units Pa^{-n} yr^{-1}).
    ! Also multiply by fudge factor

    default_flwa = flow_factor * default_flwa_arg / (vis0*scyr) 
    !write(*,*)"Default flwa = ",default_flwa

    select case(flag)

    case(FLWA_PATERSON_BUDD)

!LOOP TODO - Loop over locally owned cells?  Alternatively, just compute over all cells,
!        but make sure temp and thck are up to date in halo cells.
!       Might be cheaper to compute in all cells and avoid a halo call.

      ! This is the Paterson and Budd relationship
      ! BDM add waterfrac relationship for whichtemp=TEMP_ENTHALPY case

      do ns = 1,nsn
         do ew = 1,ewn

            if (thck(ew,ns) > thklim) then
            
               do up = 1, uflwa   ! uflwa = upn - 1 (values at layer midpoints)

                  ! Calculate the corrected temperature

                  tempcor = min(0.0d0, temp(up,ew,ns) + thck(ew,ns)*fact*stagsigma(up))
                  tempcor = max(-50.0d0, tempcor)

                  ! Calculate Glen's A (including flow fudge factor)

                  if (tempcor >= -10.d0) then
                     flwa(up,ew,ns) = flow_factor * arrfact(1) * exp(arrfact(3)/(tempcor + trpt))
                  else
                     flwa(up,ew,ns) = flow_factor * arrfact(2) * exp(arrfact(4)/(tempcor + trpt))
                  endif

                  ! BDM add correction for a liquid water fraction 
                  ! Using Greve and Blatter, 2009 formulation for Glen's A flow rate factor:
                  !    A = A(theta_PMP) * (1 + 181.25 * waterfrac)
                  if (present(waterfrac)) then
                     if (waterfrac(up,ew,ns) > 0.0d0) then
                        flwa(up,ew,ns) = flwa(up,ew,ns) * (1.d0 + 181.25d0 * waterfrac(up,ew,ns))      
                     endif
                  endif
               enddo

            else   ! thck < thklim

               flwa(:,ew,ns) = default_flwa

            end if

         end do
      end do

    case(FLWA_PATERSON_BUDD_CONST_TEMP)

      ! This is the Paterson and Budd relationship, but with the temperature held constant
      ! at -5 deg C

      do ns = 1,nsn
         do ew = 1,ewn

            if (thck(ew,ns) > thklim) then

               ! Calculate Glen's A with a fixed temperature (including flow fudge factor)

               if (const_temp >= -10.d0) then
                  flwa(:,ew,ns) = flow_factor * arrfact(1) * exp(arrfact(3)/(const_temp + trpt))
               else
                  flwa(:,ew,ns) = flow_factor * arrfact(2) * exp(arrfact(4)/(const_temp + trpt))
               endif

            else

               flwa(:,ew,ns) = default_flwa

            end if

         end do
      end do

    case(FLWA_CONST_FLWA) 

      flwa(:,:,:) = default_flwa
  
    end select

  end subroutine glissade_calcflwa 

!------------------------------------------------------------------------------------

end module glissade_temp

!------------------------------------------------------------------------------------
