!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_temp.F90 - part of the Community Ice Sheet Model (CISM)  
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

  ! This module is based on glide_temp.F90, but has been modified for Glissade
  !  by William Lipscomb (LANL).
  ! It computes temperature diffusion and strain heating in a local column 
  !  without doing horizontal or vertical advection.
  ! Temperature advection is done separately, e.g. using the incremental 
  !  remapping transport scheme.
  ! It is assumed here that temperature values are staggered in the
  !  vertical compared to the velocity.  That is, the temperature lives
  !  at the midpoint of each layer instead of at layer interfaces.
  ! The temperature is also defined at the upper and lower surfaces with 
  !  appropriate boundary conditions.   

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade_temp

    use glimmer_global, only : dp 
    use glide_types
    use glimmer_log
    use parallel, only: this_rank

    implicit none

    private
    public :: glissade_init_temp, glissade_temp_driver, glissade_calcflwa, glissade_calcbpmp

    ! time stepping scheme

    !NOTE:  For the dome test case, the Crank-Nicolson scheme can give unstable 
    !        temperature fluctuations for thin ice immediately after the ice 
    !        becomes thick enough for the temperature calculation.
    !       The fully implicit scheme has been stable for all cases (but is only
    !        first-order accurate in time). 

    logical, parameter::   &
         crank_nicolson = .false.  ! if true, use Crank-Nicolson time-stepping
                                   ! if false, use fully implicit

    ! max and min allowed temperatures
    ! Temperatures sometimes go below -100 for cases where Crank-Nicholson is unstable
    real(dp), parameter ::   &
       maxtemp_threshold = 1.d11,   &
       mintemp_threshold = -100.d0

    !WHL - debug
    integer :: itest, jtest, rtest

contains

!****************************************************    

  subroutine glissade_init_temp (model)

    ! initialization subroutine for the case that temperature lives on the
    ! vertically staggered grid (i.e., at layer centers)

    use glimmer_physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow, trpt
    use glimmer_paramets, only : tim0, thk0, len0, vis0, vel0, tau0, unphys_val
    use parallel, only: lhalo, uhalo
    use glissade_enthalpy, only: glissade_init_enthalpy

    type(glide_global_type),intent(inout) :: model       !> Ice model parameters.

    integer, parameter :: p1 = gn + 1  
    integer up, ns, ew

    !TODO - Should these allocations be done in glide_allocarr?

    ! Note vertical dimensions here.  Dissipation is computed for each of (upn-1) layers.
    ! Temperature is defined at midpoint of each layer, plus upper and lower surfaces.
    !TODO - Allocate dissip in glide_types?
    allocate(model%tempwk%dups(model%general%upn+1,2))
    allocate(model%tempwk%inittemp(model%general%upn+1,model%general%ewn,model%general%nsn))
    !WHL - Moved dissip to model%temper and allocated in glide_types
!!    allocate(model%tempwk%dissip  (model%general%upn-1,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%compheat(model%general%upn-1,model%general%ewn,model%general%nsn))
    model%tempwk%compheat = 0.0d0

    allocate(model%tempwk%c1(model%general%upn-1))   ! upn-1 for staggered grid

    allocate(model%tempwk%dupa(model%general%upn))
    allocate(model%tempwk%dupb(model%general%upn))
    allocate(model%tempwk%dupc(model%general%upn))

    model%tempwk%dups = 0.0d0

    !Note: The 'dups' grid coefficients are not the same as for unstaggered temperatures.
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

!WHL - The factor of 2 in the numerator in the original code can be traced to a missing factor of 0.5
!      in the denominator of the dups coefficients.  On the vertically staggered grid, there is no
!      factor of 0.5 in the dups coefficients, so there is no factor of 2 here.
!WHL - The factor of 2 in the denominator is a Crank-Nicolson averaging factor.
!      If doing fully implicit timestepping, model%tempwk%cons(1) is multiplied by 2 below.
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
    !TODO - When reading temperature from restart or input file, make sure that halo values are correct.

    if (model%options%is_restart == RESTART_TRUE) then

       ! Temperature has already been initialized from a restart file. 
       ! (Temperature is always a restart variable.)

       call write_log('Initializing ice temperature from the restart file')

    elseif ( minval(model%temper%temp(1:model%general%upn, &
                    1+lhalo:model%general%ewn-lhalo, 1+uhalo:model%general%nsn-uhalo)) > &
                    (-1.0d0 * trpt) ) then    ! trpt = 273.15 K
                                              ! Default initial temps in glide_types are unphys_val = -999

       !TODO - Verify vertical dimension here (should be 0 or upn?)
       if ( (maxval(model%temper%temp(model%general%upn, &
                    1+lhalo:model%general%ewn-lhalo, 1+uhalo:model%general%nsn-uhalo)) == unphys_val) .and.    &
                    (model%options%whichdycore /= DYCORE_BISICLES) )then
          ! Throw a fatal error if we think the user has supplied temp instead of tempstag
          ! (We don't want to implicitly shift the vertical layers from one coordinate system
          !  to another without the user knowing.)
          ! This case will look like good data in all the layers except the top layer.
          ! MJH: Letting BISICLES run with this situation as per Dan Martin's request.
          call write_log("The variable 'temp' has been read from an input file, but it only is appropriate " &
             // "for the Glide dycore.  Use the 'tempstag' variable with higher-order dycores instead.", GM_FATAL)
       else
          ! Temperature has already been initialized from an input file.
          ! (We know this because the default initial temps of unphys_val = -999 have been overwritten.)

          call write_log('Initializing ice temperature from an input file')
       endif

    else   ! not reading temperature from restart or input file
           ! initialize it here based on model%options%temp_init

       ! First set T = 0 C everywhere

       model%temper%temp(:,:,:) = 0.0d0                                              
                                                    
       if (model%options%temp_init == TEMP_INIT_ZERO) then

          call write_log('Initializing ice temperature to 0 deg C')

          ! No call is needed to glissade_init_temp_column because the
          ! ice temperature has been set to zero above

       elseif (model%options%temp_init == TEMP_INIT_ARTM) then

          ! Initialize ice column temperature to min(artm, 0 C).

          !Note: Old glide sets temp = artm everywhere without regard to whether ice exists in a column.

          call write_log('Initializing ice temperature to the surface air temperature')

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

          call write_log('Initializing ice temperature to a linear profile in each column')

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
  ! Three possibilities:
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

  !TODO - Define pmpt_offset elsewhere?
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
        temp(:) = min(0.0d0, artm)
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

    ! Calculates the new ice temperature 

    use glimmer_utils,  only : tridiag
    use glimmer_paramets, only : thk0, tim0
    use glimmer_physcon, only: shci, coni, rhoi
    use glide_mask
    use glissade_enthalpy
    use glissade_grid_operators, only: glissade_stagger
    use glissade_masks, only: glissade_get_masks

    !TODO - Modify glissade_temp_driver to compute over locally owned cells only?
    !       This would make the module a bit cheaper but would require halo updates at the end.

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_global_type),intent(inout) :: model       ! Ice model parameters
    integer,                intent(in)    :: whichtemp   ! Flag to choose method

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
    real(dp) :: einit, efinal, delta_e, dTtop, dTbot, denth_top, denth_bot

    real(dp) :: maxtemp, mintemp   ! max and min temps in column

    real(dp), dimension(0:model%general%upn) :: pmptemp   ! pressure melting pt temperature

    real(dp), dimension(1:model%general%upn) :: alpha_enth   ! diffusivity at interfaces (m2/s) for enthalpy solver
                                                             ! = coni / (rhoi*shci) for cold ice

    integer, dimension(model%general%ewn,model%general%nsn) ::  &
         ice_mask   ! = 1 where thck > thklim_temp, else = 0

!!    logical, parameter:: verbose_temp = .false.
    logical, parameter:: verbose_temp = .true.
    integer :: k

    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    upn = model%general%upn

    select case(whichtemp)

    case(TEMP_SURFACE_AIR_TEMP)  ! Set column to surface air temperature ------------------

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
       
       !TODO - Change Tstagsigma to stagwbndsigma
       !       Change model%general%upn to upn
       Tstagsigma(0) = 0.d0
       Tstagsigma(1:model%general%upn-1) = model%numerics%stagsigma(1:model%general%upn-1)
       Tstagsigma(model%general%upn) = 1.d0

       model%tempwk%inittemp = 0.0d0

       ! Calculate interior heat dissipation -------------------------------------

       call glissade_finddisp( model,                      &
                               model%geometry%thck,        &
                               model%options%which_ho_disp,&
                               model%stress%efvs,          &
                               model%geomderv%stagthck,    &
                               model%geomderv%dusrfdew,    &
                               model%geomderv%dusrfdns,    &
                               model%temper%flwa)

       ! Calculate heating from basal friction (if not already computed by the Glissade velocity solver)

       call glissade_calcbfric( model,                        &
                                model%options%whichdycore,    &
                                model%geometry%thck,          &
                                model%velocity%btraction,     &
                                model%velocity%ubas,          &
                                model%velocity%vbas,          &
                                GLIDE_IS_FLOAT(model%geometry%thkmask), &
                                model%temper%bfricflx )

       ! Note: No iteration is needed here since we are doing a local tridiagonal solve without advection.

       do ns = 2,model%general%nsn-1
       do ew = 2,model%general%ewn-1

          if(model%geometry%thck(ew,ns) > model%numerics%thklim_temp) then

             if (verbose_temp .and. this_rank==rtest .and. ew==itest .and. ns==jtest) then
                print*, ' '
                print*, 'Before prognostic temp, i, j =', ew, ns
                print*, 'thck =', model%geometry%thck(ew,ns)*thk0
                print*, 'Temp:'
                do k = 0, upn
                   print*, k, model%temper%temp(k,ew,ns)
                enddo
             endif

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
                
             if (verbose_temp .and. this_rank==rtest .and. ew==itest .and. ns==jtest) then
                print*, 'After glissade_findvtri, i, j =', ew,ns
                print*, 'k, subd, diag, supd, rhsd:'
                do k = 1, upn+1
                   print*, k, subd(k), diag(k), supd(k), rhsd(k)
                enddo
             endif

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

             if (verbose_temp .and. this_rank==rtest .and. ew==itest .and. ns==jtest) then
                print*, ' '
                print*, 'After prognostic temp, i, j =', ew, ns
                print*, 'Temp:'
                do k = 0, upn
                   print*, k, model%temper%temp(k,ew,ns)
                enddo
             endif

             ! Check that the net input of energy to the column is equal to the difference
             !  between the initial and final internal energy.
             !TODO - Make this energy check optional and/or move it to a subroutine?

             ! compute the final internal energy

             efinal = 0.0d0
             do up = 1, upn-1
                efinal = efinal + model%temper%temp(up,ew,ns) *  &
                                 (model%numerics%sigma(up+1) - model%numerics%sigma(up))
             enddo
             efinal = efinal * rhoi*shci * model%geometry%thck(ew,ns)*thk0

             ! compute net heat flux to the column

             ! conductive flux = (k/H * dT/dsigma) at upper and lower surfaces; positive down

             if (crank_nicolson) then
                ! average temperatures between start and end of timestep
                dTtop = 0.5d0 * ( model%temper%temp(1,ew,ns) - model%temper%temp(0,ew,ns) &
                                +     prevtemp_stag(1)       -     prevtemp_stag(0) )
                dTbot = 0.5d0 * ( model%temper%temp(upn,ew,ns) - model%temper%temp(upn-1,ew,ns) &
                                +     prevtemp_stag(upn)       -   prevtemp_stag(upn-1) )
             else    ! fully implicit
                ! use temperatures at end of timestep
                dTtop = model%temper%temp(1,ew,ns) - model%temper%temp(0,ew,ns)
                dTbot = model%temper%temp(upn,ew,ns) - model%temper%temp(upn-1,ew,ns)
             endif

             model%temper%ucondflx(ew,ns) = (-coni / (model%geometry%thck(ew,ns)*thk0) )         &
                                           * dTtop / (Tstagsigma(1) - Tstagsigma(0))

             model%temper%lcondflx(ew,ns) = (-coni / (model%geometry%thck(ew,ns)*thk0) )         &
                                           * dTbot / (Tstagsigma(upn) - Tstagsigma(upn-1))

             ! total dissipation in column (W/m^2)

             model%temper%dissipcol(ew,ns) = 0.0d0
             do up = 1, upn-1
                model%temper%dissipcol(ew,ns) = model%temper%dissipcol(ew,ns) + &
                                              model%temper%dissip(up,ew,ns)  &
                                           * (model%numerics%sigma(up+1) - model%numerics%sigma(up))  
             enddo 
             model%temper%dissipcol(ew,ns) = model%temper%dissipcol(ew, ns)     &
                                     * thk0*model%geometry%thck(ew,ns)*rhoi*shci / (tim0*model%numerics%dttem)  

             ! Verify that the net input of energy into the column is equal to the change in
             ! internal energy.  

             delta_e = (model%temper%ucondflx(ew,ns) - model%temper%lcondflx(ew,ns)  &
                      + model%temper%dissipcol(ew,ns)) * tim0*model%numerics%dttem

             if ( abs((efinal-einit-delta_e)/(tim0*model%numerics%dttem)) > 1.0d-8 ) then

                if (verbose_temp) then
                   print*, 'Ice thickness:', thk0*model%geometry%thck(ew,ns)
                   print*, 'thklim_temp:', thk0*model%numerics%thklim_temp
                   print*, ' '
                   print*, 'Interior fluxes:'
                   print*, 'ftop (pos up)=', -model%temper%ucondflx(ew,ns) 
                   print*, 'fbot (pos up)=', -model%temper%lcondflx(ew,ns)
                   print*, 'fdissip =',       model%temper%dissipcol(ew,ns)
                   print*, 'Net flux =', delta_e/(tim0*model%numerics%dttem)
                   print*, ' '
                   print*, 'delta_e =', delta_e
                   print*, 'einit =',  einit
                   print*, 'efinal =', efinal
                   print*, 'einit + delta_e =', einit + delta_e
                   print*, ' '
                   print*, 'Energy imbalance =', efinal - einit - delta_e
                   print*, ' '
                   print*, 'Basal fluxes:'
                   print*, 'ffric =', model%temper%bfricflx(ew,ns)
                   print*, 'fgeo =', -model%temper%bheatflx(ew,ns)
                   print*, 'flux for bottom melting =', model%temper%bfricflx(ew,ns)   &
                                                          - model%temper%bheatflx(ew,ns)   &
                                                          + model%temper%lcondflx(ew,ns)
                endif   ! verbose_temp

                write(message,*) 'WARNING: Energy conservation error, ew, ns =', ew, ns
                call write_log(message,GM_FATAL)
             endif

             !WHL - No call here to corrpmpt.  Temperatures above pmpt are set to pmpt 
             !      in glissade_calcbmlt (conserving energy).

          endif  ! thck > thklim_temp
       end do    ! ew
       end do    ! ns

       ! set temperature of thin ice to the air temperature and set ice-free nodes to zero

       do ns = 1, model%general%nsn
          do ew = 1, model%general%ewn

!             if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
!                model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
!             else if (model%geometry%thkmask(ew,ns) < 0) then
!                model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
!             !else if (model%geometry%thkmask(ew,ns) < -1) then
!             !   model%temper%temp(:,ew,ns) = 0.0d0
!             end if

             !WHL - Changed threshold from thklim to thklim_temp
              if (model%geometry%thck(ew,ns) <= model%numerics%thklim_temp) then
                 model%temper%temp(:,ew,ns) = min(0.d0, model%climate%artm(ew,ns))
              endif

              !TODO - Maybe it should be done in the following way, so that the temperature profile for thin ice
              !       is consistent with the temp_init option, with T = 0 for ice-free cells.

             ! NOTE: Calling this subroutine will maintain a sensible temperature profile
             !        for thin ice, but in general does *not* conserve energy.
             !       To conserve energy, we need either thklim_temp = 0, or some additional
             !        energy accounting and correction.
 
!             if (model%geometry%thck(ew,ns) <= model%numerics%thklim_temp) then
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

       call glissade_calcbmlt( model,                     &
                               model%temper%temp,         &
                               Tstagsigma,                &
                               model%geometry%thck,       &
                               model%temper%bmlt_ground,  &
                               GLIDE_IS_FLOAT(model%geometry%thkmask))

       ! Interpolate basal temperature and pressure melting point onto velocity grid

       call glissade_get_masks(model%general%ewn,    model%general%nsn,           &
                               model%geometry%thck,  model%geometry%topg,         &
                               model%climate%eus,    model%numerics%thklim_temp,  &
                               ice_mask)

       call glissade_stagger(model%general%ewn,     model%general%nsn,  &
                             model%temper%temp(model%general%upn,:,:),  &
                             model%temper%stagbtemp(:,:),               &
                             ice_mask(:,:),                             &
                             stagger_margin_in = 1)

       call glissade_stagger(model%general%ewn,     model%general%nsn,  &
                             model%temper%bpmp(:,:),                    &
                             model%temper%stagbpmp(:,:),                &
                             ice_mask(:,:),                             &
                             stagger_margin_in = 1)

   case(TEMP_STEADY)! do nothing

   case(TEMP_ENTHALPY)! BDM Local column calculation (with advection done elsewhere)

      !WHL - debug
      print*, 'Starting enthalpy calculation'

      !TODO - Modify enthalpy code to use a backward Euler timestep.
      !       With a Crank-Nicolson timestep, the ice in thin cells can become excessively cold, and the code aborts.
      !       To allow the code to keep running, set thklim_temp to a larger value (e.g., 100 m)
      !TODO - Rearrange prognostic/enthalpy code to avoid duplication of many calls.

      ! No horizontal or vertical advection; vertical diffusion and strain heating only.
      ! Enthalpy is vertically staggered relative to velocities.  
      ! That is, enthalpy is defined at the midpoint of each layer 
      ! (and at the top and bottom surfaces).

      ! BDM Enthalpy Gradient Method is used here to solve for temp. and water content.

      !TODO - If I'm using TEMP_ENTHALPY, should I make a call to glissade_init_enthalpy here?

      !TODO - Change Tstagsigma to stagwbndsigma
      Tstagsigma(0) = 0.d0
      Tstagsigma(1:model%general%upn-1) = model%numerics%stagsigma(1:model%general%upn-1)
      Tstagsigma(model%general%upn) = 1.d0

      !WHL - Commenting out the next two calls because dissip and bfricflx are now
      !      computed at the end of the previous time step

      ! Calculate interior heat dissipation -------------------------------------

!      call glissade_finddisp(  model,                      &
!                               model%geometry%thck,        &
!                               model%options%which_ho_disp,&
!                               model%stress%efvs,          &
!                               model%geomderv%stagthck,    &
!                               model%geomderv%dusrfdew,    &
!                               model%geomderv%dusrfdns,    &
!                               model%temper%flwa)

      ! Calculate heating from basal friction -----------------------------------

!      call glissade_calcbfric( model,                        &
!                               model%options%whichdycore,    &
!                               model%geometry%thck,          &
!                               model%velocity%btraction,     &
!                               model%velocity%ubas,          &
!                               model%velocity%vbas,          &
!                               GLIDE_IS_FLOAT(model%geometry%thkmask), &
!                               model%temper%bfricflx )

      ! Note: No iteration is needed here since we are doing a local tridiagonal solve without advection.

      do ns = 2,model%general%nsn-1
         do ew = 2,model%general%ewn-1
            if (model%geometry%thck(ew,ns) > model%numerics%thklim_temp) then

               ! Convert model%temper%temp and model%temper%waterfrac to enthalpy (dimension 0:upn).
               ! For interior and boundary nodes; assume waterfrac = 0 at boundaries.
               ! BDM enthalpy will be size 0:upn                                                                                                                                
               call temp2enth(model%temper%enthalpy(0:model%general%upn,ew,ns),    &
                              model%temper%temp(0:model%general%upn,ew,ns),        &
                              model%temper%waterfrac(1:model%general%upn-1,ew,ns), &
                              model%geometry%thck(ew,ns),                          &
                              model%numerics%stagsigma(1:model%general%upn-1))

               if (verbose_temp .and. this_rank==rtest .and. ew==itest .and. ns==jtest) then
                  print*, ' '
                  print*, 'Before prognostic enthalpy, i, j =', ew, ns
                  print*, 'thck =', model%geometry%thck(ew,ns)*thk0
                  print*, 'Temp, waterfrac, enthalpy:'
                  k = 0
                  print*, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)
                  do k = 1, upn-1
                     print*, k, model%temper%temp(k,ew,ns), model%temper%waterfrac(k,ew,ns), model%temper%enthalpy(k,ew,ns)
                  enddo
                  k = upn
                  print*, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)
               endif

               ! compute initial internal energy in column (for energy conservation check)
               einit = 0.0d0
               do up = 1, upn-1
                  einit = einit + model%temper%enthalpy(up,ew,ns) *  &
                                  (model%numerics%sigma(up+1) - model%numerics%sigma(up) )
               enddo
               einit = einit * model%geometry%thck(ew,ns)*thk0

               ! BDM compute matrix elements using Enthalpy Gradient Method

               call glissade_enthalpy_findvtri(model, ew,   ns,         &
                                               subd,  diag, supd, rhsd, &
                                               GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)),  &
                                               alpha_enth)
					
               ! BDM leave as prevtemp because it's only used for dT/dsigma at top and bottom boundaries,
               ! we don't want this as enthalpy
               prevtemp_stag(:) = model%temper%temp(:,ew,ns)

               ! solve the tridiagonal system
               ! Note: Temperature is indexed from 0 to upn, with indices 1 to upn-1 colocated
               ! with stagsigma values of the same index.
               ! However, the matrix elements are indexed 1 to upn+1, with the first row
               ! corresponding to the surface temperature, temp(0,:,:).

               !WHL - debug
               if (ew==itest .and. ns==jtest) then
                  print*, ' '
                  print*, 'After vtri, i, j =', ew, ns
                  print*, 'k, subd, diag, supd, rhs/(rhoi*ci):'
                  do k = 1, upn+1
                     print*, k-1, subd(k), diag(k), supd(k), rhsd(k)/(rhoi*shci)
                  enddo
               endif

               call tridiag(subd(1:upn+1),   &
                            diag(1:upn+1),   &
                            supd(1:upn+1),   &
                            enthalpy(0:upn), &
                            rhsd(1:upn+1))
               
               ! Copy the local enthalpy array into the global derived type
               model%temper%enthalpy(:,ew,ns) = enthalpy(:)

               ! BDM convert back to temperature and water content
               call enth2temp(model%temper%enthalpy(0:upn,ew,ns),    &
                              model%temper%temp(0:upn,ew,ns),        &
                              model%temper%waterfrac(1:upn-1,ew,ns), &
                              model%geometry%thck(ew,ns),            &
                              model%numerics%stagsigma(1:upn-1))

               if (verbose_temp .and. this_rank==rtest .and. ew==itest .and. ns==jtest) then
                  print*, ' '
                  print*, 'After prognostic enthalpy, i, j =', ew, ns
                  print*, 'thck =', model%geometry%thck(ew,ns)*thk0
                  print*, 'Temp, waterfrac, enthalpy:'
                  k = 0
                  print*, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)
                  do k = 1, upn-1
                     print*, k, model%temper%temp(k,ew,ns), model%temper%waterfrac(k,ew,ns), model%temper%enthalpy(k,ew,ns)
                  enddo
                  k = upn
                  print*, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)
               endif

               ! Check that the net input of energy to the column is equal to the difference
               !  between the initial and final internal energy.

               ! compute the final internal energy

               efinal = 0.0d0
               do up = 1, upn-1
                  efinal = efinal + model%temper%enthalpy(up,ew,ns) *  &
                                   (model%numerics%sigma(up+1) - model%numerics%sigma(up) )
               enddo
               efinal = efinal * model%geometry%thck(ew,ns)*thk0

               ! compute net heat flux to the column

               ! conductive flux = (alpha/H * denth/dsigma) at upper and lower surfaces; positive down.
               ! Here alpha = coni / (rhoi*shci) for cold ice, with a smaller value for temperate ice.
               ! Assume fully implicit backward Euler time step

               denth_top = enthalpy(1) - enthalpy(0)
               denth_bot = enthalpy(upn) - enthalpy(upn-1)

               model%temper%ucondflx(ew,ns) = -alpha_enth(1) / (model%geometry%thck(ew,ns)*thk0)         &
                                             * denth_top / (Tstagsigma(1) - Tstagsigma(0))

               model%temper%lcondflx(ew,ns) = -alpha_enth(upn) / (model%geometry%thck(ew,ns)*thk0)         &
                                             * denth_bot / (Tstagsigma(upn) - Tstagsigma(upn-1))

               !TODO - From here on, the energy conservation check is the same as for temperature.
               ! total dissipation in column (W/m^2)

               model%temper%dissipcol(ew,ns) = 0.0d0
               do up = 1, upn-1
                  model%temper%dissipcol(ew,ns) = model%temper%dissipcol(ew,ns) + &
                                                  model%temper%dissip(up,ew,ns)  &
                                               * (model%numerics%sigma(up+1) - model%numerics%sigma(up))  
               enddo
               model%temper%dissipcol(ew,ns) = model%temper%dissipcol(ew,ns)     &
                                             * thk0*model%geometry%thck(ew,ns)*rhoi*shci / (tim0*model%numerics%dttem)  

               ! Verify that the net input of energy into the column is equal to the change in internal energy.  

               delta_e = (model%temper%ucondflx(ew,ns) - model%temper%lcondflx(ew,ns)  &
                        + model%temper%dissipcol(ew,ns)) * tim0*model%numerics%dttem

               if ( abs((efinal-einit-delta_e)/(tim0*model%numerics%dttem)) > 1.0d-8 ) then

                  if (verbose_temp) then
                     print*, 'Ice thickness:', thk0*model%geometry%thck(ew,ns)
                     print*, 'thklim_temp:', thk0*model%numerics%thklim_temp
                     print*, ' '
                     print*, 'Interior fluxes:'
                     print*, 'ftop (pos up)=', -model%temper%ucondflx(ew,ns) 
                     print*, 'fbot (pos up)=', -model%temper%lcondflx(ew,ns)
                     print*, 'fdissip =',       model%temper%dissipcol(ew,ns)
                     print*, 'Net flux =', delta_e/(tim0*model%numerics%dttem)
                     print*, ' '
                     print*, 'delta_e =', delta_e
                     print*, 'einit =',  einit
                     print*, 'efinal =', efinal
                     print*, 'einit + delta_e =', einit + delta_e
                     print*, ' '
                     print*, 'Energy imbalance =', efinal - einit - delta_e
                     print*, ' '
                     print*, 'Basal fluxes:'
                     print*, 'ffric =', model%temper%bfricflx(ew,ns)
                     print*, 'fgeo =', -model%temper%bheatflx(ew,ns)
                     print*, 'flux for bottom melting =', model%temper%bfricflx(ew,ns)   &
                                                        - model%temper%bheatflx(ew,ns)   &
                                                        + model%temper%lcondflx(ew,ns)
                  endif   ! verbose_temp
                  
                  write(message,*) 'WARNING: Energy conservation error, ew, ns =', ew, ns
                  call write_log(message,GM_FATAL)
               endif
               
               !WHL - No call here to corrpmpt.  Temperatures above pmpt are set to pmpt 
               !      in glissade_calcbmlt (conserving energy).

               !WHL - debug
               if (ew==itest .and. ns==jtest) then
                  print*, 'k/(rho*c) =', coni/(rhoi*shci)
                  print*, 'alpha_enth(upn) =', alpha_enth(upn)
                  print*, ' '
                  print*, 'After enthalpy calc, i, j =', ew, ns
                  print*, 'k, temp, wfrac, enthalpy/(rhoi*ci):'
                  k = 0
                  print*, k, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)/(rhoi*shci)
                  do k = 1, model%general%upn-1
                     print*, k, model%temper%temp(k,ew,ns), model%temper%waterfrac(k,ew,ns), &
                                model%temper%enthalpy(k,ew,ns)/(rhoi*shci)
                  enddo
                  k = model%general%upn
                  print*, k, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)/(rhoi*shci)
                  print*, ' '
                  print*, 'bheatflx, bfricflx, lcondflx, sum:', &
                      -model%temper%bheatflx(ew,ns), model%temper%bfricflx(ew,ns), model%temper%lcondflx(ew,ns), &
                      -model%temper%bheatflx(ew,ns)+ model%temper%bfricflx(ew,ns)+ model%temper%lcondflx(ew,ns)
               endif

            endif  ! thck > thklim_temp
         end do    ! ew
      end do    ! ns

      ! set temperature of thin ice to the air temperature and set ice-free nodes to zero

      do ns = 1, model%general%nsn
         do ew = 1, model%general%ewn

!            if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
!               model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
!            else if (model%geometry%thkmask(ew,ns) < 0) then
!               model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
!            !else if (model%geometry%thkmask(ew,ns) < -1) then
!            !   model%temper%temp(:,ew,ns) = 0.0d0
!            end if

            !WHL - Changed threshold to thklim_temp
            if (model%geometry%thck(ew,ns) <= model%numerics%thklim_temp) then

               !WHL - Make sure T <= T_pmp
               !TODO - Impose this condition for standard temperature calculation too?
               pmptemp(0) = 0.0d0
               call glissade_calcpmpt(pmptemp(1:upn-1), model%geometry%thck(ew,ns), &
                                                        model%numerics%stagsigma(1:upn-1))
               call glissade_calcpmpt_bed(pmptemp(upn), model%geometry%thck(ew,ns))
               model%temper%temp(:,ew,ns) = min(pmptemp(:), dble(model%climate%artm(ew,ns)))
            endif

            !NOTE - See comments above about setting temperature in thin ice

         end do !ew
      end do !ns

      ! BDM since there will be no temps above PMP, need new subroutine to calculate basal melt
      !     and basal water depth
      call glissade_enthalpy_calcbmlt(model,                     &
                                      model%temper%temp,         &
                                      model%temper%waterfrac,    &
                                      Tstagsigma,                &
                                      model%geometry%thck,       &
                                      model%temper%bmlt_ground,  &
                                      GLIDE_IS_FLOAT(model%geometry%thkmask))


      !WHL - debug
      ew = itest
      ns = jtest
      k = model%general%upn
      print*, ' '
      print*, 'After calcbmlt, i, j, basal temp =', ew, ns, model%temper%temp(k,ew,ns)

      ! Interpolate basal temperature and pressure melting point onto velocity grid
       !WHL - Replaced calls to stagvarb (an old Glide routine) with calls to glissade_stagger.
       !       stagger_margin_in = 1 implies that ice-free cells (where the basal temperature has
       !       no physical meaning) are not included in the average.
       !      With stagvarb, values in ice-free cells are (erroneously) included.

       call glissade_get_masks(model%general%ewn,    model%general%nsn,           &
                               model%geometry%thck,  model%geometry%topg,         &
                               model%climate%eus,    model%numerics%thklim_temp,  &
                               ice_mask)

       call glissade_stagger(model%general%ewn,     model%general%nsn,  &
                             model%temper%temp(model%general%upn,:,:),  &
                             model%temper%stagbtemp(:,:),               &
                             ice_mask(:,:),                             &
                             stagger_margin_in = 1)

       call glissade_stagger(model%general%ewn,     model%general%nsn,  &
                             model%temper%bpmp(:,:),                    &
                             model%temper%stagbpmp(:,:),                &
                             ice_mask(:,:),                             &
                             stagger_margin_in = 1)

    end select   ! whichtemp

    ! Check for temperatures that are physically unrealistic.
    ! Thresholds are set at the top of this module.

    do ns = 1, model%general%nsn
       do ew = 1, model%general%ewn

          maxtemp = maxval(model%temper%temp(:,ew,ns))
          mintemp = minval(model%temper%temp(:,ew,ns))
          
          if (maxtemp > maxtemp_threshold) then
             write(message,*) 'maxtemp > 0: i, j, maxtemp =', ew, ns, maxtemp
             call write_log(message,GM_FATAL)
          endif
          
          if (mintemp < mintemp_threshold) then
             !uncommment these line to get more info
!             print*, 'thck =', model%geometry%thck(ew,ns) * thk0
!             print*, 'temp:'
!             do k = 1, model%general%upn
!                print*, k, model%temper%temp(k,ew,ns)
!             enddo
             write(message,*) 'mintemp < mintemp_threshold: i, j, mintemp =', ew, ns, mintemp
             call write_log(message,GM_FATAL)
          endif
          
       enddo
    enddo

    ! Rescale dissipation term to deg C/s (instead of deg C)
    !WHL - Treat dissip above as a rate (deg C/s) instead of deg C
    model%temper%dissip(:,:,:) =  model%temper%dissip(:,:,:) /  (model%numerics%dttem*tim0)

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

    ! model%tempwk%cons(1) = tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)

    if (crank_nicolson) then

       fact = model%tempwk%cons(1) / model%geometry%thck(ew,ns)**2
       subd(2:model%general%upn) = -fact * model%tempwk%dups(1:model%general%upn-1,1)
       supd(2:model%general%upn) = -fact * model%tempwk%dups(1:model%general%upn-1,2)
       diag(2:model%general%upn) = 1.0d0 - subd(2:model%general%upn)     &
                                         - supd(2:model%general%upn)

       model%tempwk%inittemp(1:model%general%upn-1,ew,ns) =   &
                model%temper%temp(1:model%general%upn-1,ew,ns) * (2.0d0 - diag(2:model%general%upn)) &
              - model%temper%temp(0:model%general%upn-2,ew,ns) * subd(2:model%general%upn) &
              - model%temper%temp(2:model%general%upn,  ew,ns) * supd(2:model%general%upn) & 
              + model%temper%dissip(1:model%general%upn-1,ew,ns)
    
       rhsd(2:model%general%upn) = model%tempwk%inittemp(1:model%general%upn-1,ew,ns)

    else   ! fully implicit

       fact = 2.d0 * model%tempwk%cons(1) / model%geometry%thck(ew,ns)**2  ! Remove factor of 2 in denominator
       subd(2:model%general%upn) = -fact * model%tempwk%dups(1:model%general%upn-1,1)
       supd(2:model%general%upn) = -fact * model%tempwk%dups(1:model%general%upn-1,2)
       diag(2:model%general%upn) = 1.0d0 - subd(2:model%general%upn)     &
                                         - supd(2:model%general%upn)
       
       model%tempwk%inittemp(1:model%general%upn-1,ew,ns) =   &
                model%temper%temp(1:model%general%upn-1,ew,ns)  &
              + model%temper%dissip(1:model%general%upn-1,ew,ns)
    
       rhsd(2:model%general%upn) = model%tempwk%inittemp(1:model%general%upn-1,ew,ns)

    endif    ! crank_nicolson

    ! basal boundary:
    ! for grounded ice, a heat flux is applied
    ! for floating ice, the basal temperature is held constant

    !NOTE: This lower BC is different from the one in glide_temp.
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
             (model%temper%bfricflx(ew,ns) - model%temper%bheatflx(ew,ns)) &
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

  subroutine glissade_calcbfric (model,    whichdycore, &
                                 thck,     btraction,   &
                                 ubas,     vbas,        &
                                 float,    bfricflx)

    ! compute frictional heat source due to sliding at the bed

    use glimmer_physcon,  only: rhoi, grav
    use glimmer_paramets, only: thk0, vel0, vel_scale

    type(glide_global_type) :: model
    integer, intent(in) :: whichdycore   ! 1 = Glam, 2 = Glissade
    real(dp), dimension(:,:), intent(in) :: thck
    real(dp), dimension(:,:), intent(in) :: ubas, vbas
    real(dp), dimension(:,:,:), intent(in) :: btraction
    logical, dimension(:,:), intent(in) :: float
    ! Note: bfricflx needs to have intent (inout) in case it has already been computed by Glissade.
    real(dp), dimension(:,:), intent(inout) :: bfricflx

    real(dp) :: slterm       ! sliding friction
 
    integer :: ewp, nsp, ew, ns
    integer :: slide_count   ! number of neighbor cells with nonzero sliding

    if (whichdycore == DYCORE_GLISSADE) then

       ! basal friction heat flux (model%temper%bfricflx) already computed in velocity solver
       ! do nothing and return

    else   ! Glam dycore

       ! compute heat source due to basal friction
       ! Note: slterm and bfricflx are defined to be >= 0

       do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1

          slterm = 0.d0
          slide_count = 0

             !WHL - copied Steve Price's formulation from calcbmlt
             ! btraction is computed in glam_strs2.F90

             !WHL - Using thklim instead of thklim_temp because ice thinner than thklim
             !      is assumed to be at rest.

             if (thck(ew,ns) > model%numerics%thklim .and. .not. float(ew,ns)) then
                do nsp = ns-1,ns
                do ewp = ew-1,ew

                   !SCALING - WHL: Multiplied ubas by vel0/vel_scale so we get the same result in these two cases:
                   !           (1) With scaling:     vel0 = vel_scale = 500/scyr, and ubas is non-dimensional
                   !           (2) Without scaling:  vel0 = 1, vel_scale = 500/scyr, and ubas is in m/s.

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
          !NOTE - The following logic may result in non-conservation of energy.  Include all nonzero terms? 

          if (slide_count >= 4) then
             slterm = 0.25d0 * slterm
          else
             slterm = 0.0d0
          end if

          bfricflx(ew,ns) = slterm

       enddo    ! ns
       enddo    ! ew

    endif       ! whichdycore

  end subroutine glissade_calcbfric

  !-----------------------------------------------------------------------------------

  subroutine glissade_calcbmlt( model,                   &
                                temp,        stagsigma,     &
                                thck,                    &
                                bmlt_ground, floater)

    ! Compute the amount of basal melting.
    ! The basal melting computed here is applied to the ice thickness
    !  by glissade_transport_driver, conserving mass and energy.
    !
    ! Any internal temperatures above the pressure melting point are reset to the
    !  pmp temperature, with excess energy applied toward basal melting.
    !  Hopefully this is rare.
    ! TODO: Moving all internal melting to the basal surface is not very realistic 
    !       and should be revisited.
    ! Note: Since this module is deprecated, the calculation has not been updated
    !       to include bmlt_float.

    use glimmer_physcon, only: shci, rhoi, lhci
    use glimmer_paramets, only : thk0, tim0

    type(glide_global_type) :: model

    real(dp), dimension(0:,:,:), intent(inout) :: temp
    real(dp), dimension(0:),     intent(in) :: stagsigma   ! This is Tstagsigma, (0:upn)
    real(dp), dimension(:,:),    intent(in) :: thck
    real(dp), dimension(:,:),    intent(out):: bmlt_ground ! scaled melt rate (m/s * tim0/thk0)
                                                           ! > 0 for melting, < 0 for freeze-on
    logical,  dimension(:,:),    intent(in) :: floater

    real(dp), dimension(size(stagsigma))    :: pmptemp   ! pressure melting point temperature
    real(dp) :: bflx    ! heat flux available for basal melting (W/m^2)
    integer :: up, ew, ns

    real(dp) :: layer_thck    ! layer thickness (m)
    real(dp) :: melt_energy   ! energy available for internal melting (J/m^2)
    real(dp) :: internal_melt_rate   ! internal melt rate, transferred to bed (m/s)

    real(dp), parameter :: eps11 = 1.d-11       ! small number

    bmlt_ground(:,:) = 0.0d0

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1

          if (thck(ew,ns) > model%numerics%thklim_temp .and. .not. floater(ew,ns)) then

             ! Basal friction term is computed above in subroutine glissade_calcbfric,
             !  or in the Glissade velocity solver.
             !
             ! Compute basal melting
             ! Note: bmlt > 0 for melting, < 0 for freeze-on
             !       bfricflx >= 0 by definition
             !       bheatflx is positive down, so usually bheatflx < 0 (with negative values contributing to melt)
             !       lcondflx is positive down, so lcondflx < 0 for heat is flowing from the bed toward the surface

             !TODO - This equation allows for freeze-on (bmlt < 0) if the conductive term 
             !       (lcondflx, positive down) is carrying enough heat away from the boundary.  
             !       But freeze-on requires a local water supply, bwat > 0.
             !       What should we do if bwat = 0?

             bflx = model%temper%bfricflx(ew,ns) + model%temper%lcondflx(ew,ns) - model%temper%bheatflx(ew,ns)  ! W/m^2

             ! bflx might be slightly different from zero because of rounding errors; if so, then set bflx = 0
             if (abs(bflx) < eps11) bflx = 0.d0

             bmlt_ground(ew,ns) = bflx * model%tempwk%f(2)   ! f(2) = tim0 / (thk0 * lhci * rhoi)

            ! Add internal melting associated with temp > pmptemp
            ! Note: glissade_calcpmpt does not compute pmpt at the top surface or the bed.

             call glissade_calcpmpt(pmptemp(1:model%general%upn-1), thck(ew,ns), stagsigma(1:model%general%upn-1))

             do up = 1, model%general%upn-1
                 if (temp(up,ew,ns) > pmptemp(up)) then
                    ! compute excess energy available for melting
                    layer_thck = thck(ew,ns) * (model%numerics%sigma(up+1) - model%numerics%sigma(up)) * thk0  ! m
                    melt_energy = rhoi * shci * (temp(up,ew,ns) - pmptemp(up)) * layer_thck         ! J/m^2
                    ! compute melt rate 
                    internal_melt_rate = melt_energy / (rhoi * lhci * model%numerics%dttem * tim0)  ! m/s
                    ! transfer internal melting to the bed
                    bmlt_ground(ew,ns) = bmlt_ground(ew,ns) + internal_melt_rate * tim0/thk0  ! m/s * tim0/thk0
                    ! reset T to Tpmp
                    temp(up,ew,ns) = pmptemp(up)
                 endif
             enddo

             ! Cap basal temp at pmptemp, if necessary

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

  end subroutine glissade_calcbmlt

!-------------------------------------------------------------------
 
  subroutine glissade_finddisp (model,     &
                                thck,      &
                                whichdisp, &
                                efvs,      &
                                stagthck,  &
                                dusrfdew,  &
                                dusrfdns,  &
                                flwa)

    ! Compute the dissipation source term associated with strain heating.
    ! Note that the dissipation is computed in the same way on either a staggered or an
    !  unstaggered vertical grid.  
    ! Note also that dissip and flwa must have the same vertical dimension 
    !  (1:upn on an unstaggered vertical grid, or 1:upn-1 on a staggered vertical grid).
    
    use glimmer_physcon, only : gn   ! Glen's n

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(in) :: thck, stagthck, dusrfdew, dusrfdns
    real(dp), dimension(:,:,:), intent(in) :: flwa, efvs
    integer, intent(in) :: whichdisp

    integer, parameter :: p1 = gn + 1  
    integer :: ew, ns

    real(dp) :: c2 
    real(dp), dimension(model%general%upn-1) :: c5     

    model%temper%dissip(:,:,:) = 0.0d0
    
    select case( whichdisp ) 

    case(HO_DISP_NONE)

       ! do nothing; return dissip = 0 everywhere

    case(HO_DISP_SIA)   ! required for whichapprox = HO_APPROX_LOCAL_SIA

    !*sfp* 0-order SIA case only
    ! two methods of doing this: 
    ! 1. find dissipation at u-pts and then average
    ! 2. find dissipation at H-pts by averaging quantities from u-pts
    ! (2) works best for eismint divide (symmetry) but I likely to be better for full expts

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1

          !WHL - Using thklim instead of thklim_temp because ice thinner than thklim
          !      is assumed to be at rest.
          if (thck(ew,ns) > model%numerics%thklim) then
             
             c2 = (0.25d0*sum(stagthck(ew-1:ew,ns-1:ns)) * dsqrt((0.25d0*sum(dusrfdew(ew-1:ew,ns-1:ns)))**2 &
                                                               + (0.25d0*sum(dusrfdns(ew-1:ew,ns-1:ns)))**2))**p1
             
             model%temper%dissip(:,ew,ns) = c2 * model%tempwk%c1(:) * ( &
                  flwa(:,ew-1,ns-1) + flwa(:,ew-1,ns+1) + flwa(:,ew+1,ns+1) + flwa(:,ew+1,ns-1) + &
                  2*(flwa(:,ew-1,ns)+flwa(:,ew+1,ns)+flwa(:,ew,ns-1)+flwa(:,ew,ns+1)) + &
                  4*flwa(:,ew,ns))

          end if
       end do
    end do

    case(HO_DISP_FIRSTORDER)

    ! 3D, 1st-order case
    ! Note: Glissade computes efvs and tau%scalar using the strain rate terms appropriate for the approximation.
    ! E.g, the SIA quantities are computed based on (du_dz, dv_dz) only, and the SSA quantities
    !  are computed based on (du_dx, du_dy, dv_dx, dv_dy) only.
    ! So this computation should give the appropriate heating for whichapprox = HO_APPROX_SIA,
    !  HO_APPROX_SSA, HO_APPROX_L1L2 or HO_APPROX_BP.
    !
    ! NOTE (SFP): For simplicity, tau can be calculated from: tau = 2*efvs*eps_eff,
    ! where eps_eff is the eff. strain rate. Further, eps_eff can be calculated from the efvs according to a 
    ! re-arrangement of: efvs = 1/2 * ( 1 / A(T) )^(1/n) * eps_eff^((1-n)/n), in which case only the efvs and rate
    ! factor arrays need to be passed in for this calculation.

    if (size(model%temper%dissip,1) /= model%general%upn-1) then  ! staggered vertical grid
        !TODO - Write an error message and exit gracefully
    endif

    do ns = 1, model%general%nsn
       do ew = 1, model%general%ewn

          !WHL - Using thklim instead of thklim_temp because ice thinner than thklim
          !      is assumed to be at rest.
          if (thck(ew,ns) > model%numerics%thklim) then

             c5(:) = 0.0d0

             if ( sum( efvs(:,ew,ns) ) /= 0.0d0) then

                ! Use space in c5 vector to store dissip terms that apply at layer midpoints 
                ! (i.e. on staggered vertical grid).  No vertical averaging is needed, since
                ! temp and dissip are colocated with eff stress and eff viscosity.

                 c5(:) = model%stress%tau%scalar(:,ew,ns)**2 / efvs(:,ew,ns)
             endif

             !Note: model%tempwk%cons(5) = (tau0*vel0/len0) / (rhoi*shci) * (model%numerics%dttem*tim0)

             model%temper%dissip(:,ew,ns) = c5(:) * model%tempwk%cons(5)

          endif
       enddo
    enddo

    end select

  end subroutine glissade_finddisp

  !-----------------------------------------------------------------------------------
 
  !TODO - Inline glissade_calcpmpt and glissade_calcbpmp above?

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

  subroutine glissade_calcbpmp(ewn,   nsn,  &
                               thck,  bpmp)

    ! Calculate the pressure melting point at the base of the ice sheet

    integer, intent(in) ::  ewn, nsn    ! grid dimensions

    real(dp), dimension(:,:), intent(in)  :: thck  ! ice thickness (dimensionless)
    real(dp), dimension(:,:), intent(out) :: bpmp  ! bed pressure melting point (deg C)

    integer :: ew,ns

    bpmp(:,:) = 0.d0

    do ns = 1, nsn
       do ew = 1, ewn
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

    pmptemp_bed = - grav * rhoi * pmlt * thk0 * thck 

  end subroutine glissade_calcpmpt_bed

  !-------------------------------------------------------------------

  subroutine glissade_calcflwa(stagsigma,   thklim,   &
                               flwa,        temp,     &
                               thck,        flow_enhancement_factor, &
                               default_flwa_arg,      &
                               flag,        waterfrac)

    ! Calculates Glen's $A$ over the three-dimensional domain,
    ! using one of three possible methods.
    !
    ! The primary method is to use this equation from \emph{Paterson and Budd} [1982]:
    ! \[
    ! A(T^{*})=a \exp \left(\frac{-Q}{RT^{*}}\right)
    ! \]
    ! This is equation 9 in {\em Payne and Dongelmans}. $a$ is a constant of proportionality,
    ! $Q$ is the activation energy for for ice creep, and $R$ is the universal gas constant.
    ! The pressure-corrected temperature, $T^{*}$ is given by:
    ! \[
    ! T^{*}=T-T_{\mathrm{pmp}}+T_0
    ! \] 
    ! \[
    ! T_{\mathrm{pmp}}=T_0-\sigma \rho g H \Phi
    ! \]
    ! $T$ is the ice temperature, $T_{\mathrm{pmp}}$ is the pressure melting point 
    ! temperature, $T_0$ is the triple point of water, $\rho$ is the ice density, and 
    ! $\Phi$ is the (constant) rate of change of melting point temperature with pressure.

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
    real(dp)                                  :: flow_enhancement_factor ! flow enhancement factor in Arrhenius relationship
    real(dp),                   intent(in)    :: default_flwa_arg ! Glen's A to use in isothermal case 
                                                                  ! Units: Pa^{-n} yr^{-1} 
    integer,                    intent(in)    :: flag      !> Flag to select the method of calculation
    real(dp),dimension(:,:,:),  intent(out)   :: flwa      !> The calculated values of $A$
    real(dp),dimension(:,:,:),  intent(in), optional :: waterfrac !> internal water content fraction, 0 to 1

    !> \begin{description}
    !> \item[0] {\em Paterson and Budd} relationship.
    !> \item[1] {\em Paterson and Budd} relationship, with temperature set to -5$^{\circ}$C.
    !> \item[2] Set to prescribed constant value.
    !> \end{description}

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp) :: default_flwa
    integer :: ew, ns, up, ewn, nsn, uflwa
    real(dp) :: tempcor

    real(dp),dimension(4), parameter ::  &
       arrfact = (/ arrmlh / vis0,      &   ! Value of A when T* is above -263K
                    arrmll / vis0,      &   ! Value of A when T* is below -263K
                   -actenh / gascon,    &   ! Value of -Q/R when T* is above -263K
                   -actenl / gascon/)       ! Value of -Q/R when T* is below -263K
    
    real(dp), parameter :: const_temp = -5.0d0
    real(dp), parameter :: flwa_waterfrac_enhance_factor = 181.25d0

    !------------------------------------------------------------------------------------ 
   
    uflwa=size(flwa,1) ; ewn=size(flwa,2) ; nsn=size(flwa,3)

    ! Check that the temperature array has the desired vertical dimension

    if (size(temp,1) /= size(flwa,1)) then
       call write_log('glissade_calcflwa: temp and flwa must have the same vertical dimensions', GM_FATAL)
    endif

    ! Scale the default rate factor (default value has units Pa^{-n} yr^{-1}).
    ! Also multiply by flow enhancement factor

    default_flwa = flow_enhancement_factor * default_flwa_arg / (vis0*scyr)
    !write(*,*)"Default flwa = ",default_flwa

    select case(flag)

    case(FLWA_PATERSON_BUDD)

      ! This is the Paterson and Budd relationship
      ! BDM added waterfrac relationship for whichtemp=TEMP_ENTHALPY case

      do ns = 1,nsn
         do ew = 1,ewn

            if (thck(ew,ns) > thklim) then
            
               do up = 1, uflwa   ! uflwa = upn - 1 (values at layer midpoints)

                  ! Calculate the corrected temperature

                  tempcor = min(0.0d0, temp(up,ew,ns) + thck(ew,ns)*grav*rhoi*pmlt*thk0*stagsigma(up))
                  tempcor = max(-50.0d0, tempcor)

                  ! Calculate Glen's A (including flow enhancement factor)

                  if (tempcor >= -10.d0) then
                     flwa(up,ew,ns) = flow_enhancement_factor * arrfact(1) * exp(arrfact(3)/(tempcor + trpt))
                  else
                     flwa(up,ew,ns) = flow_enhancement_factor * arrfact(2) * exp(arrfact(4)/(tempcor + trpt))
                  endif

                  ! BDM added correction for a liquid water fraction 
                  ! Using Greve and Blatter, 2009 formulation for Glen's A flow rate factor:
                  !    A = A(theta_PMP) * (1 + 181.25 * waterfrac)
		  ! RJH - commenting out waterfrac correction to explore causes of
		  ! oscillations in thk and vel for EISMINT-2 test cases
                  if (present(waterfrac)) then
                     if (waterfrac(up,ew,ns) > 0.0d0) then
                        flwa(up,ew,ns) = flwa(up,ew,ns) * (1.d0 + flwa_waterfrac_enhance_factor * waterfrac(up,ew,ns))      
                     endif
                  endif

               enddo   ! up

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

               ! Calculate Glen's A with a fixed temperature (including flow enhancement factor)

               if (const_temp >= -10.d0) then
                  flwa(:,ew,ns) = flow_enhancement_factor * arrfact(1) * exp(arrfact(3)/(const_temp + trpt))
               else
                  flwa(:,ew,ns) = flow_enhancement_factor * arrfact(2) * exp(arrfact(4)/(const_temp + trpt))
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
