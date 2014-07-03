!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_velo.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_velo

  !*FD Contains routines which handle various aspects of velocity in the model,
  !*FD not only the bulk ice velocity, but also basal sliding, and vertical grid 
  !*FD velocities, etc.

  use glide_types
  use glimmer_global, only : dp
  use glimmer_physcon, only : rhoi, grav, gn
  use glimmer_paramets, only : thk0, len0, vis0, vel0

  implicit none

  private vertintg

  ! some private parameters
  integer, private, parameter :: p1 = gn+1
  integer, private, parameter :: p2 = gn-1
  integer, private, parameter :: p3 = 2*gn+1
  integer, private, parameter :: p4 = gn+2

  real(dp),private, parameter :: cflow = -2.0d0*vis0*(rhoi*grav)**gn*thk0**p3/(8.0d0*vel0*len0**gn)

contains

!TODO - Pretty sure that none of the arrays in this subroutine are needed for HO,
!        so we may not need to call this subroutine from glissade_initialise.

!       Some velowk arrays are used in wvelintg, but not hard to rewrite wvelintg without these arrays.

  subroutine init_velo(model)
    !*FD initialise velocity module
    use glimmer_physcon, only : arrmll, arrmlh, gascon, actenl, actenh, scyr, pi 
    implicit none
    type(glide_global_type) :: model

    integer ewn, nsn, upn
    integer up

    ewn=model%general%ewn
    nsn=model%general%nsn
    upn=model%general%upn

    allocate(model%velowk%fslip(ewn-1,nsn-1))

    allocate(model%velowk%depth(upn))
    allocate(model%velowk%dintflwa(ewn-1,nsn-1))

    model%velowk%depth = (/ (((model%numerics%sigma(up+1)+model%numerics%sigma(up))/2.0d0)**gn &
         *(model%numerics%sigma(up+1)-model%numerics%sigma(up)),up=1,upn-1),0.0d0 /)

    allocate(model%velowk%dups(upn)) 
    model%velowk%dups = (/ (model%numerics%sigma(up+1) - model%numerics%sigma(up), up=1,upn-1),0.0d0 /)

    allocate(model%velowk%dupsw (upn))
    allocate(model%velowk%depthw(upn))
    allocate(model%velowk%suvel (upn))
    allocate(model%velowk%svvel (upn))

    ! Calculate the differences between adjacent sigma levels -------------------------

    model%velowk%dupsw  = (/ (model%numerics%sigma(up+1)-model%numerics%sigma(up), up=1,upn-1), 0.0d0 /) 

    ! Calculate the value of sigma for the levels between the standard ones -----------

    model%velowk%depthw = (/ ((model%numerics%sigma(up+1)+model%numerics%sigma(up)) / 2.0d0, up=1,upn-1), 0.0d0 /)

    model%velowk%fact = (/ model%paramets%flow_factor* arrmlh / vis0, &   ! Value of a when T* is above -263K
                           model%paramets%flow_factor* arrmll / vis0, &   ! Value of a when T* is below -263K
                          -actenh / gascon,        &                      ! Value of -Q/R when T* is above -263K
                          -actenl / gascon/)                              ! Value of -Q/R when T* is below -263K

    model%velowk%watwd  = model%paramets%bpar(1)
    model%velowk%watct  = model%paramets%bpar(2)
    model%velowk%trcmin = model%paramets%bpar(3) / scyr
    model%velowk%trcmax = model%paramets%bpar(4) / scyr
    model%velowk%marine = model%paramets%bpar(5)
    model%velowk%trcmax = model%velowk%trcmax / model%velowk%trc0
    model%velowk%trcmin = model%velowk%trcmin / model%velowk%trc0  
    model%velowk%c(1)   = (model%velowk%trcmax + model%velowk%trcmin) / 2.0d0 
    model%velowk%c(2)   = (model%velowk%trcmax - model%velowk%trcmin) / 2.0d0
    model%velowk%c(3)   = (thk0 * pi) / model%velowk%watwd  
    model%velowk%c(4)   = pi*(model%velowk%watct / model%velowk%watwd)
   

  end subroutine init_velo

!TODO - What is the 'if 0' option?
!TODO - May be able to remove this subroutine; it is never called. 
#if 0
  subroutine velo_compute_strain_rates(strain_zx, strain_zy, 
                                       stagthck, dusrfdew, dusrfdns, sigma, 
                                       flwa, ho_stress_zx, ho_stress_zy)

    !*FD Computes the strain rates \epsilon_{zx} and \epsilon{zy} using the
    !*FD formula:
    !*FD \epsilon_{zi}(\sigma) = A(\sigma)(-\rho g H \sigma \frac{\partial s}{\partial i} - HO_i)^n
    !*FD No vertical integration is done at this point
    
    ewn = size(strain_zx, 2)
    nsn = size(strain_zx, 3)
    upn = size(strain_zx, 1)

    !TODO: Make this work with rescaling!!
    !TODO: Vectorize
    do i = 1, upn
      do j = 1, ewn
        do k = 1, nsn
          !Compute everything inside the exponentiation
          !(we factor out -rhoi*g*H*\sigma so it's only computed once
!TODO - Where is g defined as the gravitational constant?  Should be grav.
          zx = -rhoi*g*sigma(i)*stagthck(j,k)
          zy = zx*dusrfdns(j,k)
          zx = zx*dusrfdew(j,k)

          !Add higher-order stresses if they were provided
          if (present(ho_stress_zx) .and. present(ho_stress_zy)) then
            zx = zx - ho_stress_zx(i,j,k)
            zy = zy - ho_stress_zy(i,j,k)
          end if

          !Exponentiate and finish computation
          zx = A(i,j,k)*zx**gn
          zy = A(i,j,k)*zy**gn

          strain_zx(i,j,k) = zx
          strain_zy(i,j,k) = zy
        end do
      end do
    end do
  end subroutine

!TODO - This subroutine is not called.  Remove it?
  
  !*FD Integrates the strain rates to compute both the 3d velocity fields and the
  !*FD vertically averaged velocities

  subroutine velo_integrate_strain(strain_zx, strain_zy, ubas, vbas,
    !Find the 3d velocity field by vertically integrating the strain rate from
    !the base up
    call vertint_output3d(strain_zx, uvel, sigma, .false., ubas)
    call vertint_output3d(strain_zy, vvel, sigma, .false., vbas)
    
    !Normally, we have \epsilon_{xz} = \frac{1}{2}(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x})
    !However, we assume that $u_z >> w_x$ by assuming simple shear deformation.
    !This gives us $u_z$ = 2\epsilon_{xz}
    !This 2 is pulled out of the integral.  We correct for it here.
    uvel = 2*uvel
    vvel = 2*vvel
           
  end subroutine velo_integrate_strain
#endif
  
  !*****************************************************************************
  ! new velo functions come here
  !*****************************************************************************

  subroutine velo_integrate_flwa(velowk,stagthck,flwa)
    
    !*FD this routine calculates the part of the vertically averaged velocity 
    !*FD field which solely depends on the temperature
    !*FD (The integral in eq. 3.22d)

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glide_velowk),     intent(inout) :: velowk           
    real(dp),dimension(:,:),  intent(in)    :: stagthck       !*FD ice thickness on staggered grid
    real(dp),dimension(:,:,:),intent(in)    :: flwa           !*FD ice flow factor

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------
    real(dp),dimension(size(flwa,1)) :: hrzflwa, intflwa 
    integer :: ew,ns,up,ewn,nsn,upn

    upn=size(flwa,1) ; ewn=size(flwa,2) ; nsn=size(flwa,3)

    do ns = 1,nsn-1
       do ew = 1,ewn-1
          if (stagthck(ew,ns) /= 0.0d0) then
             
             hrzflwa = flwa(:,ew,ns) + flwa(:,ew,ns+1) + flwa(:,ew+1,ns) + flwa(:,ew+1,ns+1)
             intflwa(upn) = 0.0d0

             !Perform inner integration.
             do up = upn-1, 1, -1
                intflwa(up) = intflwa(up+1) + velowk%depth(up) * (hrzflwa(up)+hrzflwa(up+1))
             end do

             velowk%dintflwa(ew,ns) = cflow * vertintg(velowk,intflwa)

          else 

             velowk%dintflwa(ew,ns) = 0.0d0

          end if
       end do
    end do
  end subroutine velo_integrate_flwa

  !*****************************************************************************

  subroutine velo_calc_diffu(velowk,stagthck,dusrfdew,dusrfdns,diffu)

    !*FD calculate diffusivities

    implicit none
    
    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glide_velowk),     intent(inout) :: velowk
    real(dp),dimension(:,:),  intent(in)    :: stagthck
    real(dp),dimension(:,:),  intent(in)    :: dusrfdew
    real(dp),dimension(:,:),  intent(in)    :: dusrfdns
    real(dp),dimension(:,:),  intent(out)   :: diffu


    where (stagthck  /=  0.0d0)
       diffu = velowk%dintflwa * stagthck**p4 * sqrt(dusrfdew**2 + dusrfdns**2)**p2 
    elsewhere
       diffu = 0.0d0
    end where

  end subroutine velo_calc_diffu

  !*****************************************************************************

  subroutine velo_calc_velo(velowk,   stagthck,    &
                            dusrfdew, dusrfdns,    &
                            flwa,     diffu,       &
                            ubas,     vbas,        &
                            uvel,     vvel,        &
                            uflx,     vflx,        &
                            velnorm)

    !*FD calculate 3D horizontal velocity field and 2D flux field from diffusivity
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glide_velowk),     intent(inout) :: velowk
    real(dp),dimension(:,:),  intent(in)    :: stagthck
    real(dp),dimension(:,:),  intent(in)    :: dusrfdew
    real(dp),dimension(:,:),  intent(in)    :: dusrfdns
    real(dp),dimension(:,:,:),intent(in)    :: flwa
    real(dp),dimension(:,:),  intent(in)    :: diffu
    real(dp),dimension(:,:),  intent(in)    :: ubas
    real(dp),dimension(:,:),  intent(in)    :: vbas
    real(dp),dimension(:,:,:),intent(out)   :: uvel
    real(dp),dimension(:,:,:),intent(out)   :: vvel
    real(dp),dimension(:,:),  intent(out)   :: uflx
    real(dp),dimension(:,:),  intent(out)   :: vflx
    real(dp),dimension(:,:,:), intent(out)  :: velnorm

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------
    real(dp),dimension(size(flwa,1)) :: hrzflwa
    real(dp) :: factor
    real(dp),dimension(3)           :: const
    integer :: ew,ns,up,ewn,nsn,upn

    upn=size(flwa,1) ; ewn=size(stagthck,1) ; nsn=size(stagthck,2)
    
!      Note: Here (confusingly), nsn = size(stagthck,2) = model%general%nsn-1
!                                ewn = size(stagthck,1) = model%general%ewn-1
!TODO - Change loop limits to nsn-1, ewn-1
!       (provided ewn = size(stagthck,1)-1, etc.)

    do ns = 1,nsn
       do ew = 1,ewn

          if (stagthck(ew,ns) /= 0.0d0) then

             vflx(ew,ns) = diffu(ew,ns) * dusrfdns(ew,ns) + vbas(ew,ns) * stagthck(ew,ns)
             uflx(ew,ns) = diffu(ew,ns) * dusrfdew(ew,ns) + ubas(ew,ns) * stagthck(ew,ns)

             uvel(upn,ew,ns) = ubas(ew,ns)
             vvel(upn,ew,ns) = vbas(ew,ns)

             hrzflwa = flwa(:,ew,ns) + flwa(:,ew,ns+1) + flwa(:,ew+1,ns) + flwa(:,ew+1,ns+1)

             factor = velowk%dintflwa(ew,ns)*stagthck(ew,ns)
             if (factor /= 0.0d0) then

                const(2) = cflow * diffu(ew,ns) / factor

                const(3) = const(2) * dusrfdns(ew,ns)  
                const(2) = const(2) * dusrfdew(ew,ns) 
             else
                const(2:3) = 0.0d0
             end if

             do up = upn-1, 1, -1
                const(1) = velowk%depth(up) * (hrzflwa(up)+hrzflwa(up+1))
                uvel(up,ew,ns) = uvel(up+1,ew,ns) + const(1) * const(2)
                vvel(up,ew,ns) = vvel(up+1,ew,ns) + const(1) * const(3) 
             end do

          else 

             uvel(:,ew,ns) = 0.0d0
             vvel(:,ew,ns) = 0.0d0
             uflx(ew,ns) = 0.0d0
             vflx(ew,ns) = 0.0d0 

          end if
       end do
    end do
    
    ! horizontal ice speed (mainly for diagnostic purposes)

    velnorm(:,:,:) = sqrt((uvel(:,:,:)**2 + vvel(:,:,:)**2))

  end subroutine velo_calc_velo

  !*****************************************************************************
  ! old velo functions come here
  !*****************************************************************************

  subroutine slipvelo(model,flag1,btrc,ubas,vbas)

    !*FD Calculate the basal slip velocity and the value of $B$, the free parameter
    !*FD in the basal velocity equation (though I'm not sure that $B$ is used anywhere 
    !*FD else).

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_global_type) :: model                  !*FD model instance
    integer, intent(in)                 :: flag1      !*FD \texttt{flag1} sets the calculation
                                                      !*FD method to use for the basal velocity
                                                      !*FD (corresponded to \texttt{whichslip} in the
                                                      !*FD old model. 
    real(dp),dimension(:,:),intent(in)   :: btrc     !*FD The basal slip coefficient.
    real(dp),dimension(:,:),intent(out)   :: ubas     !*FD The $x$ basal velocity (scaled)
    real(dp),dimension(:,:),intent(out)   :: vbas     !*FD The $y$ basal velocity (scaled)

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp), parameter :: rhograv = - rhoi * grav
    integer :: nsn,ewn

    ! Get array sizes -------------------------------------------------------------------

    ewn=size(btrc,1) ; nsn=size(btrc,2)    

    !------------------------------------------------------------------------------------
    ! Main calculation starts here
    !------------------------------------------------------------------------------------

    select case(flag1)
    case(0)  
    
      ! Linear function of gravitational driving stress ---------------------------------

      where (model%numerics%thklim < model%geomderv%stagthck)
        ubas = btrc * rhograv * model%geomderv%stagthck * model%geomderv%dusrfdew
        vbas = btrc * rhograv * model%geomderv%stagthck * model%geomderv%dusrfdns
      elsewhere
        ubas = 0.0d0
        vbas = 0.0d0
      end where

    case(1)

      ! *tp* option to be used in picard iteration for thck
      ! *tp* start by find constants which dont vary in iteration

      model%velowk%fslip = rhograv * btrc

    case(2)

      ! *tp* option to be used in picard iteration for thck
      ! *tp* called once per non-linear iteration, set uvel to ub * H /(ds/dx) which is
      ! *tp* a diffusivity for the slip term (note same in x and y)

      where (model%numerics%thklim < model%geomderv%stagthck)
        ubas = model%velowk%fslip * model%geomderv%stagthck**2  
      elsewhere
        ubas = 0.0d0
      end where

    case(3)

      ! *tp* option to be used in picard iteration for thck
      ! *tp* finally calc ub and vb from diffusivities

      where (model%numerics%thklim < model%geomderv%stagthck)
        vbas = ubas *  model%geomderv%dusrfdns / model%geomderv%stagthck
        ubas = ubas *  model%geomderv%dusrfdew / model%geomderv%stagthck
      elsewhere
        ubas = 0.0d0
        vbas = 0.0d0
      end where

    case default
      ubas = 0.0d0
      vbas = 0.0d0
    end select

  end subroutine slipvelo

!------------------------------------------------------------------------------------------

  subroutine zerovelo(velowk,sigma,flag,stagthck,dusrfdew,dusrfdns,flwa,ubas,vbas,uvel,vvel,uflx,vflx,diffu)

    !*FD Performs the velocity calculation. This subroutine is called with
    !*FD different values of \texttt{flag}, depending on exactly what we want to calculate.

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_velowk),     intent(inout) :: velowk
    real(dp),dimension(:),    intent(in)    :: sigma
    integer,                  intent(in)    :: flag
    real(dp),dimension(:,:),  intent(in)    :: stagthck
    real(dp),dimension(:,:),  intent(in)    :: dusrfdew
    real(dp),dimension(:,:),  intent(in)    :: dusrfdns
    real(dp),dimension(:,:,:),intent(in)    :: flwa
    real(dp),dimension(:,:),  intent(in)    :: ubas
    real(dp),dimension(:,:),  intent(in)    :: vbas
    real(dp),dimension(:,:,:),intent(out)   :: uvel
    real(dp),dimension(:,:,:),intent(out)   :: vvel
    real(dp),dimension(:,:),  intent(out)   :: uflx
    real(dp),dimension(:,:),  intent(out)   :: vflx
    real(dp),dimension(:,:),  intent(out)   :: diffu

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    
    real(dp),dimension(size(sigma)) :: hrzflwa, intflwa 
    real(dp),dimension(3)           :: const

    integer :: ew,ns,up,ewn,nsn,upn

    !------------------------------------------------------------------------------------

    upn=size(sigma) ; ewn=size(ubas,1) ; nsn=size(ubas,2)


    !------------------------------------------------------------------------------------

    select case(flag)
    case(0)

      do ns = 1,nsn
        do ew = 1,ewn

          if (stagthck(ew,ns) /= 0.0d0) then

            ! Set velocity to zero at base of column

            uvel(upn,ew,ns) = 0.0d0
            vvel(upn,ew,ns) = 0.0d0

            ! Get column profile of Glen's A

            hrzflwa = flwa(:,ew,ns) + flwa(:,ew,ns+1) + flwa(:,ew+1,ns) + flwa(:,ew+1,ns+1)

            ! Calculate coefficient for integration

            const(1) = cflow * stagthck(ew,ns)**p1 * sqrt(dusrfdew(ew,ns)**2 + dusrfdns(ew,ns)**2)**p2  

            ! Do first step of finding u according to (8) in Payne and Dongelmans 

            do up = upn-1, 1, -1
              uvel(up,ew,ns) = uvel(up+1,ew,ns) + const(1) * velowk%depth(up) * sum(hrzflwa(up:up+1)) 
            end do

            ! Calculate u diffusivity (?)

            diffu(ew,ns) = vertintg(velowk,uvel(:,ew,ns)) * stagthck(ew,ns)

            ! Complete calculation of u and v

            vvel(:,ew,ns) = uvel(:,ew,ns) * dusrfdns(ew,ns) + vbas(ew,ns)
            uvel(:,ew,ns) = uvel(:,ew,ns) * dusrfdew(ew,ns) + ubas(ew,ns)

            ! Calculate ice fluxes

            uflx(ew,ns) = diffu(ew,ns) * dusrfdew(ew,ns) + ubas(ew,ns) * stagthck(ew,ns)
            vflx(ew,ns) = diffu(ew,ns) * dusrfdns(ew,ns) + vbas(ew,ns) * stagthck(ew,ns)

          else 

            ! Where there is no ice, set everything to zero.

            uvel(:,ew,ns) = 0.0d0
            vvel(:,ew,ns) = 0.0d0
            uflx(ew,ns)   = 0.0d0
            vflx(ew,ns)   = 0.0d0
            diffu(ew,ns)  = 0.0d0

          end if

        end do
      end do

    case(1)

      do ns = 1,nsn
        do ew = 1,ewn
          if (stagthck(ew,ns) /= 0.0d0) then

            hrzflwa = flwa(:,ew,ns) + flwa(:,ew,ns+1) + flwa(:,ew+1,ns) + flwa(:,ew+1,ns+1)
            intflwa(upn) = 0.0d0

            do up = upn-1, 1, -1
               intflwa(up) = intflwa(up+1) + velowk%depth(up) * sum(hrzflwa(up:up+1)) 
            end do

            velowk%dintflwa(ew,ns) = cflow * vertintg(velowk,intflwa)

          else 

            velowk%dintflwa(ew,ns) = 0.0d0

          end if
        end do
      end do

    case(2)

      where (stagthck /= 0.0d0)
        diffu = velowk%dintflwa * stagthck**p4 * sqrt(dusrfdew**2 + dusrfdns**2)**p2 
      elsewhere
        diffu = 0.0d0
      end where

    case(3)

      do ns = 1,nsn
        do ew = 1,ewn
          if (stagthck(ew,ns) /= 0.0d0) then

            vflx(ew,ns) = diffu(ew,ns) * dusrfdns(ew,ns) + vbas(ew,ns) * stagthck(ew,ns)
            uflx(ew,ns) = diffu(ew,ns) * dusrfdew(ew,ns) + ubas(ew,ns) * stagthck(ew,ns)

            uvel(upn,ew,ns) = ubas(ew,ns)
            vvel(upn,ew,ns) = vbas(ew,ns)

            hrzflwa = flwa(:,ew,ns) + flwa(:,ew,ns+1) + flwa(:,ew+1,ns) + flwa(:,ew+1,ns+1)

            if (velowk%dintflwa(ew,ns) /= 0.0d0) then

               const(2) = cflow * diffu(ew,ns) / velowk%dintflwa(ew,ns)/stagthck(ew,ns)

               const(3) = const(2) * dusrfdns(ew,ns)  
               const(2) = const(2) * dusrfdew(ew,ns) 
            else
               const(2:3) = 0.0d0
            end if

            do up = upn-1, 1, -1
              const(1) = velowk%depth(up) * sum(hrzflwa(up:up+1)) 
              uvel(up,ew,ns) = uvel(up+1,ew,ns) + const(1) * const(2)
              vvel(up,ew,ns) = vvel(up+1,ew,ns) + const(1) * const(3) 
            end do

          else 

            uvel(:,ew,ns) = 0.0d0
            vvel(:,ew,ns) = 0.0d0
            uflx(ew,ns) = 0.0d0
            vflx(ew,ns) = 0.0d0 

          end if
        end do
      end do

    end select

  end subroutine zerovelo

!------------------------------------------------------------------------------------------

  subroutine glide_velo_vertical(model)

     type(glide_global_type), intent(inout) :: model     ! model instance

     ! Compute the ice vertical velocity

     ! This is a new subroutine created by combining calls to several existing subroutines.

     ! Note: It is now called at the end of glide_tstep_p3, so that exact restart is easier.
     !       In older versions of Glimmer the vertical velocity was computed at the start of 
     !        the temperature calculation in glide_tstep_p1.

     ! Calculate time-derivatives of thickness and upper surface elevation ------------

     call timeders(model%thckwk,              &
                   model%geometry%thck,       &
                   model%geomderv%dthckdtm,   &
                   model%geometry%thck_index, &
                   model%numerics%time,       &
                   1)

     call timeders(model%thckwk,     &
                   model%geometry%usrf,       &
                   model%geomderv%dusrfdtm,   &
                   model%geometry%thck_index, &
                   model%numerics%time,       &
                   2)

     ! Calculate the vertical velocity of the grid ------------------------------------

     call gridwvel(model%numerics%sigma,  &
                   model%numerics%thklim, &
                   model%velocity%uvel,   &
                   model%velocity%vvel,   &
                   model%geomderv,        &
                   model%geometry%thck,   &
                   model%velocity%wgrd)

     ! Calculate the actual vertical velocity; method depends on whichwvel ------------

     select case(model%options%whichwvel)

     case(VERTINT_STANDARD)     ! Usual vertical integration

        call wvelintg(model%velocity%uvel,                        &
                      model%velocity%vvel,                        &
                      model%geomderv,                             &
                      model%numerics,                             &
                      model%velowk,                               &
                      model%velocity%wgrd(model%general%upn,:,:), &
                      model%geometry%thck,                        &
                      model%temper%bmlt,                          &
                      model%velocity%wvel)

     case(VERTINT_KINEMATIC_BC)     ! Vertical integration constrained so kinematic upper BC obeyed.

        call wvelintg(model%velocity%uvel,                        &
                      model%velocity%vvel,                        &
                      model%geomderv,                             &
                      model%numerics,                             &
                      model%velowk,                               &
                      model%velocity%wgrd(model%general%upn,:,:), &
                      model%geometry%thck,                        &
                      model%temper%  bmlt,                        &
                      model%velocity%wvel)

        call chckwvel(model%numerics,                             &
                      model%geomderv,                             &
                      model%velocity%uvel(1,:,:),                 &
                      model%velocity%vvel(1,:,:),                 &
                      model%velocity%wvel,                        &
                      model%geometry%thck,                        &
                      model%climate% acab)

     end select

     ! apply periodic ew BC

     if (model%options%periodic_ew) then
        call wvel_ew(model)
     end if

  end subroutine glide_velo_vertical

!---------------------------------------------------------------

  subroutine timeders(thckwk,ipvr,opvr,mask,time,which)

    !*FD Calculates the time-derivative of a field. This subroutine is used by 
    !*FD the Glimmer temperature solver only.

    use glimmer_global, only : dp, sp
    use glimmer_paramets, only : tim0
    use glimmer_physcon, only: scyr

    implicit none 

    type(glide_thckwk) :: thckwk    !*FD Derived-type containing work data
    real(dp), intent(out), dimension(:,:) :: opvr  !*FD output (time derivative) field
    real(dp), intent(in),  dimension(:,:) :: ipvr  !*FD input field
    real(dp), intent(in)                  :: time  !*FD current time
    integer,  intent(in),  dimension(:,:) :: mask  !*FD mask for calculation
    integer,  intent(in)                  :: which !*FD selector for stored field

    real(dp) :: factor

    factor = (time - thckwk%oldtime)
    if (factor == 0.d0) then
       opvr = 0.0d0
    else
       factor = 1.d0/factor
       where (mask /= 0)
          opvr = (tim0/scyr) * (ipvr - thckwk%olds(:,:,which)) * factor
       elsewhere
          opvr = 0.0d0
       end where
    end if

    thckwk%olds(:,:,which) = ipvr

    if (which == thckwk%nwhich) then
      thckwk%oldtime = time
    end if

  end subroutine timeders

!------------------------------------------------------------------------------------------

  subroutine gridwvel(sigma,thklim,uvel,vvel,geomderv,thck,wgrd)

    !*FD Calculates the vertical velocity of the grid, and returns it in \texttt{wgrd}. This
    !*FD is necessary because the model uses a sigma coordinate system.
    !*FD The equation for grid velocity is:
    !*FD \[
    !*FD \mathtt{wgrd}(x,y,\sigma)=\frac{\partial s}{\partial t}+\mathbf{U}\cdot\nabla s
    !*FD -\sigma\left(\frac{\partial H}{\partial t}+\mathbf{U}\cdot\nabla H\right)
    !*FD \]
    !*FD Compare this with equation A1 in {\em Payne and Dongelmans}.

    !TODO The name of this subroutine is confusing.  It is called wvel but it does not calculate wvel, only wgrd.

    use parallel
!!    use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar
    implicit none 

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    real(dp),dimension(:),    intent(in)  :: sigma     !*FD Array holding values of sigma
                                                       !*FD at each vertical level
    real(dp),                 intent(in)  :: thklim    !*FD Minimum thickness to be considered
                                                       !*FD when calculating the grid velocity.
                                                       !*FD This is in m, divided by \texttt{thk0}.
    real(dp),dimension(:,:,:),intent(in)  :: uvel      !*FD The $x$-velocity field (scaled). Velocity
                                                       !*FD is on the staggered grid
    real(dp),dimension(:,:,:),intent(in)  :: vvel      !*FD The $y$-velocity field (scaled). Velocity
                                                       !*FD is on the staggered grid
    type(glide_geomderv),   intent(in)  :: geomderv  !*FD Derived type holding temporal
                                                       !*FD and horizontal derivatives of
                                                       !*FD ice-sheet thickness and upper
                                                       !*FD surface elevation
    real(dp),dimension(:,:),  intent(in)  :: thck      !*FD Ice-sheet thickness (divided by 
                                                       !*FD \texttt{thk0})
    real(dp),dimension(:,:,:),intent(out) :: wgrd      !*FD The grid velocity at each point. This
                                                       !*FD is the output.

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    integer :: ns,ew,nsn,ewn

    !------------------------------------------------------------------------------------

    ewn=size(wgrd,2) ; nsn=size(wgrd,3)

    do ns = 2,nsn-1
      do ew = 2,ewn-1
        if (thck(ew,ns) > thklim) then
          wgrd(:,ew,ns) = geomderv%dusrfdtm(ew,ns) - sigma * geomderv%dthckdtm(ew,ns) + & 
                      ((uvel(:,ew-1,ns-1) + uvel(:,ew-1,ns) + uvel(:,ew,ns-1) + uvel(:,ew,ns)) * &
                      (sum(geomderv%dusrfdew(ew-1:ew,ns-1:ns)) - sigma * &
                       sum(geomderv%dthckdew(ew-1:ew,ns-1:ns))) + &
                       (vvel(:,ew-1,ns-1) + vvel(:,ew-1,ns) + vvel(:,ew,ns-1) + vvel(:,ew,ns)) * &
                      (sum(geomderv%dusrfdns(ew-1:ew,ns-1:ns)) - sigma * &
                       sum(geomderv%dthckdns(ew-1:ew,ns-1:ns)))) / 16.0d0
        else
          wgrd(:,ew,ns) = 0.0d0
        end if
      end do
    end do

!TODO - Remove halo call? wgrd is needed only for the old temperature code, which is not supported in parallel.
    call parallel_halo(wgrd)
!    call horiz_bcs_unstag_scalar(wgrd)

  end subroutine gridwvel

!------------------------------------------------------------------------------------------

  subroutine wvelintg(uvel,vvel,geomderv,numerics,velowk,wgrd,thck,bmlt,wvel)

    !*FD Calculates the vertical velocity field, which is returned in \texttt{wvel}.
    !*FD This is found by doing this integration:
    !*FD \[
    !*FD w(\sigma)=-\int_{1}^{\sigma}\left[\frac{\partial \mathbf{U}}{\partial \sigma}
    !*FD (\sigma) \cdot (\nabla s - \sigma \nabla H) +H\nabla \cdot \mathbf{U}(\sigma)\right]d\sigma
    !*FD + w(1)
    !*FD \]
    !*FD (This is equation 13 in {\em Payne and Dongelmans}.) Note that this is only 
    !*FD done if the thickness is greater than the threshold given by \texttt{numerics\%thklim}.
    use parallel
!!    use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    real(dp),dimension(:,:,:), intent(in)    :: uvel      !*FD The $x$-velocity on the
                                                          !*FD staggered grid (scaled)
    real(dp),dimension(:,:,:), intent(in)    :: vvel      !*FD The $y$-velocity on the
                                                          !*FD staggered grid (scaled)
    real(dp),dimension(:,:),   intent(in)    :: thck      !*FD The ice thickness, divided
                                                          !*FD by \texttt{thk0}
    type(glide_geomderv),    intent(in)    :: geomderv  !*FD Derived type holding the
                                                          !*FD horizontal and temporal derivatives
                                                          !*FD of the thickness and upper surface
                                                          !*FD elevation.
    type(glide_numerics),    intent(in)    :: numerics  !*FD Derived type holding numerical
                                                          !*FD parameters, including sigma values.
    type(glide_velowk),      intent(inout) :: velowk    !*FD Derived type holding working arrays
                                                          !*FD used by the subroutine
    real(dp),dimension(:,:),   intent(in)    :: wgrd      !*FD The grid vertical velocity at
                                                          !*FD the lowest model level.
    real(dp),dimension(:,:),   intent(in)    :: bmlt      !*FD Basal melt-rate (scaled?) This
                                                          !*FD is required in the basal boundary
                                                          !*FD condition. See {\em Payne and Dongelmans}
                                                          !*FD equation 14.
    real(dp),dimension(:,:,:), intent(out)   :: wvel      !*FD The vertical velocity field.

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp) :: dew16, dns16        ! The grid-spacings multiplied by 16
    real(dp),dimension(6) :: cons   ! Holds temporary local values of derivatives
    integer :: ns,ew,up             ! Loop indicies
    integer :: nsn,ewn,upn          ! Domain sizes

    !------------------------------------------------------------------------------------
    ! Get some values for the domain size by checking sizes of input arrays
    !------------------------------------------------------------------------------------

    upn=size(uvel,1) ; ewn=size(uvel,2) ; nsn=size(uvel,3)


    ! Multiply grid-spacings by 16 -----------------------------------------------------

    dew16 = 1.d0/(16.0d0 * numerics%dew)
    dns16 = 1.d0/(16.0d0 * numerics%dns)

    ! ----------------------------------------------------------------------------------
    ! Main loop over each grid-box
    ! ----------------------------------------------------------------------------------

    do ns = 2,nsn
      do ew = 2,ewn
        if (thck(ew,ns) > numerics%thklim) then
  
          ! Set the bottom boundary condition ------------------------------------------

          wvel(upn,ew,ns) = wgrd(ew,ns) - bmlt(ew,ns)

          ! Calculate temporary local values of thickness and surface ------------------
          ! elevation derivatives.

          cons(1) = sum(geomderv%dusrfdew(ew-1:ew,ns-1:ns)) / 16.0d0
          cons(2) = sum(geomderv%dthckdew(ew-1:ew,ns-1:ns)) / 16.0d0
          cons(3) = sum(geomderv%dusrfdns(ew-1:ew,ns-1:ns)) / 16.0d0
          cons(4) = sum(geomderv%dthckdns(ew-1:ew,ns-1:ns)) / 16.0d0
          cons(5) = sum(geomderv%stagthck(ew-1:ew,ns-1:ns))
          cons(6) = cons(5)*dns16
          cons(5) = cons(5)*dew16
          ! * better? (an alternative from TP's original code)
          !cons(5) = (thck(ew-1,ns)+2.0d0*thck(ew,ns)+thck(ew+1,ns)) * dew16
          !cons(6) = (thck(ew,ns-1)+2.0d0*thck(ew,ns)+thck(ew,ns+1)) * dns16

          velowk%suvel(:) = uvel(:,ew-1,ns-1) + uvel(:,ew-1,ns) + uvel(:,ew,ns-1) + uvel(:,ew,ns)
          velowk%svvel(:) = vvel(:,ew-1,ns-1) + vvel(:,ew-1,ns) + vvel(:,ew,ns-1) + vvel(:,ew,ns)

          ! Loop over each model level, starting from the bottom ----------------------

          do up = upn-1, 1, -1
            wvel(up,ew,ns) = wvel(up+1,ew,ns) &
                       - velowk%dupsw(up) * cons(5) * (sum(uvel(up:up+1,ew,ns-1:ns))  - sum(uvel(up:up+1,ew-1,ns-1:ns))) &
                       - velowk%dupsw(up) * cons(6) * (sum(vvel(up:up+1,ew-1:ew,ns))  - sum(vvel(up:up+1,ew-1:ew,ns-1))) &
                       - (velowk%suvel(up+1) - velowk%suvel(up)) * (cons(1) - velowk%depthw(up) * cons(2)) &
                       - (velowk%svvel(up+1) - velowk%svvel(up)) * (cons(3) - velowk%depthw(up) * cons(4)) 
          end do
        else 

          ! If there isn't enough ice, set velocities to zero ----------------------------

          wvel(:,ew,ns) = 0.0d0  

        end if
      end do
    end do

!TODO - Remove halo call? wvel is needed only for the old temperature code, which is not supported in parallel.
    call parallel_halo(wvel)
!    call horiz_bcs_unstag_scalar(wvel)

  end subroutine wvelintg

  subroutine wvel_ew(model)

    !*FD set periodic EW boundary conditions
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.

    model%velocity%wgrd(:,1,:)                  = model%velocity%wgrd(:,model%general%ewn-1,:)
    model%velocity%wgrd(:,model%general%ewn,:) = model%velocity%wgrd(:,2,:)
    model%velocity%wvel(:,1,:)                  = model%velocity%wvel(:,model%general%ewn-1,:)
    model%velocity%wvel(:,model%general%ewn,:) = model%velocity%wvel(:,2,:)

  end subroutine wvel_ew

!------------------------------------------------------------------------------------------
!TODO - Remove 'use parallel'?

  subroutine chckwvel(numerics,geomderv,uvel,vvel,wvel,thck,acab)

    !*FD Constrain the vertical velocity field to obey a kinematic upper boundary 
    !*FD condition.
    use parallel
    use glimmer_global, only : sp 
!!    use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_numerics),   intent(in)    :: numerics !*FD Numerical parameters of model
    type(glide_geomderv),   intent(in)    :: geomderv !*FD Temporal and horizontal derivatives
                                                        !*FD of thickness and upper ice surface
                                                        !*FD elevation.
    real(dp),dimension(:,:),  intent(in)    :: uvel     !*FD $x$ velocity field at top model
                                                        !*FD level (scaled, on staggered grid).
    real(dp),dimension(:,:),  intent(in)    :: vvel     !*FD $y$ velocity field at top model
                                                        !*FD level (scaled, on staggered grid).
    real(dp),dimension(:,:,:),intent(inout) :: wvel     !*FD Vertical velocity field, 
    real(dp),dimension(:,:),  intent(in)    :: thck     !*FD Ice thickness (scaled)
    real(sp),dimension(:,:),  intent(in)    :: acab     !*FD Mass-balance (scaled)

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp) :: wchk
    real(dp) :: tempcoef
    integer  :: ns,ew,nsn,ewn

    ! Get array sizes -------------------------------------------------------------------

    ewn=size(thck,1) ; nsn=size(thck,2)

    ! Allocate temporary work array -----------------------------------------------------


    ! Loop over all grid-boxes ----------------------------------------------------------

    do ns = 2,nsn-1
      do ew = 2,ewn-1
         if (thck(ew,ns) > numerics%thklim .and. wvel(1,ew,ns) /= 0) then

            wchk = geomderv%dusrfdtm(ew,ns) &
                 - acab(ew,ns) &
                 + (sum(uvel(ew-1:ew,ns-1:ns)) * sum(geomderv%dusrfdew(ew-1:ew,ns-1:ns)) &
                 +  sum(vvel(ew-1:ew,ns-1:ns)) * sum(geomderv%dusrfdns(ew-1:ew,ns-1:ns))) &
                 / 16.0d0

            
            tempcoef = wchk - wvel(1,ew,ns)

            wvel(:,ew,ns) = wvel(:,ew,ns) + tempcoef * (1.0d0 - numerics%sigma) 
         end if
      end do
    end do

!TODO - Remove halo call?  wvel is needed only for the old temperature code, which is not supported in parallel.
    call parallel_halo(wvel)
!    call horiz_bcs_unstag_scalar(wvel)

  end subroutine chckwvel

!------------------------------------------------------------------------------------------
! PRIVATE subroutines
!------------------------------------------------------------------------------------------

!TODO - Note: There is another copy of this function in glam_strs2.  
!       Maybe better to move this subroutine to another module to avoid duplication.
 
  function vertintg(velowk,in)

    !*FD Performs a depth integral using the trapezium rule.
    !*RV The value of in integrated over depth.


    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_velowk), intent(inout) :: velowk !*FD Work arrays and things for this module
    real(dp),dimension(:),intent(in)    :: in     !*FD Input array of vertical velocities (size = upn)
    real(dp) :: vertintg

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    integer :: up, upn

    ! Set up array of sigma intervals, if not done already ------------------------------

    upn=size(in)


    ! Do integration --------------------------------------------------------------------

    vertintg = 0.0d0

    do up = upn-1, 1, -1
      vertintg = vertintg + (in(up)+in(up+1)) * velowk%dups(up)                   
    end do

    vertintg = 0.5d0*vertintg

  end function vertintg

!------------------------------------------------------------------------------------------

  subroutine calc_btrc(model,flag,btrc)

    !*FD Calculate the value of $B$ used for basal sliding calculations.
    use glimmer_global, only : dp 
    use glimmer_physcon, only : rhoo, rhoi
    use glimmer_paramets, only : len0, thk0, scyr, vel0
    implicit none

    type(glide_global_type) :: model        !*FD model instance
    integer,                intent(in)    :: flag     !*FD Flag to select method of
    real(dp),dimension(:,:),intent(out)   :: btrc     !*FD Array of values of $B$.

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp) :: stagbwat, stagbmlt 
    integer :: ew,ns,nsn,ewn
    real :: Asl = 1.8d-10 !in units N^-3 yr^-1 m^8 for case(5)
    real :: Z !accounts for reduced basal traction due to pressure of
              !subglacial water for case(5)
    real :: tau !basal shear stress

    !scaling
    real :: tau_factor = 1.d-3*thk0*thk0/len0
    !real :: tau_factor = 1.0d0
    !------------------------------------------------------------------------------------

    ewn=model%general%ewn
    nsn=model%general%nsn

    !------------------------------------------------------------------------------------

    select case(flag)

    case(BTRC_CONSTANT)
       ! constant everywhere
       ! This option is used for EISMINT-2 experiment G
       btrc = model%velocity%bed_softness

    case(BTRC_CONSTANT_BWAT)
       ! constant where basal melt water is present, else = 0
       ! This option can be used for EISMINT-2 experiment H, provided that 
       ! basal water is present where T = Tpmp (e.g., BWATER_LOCAL)

       do ns = 1,nsn-1
          do ew = 1,ewn-1
             if (0.0d0 < model%temper%stagbwat(ew,ns)) then
                btrc(ew,ns) = model%velocity%bed_softness(ew,ns)
             else
                btrc(ew,ns) = 0.0d0
             end if
          end do
       end do

    case(BTRC_CONSTANT_TPMP)
       ! constant where basal temperature equal to pressure melting point, else = 0
       ! This is the actual condition for EISMINT-2 experiment H, which may not be 
       ! the same as case BTRC_CONSTANT_BWAT above, depending on the hydrology

       do ns = 1,nsn-1
          do ew = 1,ewn-1
             if (abs(model%temper%stagbpmp(ew,ns) - model%temper%stagbtemp(ew,ns))<0.001) then
                btrc(ew,ns) = model%velocity%bed_softness(ew,ns)
             else
                btrc(ew,ns) = 0.0d0
             end if
          end do
       end do

    case(BTRC_LINEAR_BMLT)
       ! linear function of basal melt rate

       do ns = 1,nsn-1
          do ew = 1,ewn-1
             stagbmlt = 0.25d0*sum(model%temper%bmlt(ew:ew+1,ns:ns+1))
             
             if (stagbmlt > 0.0d0) then
                btrc(ew,ns) = min(model%velowk%btrac_max, &
                                  model%velocity%bed_softness(ew,ns) + model%velowk%btrac_slope*stagbmlt)
             else
                btrc(ew,ns) = 0.0d0
             end if
          end do
       end do

    case(BTRC_TANH_BWAT)
       ! tanh function of basal water depth
       ! The 'velowk%c' parameters are derived above from the 5-part parameter bpar

       do ns = 1,nsn-1
          do ew = 1,ewn-1
             if (0.0d0 < model%temper%stagbwat(ew,ns)) then
               
                btrc(ew,ns) = model%velowk%c(1) + model%velowk%c(2) * tanh(model%velowk%c(3) * &
                     model%temper%stagbwat(ew,ns) - model%velowk%c(4))
                
                if (0.0d0 > sum(model%isostasy%relx(ew:ew+1,ns:ns+1))) then
                   btrc(ew,ns) = btrc(ew,ns) * model%velowk%marine  
                end if
             else
                btrc(ew,ns) = 0.0d0
             end if
          end do
       end do

!WHL - I'm not aware of anyone using this parameterization. Commented out for now.
!!    case(6)
!!       ! increases with the third power of the basal shear stress, from Huybrechts

!!       Asl = model%climate%slidconst
!!       do ns = 1, nsn-1
!!         do ew = 1, ewn-1
!TODO - Scaling looks wrong here: stagthck and thklim should have the same scaling.
!!           if ((model%geomderv%stagthck(ew,ns)*thk0) > model%numerics%thklim) then 
!!             if((model%geomderv%stagtopg(ew,ns)*thk0) > (model%climate%eus*thk0)) then
!!               Z = model%geomderv%stagthck(ew,ns)*thk0
!!             else
!!               Z = model%geomderv%stagthck(ew,ns)*thk0 + rhoi*((model%geomderv%stagtopg(ew,ns) *thk0 &
!!                   - model%climate%eus*thk0)/ rhoo)   
!!             end if 
              
!!            if(Z <= model%numerics%thklim) then !avoid division by zero
!!                Z = model%numerics%thklim
!!            end if 
            
!!             tau = ((tau_factor*model%stress%tau_x(ew,ns))**2 +&
!!             (model%stress%tau_y(ew,ns)*tau_factor)**2)**(0.5d0)
             
!!             btrc(ew,ns) = (Asl*(tau)**2)/Z !assuming that that btrc is later
!!                                             !multiplied again by the basal shear stress
       
!!           end if  
!!          end do
!!       end do

    case default   ! includes BTRC_ZERO
       ! zero everywhere
       ! This is used for EISMINT-2 experiments A to F
       btrc = 0.0d0

    end select

  end subroutine calc_btrc

!TODO - Remove this version of the subroutine?

#ifdef JEFFORIG
  subroutine calc_basal_shear(model)
    !*FD calculate basal shear stress: tau_{x,y} = -ro_i*g*H*d(H+h)/d{x,y}
    use glimmer_physcon, only : rhoi,grav
    implicit none
    type(glide_global_type) :: model        !*FD model instance


    model%velocity%tau_x = -rhoi*grav*model%geomderv%stagthck
    model%velocity%tau_y = model%velocity%tau_x * model%geomderv%dusrfdns
    model%velocity%tau_x = model%velocity%tau_x * model%geomderv%dusrfdew
  end subroutine calc_basal_shear
#endif

  subroutine calc_basal_shear(stagthck, dusrfdew, dusrfdns, tau_x, tau_y)

    ! calculate basal shear stress: tau_{x,y} = -ro_i*g*H*d(H+h)/d{x,y}
    use glimmer_physcon, only : rhoi,grav

    implicit none
    real(dp),dimension(:,:),intent(in) :: stagthck    !*FD Ice thickness (scaled)
    real(dp),dimension(:,:),intent(in) :: dusrfdew, dusrfdns
    real(dp),dimension(:,:),intent(out) :: tau_x
    real(dp),dimension(:,:),intent(out) :: tau_y

    tau_x(:,:) = -rhoi*grav*stagthck(:,:)
    tau_y(:,:) = tau_x * dusrfdns(:,:)
    tau_x(:,:) = tau_x * dusrfdew(:,:)

    !JEFF Are these replaced by the three lines above? They are not compiling.  7/28/11 
    ! model%stress%tau_x = -rhoi*grav*model%geomderv%stagthck
    ! model%stress%tau_y = model%stress%tau_x * model%geomderv%dusrfdns
    ! model%stress%tau_x = model%stress%tau_x * model%geomderv%dusrfdew

  end subroutine calc_basal_shear

!-------------------------------------------------------------------

end module glide_velo

!-------------------------------------------------------------------
