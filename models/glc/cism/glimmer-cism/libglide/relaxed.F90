!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   relaxed.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#define NCI infile%nc
#define NCO outfile%nc

!TODO - This is a program.  Move to a utilities directory?

program relaxed

  ! utility to add relaxed bedrock topography to GLIMMER netcdf input files

  use glimmer_global
  use glimmer_ncdf
  use glimmer_ncinfile
  use glimmer_ncfile
  use glimmer_types
  use glimmer_setup
  use glimmer_cfproj
  use netcdf
  use glimmer_paramets, only: len0
  implicit none

  ! File-handling stuff

  integer unit, status
  type(glimmer_nc_input),  pointer :: infile
  type(glimmer_nc_output), pointer :: outfile
  type(glimmer_global_type) :: model
 
  ! Array

  real(dp),dimension(:,:),allocatable :: depr

  ! Densities of ice, ocean and mantle

  real(dp), parameter :: rhoi = 910.0d0  
  real(dp), parameter :: rhoo = 1028.0d0 
  real(dp), parameter :: rhom = 3300.0d0

  ! Other parameters

  integer :: nx,ny,i,j,dimsize
  logical :: flex = .true.

  ! Allocate pointer variables

  allocate(infile,outfile)

  ! Get some input parameters

  Print*,'RELAXED - utility to add relaxed bedrock topography'
  print*,'to GLIMMER netcdf input files'
  Print*,'Name of input file:'
  Read*,NCI%filename
  Print*,'Name of output file:'
  Read*,NCO%filename

  ! Open an input file and retrieve parameters
  status = nf90_open(NCI%filename,NF90_NOWRITE,NCI%id)
  ! getting ids
  status = nf90_inq_dimid(NCI%id, 'x1', NCI%x1dim)
  call nc_errorhandle(__FILE__,__LINE__,status)
  status = nf90_inq_dimid(NCI%id, 'y1', NCI%y1dim)
  call nc_errorhandle(__FILE__,__LINE__,status)
  ! getting dimensions
  status = nf90_inquire_dimension(NCI%id,NCI%x1dim,len=model%general%ewn)
  call nc_errorhandle(__FILE__,__LINE__,status)
  status = nf90_inquire_dimension(NCI%id,NCI%y1dim,len=model%general%nsn)
  call nc_errorhandle(__FILE__,__LINE__,status)
  ! close file
  status = nf90_close(NCI%id)

  ! Allocate necessary arrays

  model%general%upn=1
  call allocarr(model)
  allocate(depr(model%general%ewn,model%general%nsn))

  ! Reopen file

  call glimmer_nc_openfile(infile,model)

  ! Read data

  call glimmer_nc_read(infile,model,.false.)

  ! Get projection data

  model%projection=CFproj_GetProj(NCI%id)

  ! Calculate thickness

  where (model%climate%out_mask == 1.0)
     model%geometry%thck = max(0.,real(model%climate%presusrf - model%geometry%topg))
  elsewhere
     model%geometry%thck = 0.0
  end where

  ! Do flexure of some kind

  if (flex) then

     ! calculate present-day load (ice plus ocean)
     call flextopg(0,depr,model%geometry%thck,model%geometry%topg,rhoo,rhoi,rhom,model%numerics%dew,model%numerics%dns)

     ! rebound the bedrock
     model%geometry%relx = model%geometry%topg - depr

     ! calculate load from ocean
     call flextopg(1,depr,model%geometry%thck,model%geometry%relx,rhoo,rhoi,rhom,model%numerics%dew,model%numerics%dns)

     ! depress bedrock using this load
     model%geometry%relx = model%geometry%relx + depr

     ! Flatten masked areas
     
     where (model%climate%out_mask==0.0)
        model%geometry%relx=min(-1.0,model%geometry%relx)
     end where

  else

     where (model%geometry%thck > 0.0d0)
        model%geometry%relx = model%geometry%topg + model%geometry%thck * rhoi / rhom
     elsewhere
        model%geometry%relx = model%geometry%topg
     end where

     where (model%geometry%relx < 0.0d0 .and. model%geometry%thck > 0.0d0)
        model%geometry%relx = model%geometry%relx * rhom / (rhom - rhoo)
     end where

  end if


  ! Copy some data to the output structure

  NCO%do_var= NCI%do_var

  ! Create new file
  model%numerics%dew = model%numerics%dew / len0
  model%numerics%dns = model%numerics%dns / len0
  call glimmer_nc_createfile(outfile, model)

  status = nf90_put_var(NCO%id, NCO%varids(NC_B_TOPG), model%geometry%topg, (/1,1,1/))
  call nc_errorhandle(__FILE__,__LINE__,status)

  status = nf90_put_var(NCO%id, NCO%timevar,0.0,(/1/))
  call nc_errorhandle(__FILE__,__LINE__,status)

  if (NCI%do_var(NC_B_LAT)) then
     status = nf90_put_var(NCO%id, NCO%varids(NC_B_LAT), model%climate%lati, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NCI%do_var(NC_B_RELX)) then
     status = nf90_put_var(NCO%id, NCO%varids(NC_B_RELX), model%geometry%relx, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NCI%do_var(NC_B_USURF)) then
     status = nf90_put_var(NCO%id, NCO%varids(NC_B_USURF), model%geometry%usrf, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NCI%do_var(NC_B_PRESPRCP)) then
     status = nf90_put_var(NCO%id, NCO%varids(NC_B_PRESPRCP), model%climate%presprcp, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NCI%do_var(NC_B_PRESUSRF)) then
     status = nf90_put_var(NCO%id, NCO%varids(NC_B_PRESUSRF), model%climate%presusrf, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NCI%do_var(NC_B_MASK)) then
     status = nf90_put_var(NCO%id, NCO%varids(NC_B_MASK), model%climate%out_mask, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  status = nf90_close(NCO%id)

contains

  subroutine flextopg(flag,flex,thck,topg,rhoo,rhoi,rhom,dew,dns)

    implicit none

    ! Arguments

    integer,                  intent(in)  :: flag
    real(dp), dimension(:,:), intent(out) :: flex
    real(dp), dimension(:,:), intent(in)  :: thck
    real(dp), dimension(:,:), intent(in)  :: topg
    real(dp),                 intent(in)  :: rhoo
    real(dp),                 intent(in)  :: rhoi
    real(dp),                 intent(in)  :: rhom
    real(dp),                 intent(in)  :: dew
    real(dp),                 intent(in)  :: dns

    integer :: ew, ns,nsn,ewn

    real(dp), save :: thklim = 100.0d0   
    real(dp), parameter :: grav = 9.81 
    real(dp), parameter :: pi = 3.1416

    integer :: ewpt, nspt, ewflx, nsflx, ikelv
    integer, save :: nflx

    integer, parameter :: nkelv = 82

    ! ** constants used in calc

    ! ** Young's modulus (Nm^-2)
    ! ** thickness of lithosphere (m)
    ! ** radius of lithosphere (m)
    ! ** Poisson's ratio

    real(sp), parameter :: youngs = 8.35e10
    real(sp), parameter :: thklith = 110.0e3
    real(sp), parameter :: radlith = 6.244e6
    real(sp), parameter :: poiss = 0.25
    real(sp), parameter :: dkelv = 0.1 

    ! ** zero order kelvin function (for every dkelv from 0.0 to 8.0)

    real(sp), parameter, dimension(nkelv) :: &
         kelvin0 = (/ -0.785, -0.777, -0.758, -0.733, -0.704, &
         -0.672, -0.637, -0.602, -0.566, -0.531, &
         -0.495, -0.460, -0.426, -0.393, -0.362, &
         -0.331, -0.303, -0.275, -0.249, -0.225, &
         -0.202, -0.181, -0.161, -0.143, -0.126, &
         -0.111, -0.096, -0.083, -0.072, -0.061, &
         -0.051, -0.042, -0.035, -0.028, -0.021, &
         -0.016, -0.011, -0.007, -0.003,  0.000, & 
         0.002,  0.004,  0.006,  0.008,  0.009, & 
         0.010,  0.010,  0.011,  0.011,  0.011, & 
         0.011,  0.011,  0.011,  0.011,  0.010, & 
         0.010,  0.009,  0.009,  0.008,  0.008, & 
         0.007,  0.007,  0.006,  0.006,  0.005, & 
         0.005,  0.004,  0.004,  0.003,  0.003, & 
         0.003,  0.002,  0.002,  0.002,  0.002, & 
         0.001,  0.001,  0.001,  0.001,  0.001, & 
         0.000,  0.000 /)

    real(sp), dimension(:,:), allocatable, save :: dflct

    ! ** quantities calculated

    ! ** flexural rigidity
    ! ** radius of stiffness
    ! ** multiplier for loads

    real(sp) :: rigid, alpha, multi, dist, load

    logical, save :: first = .true.

    ! ******* CODE STARTS HERE ********

    nsn=size(flex,2) ; ewn=size(flex,1)

    if (first) then                                                  

       rigid = (youngs * thklith**3) / (12.0 * (1.0 - poiss**2))

       alpha = (rigid / ((youngs * thklith / radlith**2) + rhom * grav))**0.25

       multi = grav * dew**2 * alpha**2 / (2.0 * pi * rigid)

       nflx = 7 * int(alpha / dew) + 1

       allocate(dflct(nflx,nflx))

       do nsflx = 1,nflx
          do ewflx = 1,nflx

             dist = dew * sqrt(real(ewflx-1)**2 + real(nsflx-1)**2) / alpha

             ikelv = min(nkelv-1,int(dist/dkelv) + 1)

             dflct(ewflx,nsflx) = multi * &
                  (kelvin0(ikelv) + &
                  (dist - dkelv * (ikelv-1)) * &
                  (kelvin0(ikelv+1)- kelvin0(ikelv)) / dkelv)

          end do
       end do

       first = .false.

    end if

    flex = 0.0

    ! ** now loop through all of the points on the model grid

    do ns = 1,nsn
       do ew = 1,ewn

          ! ** for each point find the load it is imposing which
          ! ** depends on whether there it is ice covered or ocean
          ! ** covered

          ! flag == 0 
          ! calculate present-day load from ice and ocean water
          ! flag == 1
          ! calculate past load from ocean alone

          if (flag == 0) then

             load = rhoi * thck(ew,ns) 

          else 

             if (thck(ew,ns) > 0.0d0 .and. topg(ew,ns) < 0.0d0) then
                load = - rhoo * topg(ew,ns) 
             else
                load = 0.0d0
             end if
             
          end if

          ! ** now apply the calculated deflection field using the 
          ! ** ice/ocean load at that point and the function of
          ! ** how it will affect its neighbours

          ! ** the effect is linear so that we can sum the deflection
          ! ** from all imposed loads
          
          ! ** only do this if there is a load to impose and 
          ! ** be careful not to extend past grid domain

          do nsflx = max(1,ns-nflx+1), min(nsn,ns+nflx-1)
             do ewflx = max(1,ew-nflx+1), min(ewn,ew+nflx-1)

                ! ** find the correct function value to use
                ! ** the array dflct is one quadrant of the a square
                ! ** centered at the point imposing the load

                nspt = abs(ns - nsflx) + 1
                ewpt = abs(ew - ewflx) + 1

                flex(ewflx,nsflx) = load * dflct(ewpt,nspt) + flex(ewflx,nsflx)

             end do
          end do

        end do
    end do

  end subroutine flextopg

end program relaxed
