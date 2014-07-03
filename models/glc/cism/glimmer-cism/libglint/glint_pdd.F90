
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_pdd.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module glint_pdd

  !*FD The Glint annual poe degree day mass-balance scheme.
  !*FD Based on the original pdd mass-balance code from Tony's model.
  !*FD 
  !*FD {\bf N.B.} the module variables in this module are used for back-door 
  !*FD message passing, to make the integration of the PDD table look more 
  !*FD comprehensible, and avoid the need to have two customised copies of 
  !*FD the integration code.
  !*FD
  !*FD Note also that this code now deals in {\it unscaled} variables.
  !*FD 
  !*FD All precip and mass-balance amounts are as meters of water equivalent. PDD
  !*FD factors are no longer converted in the code.

  use glimmer_global, only : dp,sp

  implicit none

  private :: dp, sp

  type glint_pdd_params

    !*FD Holds parameters for positive-degree-day mass-balance
    !*FD calculation. The table has two axes - the $x$ axis is the
    !*FD difference between mean annual and July temps, while the
    !*FD $y$- axis is the mean annual temp

    integer  :: dx        = 1   !*FD Spacing of values in x-direction ($^{\circ}$C)
    integer  :: dy        = 1   !*FD Spacing of values in y-direction ($^{\circ}$C)
    integer  :: ix        = 0   !*FD Lower bound of $x$-axis ($^{\circ}$C)
    integer  :: iy        = -50 !*FD Lower bound of $y$-axis ($^{\circ}$C)
    integer  :: nx        = 31  !*FD Number of values in x-direction
    integer  :: ny        = 71  !*FD Number of values in y-direction
    real(sp) :: dailytemp = 0.0 
    real(sp) :: tma       = 0.0
    real(sp) :: tmj       = 0.0
    real(sp) :: dtmj      = 0.0
    real(sp) :: dd_sigma  = 5.0 !*FD Standard deviation of daily temperature (K)
 
    ! The actual PDD table ---------------------------------------------

    real(sp),dimension(:,:),pointer :: pddtab  => null() 
    
    !*FD PDD table - must be allocated with dimensions nx,ny.

    ! Parameters for the PDD calculation

    real(sp) :: wmax        = 0.6   !*FD Fraction of melted snow that refreezes
    real(dp) :: pddfac_ice  = 0.008 !*FD PDD factor for ice (m water day$^{-1}$ $^{\circ}C$^{-1}$)
    real(dp) :: pddfac_snow = 0.003 !*FD PDD factor for snow (m water day$^{-1}$ $^{\circ}C$^{-1}$)

  end type glint_pdd_params
 	 
  ! Module parameters use for back-door message-passing

  real(sp) :: dd_sigma            !*FD The value of $\sigma$ in the PDD integral
  real(sp) :: t_a_prime           !*FD The value of $T'_{a}$ in the PDD integral
  real(sp) :: mean_annual_temp    !*FD Mean annual temperature
  real(sp) :: mean_july_temp      !*FD Mean july temperature

  private
  public :: glint_pdd_params, glint_pdd_init, glint_pdd_mbal

contains

!-------------------------------------------------------------------------------
! PUBLIC subroutines
!-------------------------------------------------------------------------------

  subroutine glint_pdd_init(params,config)

    use glimmer_config

    type(glint_pdd_params),intent(inout) :: params !*FD The positive-degree-day parameters
    type(ConfigSection), pointer         :: config !*FD structure holding sections of configuration file   

    ! Read the config file and output to log

    call pdd_readconfig(params,config)
    call pdd_printconfig(params)

    ! Deallocate arrays

    if (associated(params%pddtab)) deallocate(params%pddtab)
    
    ! Allocate pdd table

    allocate(params%pddtab(params%nx,params%ny))
    call pddtabgrn(params)  

  end subroutine glint_pdd_init

!-------------------------------------------------------------------------------

  subroutine glint_pdd_mbal(params,artm,arng,prcp,ablt,acab,landsea)

    !*FD Calculates mass-balance over the ice model domain, by the
    !*FD positive-degree-day method.

    implicit none 
 
    type(glint_pdd_params),   intent(inout) :: params  !*FD The positive-degree-day parameters
    real(sp), dimension(:,:), intent(in)    :: artm    !*FD Annual mean air-temperature 
                                                       !*FD ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: arng    !*FD Annual temperature half-range ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: prcp    !*FD Annual accumulated precipitation 
                                                       !*FD (m water equivalent)
    real(sp), dimension(:,:), intent(out)   :: ablt    !*FD Annual ablation (m water equivalent)
    real(sp), dimension(:,:), intent(out)   :: acab    !*FD Annual mass-balance (m water equivalent)
    logical,  dimension(:,:), intent(in)    :: landsea !*FD Land-sea mask (land is TRUE)

    ! Internal variables

    real(sp) :: wfrac, pablt, tx, ty, pdd
    integer  :: ns,ew,nsn,ewn,kx,ky,jx,jy

    ! Get size of arrays. All arrays should be the same size as this.

    ewn=size(artm,1) ; nsn=size(artm,2)

    !-----------------------------------------------------------------------
    ! Main loop
    !-----------------------------------------------------------------------

    do ns = 1, nsn
      do ew = 1, ewn 

         if (landsea(ew,ns)) then
          ! Only calculate mass-balance over 'land'
          ! Find the no. of pdd from the mean annual temp and its range

          ky = int((artm(ew,ns)-params%iy)/params%dy)
          kx = int((arng(ew,ns)-params%ix)/params%dx) 
    
          ! Check to see if indicies are in range

          if ( kx < 0 ) then 
            tx = 0
            jx = 2
            kx = 1
          else if ( kx > params%nx-2 ) then
            tx = 1.0
            jx = params%nx
            kx = params%nx-1
          else
            tx = arng(ew,ns) - kx * params%dx - params%ix
            jx = kx + 2
            kx = kx + 1
          end if

          if ( ky < 0 ) then 
            ty = 0.0
            jy = 2
            ky = 1
          else if ( ky > params%ny-2 ) then
            ty = 1.0
            jy = params%ny
            ky = params%ny-1
          else
            ty = artm(ew,ns) - ky * params%dy - params%iy;
            jy = ky + 2
            ky = ky + 1
          end if
            
          ! this is done using a look-up table constructed earlier

          pdd = params%pddtab(kx,ky)*(1.0-tx)*(1.0-ty) + &
                params%pddtab(jx,ky) * tx * (1.0 - ty) + &
                params%pddtab(jx,jy) * tx * ty +         &
                params%pddtab(kx,jy) * (1.0 - tx) * ty

          ! this is the depth of superimposed ice that would need to be
          ! melted before runoff can occur

          wfrac = params%wmax * prcp(ew,ns)

          ! this is the total potential ablation of SNOW
    
          pablt = pdd * params%pddfac_snow

          ! if the total snow ablation is less than the depth of 
          ! superimposed ice - no runoff occurs

          ! else if the total snow ablation is more than the depth
          ! of superimposed ice BUT less than the total amount of
          ! prcpitation - runoff occurs (at a rate equal to the
          ! total potential snowmelt minus that which forms superimposed ice)

          ! else if the total snow ablation is more than the amount
          ! of prcpitation - all snow that is not superimposed ice is lost 
          ! and the potential ablation not used on snow is used on ice
          ! (including the superimposed ice)

          ! there is a change in the pddfi term, replaced wfrac with prcp
          ! error spotted by jonathan 18-04-00

          if ( pablt <= wfrac ) then
            ablt(ew,ns) = 0.0
          else if(pablt > wfrac .and.pablt <= prcp(ew,ns)) then   
            ablt(ew,ns) = pablt - wfrac 
          else
            ablt(ew,ns) = prcp(ew,ns) - wfrac + params%pddfac_ice*(pdd-prcp(ew,ns)/params%pddfac_snow) 
          end if

          ! Finally, mass-balance is difference between accumulation and
          ! ablation.

          acab(ew,ns) = prcp(ew,ns) - ablt(ew,ns)
         else
            ablt(ew,ns)=prcp(ew,ns)
            acab(ew,ns)=0.0
         endif
      end do
    end do

  end subroutine glint_pdd_mbal

!-------------------------------------------------------------------------------
! PRIVATE subroutines and functions
!-------------------------------------------------------------------------------

  subroutine pdd_readconfig(params,config)

    !*FD Reads in configuration data for the annual PDD scheme.

    use glimmer_config

    type(glint_pdd_params),intent(inout) :: params !*FD The positive-degree-day parameters
    type(ConfigSection), pointer         :: config !*FD structure holding sections of configuration file   

    ! local variables
    type(ConfigSection), pointer :: section
    
    call GetSection(config,section,'GLINT annual pdd')
    if (associated(section)) then
       call GetValue(section,'dx',params%dx)
       call GetValue(section,'dy',params%dy)
       call GetValue(section,'ix',params%ix)
       call GetValue(section,'iy',params%iy)
       call GetValue(section,'nx',params%nx)
       call GetValue(section,'ny',params%ny)
       call GetValue(section,'wmax',params%wmax)
       call GetValue(section,'pddfac_ice',params%pddfac_ice)
       call GetValue(section,'pddfac_snow',params%pddfac_snow)
       call GetValue(section,'dd_sigma',params%dd_sigma)
    end if

  end subroutine pdd_readconfig

!-------------------------------------------------------------------------------

  subroutine pdd_printconfig(params)

    use glimmer_log

    type(glint_pdd_params),intent(inout) :: params !*FD The positive-degree-day parameters
    character(len=100) :: message

    call write_log_div

    call write_log('GLINT annual PDD scheme parameters:')
    call write_log('-----------------------------------')
    write(message,*) 'x-spacing of pdd table',params%dx,' degC'
    call write_log(message)
    write(message,*) 'y-spacing of pdd table',params%dy,' degC'
    call write_log(message)
    write(message,*) 'Lower bound of x-axis',params%ix,' degC'
    call write_log(message)
    write(message,*) 'Lower bound of y-axis',params%iy,' degC'
    call write_log(message)
    write(message,*) 'Number of points in x',params%nx
    call write_log(message)
    write(message,*) 'Number of points in y',params%ny
    call write_log(message)
    write(message,*) 'Snow refreezing fraction',params%wmax
    call write_log(message)
    write(message,*) 'PDD factor for ice',params%pddfac_ice
    call write_log(message)
    write(message,*) 'PDD factor for snow',params%pddfac_snow
    call write_log(message)
    write(message,*) 'Standard deviation of temperature cycle',params%dd_sigma,' degC'
    call write_log(message)
    call write_log('')

  end subroutine pdd_printconfig

!-------------------------------------------------------------------------------
        
  subroutine pddtabgrn(params)

    !*FD Initialises the positive-degree-day-table.

    use glimmer_global, only: sp
    use glimmer_integrate
    use glimmer_log

    implicit none

    type(glint_pdd_params),intent(inout) :: params !*FD PDD parameters

    ! Internal variables

    real(sp)           :: tma,dtmj
    real(sp),parameter :: twopi = 3.1416 * 2.0   !TODO - Change to 2.0 * pi 
    integer  :: kx,ky, i,j

    !--------------------------------------------------------------------
    ! Main loops:
    !  tma  -- the mean annual temperature (y-axis of pdd table)
    !  dtmj -- difference from mean july temperature (x-axis of table)
    !  tmj -- the actual july temperature
    !--------------------------------------------------------------------

    call write_log('Calculating PDD table...',GM_DIAGNOSTIC)

    do j=0,params%ny-1
       tma=params%iy + j*params%dy

       ky = findgrid(tma,real(params%iy),real(params%dy))

       do i=0,params%nx-1
          dtmj = params%ix + i*params%dx

          mean_july_temp = tma + dtmj   
          kx  = findgrid(dtmj,real(params%ix),real(params%dx)) 

          ! need these lines to take account of the backdoor message passing used here

          mean_annual_temp=tma
          dd_sigma=params%dd_sigma

          params%pddtab(kx,ky)=(1.0/(dd_sigma*sqrt(twopi)))*romberg_int(inner_integral,0.0,twopi)

          ! convert to days     

          params%pddtab(kx,ky) = 365.0 * params%pddtab(kx,ky) / twopi

       end do
    end do

    call write_log('   ...done.',GM_DIAGNOSTIC)

  end subroutine pddtabgrn

!-------------------------------------------------------------------------------

  real(sp) function inner_integral(day)

    !*FD Calculates the value of the inner integral, i.e.
    !*FD \begin{equation}
    !*FD \int^{T_{a}'+2.5\sigma}_{0}T_{a}\times
    !*FD \exp\left(\frac{-(T_a-T_{a}')^2}{2\sigma^2}\right)\,dT
    !*FD \end{equation}
    use glimmer_integrate

    implicit none

    real(sp), intent(in) :: day !*FD The `day', in radians, so that a year is $2\pi$ long.

    real(sp) :: upper_limit

    t_a_prime=mean_annual_temp+(mean_july_temp-mean_annual_temp)*cos(day)

    upper_limit=t_a_prime+2.5*dd_sigma

    if (upper_limit<=0.0) then
      inner_integral=0.0
    else
      inner_integral=romberg_int(pdd_integrand,0.0,upper_limit)
    endif

  end function inner_integral

!-------------------------------------------------------------------------------
        
  real(sp) function pdd_integrand(artm)

    !*FD The expression to be integrated in the calculation of the PDD table. The whole
    !*FD integral is:
    !*FD \begin{equation}
    !*FD D=\frac{1}{\sigma\sqrt{2\pi}}\int^{A}_{0}\int^{T_{a}'+2.5\sigma}_{0}T_{a}\times
    !*FD \exp\left(\frac{-(T_a-T_{a}')^2}{2\sigma^2}\right)\,dTdt
    !*FD \end{equation}

    implicit none

    real(sp), intent(in) :: artm      !*FD The annual mean air temperature (degC)

     pdd_integrand = artm *  exp(- (artm - t_a_prime)**2 / (2.0 * dd_sigma**2))

  end function pdd_integrand

!-------------------------------------------------------------------------------

  integer function findgrid(rin,init,step)

    !*FD Calculates which row or column of the pdd table corresponds
    !*FD to a given value on the appropriate axis, so that:
    !*FD \[
    !*FD \mathtt{findgrid}=\frac{\mathtt{rin}-\mathtt{init}}{\mathtt{step}+1}
    !*FD \] 
    !*RV The relevant array index.

    use glimmer_global, only : sp
    
    implicit none
    
    real(sp), intent(in) :: rin  !*FD Value of axis variable at current point.
    real(sp), intent(in) :: init !*FD Value of axis variable at first point.
    real(sp), intent(in) :: step !*FD Grid spacing.
    
    findgrid = (rin - init) / step + 1

  end function findgrid
 
end module glint_pdd
