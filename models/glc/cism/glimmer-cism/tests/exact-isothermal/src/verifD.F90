!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   verifD.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

! verifD
! implement test case D for verifying isothermal ice sheet models.
! This implementation is based on matlab code by Bueler et al 2005.
!
! Ed Bueler, Craig S. Lingle, Jed A. Kallen-Brown, David N. Covey, and Latrice
! N. Bowman, Exact solutions and the verification of numerical models for
! isothermal ice sheets, to appear in J. Glaciology

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module verifD
  
  use glimmer_global, only : rk, sp

  private
  public :: verifD_type, verifD_config, verifD_printconfig, verifD_init, verifD_update

  real(rk), parameter, private :: eps =  2.2204d-16

  type verifD_type
     real(rk) :: Cp = 200.   !*FD magnitude of perturbation (meters)
     real(rk) :: Tp = 5000.  !*FD period of perturbation (years)
     real(rk) :: H0 = 3600.  !*FD central thickness at $t_0$ (m)
     real(rk) :: R0 = 750.   !*FD margin radius (km)
     real(rk) :: tf = 25000. !*FD run from 0 to tf years
     real(rk) :: gamma
     real(rk) :: C,c1
     real(rk) :: gn
     real(rk) :: centre

     real(sp), dimension(:,:), pointer :: mb_steady => NULL() !*FD steady state mass balance
  end type verifD_type

contains
  subroutine verifD_config(section,veri)
    !*FD read configuration
    use glimmer_config
    implicit none
    type(verifD_type)            :: veri    !*FD structure holding test setup
    type(ConfigSection), pointer :: section !*FD structure holding configuration section

    call GetValue(section,'Cp',veri%Cp)
    call GetValue(section,'Tp',veri%Tp)
    call GetValue(section,'H0',veri%H0)
    call GetValue(section,'R0',veri%R0)
    call GetValue(section,'tf',veri%tf)
  end subroutine verifD_config

  subroutine verifD_printconfig(veri)
    !*FD print configuration to log
    use glimmer_log
    implicit none
    type(verifD_type)           :: veri    !*FD structure holding test setup
    ! local variables
    character(len=100) :: message
    call write_log('VerifD test')
    call write_log('------------')
    write(message,*) 'magnitude of perturbation        :',veri%Cp
    call write_log(message)
    write(message,*) 'period of perturbation           :',veri%Tp
    call write_log(message)
    write(message,*) 'initial central ice thickness, H0:',veri%H0
    call write_log(message)
    write(message,*) 'initial ice margin radius, R0    :',veri%R0
    call write_log(message)
    write(message,*) 'time to advance from 0           :',veri%tf
    call write_log(message)
    call write_log('')
  end subroutine verifD_printconfig
  
  subroutine verifD_init(model, veri)
    !*FD initialise test case D
    use glide_types
    use glimmer_physcon
    use glimmer_paramets
    use glimmer_log
    implicit none
    
    type(glide_global_type) :: model !*FD model instance
    type(verifD_type)      :: veri  !*FD structure holding test setup

    ! local variables
    integer num, i, j
    real(rk) :: x,y,r
    character(len=100) :: message

    call write_log('Setting up VerifD test')
    call write_log('-----------------------')

    ! scale some parameters
    veri%R0 = veri%R0*1000.
    veri%tf = veri%tf*scyr
    veri%Tp = veri%Tp*scyr
    veri%gn = real(gn,kind=rk)

    ! calculate Gamma
    veri%gamma = 2.*(rhoi*grav)**veri%gn*vis0/(veri%gn+2.)
    veri%C = (veri%gamma*veri%H0**(2.*veri%gn+2.))/(2.*veri%R0*(1.-1./veri%gn))**veri%gn

    veri%c1 = (3.*veri%H0)/(8.*(2./3.)**(3./8.))
    
    ! set model start and end time
    model%numerics%tstart = 0.
    model%numerics%tend = veri%tf/scyr
    write(message,*) 'Setting starting time to : ',model%numerics%tstart
    call write_log(trim(message))
    write(message,*) 'Setting end time to      : ',model%numerics%tend
    call write_log(trim(message))

    ! check delta x and y
    if (model%numerics%dew.ne.model%numerics%dns) then
       model%numerics%dew = max(model%numerics%dew,model%numerics%dns)
       model%numerics%dns = model%numerics%dew
       write(message,*) 'dx and dy do not match, setting dx=dy to :',model%numerics%dew
       call write_log(trim(message),GM_WARNING)
    end if

    ! set model dimensions
    num = int(1.2*veri%R0/model%numerics%dns+0.5)
    veri%centre = num*model%numerics%dns
    num = 2*num+1

    model%general%ewn = num
    model%general%nsn = num
    write(message,*) 'Setting grid size to :',model%general%ewn
    call write_log(trim(message))

    ! allocate memory
    allocate(veri%mb_steady(num,num))
    ! calculate steady state mb
    do j=1,model%general%nsn
       y = (j-1)*model%numerics%dns - veri%centre
       do i=1,model%general%ewn
          x = (i-1)*model%numerics%dew - veri%centre
          r = sqrt(x**2+y**2)
          veri%mb_steady(i,j) = calc_ms(veri,r)*scyr
       end do
    end do
    call write_log('')
  end subroutine verifD_init

  subroutine verifD_update(model, veri, time, exact_h, mb)
    !*FD calculate exact ice thickness and mass balance
    use glide_types
    use glimmer_paramets, only: scyr, len0
    implicit none
    type(glide_global_type)   :: model !*FD model instance
    type(verifD_type)        :: veri  !*FD structure holding test setup
    real(kind=rk), intent(in) :: time  !*FD current time
    real(kind=rk), dimension(:,:), intent(out) :: exact_h
    real(sp), dimension(:,:), intent(out) :: mb

    ! local variables
    real(rk) t,r,x,y
    integer i,j

    t = time*scyr

    do j=1,model%general%nsn
       y = (j-1)*model%numerics%dns*len0 - veri%centre
       do i=1,model%general%ewn
          x = (i-1)*model%numerics%dew*len0 - veri%centre
          r = sqrt(x**2+y**2)
          ! calculate ice thickness
          exact_h(i,j) = calc_hp(veri,r,t)

          ! calculate mass balance
          mb(i,j) = veri%mb_steady(i,j) + calc_mc(veri,r,t)*scyr
       end do
    end do
  end subroutine verifD_update

  function calc_hp(veri,r,t)
    !*FD calculate perturbed ice thickness
    use glimmer_physcon, only : pi
    implicit none
    type(verifD_type)    :: veri !*FD structure holding test setup
    real(rk), intent(in) :: r    !*FD radius (m)
    real(rk), intent(in) :: t    !*FD time

    real(rk) :: calc_hp

    calc_hp = calc_hs(veri,r) + veri%cp*sin(2.*pi*t/veri%tp)*getgp(veri,r)
  end function calc_hp

  function calc_hs(veri,r)
    !*FD calculate steady state ice thickness
    implicit none
    type(verifD_type)    :: veri !*FD structure holding test setup
    real(rk), intent(in) :: r    !*FD radius (m)

    real(rk) :: calc_hs

    ! local variables
    real(rk) chi, sr

    if (r.lt.veri%R0) then
       sr = r/veri%R0
       chi = (1.+1./veri%gn)*sr - 1./veri%gn + (1.-sr)**(1.+1./veri%gn) - sr**(1.+1./veri%gn)
       calc_hs = ( veri%H0/(1.-1./veri%gn)**(veri%gn/(2.*veri%gn+2.)) ) * chi**(veri%gn/(2.*veri%gn+2.))
    else
       calc_hs = 0.
    end if
  end function calc_hs

  function calc_mc(veri,r,t)
    !*FD massbalance 1
    use glimmer_physcon, only : pi
    implicit none
    type(verifD_type)    :: veri !*FD structure holding test setup
    real(rk), intent(in) :: r    !*FD radius (m)
    real(rk), intent(in) :: t    !*FD time   (s)

    real(rk) :: calc_mc

    real(rk) :: hp, divterms
    real(rk), dimension(2) :: ddhp

    if (r.ge.0.3*veri%R0 .and. r.le.0.9*veri%R0) then
       hp = calc_hp(veri,r,t)
       ddhp = getddr(veri,r,t)
       divterms = veri%gamma*hp**4.*ddhp(1)**2. * ((1./r)*hp*ddhp(1) + 5.*ddhp(1)**2. +3.*hp*ddhp(2) )
       calc_mc = (2.*pi*veri%cp/veri%tp)*cos(2.*pi*t/veri%tp)*getgp(veri,r)-calc_ms(veri,r)-divterms
    else
       calc_mc = 0.
    end if
  end function calc_mc

  function calc_ms(veri,r)
    !*FD steady state mass balance
    use glimmer_physcon, only : scyr
    implicit none
    type(verifD_type)    :: veri !*FD structure holding test setup
    real(rk), intent(in) :: r    !*FD radius (m)

    real(rk) :: calc_ms

    real temp, sr

    if (r.lt.veri%R0) then
       if (r.ge.5.*eps) then
          sr = r/veri%R0
          temp = sr**(1./veri%gn) + (1.-sr)**(1./veri%gn) - 1.
          calc_ms = veri%C/r * temp**(veri%gn-1.) * (2.*sr**(1./veri%gn)+(1.-sr)**(1./veri%gn-1.)*(1.-2.*sr) -1.)
       else
          calc_ms = 2.*veri%C/veri%R0
       end if
    else
       calc_ms = -0.1/scyr
    end if
  end function calc_ms

  function getgp(veri,r)
    !*FD helper routine
    use glimmer_physcon, only : pi
    implicit none
    type(verifD_type)    :: veri !*FD structure holding test setup
    real(rk), intent(in) :: r    !*FD radius (m)
    
    real(rk) :: getgp

    if (r.ge.0.3*veri%R0 .and. r.le.0.9*veri%R0) then
       getgp = .5*cos(pi*(r-0.6*veri%R0)/(0.3*veri%R0))+0.5
    else
       getgp = 0.
    end if
  end function getgp

  function getddr(veri,r,t)
    !*FD calculate derivatives
    use glimmer_physcon, only : pi
    implicit none
    type(verifD_type)    :: veri !*FD structure holding test setup
    real(rk), intent(in) :: r    !*FD radius (m)
    real(rk), intent(in) :: t    !*FD time   (s)

    real(rk), dimension(2) :: getddr

    real(rk) :: dgp, ddgp, chi, dchi, ddchi, dHs, ddHS
    real(rk) :: sr

    sr = r/veri%R0

    !gp   = getgp(veri,r)
    dgp  = (-pi/(0.6*veri%R0))*sin(pi*(r-0.6*veri%R0)/(0.3*veri%R0))
    ddgp = (-pi**2./(.18*veri%R0**2.))*cos(pi*(r-0.6*veri%R0)/(0.3*veri%R0))

    chi  = 4.*r/(3.*veri%R0) - 1./3. + (1.-sr)**(4./3.) - sr**(4./3.)
    dchi = (-4./(3.*veri%R0))*( sr**(1./3.) + (1.-sr)**(1./3) - 1.)
    ddchi= (-4./(9.*veri%R0**2.))*( sr**(-2./3.) - (1.-sr)**(-2./3.) )
    
    dHs  = veri%c1*chi**(-5./8.) * dchi
    ddHS = veri%c1*((-5./8.)*chi**(-13./8.)*dchi**2. + chi**(-5./8.)*ddchi )

    getddr(1) = dHs + veri%Cp*sin(2.*pi*t/veri%tp)*dgp
    getddr(2) = ddHs + veri%Cp*sin(2.*pi*t/veri%tp)*ddgp
  end function getddr
end module verifD
