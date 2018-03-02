!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   verifBC.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

! verifBC
! implement test case B and C for verifying isothermal ice sheet models.
! This implementation is based on matlab code by Bueler et al 2005.
!
! Ed Bueler, Craig S. Lingle, Jed A. Kallen-Brown, David N. Covey, and Latrice
! N. Bowman, Exact solutions and the verification of numerical models for
! isothermal ice sheets, to appear in J. Glaciology

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module verifBC
  
  use glimmer_global, only : dp

  private :: dp, calc_h

  type verifBC_type
     real(dp) :: gamma
     real(dp) :: lambda = 0.   !*FD accumulation paramter ($M_\lambda = \lambda t^{-1} H$)
     real(dp) :: H0 = 3600.    !*FD central thickness at $t_0$ (m)
     real(dp) :: R0 = 750.     !*FD radius of magin at t0 (m)
     real(dp) :: tDel = 25000. !*FD time to advance from $t_0$ ($t_f = t_0+\Delta t$) (secs) if tDel<0 run from t=0 to $t_0$
     real(dp) :: t0
     real(dp) :: centre
     real(dp) :: alf,bet
     real(dp) :: gn
  end type verifBC_type

contains
  subroutine verifBC_config(section,veri)
    !*FD read configuration
    use glimmer_config
    implicit none
    type(verifBC_type)           :: veri    !*FD structure holding test setup
    type(ConfigSection), pointer :: section !*FD structure holding configuration section

    call GetValue(section,'lambda',veri%lambda)
    call GetValue(section,'H0',veri%H0)
    call GetValue(section,'R0',veri%R0)
    call GetValue(section,'tDel',veri%tDel)

  end subroutine verifBC_config

  subroutine verifBC_printconfig(veri)
    !*FD print configuration to log
    use glimmer_log
    implicit none
    type(verifBC_type)           :: veri    !*FD structure holding test setup
    ! local variables
    character(len=100) :: message
    call write_log('VerifBC test')
    call write_log('------------')
    write(message,*) 'accumulation parameter, lambda   :',veri%lambda
    call write_log(message)
    write(message,*) 'initial central ice thickness, H0:',veri%H0
    call write_log(message)
    write(message,*) 'initial ice margin radius, R0    :',veri%R0
    call write_log(message)
    write(message,*) 'time to advance from t0          :',veri%tDel
    call write_log(message)
    call write_log('')
  end subroutine verifBC_printconfig

  subroutine verifBC_init(model, veri)
    !*FD initialise test cases B and C
    use glide_types
    use glimmer_physcon
    use glimmer_paramets
    use glimmer_log
    implicit none
    
    type(glide_global_type) :: model !*FD model instance
    type(verifBC_type)      :: veri  !*FD structure holding test setup

    ! local variables
    real(dp) :: ts,tf
    real(dp) :: rmax,l
    integer :: num
    character(len=100) :: message

    call write_log('Setting up VerifBC test')
    call write_log('-----------------------')

    ! scale some parameters
    veri%R0 = veri%R0*1000.
    veri%tDel = veri%tDel*scyr
    veri%gn = real(gn,kind=dp)

    ! calculate Gamma
    veri%gamma = 2.*(rhoi*grav)**veri%gn*vis0/(veri%gn+2.)
    ! some parameters
    veri%alf = (2.-(veri%gn+1.)*veri%lambda)/(5.*veri%gn+3.)
    veri%bet = (1.+(2.*veri%gn+1.)*veri%lambda)/(5.*veri%gn+3.)

    ! time since creation
    veri%t0 = (veri%bet/veri%gamma) * ( (2.*veri%gn+1.)/(veri%gn+1.) )**veri%gn * (veri%R0**(veri%gn+1.)/veri%H0**(2.*veri%gn+1.))
    if (veri%tDel .le. 0) then
       ts = 0.
       tf = veri%t0
    else
       ts = veri%t0
       tf = veri%t0+veri%tDel
    end if

    ! set model start and end time
    model%numerics%tstart = ts/scyr
    model%numerics%tend = tf/scyr
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
    rmax = (tf/veri%t0)**(veri%bet)*veri%R0
    l = 1.1*rmax
    num = int(l/model%numerics%dns+0.5)
    veri%centre = num*model%numerics%dns
    num = 2*num+1

    model%general%ewn = num
    model%general%nsn = num
    write(message,*) 'Setting grid size to :',model%general%ewn
    call write_log(trim(message))

    call write_log('')
  end subroutine verifBC_init

  subroutine verifBC_update(model, veri, time, exact_h, mb)
    !*FD calculate exact ice thickness and mass balance
    use glide_types
    use glimmer_paramets, only: scyr, len0
    implicit none
    type(glide_global_type)   :: model !*FD model instance
    type(verifBC_type)        :: veri  !*FD structure holding test setup
    real(dp), intent(in) :: time  !*FD current time
    real(dp), dimension(:,:), intent(out) :: exact_h
    real(dp), dimension(:,:), intent(out) :: mb

    ! local variables
    real(dp) t,r,x,y
    integer i,j

    t = time*scyr

    ! calculate ice thickness
    do j=1,model%general%nsn
       y = (j-1)*model%numerics%dns*len0 - veri%centre
       do i=1,model%general%ewn
          x = (i-1)*model%numerics%dew*len0 - veri%centre
          r = sqrt(x**2+y**2)
          exact_h(i,j) = calc_h(veri,r,t)
       end do
    end do
    ! calculate mass balance
    if (t.eq.0.) then
       mb(:,:) = 0.
    else
       mb(:,:) = veri%lambda/t*exact_h(:,:)*scyr
    end if

  end subroutine verifBC_update

  function calc_h(veri,r,t)
    !*FD calculate exact ice thickness
    implicit none
    type(verifBC_type)   :: veri !*FD structure holding test setup
    real(dp), intent(in) :: r    !*FD radius (m)
    real(dp), intent(in) :: t    !*FD time   (s)

    real(dp) :: calc_h

    ! local variables
    real(dp) :: rscl, temp

    if (t.eq.0.) then
       calc_h = 0.
    else
       rscl = (veri%t0/t)**(veri%bet)*(r/veri%R0)
       temp = max(real(0.,kind=dp),1.-(rscl**((veri%gn+1.)/veri%gn)))
       calc_h = veri%H0*(t/veri%t0)**(-veri%alf) * temp**(veri%gn/(2.*veri%gn+1.))
    end if

  end function calc_h

end module verifBC
