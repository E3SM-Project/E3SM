!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   simple_forcing.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module simple_forcing
  !*FD read configuration and generate simple massbalance and 
  !*FD temperature fields

  use glimmer_global, only : sp

  type simple_climate
     ! holds parameters for the simple climate

     integer :: eismint_type = 0
     !*FD select EISMINT experiment
     !*FD \begin{description}
     !*FD \item[{\bf 1}] EISMINT-1 fixed margin
     !*FD \item[{\bf 2}] EISMINT-1 moving margin
     !*FD \end{description}
     real(kind=sp), dimension(2) :: airt = (/ -3.150, 1.e-2 /)  
     !*FD air temperature parameterisation K, K km$^{-3}$
     real(kind=sp), dimension(3) :: nmsb = (/ 0.5, 1.05e-5, 450.0e3 /)
     !*FD mass balance parameterisation m yr$^{-1}$, yr$^{-1}$, m
     real(kind=sp) :: period = 0.
     !*FD EISMINT time-dep climate forcing period, switched off when set to 0
     real(kind=sp) :: mb_amplitude = 0.2
     !*FD EISMINT amplitude of mass balance time-dep climate forcing
  end type simple_climate

  !MAKE_RESTART
#ifdef RESTARTS
#define RST_SIMPLE_FORCING
!JCC - no restarts yet
!#include "glimmer_rst_head.inc"
#undef RST_SIMPLE_FORCING
#endif

contains

#ifdef RESTARTS
#define RST_SIMPLE_FORCING
!JCC - no restarts yet
!#include "glimmer_rst_body.inc"
#undef RST_SIMPLE_FORCING
#endif

  subroutine simple_initialise(climate,config)
    !*FD initialise simple climate model
    use glimmer_paramets, only: thk0, acc0, scyr
    use glimmer_config
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   

  
    call simple_readconfig(climate,config)
    call simple_printconfig(climate)

    ! scale parameters
    select case(climate%eismint_type)
    case(1)
       climate%nmsb(1) = climate%nmsb(1) / (acc0 * scyr)
    case(2)
       climate%airt(2) = climate%airt(2) * thk0
       climate%nmsb(1) = climate%nmsb(1) / (acc0 * scyr)
       climate%nmsb(2) = climate%nmsb(2) / (acc0 * scyr)
    case(3)
       climate%nmsb(1) = climate%nmsb(1) / (acc0 * scyr)
       climate%nmsb(2) = climate%nmsb(2) / (acc0 * scyr)
    case(4)
       climate%nmsb(1) = climate%nmsb(1) / (acc0 * scyr)
    end select
       
  end subroutine simple_initialise

  subroutine simple_readconfig(climate, config)
    !*FD read configuration
    use glimmer_log
    use glimmer_config
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   

    ! local variables
    type(ConfigSection), pointer :: section
    real(kind=sp), dimension(:), pointer :: dummy

    call GetSection(config,section,'EISMINT-1 fixed margin')
    if (associated(section)) then
       climate%eismint_type = 1
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       climate%airt = (/-34.15, 8.e-8/)
       if (associated(dummy)) then
          climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,1)
       climate%nmsb = (/0.3, 0.0, 0.0/)
       if (associated(dummy)) then
          climate%nmsb(1) = dummy(1)
       end if
       call GetValue(section,'period',climate%period)
       call GetValue(section,'mb_amplitude',climate%mb_amplitude)
       return       
    end if
    call GetSection(config,section,'EISMINT-1 moving margin')
    if (associated(section)) then
       climate%eismint_type = 2
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
          climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
          climate%nmsb = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'period',climate%period)
       climate%mb_amplitude = 100000.
       call GetValue(section,'mb_amplitude',climate%mb_amplitude)
       return
    end if
    call GetSection(config,section,'EISMINT-2')
    if (associated(section)) then
       climate%eismint_type = 3
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
          climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       else
          climate%airt = (/-35., 1.67e-5/)
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
          climate%nmsb = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       return
    end if
    
    ! *SFP* standard higher-order tests
    call GetSection(config,section,'DOME-TEST')
    if (associated(section)) then
        return
    end if

    call GetSection(config,section,'ISMIP-HOM-TEST')
    if (associated(section)) then
        return
    end if

    call GetSection(config,section,'SHELF-TEST') 
    if (associated(section)) then 
        return
    end if

    call GetSection(config,section,'STREAM-TEST') 
    if (associated(section)) then 
        return
    end if

    call GetSection(config,section,'ROSS-TEST') 
    if (associated(section)) then 
        return
    end if

    !mismip tests
    call GetSection(config,section,'MISMIP-1')
    if (associated(section)) then
       climate%eismint_type = 4
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
           climate%airt = dummy
           deallocate(dummy)
           dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
           climate%nmsb = dummy
           deallocate(dummy)
           dummy=>NULL()
       end if
       return
    end if

    !exact verification
    call GetSection(config,section,'EXACT')
    if (associated(section)) then
       climate%eismint_type = 5
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
           climate%airt = dummy
           deallocate(dummy)
           dummy=>NULL()
       end if
       return
    end if

    call write_log('No EISMINT forcing selected',GM_FATAL)
  end subroutine simple_readconfig

  subroutine simple_printconfig(climate)
    !*FD print simple climate configuration
    use glimmer_log
    implicit none
    type(simple_climate) :: climate   !*FD structure holding climate info
    character(len=100) :: message

    call write_log_div
    select case(climate%eismint_type)
    case(1)
       call write_log('EISMINT-1 fixed margin configuration')
       call write_log('------------------------------------')
       write(message,*) 'temperature  : ',climate%airt(1)
       call write_log(message)
       write(message,*) '               ',climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',climate%nmsb(1)
       call write_log(message)
       write(message,*) 'period       : ',climate%period
       call write_log(message)
       if (climate%period .gt.0) then
          write(message,*) 'mb amplitude : ',climate%mb_amplitude
          call write_log(message)
       end if
    case(2)
       call write_log('EISMINT-1 moving margin configuration')
       call write_log('-------------------------------------')
       write(message,*) 'temperature  : ',climate%airt(1)
       call write_log(message)
       write(message,*) '               ',climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',climate%nmsb(1)
       call write_log(message)
       write(message,*) '               ',climate%nmsb(2)
       call write_log(message)
       write(message,*) '               ',climate%nmsb(3)
       call write_log(message)
       write(message,*) 'period       : ',climate%period
       call write_log(message)
       if (climate%period .gt.0) then
          write(message,*) 'mb amplitude : ',climate%mb_amplitude
          call write_log(message)
       end if
    case(3)
       call write_log('EISMINT-2')
       call write_log('---------')
       write(message,*) 'temperature  : ',climate%airt(1)
       call write_log(message)
       write(message,*) '               ',climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',climate%nmsb(1)
       call write_log(message)
       write(message,*) '               ',climate%nmsb(2)
       call write_log(message)
       write(message,*) '               ',climate%nmsb(3)
       call write_log(message)
    end select
    call write_log('')
  end subroutine simple_printconfig

  subroutine simple_massbalance(climate,model,time)
    use parallel
    !*FD calculate simple mass balance
    use glimmer_global, only:rk
    use glide_types
    use glimmer_paramets, only : len0, acc0, scyr
    use glimmer_physcon, only : pi
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(glide_global_type) :: model        !*FD model instance
    real(kind=rk), intent(in) :: time                !*FD current time

    ! local variables
    integer  :: ns,ew
    real :: dist, ewct, nsct, grid, rel
    real :: periodic_bc = 1.

    ewct = real(model%general%ewn+1) / 2.0
    nsct = real(model%general%nsn+1) / 2.0
    grid = model%numerics%dew * len0

    if (model%options%periodic_ew) then
        periodic_bc = 0
    else
        periodic_bc = 1
    end if

    select case(climate%eismint_type)
    case(1)
       ! EISMINT-1 fixed margin
       model%climate%acab(:,:) = climate%nmsb(1)
       if (climate%period.ne.0) then
          model%climate%acab(:,:) = model%climate%acab(:,:) + climate%mb_amplitude * sin(2.*pi*time/climate%period)/ (acc0 * scyr)
       end if
    case(2)
       ! EISMINT-1 moving margin       
       if (climate%period.ne.0) then
          rel = climate%nmsb(3) + climate%mb_amplitude*sin(2.*pi*time/climate%period)
       else
          rel = climate%nmsb(3)
       end if

       call not_parallel(__FILE__,__LINE__)
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic_bc*(real(ew) - ewct)**2 + (real(ns) - nsct)**2)
             model%climate%acab(ew,ns) = min(climate%nmsb(1), climate%nmsb(2) * (rel - dist))
          end do
       end do
    case(3)
       ! EISMINT-2
       rel = climate%nmsb(3)

       call not_parallel(__FILE__,__LINE__)
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic_bc*(real(ew) - ewct)**2 + (real(ns) - nsct)**2)
             model%climate%acab(ew,ns) = min(climate%nmsb(1), climate%nmsb(2) * (rel - dist))
          end do
       end do
    case(4)
       !mismip 1
       model%climate%acab = climate%nmsb(1)
    case(5)
       !verification 
       call not_parallel(__FILE__,__LINE__)
       call exact_surfmass(climate,model,time,1.0,climate%airt(2))
    end select
  end subroutine simple_massbalance

  subroutine simple_surftemp(climate,model,time)
    !*FD calculate simple air surface temperature
    use parallel
    use glide_types
    use glimmer_global, only:rk
    use glimmer_paramets, only : len0
    use glimmer_physcon, only : pi
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(glide_global_type) :: model        !*FD model instance
    real(kind=rk), intent(in) :: time                !*FD current time

    ! local variables
    integer  :: ns,ew
    real :: dist, ewct, nsct, grid
    real :: periodic_bc = 1.

    ewct = real(model%general%ewn+1) / 2.0
    nsct = real(model%general%nsn+1) / 2.0
    grid = model%numerics%dew * len0

    if (model%options%periodic_ew) then
        periodic_bc = 0
    else
        periodic_bc = 1
    end if

    select case(climate%eismint_type)
    case(1)
       call not_parallel(__FILE__,__LINE__)
       ! EISMINT-1 fixed margin
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * max(periodic_bc*abs(real(ew) - ewct),abs(real(ns) - nsct))*1e-3
             model%climate%artm(ew,ns) = climate%airt(1) + climate%airt(2) * dist*dist*dist
          end do
       end do
       if (climate%period.ne.0) then
          model%climate%artm(:,:) = model%climate%artm(:,:) + 10.*sin(2.*pi*time/climate%period)
       end if
    case(2)
       ! EISMINT-1 moving margin
       model%climate%artm(:,:) = climate%airt(1) - model%geometry%thck(:,:) * climate%airt(2)
       if (climate%period.ne.0) then
          model%climate%artm(:,:) = model%climate%artm(:,:) + 10.*sin(2.*pi*time/climate%period)
       end if
    case(3)
       call not_parallel(__FILE__,__LINE__)
       ! EISMINT-2
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic_bc*(real(ew) - ewct)**2 + (real(ns) - nsct)**2)
             model%climate%artm(ew,ns) = climate%airt(1)+climate%airt(2) * dist
          end do
       end do
    case(4)
       model%climate%artm = climate%airt(1)
    case(5)
       call not_parallel(__FILE__,__LINE__)
       !call both massbalance and surftemp at the same time to save computing time. 
       call exact_surfmass(climate,model,time,0.0,climate%airt(2))
    end select
  end subroutine simple_surftemp
  
  !which_call - simple_surftemp(0)/simple_massbalance(1)/both(2)
  !which_test - test f(0)/test g(1)/exact(2)
  subroutine exact_surfmass(climate,model,time,which_call,which_test)
    use glide_types
    use testsFG
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(glide_global_type) :: model        !*FD model instance
    real(kind=rk), intent(in) :: time                !*FD current time
    real(sp), intent(in) :: which_test                !*FD  Which exact test (F=0,G=1)
    real(sp), intent(in) :: which_call                !*FD  0 = surface temp, 1= mass balance
    integer  :: ns,ew,lev,center
    !verification
    real(dp) ::  t, r, z, x, y                       !in variables
    real(dp) ::  H, TT, U, w, Sig, M, Sigc        !out variables
    real(dp) :: H_0
    center = (model%general%ewn - 1)*.5
    if (which_call .eq. 0.0 .or. which_call .eq. 2.0) then
        !point by point call to the function 
        do ns = 1,model%general%nsn
            do ew = 1,model%general%ewn
                x = (ew - center)*model%numerics%dew
                y = (ns - center)*model%numerics%dns
                r = sqrt(x**2 + y**2)
                do lev = 1, model%general%upn
                    z = model%geometry%thck(ew,ns)*model%numerics%sigma(lev)
                    !the function only returns values within the radius
                    if(r>0.0 .and. r<L) then
                        if (which_test .eq. 0.0) then
                            call testF(r,z,H,TT,U,w,Sig,M,Sigc)
                        else if (which_test .eq. 1.0) then
                            call testG(time,r,z,H,TT,U,w,Sig,M,Sigc)
                        else if (which_test .eq. 2.0) then
                            !H_0 = H0
                            H_0 = model%geometry%thck(center,center)
                            !TT = 0.0
                            TT = model%tempwk%dissip(lev,ew,ns)
                            call model_exact(time,r,z,model%geometry%thck(ew,ns),H_0,TT,U,w,Sig,M,Sigc)
                        end if
                        model%tempwk%compheat(lev,ew,ns) = -Sigc*SperA*1.0e6 ! (10^(-3) deg mK/a) (*10^3)
                    else
                        model%tempwk%compheat(lev,ew,ns) = 0.0
                    end if
                end do
            end do
        end do
    else if (which_call .eq. 1.0 .or. which_call .eq. 2.0) then
        do ns = 1,model%general%nsn
            do ew = 1,model%general%ewn
                x = (ew - center)*model%numerics%dew
                y = (ns - center)*model%numerics%dns
                r = sqrt(x**2 + y**2)
                z = model%geometry%thck(ew,ns)
                !the function only returns values within the radius
                if(r>0.0 .and. r<L) then
                    if (which_test .eq. 0.0) then
                        call testF(r,z,H,TT,U,w,Sig,M,Sigc)
                    else if (which_test .eq. 1.0) then
                        call testG(time,r,z,H,TT,U,w,Sig,M,Sigc)
                    else if (which_test .eq. 2.0) then
                        H_0 = H0 !H_0 = model%geometry%thck(center,center)
                        TT = 0.0 !TT = model%tempwk%dissip(lev,ew,ns)
                        call model_exact(time,r,z,model%geometry%thck(ew,ns),H_0,TT,U,w,Sig,M,Sigc)
                    end if
                    model%climate%acab(ew,ns) = M*SperA  !m/a
                else
                    if(r .eq. 0.0) then
                        model%climate%acab(ew,ns) = 0.0
                        !model%climate%acab(ew,ns) = H0 - model%geometry%thck(center,center) !set it back to H0
                    else
                        model%climate%acab(ew,ns) = -0.1  !outside the glacier we use .1 m ablation as specified in the paper
                    end if
                end if
            end do
        end do
    end if
  end subroutine exact_surfmass
  
end module simple_forcing
