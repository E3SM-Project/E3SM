!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   eismint_forcing.F90 - part of the Community Ice Sheet Model (CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module eismint_forcing

  ! read configuration and generate eismint mass balance and 
  ! temperature fields

  use glimmer_global, only : dp
  use glide_types, only : eismint_climate_type

  !MAKE_RESTART
#ifdef RESTARTS
#define RST_EISMINT_FORCING
!JCC - no restarts yet
!#include "glimmer_rst_head.inc"
#undef RST_EISMINT_FORCING
#endif

contains

#ifdef RESTARTS
#define RST_EISMINT_FORCING
!JCC - no restarts yet
!#include "glimmer_rst_body.inc"
#undef RST_EISMINT_FORCING
#endif

  subroutine eismint_initialise(eismint_climate,config)

    ! initialise eismint_climate model

    use glimmer_global, only: dp
    use glimmer_paramets, only: thk0, scyr, tim0
    use glimmer_physcon, only: scyr
    use glimmer_config
    use glide_types
    implicit none

    type(eismint_climate_type) :: eismint_climate         ! structure holding climate info
    type(ConfigSection), pointer :: config  ! structure holding sections of configuration file   

!WHL - The old scaling looked like this: eismint_climate%nmsb(1) = eismint_climate%nmsb(1) / (acc0 * scyr)
!       where acc0 = thk0*vel0/len0.
!      I replaced (acc0 * scyr) with acab_scale = scyr*thk0/tim0, where tim0 = len0/vel0.  
!      This is the scaling used in other parts of the code, including Glint.
!      It can be shown (but is not immediately obvious) that acab_scale = acc0 * scyr.
!      This scale factor assumes that the input mass balance has units of m/yr.
!
!      Note: We should not use the parameter scale_acab in glimmer_scales because
!            it may not have been initialized yet.

    real(dp), parameter :: acab_scale = scyr*thk0/tim0

    call eismint_readconfig(eismint_climate,config)
    call eismint_printconfig(eismint_climate)

    ! scale parameters
    ! assumes that eismint_climate%nmsb starts with units of m/yr
 
    select case(eismint_climate%eismint_type)

    case(1)   ! EISMINT-1 fixed margin
       eismint_climate%nmsb(1) = eismint_climate%nmsb(1) / acab_scale

    case(2)   ! EISMINT-1 moving margin
       eismint_climate%airt(2) = eismint_climate%airt(2) * thk0
       eismint_climate%nmsb(1) = eismint_climate%nmsb(1) / acab_scale
       eismint_climate%nmsb(2) = eismint_climate%nmsb(2) / acab_scale

    case(3)   ! EISMINT-2
       eismint_climate%nmsb(1) = eismint_climate%nmsb(1) / acab_scale
       eismint_climate%nmsb(2) = eismint_climate%nmsb(2) / acab_scale

    case(4)   ! MISMIP-1
       eismint_climate%nmsb(1) = eismint_climate%nmsb(1) / acab_scale

    end select
       
  end subroutine eismint_initialise

  subroutine eismint_readconfig(eismint_climate, config)

    ! read configuration

    use glimmer_log
    use glimmer_config
    implicit none

    type(eismint_climate_type) :: eismint_climate         ! structure holding climate info
    type(ConfigSection), pointer :: config  ! structure holding sections of configuration file   

    ! local variables
    type(ConfigSection), pointer :: section
    real(kind=dp), dimension(:), pointer :: dummy

    call GetSection(config,section,'EISMINT-1 fixed margin')
    if (associated(section)) then
       eismint_climate%eismint_type = 1
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       eismint_climate%airt = (/-34.15d0, 8.d-8/)
       if (associated(dummy)) then
          eismint_climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,1)
       eismint_climate%nmsb = (/0.3d0, 0.d0, 0.d0/)
       if (associated(dummy)) then
          eismint_climate%nmsb(1) = dummy(1)
       end if
       call GetValue(section,'period',eismint_climate%period)
       call GetValue(section,'mb_amplitude',eismint_climate%mb_amplitude)
       return       
    end if

    !TODO - I think the default airt values declared above are appropriate for this case.
    !       Set them here instead.

    call GetSection(config,section,'EISMINT-1 moving margin')
    if (associated(section)) then
       eismint_climate%eismint_type = 2
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
          eismint_climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
          eismint_climate%nmsb = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'period',eismint_climate%period)
       eismint_climate%mb_amplitude = 100000.d0
       call GetValue(section,'mb_amplitude',eismint_climate%mb_amplitude)
       return
    end if

    call GetSection(config,section,'EISMINT-2')
    if (associated(section)) then
       eismint_climate%eismint_type = 3
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
          eismint_climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       else
          eismint_climate%airt = (/-35.d0, 1.67d-5/)
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
          eismint_climate%nmsb = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       return
    end if
    
    !mismip tests

    !TODO - Assign reasonable default values if not present in config file

    call GetSection(config,section,'MISMIP-1')
    if (associated(section)) then
       eismint_climate%eismint_type = 4
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
           eismint_climate%airt = dummy
           deallocate(dummy)
           dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
           eismint_climate%nmsb = dummy
           deallocate(dummy)
           dummy=>NULL()
       end if
       return
    end if

    !exact verification
    !TODO - Is this test currently supported?

    call GetSection(config,section,'EXACT')
    if (associated(section)) then
       eismint_climate%eismint_type = 5
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
           eismint_climate%airt = dummy
           deallocate(dummy)
           dummy=>NULL()
       end if
       return
    end if

    ! Standard higher-order tests
    ! These do not require EISMINT-type input parameters.

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

    call GetSection(config,section,'GIS-TEST')
    if (associated(section)) then 
        return
    end if

    !TODO - Any other allowed tests to add here?

    ! Abort if one of the above cases has not been specified.
    call write_log('No EISMINT forcing selected',GM_FATAL)

  end subroutine eismint_readconfig

  subroutine eismint_printconfig(eismint_climate)

    ! print eismint_climate configuration

    use glimmer_log
    use parallel, only: tasks
    implicit none

    type(eismint_climate_type) :: eismint_climate   ! structure holding climate info
    character(len=100) :: message

    call write_log_div

    select case(eismint_climate%eismint_type)

    case(1)
       call write_log('EISMINT-1 fixed margin configuration')
       call write_log('------------------------------------')
       write(message,*) 'temperature  : ',eismint_climate%airt(1)
       call write_log(message)
       write(message,*) '               ',eismint_climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',eismint_climate%nmsb(1)
       call write_log(message)
       write(message,*) 'period       : ',eismint_climate%period
       call write_log(message)
       if (eismint_climate%period .gt. 0.d0) then
          write(message,*) 'mb amplitude : ',eismint_climate%mb_amplitude
          call write_log(message)
       end if

    case(2)
       call write_log('EISMINT-1 moving margin configuration')
       call write_log('-------------------------------------')
       write(message,*) 'temperature  : ',eismint_climate%airt(1)
       call write_log(message)
       write(message,*) '               ',eismint_climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',eismint_climate%nmsb(1)
       call write_log(message)
       write(message,*) '               ',eismint_climate%nmsb(2)
       call write_log(message)
       write(message,*) '               ',eismint_climate%nmsb(3)
       call write_log(message)
       write(message,*) 'period       : ',eismint_climate%period
       call write_log(message)
       if (eismint_climate%period .gt. 0.d0) then
          write(message,*) 'mb amplitude : ',eismint_climate%mb_amplitude
          call write_log(message)
       end if

    case(3)
       call write_log('EISMINT-2')
       call write_log('---------')
       write(message,*) 'temperature  : ',eismint_climate%airt(1)
       call write_log(message)
       write(message,*) '               ',eismint_climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',eismint_climate%nmsb(1)
       call write_log(message)
       write(message,*) '               ',eismint_climate%nmsb(2)
       call write_log(message)
       write(message,*) '               ',eismint_climate%nmsb(3)
       call write_log(message)
    end select

    if ( (eismint_climate%eismint_type > 0) .and. (tasks > 1) ) then
       call write_log('EISMINT tests are not supported for more than one processor', GM_FATAL)
    end if

    call write_log('')

  end subroutine eismint_printconfig

  subroutine eismint_massbalance(eismint_climate,model,time)

    ! calculate eismint mass balance

!TODO - Remove acc0

    use glimmer_global, only : dp
    use glide_types
    use glimmer_paramets, only : len0, acc0, scyr
    use glimmer_physcon, only : pi
    use glimmer_scales, only : scale_acab
    implicit none

    type(eismint_climate_type) :: eismint_climate         ! structure holding climate info
    type(glide_global_type) :: model        ! model instance
    real(dp), intent(in) :: time            ! current time

    !WHL - Changed 'periodic_bc' to 'periodic' to avoid a name conflict with parallel modules
    ! local variables
    integer  :: ns,ew
    real(dp) :: dist, ewct, nsct, grid, rel
    real(dp) :: periodic = 1.d0  !TODO - Make this an integer?

    ewct = (real(model%general%ewn,dp) + 1.d0) / 2.d0
    nsct = (real(model%general%nsn,dp) + 1.d0) / 2.d0
    grid = real(model%numerics%dew,dp) * len0

    if (model%options%periodic_ew) then
        periodic = 0.d0
    else
        periodic = 1.d0
    end if

    select case(eismint_climate%eismint_type)

    case(1)
       ! EISMINT-1 fixed margin
       model%climate%acab(:,:) = eismint_climate%nmsb(1)
       if (eismint_climate%period .ne. 0.d0) then
          model%climate%acab(:,:) = model%climate%acab(:,:) + eismint_climate%mb_amplitude * &
               sin(2.d0*pi*time/eismint_climate%period)/ (acc0 * scyr)
!          model%climate%acab(:,:) = model%climate%acab(:,:) + climate%mb_amplitude * sin(2.d0*pi*time/climate%period) / scale_acab
       end if

    case(2)
       ! EISMINT-1 moving margin       
       if (eismint_climate%period .ne. 0.d0) then
          rel = eismint_climate%nmsb(3) + eismint_climate%mb_amplitude*sin(2.d0*pi*time/eismint_climate%period)
       else
          rel = eismint_climate%nmsb(3)
       end if

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic*(real(ew,kind=dp) - ewct)**2 + (real(ns,kind=dp) - nsct)**2)
             model%climate%acab(ew,ns) = min(eismint_climate%nmsb(1), eismint_climate%nmsb(2) * (rel - dist))
          end do
       end do

    case(3)
       ! EISMINT-2
       rel = eismint_climate%nmsb(3)

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic*(real(ew,kind=dp) - ewct)**2 + (real(ns,kind=dp) - nsct)**2)
             model%climate%acab(ew,ns) = min(eismint_climate%nmsb(1), eismint_climate%nmsb(2) * (rel - dist))
          end do
       end do

    case(4)
       !mismip 1
       model%climate%acab = eismint_climate%nmsb(1)

    case(5)
       !verification 
       call exact_surfmass(eismint_climate,model,time,1.d0,eismint_climate%airt(2))

    end select

  end subroutine eismint_massbalance

  subroutine eismint_surftemp(eismint_climate,model,time)

    ! calculate eismint air surface temperature

    use glide_types
    use glimmer_global, only: dp
    use glimmer_paramets, only : len0
    use glimmer_physcon, only : pi
    implicit none

    type(eismint_climate_type) :: eismint_climate         ! structure holding climate info
    type(glide_global_type) :: model        ! model instance
    real(dp), intent(in) :: time            ! current time

    ! local variables
    integer  :: ns,ew
    real(dp) :: dist, ewct, nsct, grid
    real(dp) :: periodic = 1.d0

    ewct = (real(model%general%ewn,dp)+1.d0) / 2.d0
    nsct = (real(model%general%nsn,dp)+1.d0) / 2.d0
    grid = real(model%numerics%dew,dp) * len0

    if (model%options%periodic_ew) then
        periodic = 0.d0
    else
        periodic = 1.d0
    end if

    select case(eismint_climate%eismint_type)

    case(1)
       ! EISMINT-1 fixed margin
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * max(periodic*abs(real(ew,kind=dp) - ewct),abs(real(ns,kind=dp) - nsct))*1d-3
             model%climate%artm(ew,ns) = eismint_climate%airt(1) + eismint_climate%airt(2) * dist*dist*dist
          end do
       end do
       if (eismint_climate%period .ne. 0.d0) then
          model%climate%artm(:,:) = model%climate%artm(:,:) + 10.d0*sin(2.d0*pi*time/eismint_climate%period)
       end if

    case(2)
       ! EISMINT-1 moving margin
       model%climate%artm(:,:) = eismint_climate%airt(1) - model%geometry%thck(:,:) * eismint_climate%airt(2)
       if (eismint_climate%period .ne. 0.d0) then
          model%climate%artm(:,:) = model%climate%artm(:,:) + 10.d0*sin(2.d0*pi*time/eismint_climate%period)
       end if

    case(3)
       ! EISMINT-2
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic*(real(ew,kind=dp) - ewct)**2 + (real(ns,kind=dp) - nsct)**2)
             model%climate%artm(ew,ns) = eismint_climate%airt(1)+eismint_climate%airt(2) * dist
          end do
       end do

    case(4)
       model%climate%artm = eismint_climate%airt(1)

    case(5)
       !call both massbalance and surftemp at the same time to save computing time. 
       call exact_surfmass(eismint_climate,model,time,0.d0,eismint_climate%airt(2))
    end select

  end subroutine eismint_surftemp
  
  !which_call - eismint_surftemp(0)/eismint_massbalance(1)/both(2)
  !which_test - test f(0)/test g(1)/exact(2)

  subroutine exact_surfmass(eismint_climate,model,time,which_call,which_test)

    use glide_types
    use testsFG
    implicit none

    type(eismint_climate_type) :: eismint_climate             ! structure holding climate info
    type(glide_global_type) :: model            ! model instance
    real(dp), intent(in) :: time                ! current time
    real(dp), intent(in) :: which_test          !  Which exact test (F=0, G=1)
    real(dp), intent(in) :: which_call          !  0 = surface temp, 1 = mass balance
    integer  :: ns,ew,lev,center

    !verification
    real(dp) ::  r, z, x, y                    !in variables
    real(dp) ::  H, TT, U, w, Sig, M, Sigc        !out variables
    real(dp) :: H_0

    center = (model%general%ewn - 1) * 0.5

    !TODO - Change which_call to an integer?
    !       Modify for Glissade? (dissip has smaller vertical dimension)
    if (which_call .eq. 0.d0 .or. which_call .eq. 2.d0) then

        !point by point call to the function 
        do ns = 1,model%general%nsn
            do ew = 1,model%general%ewn
                x = (ew - center)*model%numerics%dew
                y = (ns - center)*model%numerics%dns
                r = sqrt(x**2 + y**2)
                do lev = 1, model%general%upn
                    z = model%geometry%thck(ew,ns)*model%numerics%sigma(lev)
                    !the function only returns values within the radius
                    if(r>0.d0 .and. r<L) then
                        if (which_test .eq. 0.d0) then
                            call testF(r,z,H,TT,U,w,Sig,M,Sigc)
                        else if (which_test .eq. 1.d0) then
                            call testG(time,r,z,H,TT,U,w,Sig,M,Sigc)
                        else if (which_test .eq. 2.d0) then
                            !H_0 = H0
                            H_0 = model%geometry%thck(center,center)
                            !TT = 0.d0
                            !WHL - For Glissade, dissip is defined only at levels 1:upn-1
                            TT = model%temper%dissip(lev,ew,ns)
                            call model_exact(time,r,z,model%geometry%thck(ew,ns),H_0,TT,U,w,Sig,M,Sigc)
                        end if
                        model%tempwk%compheat(lev,ew,ns) = -Sigc*SperA*1.0d6 ! (10^(-3) deg mK/a) (*10^3)
                    else
                        model%tempwk%compheat(lev,ew,ns) = 0.d0
                    end if
                end do
            end do
        end do

    else if (which_call .eq. 1.d0 .or. which_call .eq. 2.d0) then

        do ns = 1,model%general%nsn
            do ew = 1,model%general%ewn
                x = (ew - center)*model%numerics%dew
                y = (ns - center)*model%numerics%dns
                r = sqrt(x**2 + y**2)
                z = model%geometry%thck(ew,ns)
                !the function only returns values within the radius
                if(r>0.d0 .and. r<L) then
                    if (which_test .eq. 0.d0) then
                        call testF(r,z,H,TT,U,w,Sig,M,Sigc)
                    else if (which_test .eq. 1.d0) then
                        call testG(time,r,z,H,TT,U,w,Sig,M,Sigc)
                    else if (which_test .eq. 2.d0) then
                        H_0 = H0 !H_0 = model%geometry%thck(center,center)
                        TT = 0.d0 !TT = model%temper%dissip(lev,ew,ns)
                        call model_exact(time,r,z,model%geometry%thck(ew,ns),H_0,TT,U,w,Sig,M,Sigc)
                    end if
                    model%climate%acab(ew,ns) = M*SperA  !m/a
                else
                    if(r .eq. 0.d0) then
                        model%climate%acab(ew,ns) = 0.d0
                        !model%climate%acab(ew,ns) = H0 - model%geometry%thck(center,center) !set it back to H0
                    else
                        model%climate%acab(ew,ns) = -0.1d0  !outside the glacier we use .1 m ablation as specified in the paper
                    end if
                end if
            end do
        end do
    end if

  end subroutine exact_surfmass
  
end module eismint_forcing
