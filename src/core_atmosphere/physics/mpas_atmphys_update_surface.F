! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_update_surface
 use mpas_configure, only: config_frac_seaice,config_sfc_albedo
 use mpas_dmpar
 use mpas_grid_types

 use mpas_atmphys_date_time
 use mpas_atmphys_constants,only: stbolt
 use mpas_atmphys_landuse, only : isice,iswater 
 use mpas_atmphys_vars

 implicit none
 private
 public:: physics_update_sst,         &
          physics_update_sstskin,     &
          physics_update_surface,     &
          physics_update_deepsoiltemp


!>\brief update surface boundary conditions.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_update_surface:
!> -------------------------------------------
!> physics_update_surface     : update the surface albedo and greeness fraction.
!> physics_update_sst         : update the sea-surface temperatures.
!> physics_update_sstskin     : add a diurnal cycle to the sea-surface temperatures.
!> physics_update_deepsoiltemp: update the deep soil temperatures.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * revised subroutine physics_update_sst.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-08-24.


 contains

!==================================================================================================
 subroutine physics_update_surface(current_date,mesh,sfc_input)
!==================================================================================================

!input variables:
 type(mesh_type),intent(in) :: mesh
 character(len=*),intent(in):: current_date

!inout variables:
 type(sfc_input_type),intent(inout):: sfc_input

!local variables:
 integer:: iCell

 integer:: nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:)  ,pointer:: sfc_albbck
 real(kind=RKIND),dimension(:,:),pointer:: albedo12m

 real(kind=RKIND),dimension(:)  ,pointer:: vegfra,shdmin,shdmax
 real(kind=RKIND),dimension(:,:),pointer:: greenfrac
 
!--------------------------------------------------------------------------------------------------

 nCellsSolve = mesh % nCellsSolve

 landmask   => sfc_input % landmask   % array
 albedo12m  => sfc_input % albedo12m  % array
 sfc_albbck => sfc_input % sfc_albbck % array

 greenfrac  => sfc_input % greenfrac  % array
 vegfra     => sfc_input % vegfra     % array
 shdmin     => sfc_input % shdmin     % array
 shdmax     => sfc_input % shdmax     % array

!updates the surface background albedo for the current date as a function of the monthly-mean
!surface background albedo valid on the 15th day of the month, if config_sfc_albedo is true:
 if(config_sfc_albedo) then

    call monthly_interp_to_date(nCellsSolve,current_date,albedo12m,sfc_albbck)

    do iCell = 1, nCellsSolve
       sfc_albbck(iCell) = sfc_albbck(iCell) / 100.
       if(landmask(iCell) .eq. 0) sfc_albbck(iCell) = 0.08
    enddo

 endif

!updates the green-ness fraction for the current date as a function of the monthly-mean green-
!ness valid on the 15th day of the month. get the min/max for each cell for the monthly green-
!ness fraction:
 call monthly_interp_to_date(nCellsSolve,current_date,greenfrac,vegfra)
 call monthly_min_max(nCellsSolve,greenfrac,shdmin,shdmax)

 end subroutine physics_update_surface

!==================================================================================================
 subroutine physics_update_sst(dminfo,mesh,sfc_input,diag_physics)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(sfc_input_type),intent(inout)   :: sfc_input
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: icheck
 integer:: iCell,iSoil,nCellsSolve,nSoilLevels
 integer:: nb_to_land,nb_to_ocean,nb_removed
 integer,dimension(:),pointer:: isltyp,ivgtyp,landmask

 real(kind=RKIND):: local_min,local_max
 real(kind=RKIND):: global_sst_min,global_sst_max
 real(kind=RKIND):: global_xice_min,global_xice_max

 real(kind=RKIND),dimension(:),pointer  :: sfc_albbck,sst,snow,tmn,tsk,vegfra,xice,seaice
 real(kind=RKIND),dimension(:),pointer  :: snowc,snowh
 real(kind=RKIND),dimension(:,:),pointer:: tslb,sh2o,smois

 real(kind=RKIND),dimension(:),pointer:: sfc_albedo,sfc_emiss,sfc_emibck
 real(kind=RKIND),dimension(:),pointer:: xicem,xland

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_update_sst:'
 write(0,*) '--- config_frac_seaice =', config_frac_seaice
 write(0,*) '--- xice_threshold     =', xice_threshold
 write(0,*) '--- isice  =', isice
 write(0,*) '--- iswater=', iswater

 nCellsSolve = mesh % nCellsSolve
 nSoilLevels = mesh % nSoilLevels

 isltyp     => sfc_input    % isltyp   % array    
 ivgtyp     => sfc_input    % ivgtyp   % array
 landmask   => sfc_input    % landmask % array
 vegfra     => sfc_input    % vegfra   % array
 sfc_albbck => sfc_input    % sfc_albbck % array
 sst        => sfc_input    % sst        % array
 tmn        => sfc_input    % tmn        % array
 tsk        => sfc_input    % skintemp   % array
 tslb       => sfc_input    % tslb       % array
 sh2o       => sfc_input    % sh2o       % array
 smois      => sfc_input    % smois      % array
 snow       => sfc_input    % snow       % array
 snowc      => sfc_input    % snowc      % array
 snowh      => sfc_input    % snowh      % array
 seaice     => sfc_input    % seaice     % array
 xice       => sfc_input    % xice       % array
 xland      => sfc_input    % xland      % array

 sfc_albedo => diag_physics % sfc_albedo % array
 sfc_emiss  => diag_physics % sfc_emiss  % array
 sfc_emibck => diag_physics % sfc_emibck % array
 xicem      => diag_physics % xicem      % array

 if(config_frac_seaice) then
    do iCell = 1,nCellsSolve
       if(xice(iCell) < xice_threshold) xice(iCell) = 0._RKIND
    enddo
 elseif(.not.config_frac_seaice) then
    do iCell = 1,nCellsSolve
       if(xice(iCell) >= xice_threshold) then
          xice(iCell) = 1._RKIND
       else
          xice(iCell) = 0._RKIND
       endif
    enddo
 endif

!update the surface albedo and surface emissivity. before updating xice, sfc_albedo and sfc_emiss
!are valid according to the earlier value of xice, xicem. now that xice has been updated, we also
!update sfc_albedo and sfc_emiss accordingly:
 if(config_frac_seaice) then
    do iCell = 1, nCellsSolve
       if(xice(iCell) /= xicem(iCell) .and. xicem(iCell) >= xice_threshold) then
          sfc_albedo(iCell) = 0.08 + (sfc_albedo(iCell) -0.08) * xice(iCell)/xicem(iCell)
          sfc_emiss(iCell)  = 0.98 + (sfc_emiss(iCell)-0.98  ) * xice(iCell)/xicem(iCell)
       endif
    enddo
 endif

 nb_to_land  = 0
 nb_to_ocean = 0
 nb_removed  = 0
 do iCell = 1, nCellsSolve
    if(xland(iCell) >= 1.5_RKIND .and. xice(iCell) >=  xice_threshold .and. &
       xicem(iCell) < xice_threshold) then

       nb_to_land = nb_to_land + 1
    !... sea-ice points are converted to land points:
       ivgtyp(iCell) = isice
       isltyp(iCell) = 16
       vegfra(iCell) = 0._RKIND
       xland(iCell)  = 1._RKIND
       tmn(iCell)    = 271.4_RKIND

       do iSoil = 1, nSoilLevels
          tslb(iSoil,iCell)  = tsk(iCell)
          smois(iSoil,iCell) = 1.0_RKIND
          sh2o(iSoil,iCell)  = 0.0_RKIND
       enddo

       !... over newly formed ice, initial guesses for the albedo and emissivity are based on
       !... default values over weater and ice. The surface albedo and emissivity can be upda
       !... ted later with the land-surface scheme.
       sfc_albedo(iCell) = 0.80 * xice(iCell) + 0.08 * (1.-xice(iCell))
       sfc_emiss(iCell)  = 0.98 * xice(iCell) + 0.98 * (1.-xice(iCell))
       sfc_albbck(iCell) = 0.80
       sfc_emibck(iCell) = 0.98

    elseif(xland(iCell) < 1.5_RKIND .and. xice(iCell) < xice_threshold .and. &
       xicem(iCell) >= xice_threshold) then

       nb_to_ocean = nb_to_ocean + 1
       !land points turn to water points:
       ivgtyp(iCell) = iswater
       isltyp(iCell) = 14
       vegfra(iCell) = 0._RKIND
       xland(iCell)  = 2._RKIND
       tmn(iCell)    = sst(iCell)

       snowc(iCell) = 0
       snow(iCell)  = 0.0_RKIND
       snowh(iCell) = 0.0_RKIND

       sfc_albedo(iCell) = 0.08_RKIND
       sfc_albbck(iCell) = 0.08_RKIND
       sfc_emiss(iCell)  = 0.98_RKIND
       sfc_emibck(iCell) = 0.98_RKIND

       do iSoil = 1, nSoilLevels
          tslb(iSoil,iCell)  = sst(iCell)
          smois(iSoil,iCell) = 1.0_RKIND
          sh2o(iSoil,iCell)  = 1.0_RKIND
       enddo

    elseif(xice(iCell) < xice_threshold .and. xicem(iCell) < xice_threshold) then
       nb_removed  = nb_removed + 1
       xice(iCell) = 0._RKIND

    endif

    if(xland(iCell) >= 1.5_RKIND) then
       tsk(iCell)    = sst(iCell)
       tslb(1,iCell) = sst(iCell)
    endif
 enddo
 write(0,*)
 write(0,*) '--- nb of seaice points converted to land points  = ',nb_to_land
 write(0,*) '--- nb of seaice points converted to ocean points = ',nb_to_ocean
 write(0,*) '--- nb of seaice points less than xice threshold  = ',nb_removed

!finally, update the sea-ice flag. save xice prior to next update:
 do iCell = 1, nCellsSolve
    xicem(iCell)  = xice(iCell)
    seaice(iCell) = 0._RKIND
    if(xice(iCell) > 0._RKIND) then
       seaice(iCell) = 1._RKIND
    endif
 enddo

!local and global max and min sea-surface temperatures and fractional sea-ice:
 local_min =  999._RKIND
 local_max = -999._RKIND
 do iCell = 1,nCellsSolve
    if(xland(iCell) == 2._RKIND .and. sst(iCell) <= local_min) local_min = sst(iCell)
    if(xland(iCell) == 2._RKIND .and. sst(iCell) >= local_max) local_max = sst(iCell)
 enddo
 call mpas_dmpar_min_real(dminfo,local_min,global_sst_min)
 call mpas_dmpar_max_real(dminfo,local_max,global_sst_max)
 write(0,*)
 write(0,*) '--- min local SST   = ',local_min
 write(0,*) '--- max local SST   = ',local_max
 write(0,*) '--- min global SST  = ',global_sst_min
 write(0,*) '--- max global SST  = ',global_sst_max

 local_min =  999._RKIND
 local_max = -999._RKIND 
 do iCell = 1,nCellsSolve
    if(xland(iCell) == 1._RKIND .and. xice(iCell) <= local_min) local_min = xice(iCell)
    if(xland(iCell) == 1._RKIND .and. xice(iCell) >= local_max) local_max = xice(iCell) 
 enddo
 call mpas_dmpar_min_real(dminfo,local_min,global_xice_min)
 call mpas_dmpar_max_real(dminfo,local_max,global_xice_max)
 if(local_min .eq.  999._RKIND) local_min = 0._RKIND
 if(local_max .eq. -999._RKIND) local_max = 0._RKIND
 write(0,*)
 write(0,*) '--- min local XICE  = ',local_min
 write(0,*) '--- max local XICE  = ',local_max
 write(0,*) '--- min global XICE = ',global_xice_min
 write(0,*) '--- max global XICE = ',global_xice_max

 write(0,*) '--- end subroutine physics_update_sst'
 write(0,*)

 end subroutine physics_update_sst

!==================================================================================================
 subroutine physics_update_sstskin(dt,mesh,diag_physics,sfc_input)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 real(kind=RKIND),intent(in):: dt

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(sfc_input_type),intent(inout)   :: sfc_input

!local parameters:
 integer, parameter:: n=1152
 real(kind=RKIND),parameter:: z1=3.,an=.3,zk=.4,rho=1.2,rhow=1025.,cw=4190.
 real(kind=RKIND),parameter:: g=9.8,znuw=1.e-6,zkw=1.4e-7,sdate=1201.6667

!local variables:
 integer:: iCell,nCellsSolve

 real(kind=RKIND):: lw, sw, q, qn, zeta, dep, dtw3, skinmax, skinmin
 real(kind=RKIND):: fs, con1, con2, con3, con4, con5, zlan, q2, ts, phi, qn1
 real(kind=RKIND):: usw, qo, swo, us, tb, dtc, dtw, alw, dtwo, delt, f1

 real(kind=RKIND),dimension(:),pointer:: sst,tsk,xland
 real(kind=RKIND),dimension(:),pointer:: glw,gsw
 real(kind=RKIND),dimension(:),pointer:: hfx,qfx
 real(kind=RKIND),dimension(:),pointer:: emiss,ust
 real(kind=RKIND),dimension(:),pointer:: sstsk,dtc1,dtw1

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_update_sstskin:'

 nCellsSolve = mesh % nCellsSolve

 tsk   => sfc_input % skintemp % array
 sst   => sfc_input % sst      % array
 xland => sfc_input % xland    % array

 sstsk => diag_physics % sstsk      % array
 dtc1  => diag_physics % sstsk_dtc  % array
 dtw1  => diag_physics % sstsk_dtw  % array
 emiss => diag_physics % sfc_emiss  % array
 glw   => diag_physics % glw        % array
 gsw   => diag_physics % gsw        % array
 hfx   => diag_physics % hfx        % array
 qfx   => diag_physics % qfx        % array
 ust   => diag_physics % ust        % array

 skinmax = -9999.
 skinmin =  9999.

!first, restore the surface temperature to the sea-surface temperature:
 do iCell = 1, nCellsSolve
    if(xland(iCell) .ge. 1.5) tsk(iCell) = sst(iCell)
 enddo

!calculate the skin sea-surface temperature: 
 do iCell = 1, nCellsSolve

    if(xland(iCell) .ge. 1.5) then

       qo   = glw(iCell)-emiss(iCell)*stbolt*(sstsk(iCell)**4)-2.5e6*qfx(iCell)-hfx(iCell)
       swo  = gsw(iCell)
       us   = max(ust(iCell),0.01)
       tb   = tsk(iCell)-273.15
       dtwo = dtw1(iCell)
       delt = dt

       q  = qo  / (rhow*cw)
       sw = swo / (rhow*cw)
!TEMPORARY KLUDGE
!      f1 = 1.-0.28*exp(-71.5*z1)-0.27*exp(-2.8*z1)-0.45*exp(-0.07*z1)
       f1 = 1.                   -0.27*exp(-2.8*z1)-0.45*exp(-0.07*z1)
!cool skin
       dtc = 0.0
!tb in C
       alw  = 1.e-5*max(tb,1.)
       con4 = 16.*g*alw*znuw**3/zkw**2
       usw  = sqrt(rho/rhow)*us
       con5 = con4/usw**4
!otherwise, iterations would be needed for the computation of fs
!iteration impact is less than 0.03C
       q2   = max(1./(rhow*cw),-q)
       zlan = 6./(1.+(con5*q2)**0.75)**0.333
       dep  = zlan*znuw/usw                    ! skin layer depth (m)
       fs   = 0.065+11.*dep-(6.6e-5/dep)*(1.-exp(-dep/8.e-4))
       fs   = max(fs,0.01)                     ! fract. of solar rad. absorbed in sublayer
       dtc  = dep*(q+sw*fs)/zkw                ! cool skin temp. diff (deg C)
       dtc  = min(dtc,0.)
!warm layer (X. Zeng)
       dtw  = 0.0
!tb in C
       alw  = 1.e-5*max(tb,1.)
       con1 = sqrt(5.*z1*g*alw/an)
       con2 = zk*g*alw
       qn   = q+sw*f1
       usw  = sqrt(rho/rhow)*us
!does not change when qn is positive
       if(dtwo.gt.0. .and. qn.lt.0.) then
          qn1 = sqrt(dtwo)*usw**2/con1
          qn  = max(qn,qn1)
       endif
       zeta = z1*con2*qn/usw**3
       if(zeta .gt. 0.) then
          phi = 1.+5.*zeta
       else
          phi = 1./sqrt(1.-16.*zeta)
       endif
       con3 = zk*usw/(z1*phi)
!use all SW flux
       dtw  = (dtwo+(an+1.) / an*(q+sw*f1)*    &
               delt/z1)/(1.+(an+1.)*con3*delt)
       dtw  = max(0.,dtw)
       dtwo = dtw
       ts   = tb + dtw + dtc

       skinmax = amax1(skinmax,ts-tb)
       skinmin = amin1(skinmin,ts-tb)
       sstsk(iCell) = ts+273.15                ! convert ts (in C) to sstsk (in K)
       dtc1(iCell)  = dtc                      ! dtc always in C
       dtw1(iCell)  = dtw                      ! dtw always in C

    endif

 enddo

!update the surface temperature over the oceans:
 do iCell = 1, nCellsSolve
    if(xland(iCell) .ge. 1.5) tsk(iCell) = sstsk(iCell)
 enddo

 end subroutine physics_update_sstskin

!==================================================================================================
 subroutine physics_update_deepsoiltemp(LeapYear,dt,julian_in,mesh,sfc_input,diag_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh
 logical,intent(in):: LeapYear
 real(kind=RKIND),intent(in):: dt,julian_in

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(sfc_input_type),intent(inout)   :: sfc_input

!local variables:
 integer:: iCell,iLag,n,nCellsSolve,nLags

 real(kind=RKIND),parameter:: tconst = 0.6
 real(kind=RKIND):: deltat,julian,tprior,yrday
 real(kind=RKIND),dimension(:),pointer:: nsteps_accum,ndays_accum
 real(kind=RKIND),dimension(:),pointer  :: tday_accum,tmn,tsk,tyear_accum,tyear_mean
 real(kind=RKIND),dimension(:,:),pointer:: tlag 

!--------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine physics_update_deepsoiltemp:' 

 nCellsSolve = mesh % nCellsSolve
 nLags       = mesh % nLags

 nsteps_accum => diag_physics % nsteps_accum % array
 ndays_accum  => diag_physics % ndays_accum  % array

 tmn  => sfc_input    % tmn       % array
 tsk  => sfc_input    % skintemp  % array
 tlag => diag_physics % tlag      % array
 tday_accum  => diag_physics % tday_accum  % array
 tyear_accum => diag_physics % tyear_accum % array
 tyear_mean  => diag_physics % tyear_mean  % array

!... defines the number of days in the year:
 if(LeapYear) then
    yrday = 366.
 else
    yrday = 365.
 endif

!... accumulate the skin temperature for current day:
 do iCell = 1, nCellsSolve
    tday_accum(iCell)  = tday_accum(iCell)  + tsk(iCell)*dt
!   tday_accum(iCell)  = tday_accum(iCell)  + tsk(iCell)
    nsteps_accum(iCell) = nsteps_accum(iCell) + dt
!   nsteps_accum(iCell) = nsteps_accum(iCell) + 1
 enddo

!... update the deep soil temperature at the end of the day:
 deltat = (julian_in-nint(julian_in))*24.*3600.

!write(0,*) '--- yrday          = ',yrday
!write(0,*) '--- julian_in      = ',julian_in
!write(0,*) '--- nint(julian_in)= ',nint(julian_in)
!write(0,*) '--- deltat         = ',deltat
!write(0,*) '--- nint(deltat)-dt= ',nint(deltat) .lt. dt

 if(abs(deltat) .le. dt/2) then
    julian = julian_in - 1. + dt/(3600.*24.)

    do iCell = 1, nCellsSolve

!--- update tmn:
       tprior = 0.
       do iLag = 1, nLags
          tprior = tprior + tlag(iLag,iCell)
       enddo
       tprior = tprior / nLags
       tmn(iCell) = tconst*tyear_mean(iCell) + (1-tconst)*tprior 

!--- update tlag:
       do iLag = 1, nLags-1
          tlag(iLag,iCell) = tlag(iLag+1,iCell)
       enddo
       tlag(nLags,iCell)   = tday_accum(iCell) / nsteps_accum(iCell)
       tday_accum(iCell)   = 0.0
       nsteps_accum(iCell) = 0.0

       !... end of year:
       if(yrday-julian .le. 1.) then
          tyear_mean(iCell)  = tyear_accum(iCell) / ndays_accum(iCell)
          tyear_accum(iCell) = 0.
          ndays_accum(iCell) = 0.0
       else
          tyear_accum(iCell) = tyear_accum(iCell) + tlag(nLags,iCell)
          ndays_accum(iCell) = ndays_accum(iCell) + 1.
       endif
       
    enddo

 endif !end of day

 end subroutine physics_update_deepsoiltemp

!==================================================================================================
 end module mpas_atmphys_update_surface
!==================================================================================================


