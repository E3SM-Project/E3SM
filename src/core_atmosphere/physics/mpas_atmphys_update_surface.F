! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_update_surface
 use mpas_dmpar
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_date_time
 use mpas_atmphys_constants,only: stbolt
 use mpas_atmphys_vars

 implicit none
 private
 public:: physics_update_sst,         &
          physics_update_sstskin,     &
          physics_update_surface,     &
          physics_update_deepsoiltemp


!Update surface boundary conditions.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_update_surface:
! -------------------------------------------
! physics_update_surface     : update the surface albedo and greeness fraction.
! physics_update_sst         : update the sea-surface temperatures.
! physics_update_sstskin     : add a diurnal cycle to the sea-surface temperatures.
! physics_update_deepsoiltemp: update the deep soil temperatures.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * revised subroutine physics_update_sst.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-24.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * now use isice and iswater initialized in the init file instead of initialized in mpas_atmphys_landuse.F.
!   Laura D. Fowler (laura@ucar.edu) / 2017-01-13.
! * corrected the initialization of the soil temperature tslb over ocean points for exact restartability, and
!   for consistency with module_sf_noahdrv.F when itimestep = 1.
!   Laura D. Fowler (laura@ucar.edu) / 2017-08-29.


 contains


!=================================================================================================================
 subroutine physics_update_surface(current_date,config_sfc_albedo,mesh,sfc_input)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: mesh
 character(len=*),intent(in):: current_date
 logical,intent(in):: config_sfc_albedo

!inout variables:
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
!logical,pointer:: config_sfc_albedo

 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:)  ,pointer:: sfc_albbck
 real(kind=RKIND),dimension(:,:),pointer:: albedo12m
 real(kind=RKIND),dimension(:)  ,pointer:: vegfra,shdmin,shdmax
 real(kind=RKIND),dimension(:,:),pointer:: greenfrac
 
!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(sfc_input,'landmask'  , landmask  )
 call mpas_pool_get_array(sfc_input,'albedo12m' , albedo12m )
 call mpas_pool_get_array(sfc_input,'sfc_albbck', sfc_albbck)

 call mpas_pool_get_array(sfc_input,'greenfrac' , greenfrac )
 call mpas_pool_get_array(sfc_input,'vegfra'    , vegfra    )
 call mpas_pool_get_array(sfc_input,'shdmin'    , shdmin    )
 call mpas_pool_get_array(sfc_input,'shdmax'    , shdmax    )

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

!=================================================================================================================
 subroutine physics_update_sst(dminfo,config_frac_seaice,mesh,sfc_input,diag_physics)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 logical,intent(in):: config_frac_seaice

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics


 integer,pointer:: nCellsSolve,nSoilLevels
 integer,pointer:: isice,iswater

 real(kind=RKIND),dimension(:),pointer  :: sfc_albbck,sst,snow,tmn,tsk,vegfra,xice,seaice
 real(kind=RKIND),dimension(:),pointer  :: snowc,snowh
 real(kind=RKIND),dimension(:,:),pointer:: tslb,sh2o,smois

 real(kind=RKIND),dimension(:),pointer:: sfc_albedo,sfc_emiss,sfc_emibck
 real(kind=RKIND),dimension(:),pointer:: xicem,xland

!local variables:
 integer:: icheck
 integer:: iCell,iSoil
 integer:: nb_to_land,nb_to_ocean,nb_removed
 integer,dimension(:),pointer:: isltyp,ivgtyp,landmask

 real(kind=RKIND):: local_min,local_max
 real(kind=RKIND):: global_sst_min,global_sst_max
 real(kind=RKIND):: global_xice_min,global_xice_max

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)

 call mpas_pool_get_array(sfc_input,'isice'     ,isice     )
 call mpas_pool_get_array(sfc_input,'iswater'   ,iswater   )
 call mpas_pool_get_array(sfc_input,'isltyp'    ,isltyp    )    
 call mpas_pool_get_array(sfc_input,'ivgtyp'    ,ivgtyp    )
 call mpas_pool_get_array(sfc_input,'landmask'  ,landmask  )
 call mpas_pool_get_array(sfc_input,'vegfra'    ,vegfra    )
 call mpas_pool_get_array(sfc_input,'sfc_albbck',sfc_albbck)
 call mpas_pool_get_array(sfc_input,'sst'       ,sst       )
 call mpas_pool_get_array(sfc_input,'tmn'       ,tmn       )
 call mpas_pool_get_array(sfc_input,'skintemp'  ,tsk       )
 call mpas_pool_get_array(sfc_input,'tslb'      ,tslb      )
 call mpas_pool_get_array(sfc_input,'sh2o'      ,sh2o      )
 call mpas_pool_get_array(sfc_input,'smois'     ,smois     )
 call mpas_pool_get_array(sfc_input,'snow'      ,snow      )
 call mpas_pool_get_array(sfc_input,'snowc'     ,snowc     )
 call mpas_pool_get_array(sfc_input,'snowh'     ,snowh     )
 call mpas_pool_get_array(sfc_input,'seaice'    ,seaice    )
 call mpas_pool_get_array(sfc_input,'xice'      ,xice      )
 call mpas_pool_get_array(sfc_input,'xland'     ,xland     )

 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'sfc_emibck',sfc_emibck)
 call mpas_pool_get_array(diag_physics,'xicem'     ,xicem     )

!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine physics_update_sst:')
!call mpas_log_write('--- config_frac_seaice =$l', logicArgs=(/config_frac_seaice/))
!call mpas_log_write('--- xice_threshold     =$r', realArgs=(/xice_threshold/))
!call mpas_log_write('--- isice  =$i', intArgs=(/isice/))
!call mpas_log_write('--- iswater=$i', intArgs=(/iswater/))

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
       do iSoil = 1, nSoilLevels
          tslb(iSoil,iCell) = 273.16
       enddo 
    endif
 enddo
!call mpas_log_write('')
!call mpas_log_write('--- nb of seaice points converted to land points  = $i',intArgs=(/nb_to_land/))
!call mpas_log_write('--- nb of seaice points converted to ocean points = $i',intArgs=(/nb_to_ocean/))
!call mpas_log_write('--- nb of seaice points less than xice threshold  = $i',intArgs=(/nb_removed/))

!finally, update the sea-ice flag. save xice prior to next update:
 do iCell = 1, nCellsSolve
    xicem(iCell)  = xice(iCell)
    seaice(iCell) = 0._RKIND
    if(xice(iCell) > 0._RKIND) then
       seaice(iCell) = 1._RKIND
    endif
 enddo

!local and global max and min sea-surface temperatures and fractional sea-ice:
!local_min =  999._RKIND
!local_max = -999._RKIND
!do iCell = 1,nCellsSolve
!   if(xland(iCell) == 2._RKIND .and. sst(iCell) <= local_min) local_min = sst(iCell)
!   if(xland(iCell) == 2._RKIND .and. sst(iCell) >= local_max) local_max = sst(iCell)
!enddo
!call mpas_dmpar_min_real(dminfo,local_min,global_sst_min)
!call mpas_dmpar_max_real(dminfo,local_max,global_sst_max)
!call mpas_log_write('')
!call mpas_log_write('--- min local SST   = $r',realArgs=(/local_min/))
!call mpas_log_write('--- max local SST   = $r',realArgs=(/local_max/))
!call mpas_log_write('--- min global SST  = $r',realArgs=(/global_sst_min/))
!call mpas_log_write('--- max global SST  = $r',realArgs=(/global_sst_max/))

!local_min =  999._RKIND
!local_max = -999._RKIND
!do iCell = 1,nCellsSolve
!   if(xland(iCell) == 1._RKIND .and. xice(iCell) <= local_min) local_min = xice(iCell)
!   if(xland(iCell) == 1._RKIND .and. xice(iCell) >= local_max) local_max = xice(iCell)
!enddo
!call mpas_dmpar_min_real(dminfo,local_min,global_xice_min)
!call mpas_dmpar_max_real(dminfo,local_max,global_xice_max)
!if(local_min .eq.  999._RKIND) local_min = 0._RKIND
!if(local_max .eq. -999._RKIND) local_max = 0._RKIND
!call mpas_log_write('')
!call mpas_log_write('--- min local XICE  = $r',realArgs=(/local_min/))
!call mpas_log_write('--- max local XICE  = $r',realArgs=(/local_max/))
!call mpas_log_write('--- min global XICE = $r',realArgs=(/global_xice_min/))
!call mpas_log_write('--- max global XICE = $r',realArgs=(/global_xice_max/))

!call mpas_log_write('--- end subroutine physics_update_sst')
!call mpas_log_write('')

 end subroutine physics_update_sst

!=================================================================================================================
 subroutine physics_update_sstskin(dt,mesh,diag_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 real(kind=RKIND),intent(in):: dt

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 integer,pointer:: nCellsSolve

 real(kind=RKIND),dimension(:),pointer:: sst,tsk,xland
 real(kind=RKIND),dimension(:),pointer:: glw,gsw
 real(kind=RKIND),dimension(:),pointer:: hfx,qfx
 real(kind=RKIND),dimension(:),pointer:: emiss,ust
 real(kind=RKIND),dimension(:),pointer:: sstsk,dtc1,dtw1

!local parameters:
 integer, parameter:: n=1152
 real(kind=RKIND),parameter:: z1=3.,an=.3,zk=.4,rho=1.2,rhow=1025.,cw=4190.
 real(kind=RKIND),parameter:: g=9.8,znuw=1.e-6,zkw=1.4e-7,sdate=1201.6667

!local variables:
 integer:: iCell

 real(kind=RKIND):: lw, sw, q, qn, zeta, dep, dtw3, skinmax, skinmin
 real(kind=RKIND):: fs, con1, con2, con3, con4, con5, zlan, q2, ts, phi, qn1
 real(kind=RKIND):: usw, qo, swo, us, tb, dtc, dtw, alw, dtwo, delt, f1

!-----------------------------------------------------------------------------------------------------------------
! call mpas_log_write('')
! call mpas_log_write('--- enter subroutine physics_update_sstskin:')

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(sfc_input,'skintemp',tsk  )
 call mpas_pool_get_array(sfc_input,'sst'     ,sst  )
 call mpas_pool_get_array(sfc_input,'xland'   ,xland)

 call mpas_pool_get_array(diag_physics,'sstsk'    ,sstsk)
 call mpas_pool_get_array(diag_physics,'sstsk_dtc',dtc1 )
 call mpas_pool_get_array(diag_physics,'sstsk_dtw',dtw1 )
 call mpas_pool_get_array(diag_physics,'sfc_emiss',emiss)
 call mpas_pool_get_array(diag_physics,'glw'      ,glw  )
 call mpas_pool_get_array(diag_physics,'gsw'      ,gsw  )
 call mpas_pool_get_array(diag_physics,'hfx'      ,hfx  )
 call mpas_pool_get_array(diag_physics,'qfx'      ,qfx  )
 call mpas_pool_get_array(diag_physics,'ust'      ,ust  )

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

!=================================================================================================================
 subroutine physics_update_deepsoiltemp(LeapYear,dt,julian_in,mesh,sfc_input,diag_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in) :: mesh
 logical,intent(in):: LeapYear
 real(kind=RKIND),intent(in):: dt,julian_in

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 integer,pointer:: nCellsSolve,nLags

 real(kind=RKIND),dimension(:),pointer:: nsteps_accum,ndays_accum
 real(kind=RKIND),dimension(:),pointer  :: tday_accum,tmn,tsk,tyear_accum,tyear_mean
 real(kind=RKIND),dimension(:,:),pointer:: tlag 

!local variables:
 integer:: iCell,iLag,n

 real(kind=RKIND),parameter:: tconst = 0.6
 real(kind=RKIND):: deltat,julian,tprior,yrday

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine physics_update_deepsoiltemp:')

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nLags'      ,nLags      )

 call mpas_pool_get_array(diag_physics,'nsteps_accum',nsteps_accum)
 call mpas_pool_get_array(diag_physics,'ndays_accum' ,ndays_accum )
 call mpas_pool_get_array(diag_physics,'tday_accum'  ,tday_accum  )
 call mpas_pool_get_array(diag_physics,'tyear_accum' ,tyear_accum )
 call mpas_pool_get_array(diag_physics,'tyear_mean'  ,tyear_mean  )
 call mpas_pool_get_array(diag_physics,'tlag'        ,tlag        )

 call mpas_pool_get_array(sfc_input,'tmn'     ,tmn )
 call mpas_pool_get_array(sfc_input,'skintemp',tsk )

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

!call mpas_log_write('--- yrday          = $r',realArgs=(/yrday/))
!call mpas_log_write('--- julian_in      = $r',realArgs=(/julian_in/))
!call mpas_log_write('--- nint(julian_in)= $i',intArgs=(/nint(julian_in)/))
!call mpas_log_write('--- deltat         = $r',realArgs=(/deltat/))
!call mpas_log_write('--- nint(deltat)-dt= $l',logicArgs=(/nint(deltat) .lt. dt/))

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

!=================================================================================================================
 end module mpas_atmphys_update_surface
!=================================================================================================================


