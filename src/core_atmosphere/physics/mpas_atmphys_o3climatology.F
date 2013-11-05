! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_o3climatology
 use mpas_kind_types
 use mpas_grid_types
 use mpas_atmphys_date_time
 use mpas_atmphys_constants
 use mpas_atmphys_utilities

!wrf physics:
 use module_ra_cam_support, only: r8, getfactors

 implicit none
 private
 public:: init_o3climatology,          &
          update_o3climatology,        &
          o3climatology_from_MPAS,     &
          vinterp_ozn

 integer,parameter:: latsiz = 64
 integer,parameter:: lonsiz = 1

!>\brief
!> mpas_atmphys_o3climatology contains the subroutines needed to initialize,interpolate,and update
!> the climatological monthly-mean ozone volume mixing ratios o3clim to the MPAS grid. Input data
!> files are the same as the ones used in the CAM long- and short-wave radiation codes.
!> when namelist parameter config_o3climatology is set to true, o3clim is used in the RRTMG long-
!> wave and short-wave radiation codes,and replaces the annual-mean ozone sounding used by default.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-07-03.

!>\details
!> subroutines in mpas_atmphys_o3climatology:
!> ------------------------------------------
!> init_o3climatology     : read the CAM ozone data files.
!> update_o3climatology   : interpolates the ozone volume mixing ratio to the current Julian day
!>                          as done for the greeness fraction in the MPAS time manager.
!> o3climatology_from_MPAS: interpolates the ozone volume mixing ratio to the current Julian day
!>                          as in the CAM radiation codes.      
!> vinterp_ozn            : vertical interpolation of the ozone volume mixing ratios from fixed
!>                          ozone pressure levels to the MPAS pressure levels.

!> modifications:
!>    none.


 contains

!==================================================================================================
 subroutine init_o3climatology(mesh,atm_input)
!==================================================================================================

!This subroutine assumes a uniform distribution of ozone concentration. It should be replaced
!with monthly climatology varying ozone distribution.

!input arguments:
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(atm_input_type), intent(inout):: atm_input

!local variables:
 integer,parameter:: pin_unit = 27
 integer,parameter:: lat_unit = 28
 integer,parameter:: oz_unit  = 29
 integer,parameter:: open_ok  = 0

 integer:: i,i1,i2,istat,k,j,m
 integer:: iCell,nCells,num_months,levsiz
 
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: pin
 real(kind=RKIND),dimension(:,:,:),pointer:: ozmixm

 real(kind=RKIND):: lat,lon,dlat,dlatCell
 real(kind=RKIND),dimension(latsiz):: lat_ozone
!real(Kind=RKIND),dimension(lonsiz,levsiz,latsiz,num_months):: ozmixin
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: ozmixin

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_init_o3:'

 nCells     = mesh % nCells
 num_months = mesh % nMonths
 levsiz     = mesh % nOznLevels
 pin     => atm_input % pin % array
 ozmixm  => atm_input % ozmixm % array
 latCell => mesh % latCell % array
 lonCell => mesh % lonCell % array

!-- read in ozone pressure data:
 open(pin_unit,file='OZONE_PLEV.TBL',action='READ',status='OLD',iostat=istat)
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                             'failure opening OZONE_PLEV.TBL')
 do k = 1,levsiz
    read(pin_unit,*) pin(k)
 enddo
 close(pin_unit)

!-- read in ozone lat data:
 open(lat_unit, file='OZONE_LAT.TBL',action='READ',status='OLD',iostat=istat) 
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                             'failure opening OZONE_LAT.TBL')
 do j = 1, latsiz
    read(lat_unit,*) lat_ozone(j)
!   write(0,101) j,lat_ozone(j)
 enddo
 close(lat_unit)

!-- read in ozone data:
 open(oz_unit,file='OZONE_DAT.TBL',action='READ',status='OLD',iostat=istat)
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                                'failure opening OZONE_DAT.TBL')

 allocate(ozmixin(lonsiz,levsiz,latsiz,num_months))
 do m=1,num_months
 do j=1,latsiz ! latsiz=64
 do k=1,levsiz ! levsiz=59
 do i=1,lonsiz ! lonsiz=1
    read(oz_unit,*) ozmixin(i,k,j,m)
 enddo
 enddo
 enddo
 enddo
 close(oz_unit)

!INTERPOLATION OF INPUT OZONE DATA TO MPAS GRID:
!write(0,*) 'max latCell=', maxval(latCell)/degrad
!write(0,*) 'min latCell=', minval(latCell)/degrad
!write(0,*) 'max lonCell=', maxval(lonCell)/degrad
!write(0,*) 'min lonCell=', minval(lonCell)/degrad
!write(0,*)
!write(0,*) 'max lat_ozone=',maxval(lat_ozone)
!write(0,*) 'min lat_ozone=',minval(lat_ozone)
 do iCell = 1,nCells
    lat = latCell(iCell)/degrad
    lon = lonCell(iCell)/degrad
    if(lat .gt. lat_ozone(latsiz)) then
     i1 = latsiz
     i2 = latsiz
    elseif(lat .lt. lat_ozone(1)) then
       i1 = 1
       i2 = 1
   else
       do i = 1, latsiz
          if(lat.ge.lat_ozone(i) .and. lat.lt.lat_ozone(i+1)) exit
       enddo
       i1 = i
       i2 = i+1
    endif

    do m = 1,num_months
    do k = 1,levsiz
    do j = 1,lonsiz
       dlat     = lat_ozone(i2)-lat_ozone(i1)
       dlatCell = lat-lat_ozone(i1)
       if(dlat == 0.) then
          ozmixm(m,k,iCell) = ozmixin(j,k,i1,m)
       else
          ozmixm(m,k,iCell) = ozmixin(j,k,i1,m) &
                     + (ozmixin(j,k,i2,m)-ozmixin(j,k,i1,m))*dlatCell/dlat
       endif
    enddo 
    enddo       
    enddo
!   do k = 1, levsiz
!      write(0,102) iCell,i1,i2,lat_ozone(i1),lat,lat_ozone(i2),ozmixin(1,k,i1,1), &
!                   ozmixm(1,k,iCell),ozmixin(1,k,i2,1)
!   enddo
 enddo
 deallocate(ozmixin)

 write(0,*) '--- end subroutine physics_init_o3.'

!formats:
 101 format(i3,12(1x,e15.8))
 102 format(i6,i6,i6,6(1x,e15.8))

 end subroutine init_o3climatology

!==================================================================================================
 subroutine update_o3climatology(current_date,mesh,atm_input,diag_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh
 character(len=*),intent(in):: current_date

!inout arguments:
 type(atm_input_type),intent(inout):: atm_input
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: iCell,iLev,nCellsSolve,nOznLevels
 real(kind=RKIND),dimension(:,:),pointer  :: o3clim
 real(kind=RKIND),dimension(:,:,:),pointer:: ozmixm

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_update_o3:'

 nCellsSolve = mesh % nCellsSolve
 nOznLevels  = mesh % nOznLevels

 ozmixm => atm_input % ozmixm % array
 o3clim => diag_physics % o3clim % array

 do iLev = 1,nOznLevels 
    call monthly_interp_to_date(nCellsSolve,current_date,ozmixm(:,iLev,:),o3clim(iLev,:))
 enddo
 write(0,*) '--- end subroutine physics_update_o3:'

 end subroutine update_o3climatology

!==================================================================================================
 subroutine o3climatology_from_MPAS(julian,mesh,atm_input,diag_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 real(kind=RKIND),intent(in):: julian
 type(atm_input_type),intent(in):: atm_input 

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 logical:: finddate
 logical:: ozncyc

 integer:: iCell,k,nCellsSolve,nVertLevels,nOznLevels,nMonths
 integer:: ijul,m,nm,np,np1
 integer, dimension(12) ::  date_oz
 data date_oz/16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350/

 real(kind=r8):: intjulian_r8
 real(kind=r8):: cdayozp_r8, cdayozm_r8
 real(kind=r8):: fact1_r8, fact2_r8
 real(kind=RKIND):: fact1,fact2

!--------------------------------------------------------------------------------------------------

!write(0,*)
!write(0,*) '--- enter subroutine o3climatology_from_MPAS:'

 nOznLevels  = mesh % nOznLevels
 nVertLevels = mesh % nVertLevels
 nCellsSolve = mesh % nCellsSolve
 nMonths     = mesh % nMonths

 ozncyc = .true.

!julian starts from 0.0 at 0Z on 1 Jan.
 intjulian_r8 = real(julian + 1.0_RKIND, r8)    ! offset by one day

!jan 1st 00z is julian=1.0 here
 ijul=int(intjulian_r8)
!note that following will drift. need to use actual month/day info to compute julian.
 intjulian_r8 = intjulian_r8 - real(ijul, r8)
 ijul = mod(ijul,365)
 if(ijul .eq. 0) ijul=365
 intjulian_r8 = intjulian_r8 + real(ijul, r8)
 np1=1

 finddate=.false.
 do m = 1, nMonths
    if(date_oz(m).gt.intjulian_r8 .and. .not.finddate) then
       np1 = m
       finddate = .true.
    endif
 enddo
 cdayozp_r8=date_oz(np1)

 if(np1 .gt. 1) then
    cdayozm_r8 = date_oz(np1-1)
    np = np1
    nm = np-1
 else
    cdayozm_r8 = date_oz(12)
    np = np1
    nm = nMonths
 endif
 call getfactors(ozncyc,np1,cdayozm_r8,cdayozp_r8,intjulian_r8,fact1_r8,fact2_r8) 
 fact1 = real(fact1_r8, RKIND)
 fact2 = real(fact2_r8, RKIND)
 
!write(0,*) '--- end subroutine getfactors:'
!write(0,*) 'fact1 =', fact1
!write(0,*) 'fact2 =', fact2 

!Time interpolation.
 do k = 1, nOznLevels
 do iCell = 1, nCellsSolve
    diag_physics % o3clim % array(k,iCell) = fact1 * atm_input % ozmixm % array(nm,k,iCell) &
                                           + fact2 * atm_input % ozmixm % array(np,k,iCell)
 end do
!write(0,101) k,diag_physics%o3clim%array(k,1),atm_input%ozmixm%array(nm,k,1), &
!             atm_input%ozmixm%array(np,k,1),diag_physics%o3clim%array(k,nCellsSolve), &
!             atm_input%ozmixm%array(nm,k,nCellsSolve),atm_input%ozmixm%array(np,k,nCellsSolve)
 end do
 101 format(i4,6(1x,e15.8))

!write(0,*) '--- end subroutine o3climatology_from_MPAS'

 end subroutine o3climatology_from_MPAS

!==================================================================================================
 subroutine vinterp_ozn (lchnk, ncol, pcols, pver, pmid_in, pin_in, levsiz, ozmix_in, o3vmr_out)
!--------------------------------------------------------------------------------------------------
! 
! Purpose: Interpolate ozone from current time-interpolated values to model levels
! 
! Method: Use pressure values to determine interpolation levels
! 
! Author: Bruce Briegleb
! 
!--------------------------------------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use phys_grid,     only: get_lat_all_p, get_lon_all_p
!  use comozp
!  use abortutils, only: endrun
!--------------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------------
!
! Arguments
!
   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: pcols, pver
   integer, intent(in) :: ncol                ! number of atmospheric columns
   integer, intent(in) :: levsiz              ! number of ozone layers

!ldf begin:
!  real(kind=RKIND), intent(in) :: pmid(pcols,pver)   ! level pressures (mks)
!  real(kind=RKIND), intent(in) :: pin(levsiz)        ! ozone data level pressures (mks)
!  real(kind=RKIND), intent(in) :: ozmix(pcols,levsiz)! ozone mixing ratio

!  real(kind=RKIND), intent(out) :: o3vmr(pcols,pver) ! ozone volume mixing ratio

   real(kind=RKIND), intent(in) :: pmid_in(pcols,pver)   ! level pressures (mks)
   real(kind=RKIND), intent(in) :: pin_in(levsiz)        ! ozone data level pressures (mks)
   real(kind=RKIND), intent(in) :: ozmix_in(pcols,levsiz)! ozone mixing ratio

   real(kind=RKIND), intent(out) :: o3vmr_out(pcols,pver)! ozone volume mixing ratio
!ldf end.
!
! local storage
!
   integer i                   ! longitude index
   integer k, kk, kkstart      ! level indices
   integer kupper(pcols)       ! Level indices for interpolation
   integer kount               ! Counter
   integer lats(pcols)         ! latitude indices
   integer lons(pcols)         ! latitude indices

   real(kind=RKIND) dpu        ! upper level pressure difference
   real(kind=RKIND) dpl        ! lower level pressure difference

!ldf begin:
   real(kind=RKIND):: pmid(pcols,pver)                ! level pressures (mks)
   real(kind=RKIND):: pin(levsiz)                     ! ozone data level pressures (mks)
   real(kind=RKIND):: ozmix(pcols,levsiz)             ! ozone mixing ratio
   real(kind=RKIND):: o3vmr(pcols,pver) ! ozone volume mixing ratio
!
! Initialize latitude indices
!
!  call get_lat_all_p(lchnk, ncol, lats)
!  call get_lon_all_p(lchnk, ncol, lons)
!

!ldf begin:
   do k = 1,levsiz
      pin(k) = pin_in(k)
   enddo
   do i = 1,pcols
      do k = 1,levsiz
         ozmix(i,k) = ozmix_in(i,k)
      enddo
   enddo
   do i = 1,pcols
      do k = 1,pver
         kk = pver-k+1
         pmid(i,kk) = pmid_in(i,k)
      enddo
   enddo
!ldf end.

! Initialize index array
!
   do i=1,ncol
      kupper(i) = 1
   end do

   do k=1,pver
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = levsiz
      do i=1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
      do kk=kkstart,levsiz-1
         do i=1,ncol
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
         if (kount.eq.ncol) then
            do i=1,ncol
               dpu = pmid(i,k) - pin(kupper(i))
               dpl = pin(kupper(i)+1) - pmid(i,k)
               o3vmr(i,k) = (ozmix(i,kupper(i))*dpl + &
                             ozmix(i,kupper(i)+1)*dpu)/(dpl + dpu)
            end do
            goto 35
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top ozone data level for at least some
! of the longitude points.
!
      do i=1,ncol
         if (pmid(i,k) .lt. pin(1)) then
            o3vmr(i,k) = ozmix(i,1)*pmid(i,k)/pin(1)
         else if (pmid(i,k) .gt. pin(levsiz)) then
            o3vmr(i,k) = ozmix(i,levsiz)
         else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            o3vmr(i,k) = (ozmix(i,kupper(i))*dpl + &
                          ozmix(i,kupper(i)+1)*dpu)/(dpl + dpu)
         end if
      end do

      if (kount.gt.ncol) then
!        call endrun ('VINTERP_OZN: Bad ozone data: non-monotonicity suspected')
      end if
35    continue
   end do

!ldf begin:
   do i = 1,pcols
   do k = 1,pver
      kk = pver-k+1
      o3vmr_out(i,kk) = o3vmr(i,k)
   enddo
   enddo
!ldf end.

   return
end subroutine vinterp_ozn

!==================================================================================================
 end module mpas_atmphys_o3climatology
!==================================================================================================
