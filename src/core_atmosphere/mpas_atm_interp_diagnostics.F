! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atm_interp_diagnostics
 use mpas_kind_types
 use mpas_grid_types
 use mpas_constants

 implicit none
 private
 public:: interp_diagnostics

 contains

!==================================================================================================
 subroutine interp_diagnostics(mesh, state, time_lev, diag, diag_physics)
!==================================================================================================

!input arguments:
 type (mpas_pool_type), intent(in)  :: mesh
 type (mpas_pool_type), intent(in) :: state
 integer, intent(in) :: time_lev              ! which time level to use from state

!inout arguments:
 type (mpas_pool_type), intent(inout) :: diag
 type (mpas_pool_type), intent(inout) :: diag_physics

!local variables:
 integer :: iCell,iVert,iVertD,k,kk
 integer, pointer :: nCells, nVertLevels, nVertices, vertexDegree
 integer :: nVertLevelsP1
 integer, pointer :: index_qv
 integer, dimension(:,:), pointer :: cellsOnVertex

 real (kind=RKIND), dimension(:), pointer :: areaTriangle
 real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
 
 real (kind=RKIND), dimension(:,:), pointer :: exner, height
 real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p 
 real (kind=RKIND), dimension(:,:), pointer :: qvapor, relhum, theta_m, vorticity
 real (kind=RKIND), dimension(:,:), pointer :: umeridional, uzonal, vvel
 real (kind=RKIND), dimension(:,:,:), pointer :: scalars

 real (kind=RKIND), dimension(:), pointer :: temperature_200hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_500hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_700hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_850hPa

 real (kind=RKIND), dimension(:), pointer :: relhum_200hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_500hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_700hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_850hPa

 real (kind=RKIND), dimension(:), pointer :: uzonal_200hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_500hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_700hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_850hPa

 real (kind=RKIND), dimension(:), pointer :: umeridional_200hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_500hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_700hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_850hPa

 real (kind=RKIND), dimension(:), pointer :: height_200hPa
 real (kind=RKIND), dimension(:), pointer :: height_500hPa
 real (kind=RKIND), dimension(:), pointer :: height_700hPa
 real (kind=RKIND), dimension(:), pointer :: height_850hPa

 real (kind=RKIND), dimension(:), pointer :: w_200hPa
 real (kind=RKIND), dimension(:), pointer :: w_500hPa
 real (kind=RKIND), dimension(:), pointer :: w_700hPa
 real (kind=RKIND), dimension(:), pointer :: w_850hPa

 real (kind=RKIND), dimension(:), pointer :: vorticity_200hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_500hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_700hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_850hPa

 real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressureCp1, pressure2, pressure_v, temperature

!local interpolated fields:
 integer :: nIntP
 real (kind=RKIND) :: w1,w2,z0,z1,z2
 real (kind=RKIND), dimension(:,:), allocatable :: field_in,press_in
 real (kind=RKIND), dimension(:,:), allocatable :: field_interp,press_interp
 
!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine interp_diagnostics:'

 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
 call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
 call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
 call mpas_pool_get_dimension(state, 'index_qv', index_qv)

 nVertLevelsP1 = nVertLevels + 1

 call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
 call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
 call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)

 call mpas_pool_get_array(mesh, 'zgrid', height)
 call mpas_pool_get_array(state, 'w', vvel, time_lev)
 call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
 call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

 qvapor  => scalars(index_qv,:,:)    !MGD does this actually work?
 
 call mpas_pool_get_array(diag, 'exner', exner)
 call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
 call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
 call mpas_pool_get_array(diag, 'vorticity', vorticity)
 call mpas_pool_get_array(diag, 'uReconstructMeridional', umeridional)
 call mpas_pool_get_array(diag, 'uReconstructZonal', uzonal)
 call mpas_pool_get_array(diag_physics, 'relhum', relhum)

 call mpas_pool_get_array(diag, 'temperature_200hPa', temperature_200hPa)
 call mpas_pool_get_array(diag, 'temperature_500hPa', temperature_500hPa)
 call mpas_pool_get_array(diag, 'temperature_700hPa', temperature_700hPa)
 call mpas_pool_get_array(diag, 'temperature_850hPa', temperature_850hPa)

 call mpas_pool_get_array(diag, 'relhum_200hPa', relhum_200hPa)
 call mpas_pool_get_array(diag, 'relhum_500hPa', relhum_500hPa)
 call mpas_pool_get_array(diag, 'relhum_700hPa', relhum_700hPa)
 call mpas_pool_get_array(diag, 'relhum_850hPa', relhum_850hPa)

 call mpas_pool_get_array(diag, 'uzonal_200hPa', uzonal_200hPa)
 call mpas_pool_get_array(diag, 'uzonal_500hPa', uzonal_500hPa)
 call mpas_pool_get_array(diag, 'uzonal_700hPa', uzonal_700hPa)
 call mpas_pool_get_array(diag, 'uzonal_850hPa', uzonal_850hPa)

 call mpas_pool_get_array(diag, 'umeridional_200hPa', umeridional_200hPa)
 call mpas_pool_get_array(diag, 'umeridional_500hPa', umeridional_500hPa)
 call mpas_pool_get_array(diag, 'umeridional_700hPa', umeridional_700hPa)
 call mpas_pool_get_array(diag, 'umeridional_850hPa', umeridional_850hPa)

 call mpas_pool_get_array(diag, 'height_200hPa', height_200hPa)
 call mpas_pool_get_array(diag, 'height_500hPa', height_500hPa)
 call mpas_pool_get_array(diag, 'height_700hPa', height_700hPa)
 call mpas_pool_get_array(diag, 'height_850hPa', height_850hPa)

 call mpas_pool_get_array(diag, 'w_200hPa', w_200hPa)
 call mpas_pool_get_array(diag, 'w_500hPa', w_500hPa)
 call mpas_pool_get_array(diag, 'w_700hPa', w_700hPa)
 call mpas_pool_get_array(diag, 'w_850hPa', w_850hPa)

 call mpas_pool_get_array(diag, 'vorticity_200hPa', vorticity_200hPa)
 call mpas_pool_get_array(diag, 'vorticity_500hPa', vorticity_500hPa)
 call mpas_pool_get_array(diag, 'vorticity_700hPa', vorticity_700hPa)
 call mpas_pool_get_array(diag, 'vorticity_850hPa', vorticity_850hPa)

 if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCells)      )
 if(.not.allocated(pressureCp1) ) allocate(pressureCp1(nVertLevels,nCells+1) )
 if(.not.allocated(pressure2)   ) allocate(pressure2(nVertLevelsP1,nCells)   )
 if(.not.allocated(pressure_v)  ) allocate(pressure_v(nVertLevels,nVertices) )
 if(.not.allocated(temperature) ) allocate(temperature(nVertLevels,nCells)   )

!calculation of total pressure at cell centers (at mass points):
 do iCell = 1, nCells
 do k = 1, nVertLevels
    pressure(k,iCell)    = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
    pressureCp1(k,iCell) = pressure(k,iCell)
 enddo
 enddo
 do iCell = nCells+1, nCells+1
 do k = 1, nVertLevels
    pressureCp1(k,iCell)   = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
 enddo
 enddo

!calculation of total pressure at cell centers (at vertical velocity points):
 k = nVertLevelsP1
 do iCell = 1, nCells
    z0 = height(k,iCell)
    z1 = 0.5*(height(k,iCell)+height(k-1,iCell)) 
    z2 = 0.5*(height(k-1,iCell)+height(k-2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pressure2(k,iCell) = exp(w1*log(pressure(k-1,iCell))+w2*log(pressure(k-2,iCell)))
 enddo
 do k = 2, nVertLevels
 do iCell = 1, nCells
    w1 = (height(k,iCell)-height(k-1,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
    w2 = (height(k+1,iCell)-height(k,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
    pressure2(k,iCell) = w1*pressure(k,iCell) + w2*pressure(k-1,iCell)
 enddo
 enddo
 k = 1
 do iCell = 1, nCells
    z0 = height(k,iCell)
    z1 = 0.5*(height(k,iCell)+height(k+1,iCell)) 
    z2 = 0.5*(height(k+1,iCell)+height(k+2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    pressure2(k,iCell) = w1*pressure(k,iCell)+w2*pressure(k+1,iCell)
 enddo

!calculation of total pressure at cell vertices (at mass points):
 do iVert = 1, nVertices
    pressure_v(:,iVert) = 0._RKIND

    do k = 1, nVertLevels
    do iVertD = 1, vertexDegree
       pressure_v(k,iVert) = pressure_v(k,iVert) &
               + kiteAreasOnVertex(iVertD,iVert)*pressureCp1(k,cellsOnVertex(iVertD,iVert))
    enddo
    pressure_v(k,iVert) = pressure_v(k,iVert) / areaTriangle(iVert)
    enddo
 enddo

!calculation of temperature at cell centers:
 do iCell = 1,nCells
 do k = 1,nVertLevels
    temperature(k,iCell) = (theta_m(k,iCell)/(1._RKIND+rvord*qvapor(k,iCell)))*exner(k,iCell) 
 enddo
 enddo

!interpolation to fixed pressure levels for fields located at cells centers and at mass points:
 nIntP = 4
 if(.not.allocated(field_interp)) allocate(field_interp(nCells,nIntP) )
 if(.not.allocated(press_interp)) allocate(press_interp(nCells,nIntP) )
 do iCell = 1, nCells
    press_interp(iCell,1) = 200.0_RKIND
    press_interp(iCell,2) = 500.0_RKIND
    press_interp(iCell,3) = 700.0_RKIND
    press_interp(iCell,4) = 850.0_RKIND
 enddo

 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    press_in(iCell,kk) = pressure(k,iCell)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
!... temperature:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = temperature(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 temperature_200hPa(1:nCells) = field_interp(1:nCells,1)
 temperature_500hPa(1:nCells) = field_interp(1:nCells,2)
 temperature_700hPa(1:nCells) = field_interp(1:nCells,3)
 temperature_850hPa(1:nCells) = field_interp(1:nCells,4)
 write(0,*) '--- end interpolate temperature:'

!... relative humidity:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = relhum(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 relhum_200hPa(1:nCells) = field_interp(1:nCells,1)
 relhum_500hPa(1:nCells) = field_interp(1:nCells,2)
 relhum_700hPa(1:nCells) = field_interp(1:nCells,3)
 relhum_850hPa(1:nCells) = field_interp(1:nCells,4)
 write(0,*) '--- end interpolate relative humidity:'

!... u zonal wind:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = uzonal(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 uzonal_200hPa(1:nCells) = field_interp(1:nCells,1)
 uzonal_500hPa(1:nCells) = field_interp(1:nCells,2)
 uzonal_700hPa(1:nCells) = field_interp(1:nCells,3)
 uzonal_850hPa(1:nCells) = field_interp(1:nCells,4)
 write(0,*) '--- end interpolate zonal wind:'

!... u meridional wind:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = umeridional(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 umeridional_200hPa(1:nCells) = field_interp(1:nCells,1)
 umeridional_500hPa(1:nCells) = field_interp(1:nCells,2)
 umeridional_700hPa(1:nCells) = field_interp(1:nCells,3)
 umeridional_850hPa(1:nCells) = field_interp(1:nCells,4)
 write(0,*) '--- end interpolate meridional wind:'

 if(allocated(field_in)) deallocate(field_in)
 if(allocated(press_in)) deallocate(press_in)

!interpolation to fixed pressure levels for fields located at cells centers and at vertical
!velocity points:
 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevelsP1))
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    press_in(iCell,kk) = pressure2(k,iCell)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevelsP1))
 !... height:
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    field_in(iCell,kk) = height(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
 height_200hPa(1:nCells) = field_interp(1:nCells,1)
 height_500hPa(1:nCells) = field_interp(1:nCells,2)
 height_700hPa(1:nCells) = field_interp(1:nCells,3)
 height_850hPa(1:nCells) = field_interp(1:nCells,4)
 write(0,*) '--- end interpolate height:'
 
!... vertical velocity
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    field_in(iCell,kk) = vvel(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
 w_200hPa(1:nCells) = field_interp(1:nCells,1)
 w_500hPa(1:nCells) = field_interp(1:nCells,2)
 w_700hPa(1:nCells) = field_interp(1:nCells,3)
 w_850hPa(1:nCells) = field_interp(1:nCells,4)
 write(0,*) '--- end interpolate vertical velocity:'
 
 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)

!interpolation to fixed pressure levels for fields located at cell vertices and at mass points:
 nIntP = 4
 if(.not.allocated(field_interp)) allocate(field_interp(nVertices,nIntP) )
 if(.not.allocated(press_interp)) allocate(press_interp(nVertices,nIntP) )
 do iVert = 1, nVertices
    press_interp(iVert,1) = 200.0_RKIND
    press_interp(iVert,2) = 500.0_RKIND
    press_interp(iVert,3) = 700.0_RKIND
    press_interp(iVert,4) = 850.0_RKIND
 enddo

 if(allocated(field_in)) deallocate(field_in)
 if(allocated(press_in)) deallocate(press_in)

 if(.not.allocated(press_in)) allocate(press_in(nVertices,nVertLevels))
 do iVert = 1, nVertices
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    press_in(iVert,kk) = pressure_v(k,iVert)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nVertices,nVertLevels))
!... relative vorticity:
 do iVert = 1, nVertices
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iVert,kk) = vorticity(k,iVert)
 enddo
 enddo
 call interp_tofixed_pressure(nVertices,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 vorticity_200hPa(1:nVertices) = field_interp(1:nVertices,1)
 vorticity_500hPa(1:nVertices) = field_interp(1:nVertices,2)
 vorticity_700hPa(1:nVertices) = field_interp(1:nVertices,3)
 vorticity_850hPa(1:nVertices) = field_interp(1:nVertices,4)
 write(0,*) '--- end interpolate relative vorticity:'

 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)
 if(allocated(pressure)    ) deallocate(pressure    )
 if(allocated(pressureCp1) ) deallocate(pressureCp1 )
 if(allocated(pressure2)   ) deallocate(pressure2   )
 if(allocated(pressure_v)  ) deallocate(pressure_v  )
 if(allocated(temperature) ) deallocate(temperature )

!formats:
! 201 format(i5,4(1x,e15.8))

 end subroutine interp_diagnostics

!==================================================================================================
 subroutine interp_tofixed_pressure(ncol,nlev_in,nlev_out,pres_in,field_in,pres_out,field_out)
!==================================================================================================

!input arguments:
 integer,intent(in):: ncol,nlev_in,nlev_out

 real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: pres_in,field_in
 real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: pres_out

!output arguments:
 real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out

!local variables:
! integer:: i1,i2,icol,k,kk
 integer:: icol,k,kk
 integer:: kkstart,kount
 integer,dimension(ncol):: kupper

 real(kind=RKIND):: dpl,dpu

!--------------------------------------------------------------------------------------------------

!formats:
! 201 format(i5,8(1x,e15.8))

!write(0,*)
!write(0,*) '--- enter subroutine interp_tofixed_pressure:'
!write(0,*) '... ncol     = ',ncol
!write(0,*) '... nlev_in  = ',nlev_in
!write(0,*) '... nlev_out = ',nlev_out
!i1=1 ; i2=ncol
!do k = 1, nlev_in
!   write(0,201) k,pres_in(i1,k),field_in(i1,k),pres_in(i2,k),field_in(i2,k)
!enddo
!write(0,*)

 do icol = 1, ncol
    kupper(icol) = 1
 enddo

 do k = 1, nlev_out 

    kkstart = nlev_in
    do icol = 1, ncol
       kkstart = min0(kkstart,kupper(icol))
    enddo
    kount = 0

    do kk = kkstart, nlev_in-1
       do icol = 1, ncol
          if(pres_out(icol,k).gt.pres_in(icol,kk).and.pres_out(icol,k).le.pres_in(icol,kk+1)) then
             kupper(icol) = kk
             kount = kount + 1
!            write(0,201) kupper(icol),pres_out(icol,k),pres_in(icol,kk),pres_in(icol,kk+1)
          endif
       enddo

       if(kount.eq.ncol) then
          do icol = 1, ncol
             dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
             dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
             field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                               + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
          end do
          goto 35
        end if
    enddo

    do icol = 1, ncol
       if(pres_out(icol,k) .lt. pres_in(icol,1)) then
          field_out(icol,k) = field_in(icol,1)*pres_out(icol,k)/pres_in(icol,1)
       elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in)) then
          field_out(icol,k) = field_in(icol,nlev_in)
       else
          dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
          dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
          field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                            + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
       endif
    enddo

 35 continue
!   write(0,201) kupper(i1),pres_out(i1,k),pres_in(i1,kupper(i1)),pres_in(i1,kupper(i1)+1), &
!                field_out(i1,k),field_in(i1,kupper(i1)),field_in(i1,kupper(i1)+1)
!   write(0,201) kupper(i2),pres_out(i2,k),pres_in(i2,kupper(i2)),pres_in(i2,kupper(i2)+1), &
!                field_out(i2,k),field_in(i2,kupper(i2)),field_in(i2,kupper(i2)+1)

 enddo

 end subroutine interp_tofixed_pressure

!==================================================================================================
 end module mpas_atm_interp_diagnostics
!==================================================================================================
