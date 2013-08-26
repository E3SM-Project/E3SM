! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_init
 use mpas_grid_types
 use mpas_configure, only: config_do_restart,          &
                           config_lsm_scheme,          &
                           config_microp_scheme,       &
                           config_conv_deep_scheme,    &
                           config_sfclayer_scheme,     &
                           config_radt_lw_scheme,      &
                           config_radt_sw_scheme,      &
                           config_o3climatology
 use mpas_timekeeping

 use mpas_atmphys_driver_convection_deep, only: init_convection_deep
 use mpas_atmphys_driver_lsm,only: init_lsm
 use mpas_atmphys_driver_microphysics
 use mpas_atmphys_driver_radiation_lw, only: init_radiation_lw
 use mpas_atmphys_driver_radiation_sw, only: init_radiation_sw
 use mpas_atmphys_driver_sfclayer

 use mpas_atmphys_landuse
 use mpas_atmphys_o3climatology

 implicit none
 private
 public:: physics_init

!>\brief MPAS main initialization subroutine for all physics parameterizations.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_init:
!> ---------------------------------
!> physics_init    : call initialization of individual physics parameterizations.
!> init_dir_forphys: needed for initialization of "reconstruct" subroutines.
!> r3_normalize    : needed for initialization of "reconstruct" subroutines.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * added structure diag in calls to subroutine init_radiation_lw and init_radiation_sw.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-07-01.
!>    * added call to subroutine init_o3climatology. reads monthly-mean climatological ozone data
!>      and interpolates ozone data to the MPAS grid.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-07-03.
!>    * added the calculation of the mean distance between cell centers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-08-22.
!>    * added initialization of variable xicem.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-08-24.


 contains

!==================================================================================================
 subroutine physics_init(dminfo,clock,config_do_restart,mesh,diag,state_1,state_2, &
                         diag_physics,atm_input,sfc_input)
!==================================================================================================

!input arguments:
 logical,intent(in):: config_do_restart
 type (dm_info), intent(in):: dminfo
 type(mesh_type),intent(in):: mesh
 type(diag_type),intent(in):: diag
 type(MPAS_Clock_type),intent(in):: clock

!inout arguments:
 type(state_type),intent(inout):: state_1,state_2
 type(diag_physics_type),intent(inout):: diag_physics
 type(atm_input_type),intent(inout):: atm_input
 type(sfc_input_type),intent(inout):: sfc_input

!local variables:
 type(MPAS_Time_Type):: currTime

 logical:: init_done
 integer:: iCell,iLag,ierr,julday 
 integer:: iEdge,nEdges

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_init:'

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,DoY=julday,ierr=ierr)

!calculation of the mean distance between cell centers:
 if(.not. config_do_restart) then
    do iCell = 1, mesh % nCellsSolve
       diag_physics % dcEdge_m % array(iCell) = 0._RKIND
       nEdges = mesh%nEdgesOnCell%array(iCell)
       do iEdge = 1, nEdges
          diag_physics%dcEdge_m%array(iCell) = diag_physics%dcEdge_m%array(iCell) &
                       + mesh%dcEdge%array(mesh%edgesOnCell%array(iEdge,iCell))
       enddo
       diag_physics % dcEdge_m % array(iCell) = diag_physics % dcEdge_m % array(iCell) &
                       / nEdges
!      write(0,102) iCell,nEdges,(mesh%dcEdge%array(mesh%edgesOnCell%array(iEdge,iCell)), &
!                   iEdge=1,nEdges),diag_physics%dcEdge_m%array(iCell)
    enddo
 endif
 101 format(8i9,10(1x,e15.8))
 102 format(2i9,10(1x,e15.8))

!initialization of east-north directions to convert u-tendencies from cell centers to cell
!edges:
 call init_dirs_forphys(mesh)

!initialization of counters i_rainc and i_rainnc. i_rainc and i_rainnc track the number of
!times the accumulated convective (rainc) and grid-scale (rainnc) rain exceed the prescribed
!threshold value:
 if(.not. config_do_restart) then
    do iCell = 1, mesh % nCellsSolve
       diag_physics % i_rainc  % array(iCell) = 0
       diag_physics % i_rainnc % array(iCell) = 0
    enddo
 endif

!initialization of counters i_acsw* and i_aclw*. i_acsw* and i_aclw* track the number of times
!the accumulated long and short-wave radiation fluxes exceed their prescribed theshold values.
 if(.not. config_do_restart) then
    do iCell = 1, mesh % nCellsSolve
       diag_physics % i_acswdnb  % array(iCell) = 0
       diag_physics % i_acswdnbc % array(iCell) = 0
       diag_physics % i_acswdnt  % array(iCell) = 0
       diag_physics % i_acswdntc % array(iCell) = 0
       diag_physics % i_acswupb  % array(iCell) = 0
       diag_physics % i_acswupbc % array(iCell) = 0
       diag_physics % i_acswupt  % array(iCell) = 0
       diag_physics % i_acswuptc % array(iCell) = 0

       diag_physics % i_aclwdnb  % array(iCell) = 0
       diag_physics % i_aclwdnbc % array(iCell) = 0
       diag_physics % i_aclwdnt  % array(iCell) = 0
       diag_physics % i_aclwdntc % array(iCell) = 0
       diag_physics % i_aclwupb  % array(iCell) = 0
       diag_physics % i_aclwupbc % array(iCell) = 0
       diag_physics % i_aclwupt  % array(iCell) = 0
       diag_physics % i_aclwuptc % array(iCell) = 0

       diag_physics % acswdnb  % array(iCell) = 0._RKIND
       diag_physics % acswdnbc % array(iCell) = 0._RKIND
       diag_physics % acswdnt  % array(iCell) = 0._RKIND
       diag_physics % acswdntc % array(iCell) = 0._RKIND
       diag_physics % acswupb  % array(iCell) = 0._RKIND
       diag_physics % acswupbc % array(iCell) = 0._RKIND
       diag_physics % acswupt  % array(iCell) = 0._RKIND
       diag_physics % acswuptc % array(iCell) = 0._RKIND

       diag_physics % aclwdnb  % array(iCell) = 0._RKIND
       diag_physics % aclwdnbc % array(iCell) = 0._RKIND
       diag_physics % aclwdnt  % array(iCell) = 0._RKIND
       diag_physics % aclwdntc % array(iCell) = 0._RKIND
       diag_physics % aclwupb  % array(iCell) = 0._RKIND
       diag_physics % aclwupbc % array(iCell) = 0._RKIND
       diag_physics % aclwupt  % array(iCell) = 0._RKIND
       diag_physics % aclwuptc % array(iCell) = 0._RKIND
    enddo
 endif

!initialization of xicem:
 if(.not.config_do_restart) then
    write(0,*) '--- initialization of xicem:'
    do iCell = 1, mesh % nCellsSolve
       diag_physics % xicem % array(iCell) = sfc_input % xice % array(iCell)
    enddo
 endif

!initialization of the local sea-surface temperature when a diurnal cycle of the
!sea-surface temperature is applied. This avoids having the array sstsk equal to
!zero over land:
 if(.not. config_do_restart) then
    write(0,*) '--- initialization of sstsk:'
    do iCell = 1, mesh % nCellsSolve
       diag_physics % sstsk % array(iCell) = sfc_input % sst % array(iCell)
    enddo
 endif

!initialization of temperatures needed for updating the deep soil temperature:
 if(.not. config_do_restart) then
    do iCell = 1, mesh % nCellsSolve
       diag_physics % nsteps_accum % array(iCell) = 0._RKIND
       diag_physics % ndays_accum  % array(iCell) = 0._RKIND
       diag_physics % tday_accum   % array(iCell) = 0._RKIND
       diag_physics % tyear_accum  % array(iCell) = 0._RKIND
       diag_physics % tyear_mean   % array(iCell) = sfc_input % tmn % array(iCell)
       do iLag = 1, mesh % nLags
          diag_physics % tlag % array(iLag,iCell) = sfc_input % tmn % array(iCell)
       enddo
    enddo
 endif

!read the input files that contain the monthly-mean ozone climatology on fixed pressure levels:
 if(config_o3climatology .and. (.not. config_do_restart)) &
    call init_o3climatology(mesh,atm_input)

!initialization of global surface properties. set here for now, but may be moved when time
!manager is implemented:
 call landuse_init_forMPAS(dminfo,julday,mesh,diag_physics,sfc_input)

!initialization of parameterized deep convective processes:
 if(config_conv_deep_scheme .ne. 'off') &
    call init_convection_deep(config_do_restart,mesh,diag_physics)

!initialization of cloud microphysics processes:
 if(config_microp_scheme .ne. 'off') call microphysics_init

!initialization of surface layer processes:
 if(config_sfclayer_scheme .ne. 'off') call init_sfclayer

!initialization of land-surface model:
!if(.not. config_do_restart) then
!   if(config_lsm_scheme .ne. 'off') call init_lsm(dminfo,mesh,diag_physics,sfc_input)
!endif
 if(config_lsm_scheme .ne. 'off') call init_lsm(dminfo,mesh,diag_physics,sfc_input)

!initialization of shortwave radiation processes:
 init_done = .false.
 if(config_radt_sw_scheme.ne.'off') then
    if(trim(config_radt_sw_scheme) .eq. 'cam_sw') then
       call init_radiation_sw(dminfo,mesh,atm_input,diag,state_1,state_2)
       init_done = .true.
    else
       call init_radiation_sw(dminfo)
    endif
 endif

!initialization of longwave radiation processes: if we run the CAM radiation codes, the initia
!lization of the longwave and shortwave parameterizations is the same, and needs to be called
!only once:
 if(config_radt_lw_scheme.ne.'off') then    
    if(trim(config_radt_lw_scheme) .eq. 'cam_lw') then
       if(.not. init_done) then
          call init_radiation_lw(dminfo,mesh,atm_input,diag,state_1,state_2)
       else
          write(0,*)
          write(0,*) '--- camrad lw initialization done above'
       endif
    else
       call init_radiation_lw(dminfo)
    endif
 endif

 write(0,*)
 write(0,*) '--- end subroutine physics_init'
 write(0,*)

 end subroutine physics_init

!==================================================================================================
 subroutine init_dirs_forphys(mesh)
!==================================================================================================

!inout arguments:
!----------------
 type(mesh_type),intent(in):: mesh

!local variables:
 integer:: iCell
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north

!---------------------------------------------------------------------------------------------

 latCell => mesh % latCell % array
 lonCell => mesh % lonCell % array
 east    => mesh % east  % array
 north   => mesh % north % array

!Compute unit vectors in east and north directions for each cell:
 do iCell = 1, mesh % nCells

    east(1,iCell) = -sin(lonCell(iCell))
    east(2,iCell) =  cos(lonCell(iCell))
    east(3,iCell) =  0.0
    call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

    north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
    north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
    north(3,iCell) =  cos(latCell(iCell))
    call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

 end do

 end subroutine init_dirs_forphys

!==================================================================================================
 subroutine r3_normalize(ax, ay, az)
!==================================================================================================
!normalizes the vector (ax, ay, az)

 real (kind=RKIND), intent(inout) :: ax, ay, az
 real (kind=RKIND) :: mi

!--------------------------------------------------------------------------------------------------

 mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
 ax = ax * mi
 ay = ay * mi
 az = az * mi

 end subroutine r3_normalize

!==================================================================================================
 end module mpas_atmphys_init
!==================================================================================================
