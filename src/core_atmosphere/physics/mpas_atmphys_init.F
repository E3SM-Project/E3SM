! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
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
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-22.


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

!local variables and arrays:
 type(MPAS_Time_Type):: currTime

 logical:: init_done
 integer:: ierr,julday 
 integer:: iCell,iLag,iEdge,nCellsSolve,nEdges,nLags

 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:),pointer  :: i_rainc,i_rainnc
 integer,dimension(:),pointer  :: i_acswdnb,i_acswdnbc,i_acswdnt,i_acswdntc, &
                                  i_acswupb,i_acswupbc,i_acswupt,i_acswuptc, &
                                  i_aclwdnb,i_aclwdnbc,i_aclwdnt,i_aclwdntc, &
                                  i_aclwupb,i_aclwupbc,i_aclwupt,i_aclwuptc
 integer,dimension(:,:),pointer:: edgesOnCell


 real(kind=RKIND),dimension(:),pointer  :: dcEdge,dcEdge_m
 real(kind=RKIND),dimension(:),pointer  :: acswdnb,acswdnbc,acswdnt,acswdntc,  &
                                           acswupb,acswupbc,acswupt,acswuptc,  &
                                           aclwdnb,aclwdnbc,aclwdnt,aclwdntc,  &
                                           aclwupb,aclwupbc,aclwupt,aclwuptc
 real(kind=RKIND),dimension(:),pointer  :: nsteps_accum,ndays_accum,tday_accum, &
                                           tyear_accum,tyear_mean
 real(kind=RKIND),dimension(:),pointer  :: sst,sstsk,tmn,xice,xicem
 real(kind=RKIND),dimension(:,:),pointer:: tlag

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_init:'

 nCellsSolve = mesh % nCellsSolve
 nLags       = mesh % nLags

 nEdgesOnCell => mesh         % nEdgesOnCell % array
 edgesOnCell  => mesh         % edgesOnCell  % array
 dcEdge       => mesh         % dcEdge       % array
 dcEdge_m     => diag_physics % dcEdge_m     % array
 
 i_rainc      => diag_physics % i_rainc      % array
 i_rainnc     => diag_physics % i_rainnc     % array

 i_acswdnb    => diag_physics % i_acswdnb    % array
 i_acswdnbc   => diag_physics % i_acswdnbc   % array
 i_acswdnt    => diag_physics % i_acswdnt    % array
 i_acswdntc   => diag_physics % i_acswdntc   % array
 i_acswupb    => diag_physics % i_acswupb    % array
 i_acswupbc   => diag_physics % i_acswupbc   % array
 i_acswupt    => diag_physics % i_acswupt    % array
 i_acswuptc   => diag_physics % i_acswuptc   % array
 i_aclwdnb    => diag_physics % i_aclwdnb    % array
 i_aclwdnbc   => diag_physics % i_aclwdnbc   % array
 i_aclwdnt    => diag_physics % i_aclwdnt    % array
 i_aclwdntc   => diag_physics % i_aclwdntc   % array
 i_aclwupb    => diag_physics % i_aclwupb    % array
 i_aclwupbc   => diag_physics % i_aclwupbc   % array
 i_aclwupt    => diag_physics % i_aclwupt    % array
 i_aclwuptc   => diag_physics % i_aclwuptc   % array

 acswdnb      => diag_physics % acswdnb      % array
 acswdnbc     => diag_physics % acswdnbc     % array
 acswdnt      => diag_physics % acswdnt      % array
 acswdntc     => diag_physics % acswdntc     % array
 acswupb      => diag_physics % acswupb      % array
 acswupbc     => diag_physics % acswupbc     % array
 acswupt      => diag_physics % acswupt      % array
 acswuptc     => diag_physics % acswuptc     % array
 aclwdnb      => diag_physics % aclwdnb      % array
 aclwdnbc     => diag_physics % aclwdnbc     % array
 aclwdnt      => diag_physics % aclwdnt      % array
 aclwdntc     => diag_physics % aclwdntc     % array
 aclwupb      => diag_physics % aclwupb      % array
 aclwupbc     => diag_physics % aclwupbc     % array
 aclwupt      => diag_physics % aclwupt      % array
 aclwuptc     => diag_physics % aclwuptc     % array

 nsteps_accum => diag_physics % nsteps_accum % array
 ndays_accum  => diag_physics % ndays_accum  % array
 tday_accum   => diag_physics % tday_accum   % array
 tyear_accum  => diag_physics % tyear_accum  % array
 tyear_mean   => diag_physics % tyear_mean   % array
 tlag         => diag_physics % tlag         % array
 sstsk        => diag_physics % sstsk        % array
 xicem        => diag_physics % xicem        % array

 sst          => sfc_input    % sst          % array
 tmn          => sfc_input    % tmn          % array
 xice         => sfc_input    % xice         % array


 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,DoY=julday,ierr=ierr)

!calculation of the mean distance between cell centers:
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       dcEdge_m(iCell) = 0._RKIND
       nEdges = nEdgesOnCell(iCell)
       do iEdge = 1, nEdges
          dcEdge_m(iCell) = dcEdge_m(iCell) + dcEdge(edgesOnCell(iEdge,iCell))
       enddo
       dcEdge_m(iCell) = dcEdge_m(iCell) / nEdges
!      write(0,102) iCell,nEdges,(dcEdge(edgesOnCell(iEdge,iCell)),iEdge=1,nEdges),dcEdge_m(iCell)
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
    do iCell = 1, nCellsSolve
       i_rainc(iCell)  = 0
       i_rainnc(iCell) = 0
    enddo
 endif

!initialization of counters i_acsw* and i_aclw*. i_acsw* and i_aclw* track the number of times
!the accumulated long and short-wave radiation fluxes exceed their prescribed theshold values.
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       i_acswdnb(iCell)  = 0
       i_acswdnbc(iCell) = 0
       i_acswdnt(iCell)  = 0
       i_acswdntc(iCell) = 0
       i_acswupb(iCell)  = 0
       i_acswupbc(iCell) = 0
       i_acswupt(iCell)  = 0
       i_acswuptc(iCell) = 0

       i_aclwdnb(iCell)  = 0
       i_aclwdnbc(iCell) = 0
       i_aclwdnt(iCell)  = 0
       i_aclwdntc(iCell) = 0
       i_aclwupb(iCell)  = 0
       i_aclwupbc(iCell) = 0
       i_aclwupt(iCell)  = 0
       i_aclwuptc(iCell) = 0

       acswdnb(iCell)  = 0._RKIND
       acswdnbc(iCell) = 0._RKIND
       acswdnt(iCell)  = 0._RKIND
       acswdntc(iCell) = 0._RKIND
       acswupb(iCell)  = 0._RKIND
       acswupbc(iCell) = 0._RKIND
       acswupt(iCell)  = 0._RKIND
       acswuptc(iCell) = 0._RKIND

       aclwdnb(iCell)  = 0._RKIND
       aclwdnbc(iCell) = 0._RKIND
       aclwdnt(iCell)  = 0._RKIND
       aclwdntc(iCell) = 0._RKIND
       aclwupb(iCell)  = 0._RKIND
       aclwupbc(iCell) = 0._RKIND
       aclwupt(iCell)  = 0._RKIND
       aclwuptc(iCell) = 0._RKIND
    enddo
 endif

!initialization of xicem:
 if(.not.config_do_restart) then
    write(0,*) '--- initialization of xicem:'
    do iCell = 1, nCellsSolve
       xicem(iCell) = xice(iCell)
    enddo
 endif

!initialization of the local sea-surface temperature when a diurnal cycle of the
!sea-surface temperature is applied. This avoids having the array sstsk equal to
!zero over land:
 if(.not. config_do_restart) then
    write(0,*) '--- initialization of sstsk:'
    do iCell = 1, nCellsSolve
       sstsk(iCell) = sst(iCell)
    enddo
 endif

!initialization of temperatures needed for updating the deep soil temperature:
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       nsteps_accum(iCell) = 0._RKIND
       ndays_accum(iCell)  = 0._RKIND
       tday_accum(iCell)   = 0._RKIND
       tyear_accum(iCell)  = 0._RKIND
       tyear_mean(iCell)   = tmn(iCell)
       do iLag = 1, nLags
          tlag(iLag,iCell) = tmn(iCell)
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
 integer:: iCell,nCells
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north

!---------------------------------------------------------------------------------------------

 nCells = mesh % nCells

 latCell => mesh % latCell % array
 lonCell => mesh % lonCell % array
 east    => mesh % east  % array
 north   => mesh % north % array

!Compute unit vectors in east and north directions for each cell:
 do iCell = 1, nCells

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
