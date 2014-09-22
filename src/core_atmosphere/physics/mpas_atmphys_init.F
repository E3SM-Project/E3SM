! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_init
 use mpas_kind_types
 use mpas_grid_types
 use mpas_timekeeping

 use mpas_atmphys_driver_convection, only: init_convection
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
!>    * modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-05-15.
!>    * added initialization of the accumulated surface pressure. Added initialization of the 
!>      tendency and accumulated tendency of the surface pressure.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-05-15.
!>    * renamed config_conv_deep_scheme to config_convection_scheme.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-09-18.


 contains


!==================================================================================================
 subroutine physics_init(dminfo,clock,configs,mesh,diag,tend,state,time_lev,diag_physics, &
                         atm_input,sfc_input)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(MPAS_Clock_type),intent(in):: clock

 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: tend
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: atm_input
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 logical,pointer:: config_do_restart,       &
                   config_o3climatology

 character(len=StrKIND),pointer::            &
                   config_lsm_scheme,        &
                   config_microp_scheme,     &
                   config_convection_scheme, &
                   config_sfclayer_scheme,   &
                   config_radt_lw_scheme,    &
                   config_radt_sw_scheme

 integer,pointer:: nCellsSolve,nLags
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

!local variables and arrays:
 type(MPAS_Time_Type):: currTime

 logical:: init_done
 integer:: ierr,julday 
 integer:: iCell,iLag,iEdge,nEdges_m

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_init:'

 call mpas_pool_get_config(configs,'config_do_restart'       ,config_do_restart       )
 call mpas_pool_get_config(configs,'config_o3climatology'    ,config_o3climatology    )
 call mpas_pool_get_config(configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call mpas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )
 call mpas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
 call mpas_pool_get_config(configs,'config_sfclayer_scheme'  ,config_sfclayer_scheme  )
 call mpas_pool_get_config(configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nLags'      ,nLags      )

 call mpas_pool_get_array(mesh,'nEdgesOnCell',nEdgesOnCell)
 call mpas_pool_get_array(mesh,'edgesOnCell' ,edgesOnCell )
 call mpas_pool_get_array(mesh,'dcEdge'      ,dcEdge      )

 call mpas_pool_get_array(diag_physics,'dcEdge_m'    ,dcEdge_m    )
 call mpas_pool_get_array(diag_physics,'i_rainc'     ,i_rainc     )
 call mpas_pool_get_array(diag_physics,'i_rainnc'    ,i_rainnc    )

 call mpas_pool_get_array(diag_physics,'i_acswdnb'   ,i_acswdnb   )
 call mpas_pool_get_array(diag_physics,'i_acswdnbc'  ,i_acswdnbc  )
 call mpas_pool_get_array(diag_physics,'i_acswdnt'   ,i_acswdnt   )
 call mpas_pool_get_array(diag_physics,'i_acswdntc'  ,i_acswdntc  )
 call mpas_pool_get_array(diag_physics,'i_acswupb'   ,i_acswupb   )
 call mpas_pool_get_array(diag_physics,'i_acswupbc'  ,i_acswupbc  )
 call mpas_pool_get_array(diag_physics,'i_acswupt'   ,i_acswupt   )
 call mpas_pool_get_array(diag_physics,'i_acswuptc'  ,i_acswuptc  )
 call mpas_pool_get_array(diag_physics,'i_aclwdnb'   ,i_aclwdnb   )
 call mpas_pool_get_array(diag_physics,'i_aclwdnbc'  ,i_aclwdnbc  )
 call mpas_pool_get_array(diag_physics,'i_aclwdnt'   ,i_aclwdnt   )
 call mpas_pool_get_array(diag_physics,'i_aclwdntc'  ,i_aclwdntc  )
 call mpas_pool_get_array(diag_physics,'i_aclwupb'   ,i_aclwupb   )
 call mpas_pool_get_array(diag_physics,'i_aclwupbc'  ,i_aclwupbc  )
 call mpas_pool_get_array(diag_physics,'i_aclwupt'   ,i_aclwupt   )
 call mpas_pool_get_array(diag_physics,'i_aclwuptc'  ,i_aclwuptc  )

 call mpas_pool_get_array(diag_physics,'acswdnb'     ,acswdnb     )
 call mpas_pool_get_array(diag_physics,'acswdnbc'    ,acswdnbc    )
 call mpas_pool_get_array(diag_physics,'acswdnt'     ,acswdnt     )
 call mpas_pool_get_array(diag_physics,'acswdntc'    ,acswdntc    )
 call mpas_pool_get_array(diag_physics,'acswupb'     ,acswupb     )
 call mpas_pool_get_array(diag_physics,'acswupbc'    ,acswupbc    )
 call mpas_pool_get_array(diag_physics,'acswupt'     ,acswupt     )
 call mpas_pool_get_array(diag_physics,'acswuptc'    ,acswuptc    )
 call mpas_pool_get_array(diag_physics,'aclwdnb'     ,aclwdnb     )
 call mpas_pool_get_array(diag_physics,'aclwdnbc'    ,aclwdnbc    )
 call mpas_pool_get_array(diag_physics,'aclwdnt'     ,aclwdnt     )
 call mpas_pool_get_array(diag_physics,'aclwdntc'    ,aclwdntc    )
 call mpas_pool_get_array(diag_physics,'aclwupb'     ,aclwupb     )
 call mpas_pool_get_array(diag_physics,'aclwupbc'    ,aclwupbc    )
 call mpas_pool_get_array(diag_physics,'aclwupt'     ,aclwupt     )
 call mpas_pool_get_array(diag_physics,'aclwuptc'    ,aclwuptc    )

 call mpas_pool_get_array(diag_physics,'nsteps_accum',nsteps_accum)
 call mpas_pool_get_array(diag_physics,'ndays_accum' ,ndays_accum )
 call mpas_pool_get_array(diag_physics,'tday_accum'  ,tday_accum  )
 call mpas_pool_get_array(diag_physics,'tyear_accum' ,tyear_accum )
 call mpas_pool_get_array(diag_physics,'tyear_mean'  ,tyear_mean  )
 call mpas_pool_get_array(diag_physics,'tlag'        ,tlag        )
 call mpas_pool_get_array(diag_physics,'sstsk'       ,sstsk       )
 call mpas_pool_get_array(diag_physics,'xicem'       ,xicem       )

 call mpas_pool_get_array(sfc_input,'sst' ,sst )
 call mpas_pool_get_array(sfc_input,'tmn' ,tmn )
 call mpas_pool_get_array(sfc_input,'xice',xice)

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,DoY=julday,ierr=ierr)

!calculation of the mean distance between cell centers:
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       dcEdge_m(iCell) = 0._RKIND
       nEdges_m = nEdgesOnCell(iCell)
       do iEdge = 1, nEdges_m
          dcEdge_m(iCell) = dcEdge_m(iCell) + dcEdge(edgesOnCell(iEdge,iCell))
       enddo
       dcEdge_m(iCell) = dcEdge_m(iCell) / nEdges_m
!      write(0,102) iCell,nEdges_m,(dcEdge(edgesOnCell(iEdge,iCell)),iEdge=1,nEdges_m),dcEdge_m(iCell)
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
 call landuse_init_forMPAS(dminfo,julday,mesh,configs,diag_physics,sfc_input)

!initialization of parameterized convective processes:
 if(config_convection_scheme .ne. 'off') &
    call init_convection(mesh,configs,diag_physics)

!initialization of cloud microphysics processes:
 if(config_microp_scheme .ne. 'off') call microphysics_init

!initialization of surface layer processes:
 if(config_sfclayer_scheme .ne. 'off') call init_sfclayer

!initialization of land-surface model:
!if(.not. config_do_restart) then
!   if(config_lsm_scheme .ne. 'off') call init_lsm(dminfo,mesh,configs,diag_physics,sfc_input)
!endif
 if(config_lsm_scheme .ne. 'off') call init_lsm(dminfo,mesh,configs,diag_physics,sfc_input)

!initialization of shortwave radiation processes:
 init_done = .false.
 if(config_radt_sw_scheme.ne.'off') then
    if(trim(config_radt_sw_scheme) .eq. 'cam_sw') then
       call init_radiation_sw(dminfo,mesh,atm_input,diag,state,time_lev)
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
          call init_radiation_lw(dminfo,mesh,atm_input,diag,state,time_lev)
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
 type(mpas_pool_type),intent(in):: mesh

!local pointers:
 integer,pointer:: nCells
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north

!local variables:
 integer:: iCell

!---------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'east'   ,east   )
 call mpas_pool_get_array(mesh,'north'  ,north  )

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
