! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define DM_BCAST_CHAR(A) call mpas_dmpar_bcast_char(dminfo,A)
#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_reals(dminfo,size(A),A)
#define DM_BCAST_INTEGER(A) call mpas_dmpar_bcast_int(dminfo,A)

!=================================================================================================================
 module mpas_atmphys_landuse
 use mpas_dmpar
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_io_units

 use mpas_atmphys_utilities
 use mpas_atmphys_vars
 
 implicit none
 private
 public:: landuse_init_forMPAS

!global variables:
 integer,public:: isurban

!This module reads the file LANDUSE.TBL which defines the land type of each cell, depending on
!the origin of the input data, as defined by the value of the variable "sfc_input_data".
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_landuse:
! ------------------------------------
! landuse_init_forMPAS: main initialization for land use types.

! The allowed values for sfc_input_data are:
! input_sfc_data = OLD.                      (13 land types / summer and winter).
! input_sfc_data = USGS.                     (33 land types / summer and winter).
! input_sfc_data = MODIFIED_IGBP_MODIS_NOAH  (33 land types / summer and winter).
! input_sfc_data = SiB                       (16 land types / summer and winter).
! input_sfc_data = LW12                      ( 3 land types / all seasons).
!
! Given the value of the input index lu_index, and the julian day julday, landuse_init_forMPAS
! initializes the variables:
! .. background roughness length     (z0).
! .. background surface albedo       (sfc_albbck).
! .. background surface emissivity   (sfc_emibck).
! .. roughness length                (znt).
! .. surface albedo                  (sfc_albedo).
! .. surface emissivity              (sfc_emiss).
! .. land mask                       (xland).
! .. thermal inertia                 (thc).
! .. surface moisture availability   (mavail).
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * deleted initialization of xicem and xland.
!   - xicem is now initialized in physics_init.
!   - xland is now initialized in physics_initialize_real and updated in physics_update_sst if needed.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-24.
! * added initialization of the background surface albedo over snow.
!   Laura D. Fowler (laura@ucar.edu) / 2013-10-19.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * in subroutine landuse_init_forMPAS, added the definition of isurban as a function of the
!   input landuse data file.
!   Dominikus Heinzeller (IMK) / 2014-07-24.
! * in subroutine landuse_int_forMPAS, added the initialization of variable ust to a very small value. this was
!   needed when the surface layer scheme was updated to that used in WRF version 3.8.1
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-27.
! * removed the parameter frac_seaice which is not used anymore and has been replaced with config_frac_seaice. 
!   Laura D. Fowler (laura@ucar.edu) / 2017-01-11.
! * now use isice and iswater initialized in the init file instead of initialized in mpas_atmphys_landuse.F.
!   Laura D. Fowler (laura@ucar.edu) / 2017-01-13.
! * added the initialization of sfc_albedo_seaice which is the surface albedo over seaice points.
!   Laura D. Fowler (laura@ucar.edu) / 2017-03-02.


 contains


!=================================================================================================================
 subroutine landuse_init_forMPAS(dminfo,julday,mesh,configs,diag_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input

 integer,intent(in):: julday

!local pointers:
 logical,pointer:: config_do_restart,  &
                   config_frac_seaice, &
                   config_sfc_albedo

 character(len=StrKIND),pointer:: mminlu

 integer,pointer:: nCells
 integer,pointer:: isice,iswater
 integer,dimension(:),pointer:: ivgtyp
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer:: latCell
 real(kind=RKIND),dimension(:),pointer:: snoalb,snowc,xice
 real(kind=RKIND),dimension(:),pointer:: albbck,embck,xicem,xland,z0
 real(kind=RKIND),dimension(:),pointer:: mavail,sfc_albedo,sfc_albedo_seaice,sfc_emiss,thc,ust,znt

!local variables:
 character(len=StrKIND) :: lutype
 character(len=StrKIND):: mess

 integer:: land_unit
 integer,parameter:: open_ok   = 0
 integer,parameter:: max_cats  = 100
 integer,parameter:: max_seas  = 12

 integer:: ierr,istat
 integer:: ic,is,isn,lucats,lumatch,luseas
 integer:: iCell
 integer:: julday_init

 real(kind=RKIND):: li
 real(kind=RKIND),dimension(max_cats,max_seas):: albd,slmo,sfem,sfz0,therin,scfx,sfhc

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine landuse_init_forMPAS:')

 call mpas_pool_get_config(configs,'config_do_restart' ,config_do_restart )
 call mpas_pool_get_config(configs,'config_frac_seaice',config_frac_seaice)
 call mpas_pool_get_config(configs,'config_sfc_albedo' ,config_sfc_albedo )

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 
 call mpas_pool_get_array(sfc_input,'mminlu'    ,mminlu  )
 call mpas_pool_get_array(sfc_input,'isice'     ,isice   )
 call mpas_pool_get_array(sfc_input,'iswater'   ,iswater )
 call mpas_pool_get_array(sfc_input,'landmask'  ,landmask)
 call mpas_pool_get_array(sfc_input,'ivgtyp'    ,ivgtyp  )
 call mpas_pool_get_array(sfc_input,'snoalb'    ,snoalb  )
 call mpas_pool_get_array(sfc_input,'snowc'     ,snowc   )
 call mpas_pool_get_array(sfc_input,'xice'      ,xice    )
 call mpas_pool_get_array(sfc_input,'xland'     ,xland   )
 call mpas_pool_get_array(sfc_input,'sfc_albbck',albbck  )

 nullify(mavail)
 nullify(ust)
 call mpas_pool_get_array(diag_physics,'sfc_emibck'       ,embck            )
 call mpas_pool_get_array(diag_physics,'mavail'           ,mavail           )
 call mpas_pool_get_array(diag_physics,'sfc_albedo'       ,sfc_albedo       )
 call mpas_pool_get_array(diag_physics,'sfc_albedo_seaice',sfc_albedo_seaice)
 call mpas_pool_get_array(diag_physics,'sfc_emiss'        ,sfc_emiss        )
 call mpas_pool_get_array(diag_physics,'thc'              ,thc              )
 call mpas_pool_get_array(diag_physics,'ust'              ,ust              )
 call mpas_pool_get_array(diag_physics,'xicem'            ,xicem            )
 call mpas_pool_get_array(diag_physics,'z0'               ,z0               )
 call mpas_pool_get_array(diag_physics,'znt'              ,znt              )

!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine landuse_init_forMPAS: julian day=$i' , intArgs=(/julday/))
!call mpas_log_write('--- config_frac_seaice = $1',logicArgs=(/config_frac_seaice/))
!call mpas_log_write('--- xice_threshold     = $r',realArgs=(/xice_threshold/))

!reads in the landuse properties from landuse.tbl:
 if(dminfo % my_proc_id == IO_NODE) then
    !get a unit to open init file:
    call mpas_new_unit(land_unit)
    if ( land_unit < 0 ) &
       call physics_error_fatal('landuse_init_forMPAS: All file units are taken.  Change maxUnits in mpas_io_units.F')

    open(land_unit,file='LANDUSE.TBL',action='READ',status='OLD',iostat=istat)
    if(istat /= open_ok) &
       call physics_error_fatal('subroutine landuse_init_forMPAS: ' // &
                                'failure opening LANDUSE.TBL')

    lumatch=0
    find_lutype : do while (lumatch == 0)
       read(unit=land_unit,fmt='(a35)') lutype
       read(unit=land_unit,fmt=*) lucats,luseas

       if(lutype .eq. mminlu)then
          write(mess,*) '   landuse type = ' // trim (lutype) // ' found', lucats, &
                        ' categories', luseas, ' seasons'
          call physics_message(mess)
          lumatch=1
       else
          write(mess,*) '   skipping over lutype = ' // trim (lutype)
          call physics_message(mess)          
          do is = 1,luseas
             read(unit=land_unit,fmt=*,iostat=ierr) 
             do ic = 1,lucats
                read(unit=land_unit,fmt=*,iostat=ierr)
             enddo
          enddo
       endif
    enddo find_lutype

    do is = 1, luseas
       read(unit=land_unit,fmt=*,iostat=ierr) 
       do ic = 1, lucats
          read(unit=land_unit,fmt=*) li,albd(ic,is),slmo(ic,is),sfem(ic,is),sfz0(ic,is), &
                                     therin(ic,is),scfx(ic,is),sfhc(ic,is)
       enddo
!      do ic = 1, lucats
!         call mpas_log_write('$i $r $r $r $r $r $r $r $r', intArgs=(/ic/), &
!                realArgs=(/albd(ic,is),slmo(ic,is),sfem(ic,is),sfz0(ic,is), &
!                therin(ic,is),scfx(ic,is),sfhc(ic,is)/))
!      enddo
!      if(is .lt. luseas) call mpas_log_write('')
    enddo

    close(land_unit)  
    call mpas_release_unit(land_unit)

!defines the index isurban as a function of sfc_input_data:
    sfc_input_select: select case(trim(lutype))
       case('OLD')
          isurban = 1
       case('USGS')
          isurban = 1
       case('MODIFIED_IGBP_MODIS_NOAH')
          isurban = 13
       case('SiB')
          isurban = 11
       case('LW12')
          isurban = 1
       case default
    end select sfc_input_select
 endif

 DM_BCAST_CHAR(lutype)
 DM_BCAST_INTEGER(luseas)
 DM_BCAST_INTEGER(lucats)
 DM_BCAST_INTEGER(iswater)
 DM_BCAST_INTEGER(isice)
 DM_BCAST_INTEGER(isurban)
 DM_BCAST_MACRO(albd)
 DM_BCAST_MACRO(slmo)
 DM_BCAST_MACRO(sfem)
 DM_BCAST_MACRO(sfz0)
 DM_BCAST_MACRO(therin)
 DM_BCAST_MACRO(sfhc)
 DM_BCAST_MACRO(scfx)
!call mpas_log_write('--- isice   =$i',intArgs=(/isice/))
!call mpas_log_write('--- iswater =$i',intArgs=(/iswater/))
!call mpas_log_write('--- isurban =$i',intArgs=(/isurban/))
 if(config_do_restart) then
    call mpas_log_write('--- config_do_restart =$l', logicArgs=(/config_do_restart/))
    call mpas_log_write('--- skip the end of landuse_init_forMPAS')
    return
 endif

!defines the surface properties over the entire domain:
 do iCell = 1, nCells

    !finds the season as function of julian day (summer=1, winter=2): summer in the Northern
    !Hemisphere is defined between March 15th and September 15th (winter otherwise).
    isn = 1
    if(julday.lt.105 .or. julday.ge.288) isn=2
    if(latCell(iCell) .lt. 0.) isn=3-isn

    is = ivgtyp(iCell)
    
    !set no data points to water:
    if(is.eq.0) is = iswater
    if(.not. config_sfc_albedo) albbck(iCell) = albd(is,isn)/100.
    sfc_albedo(iCell) = albbck(iCell)

    if(snowc(iCell) .gt. 0.5) then
       albbck(iCell) = albd(isice,isn) / 100.
       if(config_sfc_albedo) then
          sfc_albedo(iCell) = snoalb(iCell)
       else
          sfc_albedo(iCell) = albbck(iCell) / (1+scfx(is,isn))
       endif
    endif
    thc(iCell)    = therin(is,isn) / 100.
    z0(iCell)     = sfz0(is,isn) / 100.
    znt(iCell)    = z0(icell)
    embck(iCell)  = sfem(is,isn)
    sfc_emiss(iCell) = embck(iCell)

    if(associated(mavail)) mavail(iCell) = slmo(is,isn)
    if(associated(ust)) ust(iCell) = 0.0001
 
    !set sea-ice points to land with ice/snow surface properties:
    if(xice(iCell) .ge. xice_threshold) then
       albbck(iCell) = albd(isice,isn) / 100.
       embck(iCell)  = sfem(isice,isn)
       sfc_albedo_seaice(iCell) = albbck(iCell)
       if(config_frac_seaice) then
          !0.08 is the albedo over open water.
          !0.98 is the emissivity over open water.
          sfc_albedo(iCell) = xice(iCell)*albbck(iCell) + (1-xice(iCell))*0.08
          sfc_emiss(iCell)  = xice(iCell)*embck(iCell)  + (1-xice(iCell))*0.98
       else
          sfc_albedo(iCell) = albbck(iCell)
          sfc_emiss(iCell)  = embck(iCell)
       endif
       thc(iCell) = therin(isice,isn) / 100.
       z0(icell)  = sfz0(isice,isn) / 100.
       znt(iCell) = z0(iCell)

       if(associated(mavail)) mavail(iCell) = slmo(isice,isn)
    else
       !over seaice-free cells, initialize sfc_albedo_seaice with the background surface albedo:
       sfc_albedo_seaice(iCell) = albbck(iCell)
    endif

 enddo

!call mpas_log_write('--- end subroutine landuse_init_forMPAS')

 end subroutine landuse_init_forMPAS

!=================================================================================================================
 end module mpas_atmphys_landuse
!=================================================================================================================
