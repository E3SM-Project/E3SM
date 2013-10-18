! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define DM_BCAST_CHAR(A) call mpas_dmpar_bcast_char(dminfo,A)
#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_reals(dminfo,size(A),A)
#define DM_BCAST_INTEGER(A) call mpas_dmpar_bcast_int(dminfo,A)

!==================================================================================================
 module mpas_atmphys_landuse
 use mpas_configure,only: input_landuse_data, &
                          config_sfc_albedo,  &
                          config_frac_seaice, &
                          config_do_restart
 use mpas_dmpar
 use mpas_grid_types

 use mpas_atmphys_utilities
 use mpas_atmphys_vars
 
 implicit none
 private
 public:: landuse_init_forMPAS

!global variables:
 integer,public:: isice,iswater

 integer,parameter:: frac_seaice      = 0       ! = 1: treats seaice as fractional field.
                                                ! = 0: ice/no-ice flag. 

!>\brief
!> This module reads the file LANDUSE.TBL which defines the land type of each cell, depending on
!> the origin of the input data, as defined by the value of the variable "sfc_input_data".
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_landuse:
!> ------------------------------------
!> landuse_init_forMPAS: main initialization for land use types.

!> The allowed values for sfc_input_data are:
!> input_sfc_data = OLD.                      (13 land types / summer and winter).
!> input_sfc_data = USGS.                     (33 land types / summer and winter).
!> input_sfc_data = MODIFIED_IGBP_MODIS_NOAH  (33 land types / summer and winter).
!> input_sfc_data = SiB                       (16 land types / summer and winter).
!> input_sfc_data = LW12                      ( 3 land types / all seasons).
!>
!> Given the value of the input index lu_index, and the julian day julday, landuse_init_forMPAS
!> initializes the variables:
!> .. background roughness length     (z0).
!> .. background surface albedo       (sfc_albbck).        
!> .. background surface emissivity   (sfc_emibck).
!> .. roughness length                (znt).
!> .. surface albedo                  (sfc_albedo).
!> .. surface emissivity              (sfc_emiss).
!> .. land mask                       (xland).
!> .. thermal inertia                 (thc).
!> .. surface moisture availability   (mavail).
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * deleted initialization of xicem and xland.
!>      - xicem is now initialized in physics_init.
!>      - xland is now initialized in physics_initialize_real and updated in physics_update_sst if
!>        needed.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-08-24. 
!>    * added initialization of the background surface albedo over snow.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-10-19.


 contains

!==================================================================================================
 subroutine landuse_init_forMPAS(dminfo,julday,mesh,diag_physics,sfc_input)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mesh_type),intent(in):: mesh
 type(diag_physics_type),intent(in):: diag_physics
 type(sfc_input_type)   ,intent(in):: sfc_input

 integer,intent(in):: julday

!local variables:
 character(len=StrKIND) :: lutype
 character(len=StrKIND):: mess

 integer,parameter:: land_unit = 15
 integer,parameter:: open_ok   = 0
 integer,parameter:: max_cats  = 100
 integer,parameter:: max_seas  = 12

 integer:: ierr,istat
 integer:: ic,is,isn,lucats,lumatch,luseas
 integer:: iCell,nCells
 integer:: julday_init
 integer,dimension(:),pointer:: ivgtyp
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND):: li
 real(kind=RKIND),dimension(max_cats,max_seas):: albd,slmo,sfem,sfz0,therin,scfx,sfhc

 real(kind=RKIND),dimension(:),pointer:: latCell
 real(kind=RKIND),dimension(:),pointer:: snoalb,snowc,xice
 real(kind=RKIND),dimension(:),pointer:: albbck,embck,xicem,xland,z0
 real(kind=RKIND),dimension(:),pointer:: mavail,sfc_albedo,sfc_emiss,thc,znt

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine landuse_init_forMPAS: julian day=', julday
 write(0,*) '--- config_frac_seaice = ',config_frac_seaice
 write(0,*) '--- xice_threshold     = ',xice_threshold

 nCells = mesh % nCells
 latCell  => mesh % latCell % array
 
 landmask => sfc_input % landmask   % array
 ivgtyp   => sfc_input % ivgtyp     % array
 snoalb   => sfc_input % snoalb     % array
 snowc    => sfc_input % snowc      % array
 xice     => sfc_input % xice       % array
 xland    => sfc_input % xland      % array
 albbck   => sfc_input % sfc_albbck % array

 embck      => diag_physics % sfc_emibck % array
 mavail     => diag_physics % mavail     % array
 sfc_albedo => diag_physics % sfc_albedo % array
 sfc_emiss  => diag_physics % sfc_emiss  % array
 thc        => diag_physics % thc        % array
 xicem      => diag_physics % xicem      % array
 z0         => diag_physics % z0         % array
 znt        => diag_physics % znt        % array
 
!reads in the landuse properties from landuse.tbl:
 if(dminfo % my_proc_id == IO_NODE) then
    open(land_unit,file='LANDUSE.TBL',action='READ',status='OLD',iostat=istat)
    if(istat /= open_ok) &
       call physics_error_fatal('subroutine landuse_init_forMPAS: ' // &
                                'failure opening LANDUSE.TBL')

    lumatch=0
    find_lutype : do while (lumatch == 0)
       read(unit=land_unit,fmt='(a35)') lutype
       read(unit=land_unit,fmt=*) lucats,luseas

       if(lutype .eq. input_landuse_data)then
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
!         write(0,101) ic,albd(ic,is),slmo(ic,is),sfem(ic,is),sfz0(ic,is), &
!                therin(ic,is),scfx(ic,is),sfhc(ic,is)
!      enddo
!      if(is .lt. luseas) write(0,*)
    enddo

!defines the index iswater and isice as a function of sfc_input_data:
    sfc_input_select: select case(trim(lutype))
       case('OLD')
          iswater = 7
          isice   = 11
       case('USGS')
          iswater = 16
          isice   = 24
       case('MODIFIED_IGBP_MODIS_NOAH')
          iswater = 17
          isice   = 15
       case('SiB')
          iswater = 15
          isice   = 16
       case('LW12')
          iswater = 2
          isice   = 3
       case default
    end select sfc_input_select
 endif

 DM_BCAST_CHAR(lutype)
 DM_BCAST_INTEGER(luseas)
 DM_BCAST_INTEGER(lucats)
 DM_BCAST_INTEGER(iswater)
 DM_BCAST_INTEGER(isice)
 DM_BCAST_MACRO(albd)
 DM_BCAST_MACRO(slmo)
 DM_BCAST_MACRO(sfem)
 DM_BCAST_MACRO(sfz0)
 DM_BCAST_MACRO(therin)
 DM_BCAST_MACRO(sfhc)
 DM_BCAST_MACRO(scfx)
 write(0,*) '--- isice   =',isice
 write(0,*) '--- iswater =',iswater
 if(config_do_restart) then
    write(0,*) '--- config_do_restart =', config_do_restart
    write(0,*) '--- skip the end of landuse_init_forMPAS'
    return
 endif

!defines the surface properties over the entire domain:
 do iCell = 1, nCells

    !finds the season as function of julian day (summer=1, winter=2): summer in the Northern
    !Hemisphere is defined between March 15th and September 15th (winter otherwise).
    isn = 1
    if(julday.lt.105 .or. julday.ge.288) isn=2
    if(latCell(iCell) .lt. 0.) isn=3-isn

!   is = nint(ivgtyp(iCell))
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
    znt(iCell)    = z0(iCell)
    mavail(iCell) = slmo(is,isn)
    embck(iCell)  = sfem(is,isn)
    sfc_emiss(iCell) = embck(iCell)

    !set sea-ice points to land with ice/snow surface properties:
    if(xice(iCell) .ge. xice_threshold) then
       albbck(iCell) = albd(isice,isn) / 100.
       embck(iCell)  = sfem(isice,isn)
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
       mavail(iCell) = slmo(isice,isn)
    endif

 enddo

 write(0,*) '--- end subroutine landuse_init_forMPAS'

!formats:
 101 format(i6,8(1x,e15.8))
 102 format(i6,1x,8e15.8)

 end subroutine landuse_init_forMPAS

!==================================================================================================
 end module mpas_atmphys_landuse
!==================================================================================================
