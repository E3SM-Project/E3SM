! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_initialize_real
 use mpas_kind_types
 use mpas_dmpar
 use mpas_grid_types
 use mpas_init_atm_surface

 use mpas_atmphys_date_time
 use mpas_atmphys_utilities

 implicit none
 private
 public:: physics_initialize_real


!>\brief MPAS initialization of surface properties for real case initialization.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_initialize_real:
!> --------------------------------------------
!> physics_initialize_real    : main subroutine (called from subroutine init_atm_setup_test_case in
!>                              ./src/core_init_atmosphere/mpas_init_atm_test_cases.F).
!> init_soil_layers           : main subroutine for initialization of soil properties. 
!> init_soil_layers_depth     : initialize height and depth of soil layers needed in NOAH scheme.
!> init_soil_layers_properties: initialize soil temperature, soil moisture, etc.
!> adjust_input_soiltemp      : adjust the deep soil temperature to sea-level values.
!> physics_init_sst           : initialize the skin temperature to the SSTs over oceans.
!> physics_init_seaice        : correct vegetation and soil typs as function of fractional sea ice.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * revised entire module:
!>      -> changed nCells to nCellsSolve in every subroutine.
!>      -> removed modifying snoalb (surface albedo over snow) over sea-ice points.
!>      -> revised subroutine physics_init_sst.
!>      -> revised subroutine physics_init_seaice.
!>      Laura D. Fowler (laura@ucar.edu) / 2013-08-02.
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
!>    * In subroutine physics_init_seaice, assign the sea-ice land use category as a function of
!>      the land use category input file (MODIS OR USGS).
!>      Dominikus Heinzeller (IMK) / 2014-07-24.


 contains

!==================================================================================================
 subroutine physics_initialize_real(mesh, fg, dminfo, dims, configs)
!==================================================================================================
!input arguments:
 type (mpas_pool_type), intent(in) :: mesh
 type (dm_info), intent(in) :: dminfo
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!inout arguments:
 type (mpas_pool_type), intent(inout) :: fg 

!local variables:
 character(len=StrKIND):: initial_date

 integer:: iCell
 integer, pointer :: nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:)  ,pointer:: sfc_albbck
 real(kind=RKIND),dimension(:,:),pointer:: albedo12m

 real(kind=RKIND),dimension(:),pointer:: seaice,xice,xland
 real(kind=RKIND),dimension(:),pointer:: vegfra,shdmin,shdmax
 real(kind=RKIND),dimension(:),pointer:: snow,snowc,snowh
 real(kind=RKIND),dimension(:,:),pointer:: greenfrac

 real(kind=RKIND),dimension(:),pointer:: skintemp,sst

 character (len=StrKIND), pointer :: config_sfc_prefix
 character (len=StrKIND), pointer :: config_start_time
 logical, pointer :: config_input_sst
 
!temporary:
 integer:: iSoil,nSoilLevels

!--------------------------------------------------------------------------------------------------

! write(0,*)
! write(0,*) '--- enter physics_initialize_real:'

 call mpas_pool_get_config(configs, 'config_sfc_prefix', config_sfc_prefix)
 call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
 call mpas_pool_get_config(configs, 'config_input_sst', config_input_sst)

 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)

 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(mesh, 'albedo12m', albedo12m)
 call mpas_pool_get_array(mesh, 'greenfrac', greenfrac)
 call mpas_pool_get_array(mesh, 'shdmin', shdmin)
 call mpas_pool_get_array(mesh, 'shdmax', shdmax)

 call mpas_pool_get_array(fg, 'sfc_albbck', sfc_albbck)
 call mpas_pool_get_array(fg, 'vegfra', vegfra)
 call mpas_pool_get_array(fg, 'snow', snow)
 call mpas_pool_get_array(fg, 'snowc', snowc)
 call mpas_pool_get_array(fg, 'snowh', snowh)
 call mpas_pool_get_array(fg, 'skintemp', skintemp)
 call mpas_pool_get_array(fg, 'sst', sst)
 call mpas_pool_get_array(fg, 'seaice', seaice)
 call mpas_pool_get_array(fg, 'xice', xice)
 call mpas_pool_get_array(fg, 'xland', xland)

!initialization of xland:
 do iCell = 1, nCellsSolve
    xland(iCell) = 0._RKIND
    if(landmask(iCell) == 1) then
       xland(iCell) = 1._RKIND
    elseif(landmask(iCell) == 0) then
       xland(iCell) = 2._RKIND
    endif
 enddo

!initialization of the sea-surface temperature and seaice if they are read from a separate
!input file. calling this subroutine will overwrite the arrays sst and seaice already read
!in the file defined by config_input_name:
 if(config_input_sst) then
    write(0,*) '--- read sea-surface temperature from auxillary file:'
    call interp_sfc_to_MPAS(config_start_time(1:13),mesh,fg,dims,dminfo,config_sfc_prefix)
    call physics_init_sst(mesh,fg,dims,configs)
 endif

!initialization of the surface background albedo: interpolation of the monthly values to the
!initial date:
 initial_date = trim(config_start_time)
 call monthly_interp_to_date(nCellsSolve,initial_date,albedo12m,sfc_albbck)

 do iCell = 1, nCellsSolve
    sfc_albbck(iCell) = sfc_albbck(iCell) / 100._RKIND
    if(landmask(iCell) .eq. 0) sfc_albbck(iCell) = 0.08_RKIND
 enddo

!initialization of the green-ness (vegetation) fraction: interpolation of the monthly values to
!the initial date. get the min/max for each cell for the monthly green-ness fraction:
 initial_date = trim(config_start_time)
 call monthly_interp_to_date(nCellsSolve,initial_date,greenfrac,vegfra)

!calculates the maximum and minimum green-ness (vegetation) fraction:
 call monthly_min_max(nCellsSolve,greenfrac,shdmin,shdmax)

!initialization of the flag indicating the presence of snow (0 or 1) and of the snow depth
!(m) as functions of the input snow water content (kg/m2). we use a 5:1 ratio from liquid
!water equivalent to snow depth:
 do iCell = 1, nCellsSolve
    if(snow(iCell) .ge. 10._RKIND) then
       snowc(iCell) = 1._RKIND
    else
       snowc(iCell) = 0._RKIND
    endif
    snowh(iCell) = snow(iCell) * 5.0_RKIND / 1000._RKIND
 enddo

!initialization of soil layers properties:
 call init_soil_layers(mesh,fg,dminfo,dims,configs)

!initialize seaice points:
 call physics_init_seaice(mesh,fg,dims,configs)
 
! write(0,*) '--- end physics_initialize_real:'

 end subroutine physics_initialize_real

!==================================================================================================
 subroutine init_soil_layers(mesh,fg,dminfo,dims,configs)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(dm_info),intent(in)  :: dminfo
 type(mpas_pool_type),intent(in):: dims
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!--------------------------------------------------------------------------------------------------

!adjust the annual mean deep soil temperature:
 call adjust_input_soiltemps(mesh,fg,dims)

!initialize the depth of the soil layers:
 call init_soil_layers_depth(mesh,fg,dims,configs)
 
!initialize the temperature, moisture, and liquid water of the individual soil layers:
 call init_soil_layers_properties(mesh,fg,dminfo,dims,configs)

 end subroutine init_soil_layers

!==================================================================================================
 subroutine adjust_input_soiltemps(mesh, fg, dims)
!==================================================================================================

!input arguments:
 type (mpas_pool_type), intent(in) :: mesh
 type (mpas_pool_type), intent(in) :: dims

!inout arguments:
 type (mpas_pool_type), intent(inout) :: fg 

!local variables:
 integer:: iCell,ifgSoil
 integer, pointer:: nCellsSolve,nFGSoilLevels
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer  :: soilz,ter
 real(kind=RKIND),dimension(:),pointer  :: skintemp,soiltemp,tmn
 real(kind=RKIND),dimension(:,:),pointer:: st_fg

!--------------------------------------------------------------------------------------------------


 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(dims, 'nFGSoilLevels', nFGSoilLevels)

 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(mesh, 'soiltemp', soiltemp)
 call mpas_pool_get_array(mesh, 'ter', ter)

 call mpas_pool_get_array(fg, 'skintemp', skintemp)
 call mpas_pool_get_array(fg, 'tmn', tmn)
 call mpas_pool_get_array(fg, 'st_fg', st_fg)
 call mpas_pool_get_array(fg, 'soilz', soilz)


 do iCell = 1, nCellsSolve
    if(landmask(iCell) .eq. 1) then

       !adjust the annual deep mean soil temperature and skin temperatures over land: 
       tmn(iCell) = soiltemp(iCell) - 0.0065_RKIND * ter(iCell)
       skintemp(iCell) = skintemp(iCell) - 0.0065_RKIND * (ter(iCell)-soilz(iCell))

       !adjust the soil layer temperatures:
       do ifgSoil = 1, nFGSoilLevels
          st_fg(ifgSoil,iCell) = st_fg(ifgSoil,iCell) - 0.0065_RKIND * (ter(iCell)-soilz(iCell))
       end do

    elseif(landmask(iCell) .eq. 0) then

       tmn(iCell) = skintemp(iCell)

    endif
 enddo

 end subroutine adjust_input_soiltemps

!==================================================================================================
 subroutine init_soil_layers_depth(mesh, fg, dims, configs)
!==================================================================================================

!input arguments:
 type (mpas_pool_type), intent(in) :: mesh
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!inout arguments:
 type (mpas_pool_type), intent(inout) :: fg 

!local variables and arrays:
 integer :: iCell,iSoil
 integer, pointer :: nCellsSolve,nSoilLevels,nFGSoilLevels
 integer, pointer :: config_nsoillevels

 real(kind=RKIND),dimension(:,:),pointer:: dzs_fg,zs_fg
 real(kind=RKIND),dimension(:,:),pointer:: dzs,zs

!--------------------------------------------------------------------------------------------------

! write(0,*)
! write(0,*) '--- enter subroutine init_soil_layers_depth:'

 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(dims, 'nSoilLevels', nSoilLevels)
 call mpas_pool_get_dimension(dims, 'nFGSoilLevels', nFGSoilLevels)

 call mpas_pool_get_array(fg, 'zs_fg', zs_fg)
 call mpas_pool_get_array(fg, 'dzs_fg', dzs_fg)
 call mpas_pool_get_array(fg, 'zs', zs)
 call mpas_pool_get_array(fg, 'dzs', dzs)

 call mpas_pool_get_config(configs, 'config_nsoillevels', config_nsoillevels)

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')

 do iCell = 1, nCellsSolve
    iSoil = 1
    zs_fg(iSoil,iCell) = 0.5_RKIND * dzs_fg(iSoil,iCell)
    do iSoil = 2, nFGSoilLevels
       zs_fg(iSoil,iCell) = zs_fg(iSoil-1,iCell)        &
                          + 0.5_RKIND * dzs_fg(iSoil-1,iCell) &
                          + 0.5_RKIND * dzs_fg(iSoil,iCell)
    enddo
 enddo

 do iCell = 1, nCellsSolve
    dzs(1,iCell) = 0.10_RKIND
    dzs(2,iCell) = 0.30_RKIND
    dzs(3,iCell) = 0.60_RKIND
    dzs(4,iCell) = 1.00_RKIND

    iSoil = 1
    zs(iSoil,iCell)  = 0.5_RKIND * dzs(iSoil,iCell)
    do iSoil = 2, nSoilLevels
       zs(iSoil,iCell) = zs(iSoil-1,iCell)              &
                       + 0.5_RKIND * dzs(iSoil-1,iCell) &
                       + 0.5_RKIND * dzs(iSoil,iCell)
    enddo

 enddo

 end subroutine init_soil_layers_depth

!==================================================================================================
 subroutine init_soil_layers_properties(mesh, fg, dminfo, dims, configs)
!==================================================================================================

!input arguments:
 type (mpas_pool_type), intent(in) :: mesh
 type (dm_info), intent(in) :: dminfo
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!inout arguments:
 type (mpas_pool_type), intent(inout) :: fg 

!local variables:
 integer:: iCell,ifgSoil,iSoil
 integer, pointer:: nCellsSolve,nFGSoilLevels,nSoilLevels
 integer:: num_sm,num_st
 integer,dimension(:),pointer:: landmask
 
 real(kind=RKIND),dimension(:,:),allocatable:: zhave,sm_input,st_input

 real(kind=RKIND),dimension(:),pointer  :: skintemp,tmn
 real(kind=RKIND),dimension(:,:),pointer:: dzs,zs,tslb,smois,sh2o,smcrel
 real(kind=RKIND),dimension(:,:),pointer:: sm_fg,st_fg,zs_fg

 integer, pointer :: config_nsoillevels

!--------------------------------------------------------------------------------------------------

!write(0,*)
! write(0,*) '--- enter subroutine init_soil_layers_properties:'

 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(dims, 'nSoilLevels', nSoilLevels)
 call mpas_pool_get_dimension(dims, 'nFGSoilLevels', nFGSoilLevels)

 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(fg, 'zs_fg', zs_fg)
 call mpas_pool_get_array(fg, 'st_fg', st_fg)
 call mpas_pool_get_array(fg, 'sm_fg', sm_fg)
 call mpas_pool_get_array(fg, 'zs', zs)
 call mpas_pool_get_array(fg, 'dzs', dzs)
 call mpas_pool_get_array(fg, 'sh2o', sh2o)
 call mpas_pool_get_array(fg, 'smcrel', smcrel)
 call mpas_pool_get_array(fg, 'smois', smois)
 call mpas_pool_get_array(fg, 'tslb', tslb)
 call mpas_pool_get_array(fg, 'skintemp', skintemp)
 call mpas_pool_get_array(fg, 'tmn', tmn)

 call mpas_pool_get_config(configs, 'config_nsoillevels', config_nsoillevels)

 write(0,*) 'nSoilLevels   =',nSoilLevels
 write(0,*) 'nFGSoilLevels =',nFGSoilLevels


!check that interpolation of the meteorological data to the MPAS grid did not create negative
!values for the first-guess soil temperatures and soil moistures.
 num_sm = 0
 num_st = 0
 do iCell = 1, nCellsSolve
    do ifgSoil = 1, nFGSoilLevels
       if(st_fg(ifgSoil,iCell) .le. 0._RKIND) num_st = num_st + 1
       if(sm_fg(ifgSoil,iCell) .lt. 0._RKIND) num_sm = num_sm + 1
    enddo
 enddo
 if(num_st .gt. 0) then
    write(0,*) 'Error in interpolation of st_fg to MPAS grid: num_st =', num_st
    call mpas_dmpar_abort(dminfo)
 elseif(num_sm .gt. 0) then
    write(0,*) 'Error in interpolation of sm_fg to MPAS grid: num_sm =', num_sm
    call mpas_dmpar_abort(dminfo)
 endif 

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')

 if(.not.allocated(zhave)   ) allocate(zhave(nFGSoilLevels+2,nCellsSolve)   )
 if(.not.allocated(st_input)) allocate(st_input(nFGSoilLevels+2,nCellsSolve))
 if(.not.allocated(sm_input)) allocate(sm_input(nFGSoilLevels+2,nCellsSolve))

 do iCell = 1, nCellsSolve

    ifgSoil = 1
    zhave(ifgSoil,iCell)    = 0._RKIND
    st_input(ifgSoil,iCell) = skintemp(iCell)
    sm_input(ifgSoil,iCell) = sm_fg(ifgSoil+1,iCell)

    do ifgSoil = 1, nFGSoilLevels
       zhave(ifgSoil+1,iCell) = zs_fg(ifgSoil,iCell) / 100._RKIND
       st_input(ifgSoil+1,iCell) = st_fg(ifgSoil,iCell)
       sm_input(ifgSoil+1,iCell) = sm_fg(ifgSoil,iCell)
    enddo

    zhave(nFGSoilLevels+2,iCell) = 300._RKIND/100._RKIND
    st_input(nFGSoilLevels+2,iCell) = tmn(iCell)
    sm_input(nFGSoilLevels+2,iCell) = sm_input(nFGSoilLevels,iCell)

    if(iCell .eq. 1) then
       do ifgSoil = 1,nFGSoilLevels+2
          write(0,101) ifgSoil,zhave(ifgSoil,iCell)
       enddo
    endif

 enddo

!... interpolate the soil temperature, soil moisture, and soil liquid temperature to the four
!    layers used in the NOAH land surface scheme:

 do iCell = 1, nCellsSolve

    if(landmask(iCell) .eq. 1) then

       noah: do iSoil = 1 , nSoilLevels
          input: do ifgSoil = 1 , nFGSoilLevels+2-1
             if(iCell .eq. 1) write(0,102) iSoil,ifgSoil,zs(iSoil,iCell), &
                zhave(ifgSoil,iCell),zhave(ifgSoil+1,iCell)

             if(zs(iSoil,iCell).ge.zhave(ifgSoil,iCell) .and. &
                zs(iSoil,iCell).le.zhave(ifgSoil+1,iCell)) then

                tslb(iSoil,iCell) = &
                      (st_input(ifgSoil,iCell) * (zhave(ifgSoil+1,iCell)-zs(iSoil,iCell))    &
                    +  st_input(ifgSoil+1,iCell) * (zs(iSoil,iCell)-zhave(ifgSoil,iCell)))   &
                         / (zhave(ifgSoil+1,iCell)-zhave(ifgSoil,iCell))
                if(iCell .eq. 1) write(6,102) iSoil,ifgSoil,zs(iSoil,iCell), &
                   zhave(ifgSoil,iCell),zhave(ifgSoil+1,iCell)
                         
                smois(iSoil,iCell) = &
                       (sm_input(ifgSoil,iCell) * (zhave(ifgSoil+1,iCell)-zs(iSoil,iCell))   &
                    +  sm_input(ifgSoil+1,iCell) * (zs(iSoil,iCell)-zhave(ifgSoil,iCell)))   &
                    / (zhave(ifgSoil+1,iCell)-zhave(ifgSoil,iCell))

                sh2o(iSoil,iCell)   = 0._RKIND
                smcrel(iSoil,iCell) = 0._RKIND

                exit input
             endif
          enddo input
          if(iCell.eq. 1) write(0,*)
       enddo noah

    elseif(landmask(iCell) .eq. 0) then

       !fill the soil temperatures with the skin temperatures over oceans:
       do iSoil = 1, nSoilLevels
          tslb(iSoil,iCell)    = skintemp(iCell)
          smois(iSoil,iCell)   = 1._RKIND
          sh2o(iSoil,iCell)    = 1._RKIND
          smcrel(iSoil,iCell)  = 0._RKIND
       enddo

    endif

 enddo

!... final checks:

 do iCell = 1, nCellsSolve

    if(landmask(iCell).eq.1 .and. tslb(1,iCell).gt.170._RKIND .and. tslb(1,iCell).lt.400._RKIND &
       .and. smois(1,iCell).lt.0.005_RKIND) then
       do iSoil = 1, nSoilLevels
          smois(iSoil,iCell) = 0.005_RKIND
       enddo
    endif

 enddo

!formats:
 101 format(i4,4(1x,e15.8))
 102 format(2i5,5(1x,e15.8))

 if(allocated(zhave)   ) deallocate(zhave )
 if(allocated(st_input)) deallocate(st_input)
 if(allocated(sm_input)) deallocate(sm_input)

 end subroutine init_soil_layers_properties

!==================================================================================================
 subroutine physics_init_sst(mesh, input, dims, configs)
!==================================================================================================

!input arguments: 
 type (mpas_pool_type), intent(in) :: mesh
 type (mpas_pool_type), intent(inout) :: input
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!local variables:
 character(len=StrKIND):: mess

 integer, pointer:: nCellsSolve
 integer:: iCell
 integer:: num_seaice_changes
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND):: xice_threshold
 real(kind=RKIND),dimension(:),pointer  :: seaice,sst,tsk,xice

 logical, pointer :: config_frac_seaice

!--------------------------------------------------------------------------------------------------
! write(0,*)
! write(0,*) '--- enter subroutine physics_init_sst:'

 call mpas_pool_get_config(configs, 'config_frac_seaice', config_frac_seaice)

 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)

 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(input, 'sst', sst)
 call mpas_pool_get_array(input, 'seaice', seaice)
 call mpas_pool_get_array(input, 'skintemp', tsk)
 call mpas_pool_get_array(input, 'xice', xice)


 if(.not. config_frac_seaice) then
    xice_threshold = 0.5_RKIND
 elseif(config_frac_seaice) then
    xice_threshold = 0.02
 endif
 write(0,*) '--- config_frac_seaice      :', config_frac_seaice
 write(0,*) '--- xice_threshold          :', xice_threshold

 do iCell = 1, nCellsSolve
    seaice(iCell) = 0._RKIND

    !... initialize skin temperature with sea-surface temperature over ocean cells:
    if(landmask(iCell) == 0 .and. xice(iCell) < xice_threshold) tsk(iCell) = sst(iCell)

    !... make sure that cells with sea-ice fraction greater than 0 are defined as ocean cells. If
    !not, the sea-ice fraction is reset to zero:
    num_seaice_changes = 0
    if((landmask(iCell) == 1 .and. xice(iCell) > 0._RKIND) .or. xice(iCell) > 200._RKIND) then
       num_seaice_changes = num_seaice_changes + 1
       xice(iCell) = 0._RKIND
    endif
    if(xice(iCell) .ge. xice_threshold) seaice(iCell) = 1._RKIND
 enddo

 write(mess,fmt='(A,i12)') '    number of seaice cells converted to land cells 1 =', &
       num_seaice_changes
 call physics_message(mess)

! write(0,*) '--- end subroutine physics_init_sst:'

 end subroutine physics_init_sst

!==================================================================================================
 subroutine physics_init_seaice(mesh, input, dims, configs)
!==================================================================================================

!input arguments:
 type (mpas_pool_type), intent(in) :: mesh
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!inout arguments: this subroutine is called from the MPAS model side.
 type (mpas_pool_type), intent(inout) :: input

!local variables:
 character(len=StrKIND):: mess
 integer, pointer:: nCellsSolve,nSoilLevels
 integer:: iCell,iSoil
 integer:: num_seaice_changes
 integer,dimension(:),pointer:: landmask,isltyp,ivgtyp

 real(kind=RKIND):: xice_threshold
 real(kind=RKIND):: mid_point_depth
 real(kind=RKIND),dimension(:),pointer  :: vegfra
 real(kind=RKIND),dimension(:),pointer  :: seaice,xice
 real(kind=RKIND),dimension(:),pointer  :: skintemp,tmn,xland
 real(kind=RKIND),dimension(:,:),pointer:: tslb,smois,sh2o,smcrel

 logical, pointer :: config_frac_seaice
 character(len=StrKIND),pointer:: config_landuse_data
 integer:: isice_lu

!note that this threshold is also defined in module_physics_vars.F.It is defined here to avoid
!adding "use module_physics_vars" since this subroutine is only used for the initialization of
!a "real" forecast with $CORE = init_nhyd_atmos.
 real(kind=RKIND),parameter:: xice_tsk_threshold = 271.
 real(kind=RKIND),parameter:: total_depth        = 3.   ! 3-meter soil depth.

!--------------------------------------------------------------------------------------------------
! write(0,*)
! write(0,*) '--- enter physics_init_seaice:'

 call mpas_pool_get_config(configs, 'config_frac_seaice', config_frac_seaice)
 call mpas_pool_get_config(configs, 'config_landuse_data', config_landuse_data)

 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(dims, 'nSoilLevels', nSoilLevels)

 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(mesh, 'lu_index', ivgtyp)
 call mpas_pool_get_array(mesh, 'soilcat_top', isltyp)

 call mpas_pool_get_array(input, 'seaice', seaice)
 call mpas_pool_get_array(input, 'xice', xice)
 call mpas_pool_get_array(input, 'vegfra', vegfra)

 call mpas_pool_get_array(input, 'skintemp', skintemp)
 call mpas_pool_get_array(input, 'tmn', tmn)
 call mpas_pool_get_array(input, 'xland', xland)

 call mpas_pool_get_array(input, 'tslb', tslb)
 call mpas_pool_get_array(input, 'smois', smois)
 call mpas_pool_get_array(input, 'sh2o', sh2o)
 call mpas_pool_get_array(input, 'smcrel', smcrel)

!define the land use category for sea-ice as a function of the land use category input file:
 sfc_input_select1: select case(trim(config_landuse_data))
    case('OLD')
       isice_lu = 11
    case('USGS')
       isice_lu = 24
    case('MODIFIED_IGBP_MODIS_NOAH')
       isice_lu = 15
    case default
       CALL physics_error_fatal ('Invalid Land Use Dataset '//trim(config_landuse_data))
 end select sfc_input_select1
 write(0,*) '--- isice_lu   = ',isice_lu

!assign the threshold value for xice as a function of config_frac_seaice:
 if(.not. config_frac_seaice) then
    xice_threshold = 0.5_RKIND
    do iCell = 1,nCellsSolve
       if(xice(iCell) >= xice_threshold) then
          xice(iCell) = 1._RKIND
       else
          xice(iCell) = 0._RKIND
       endif
    enddo
 elseif(config_frac_seaice) then
    xice_threshold = 0.02
 endif
 write(0,*) '--- config_frac_seaice      :', config_frac_seaice
 write(0,*) '--- xice_threshold          :', xice_threshold

!convert seaice points to land points when the sea-ice fraction is greater than the
!prescribed threshold:
 num_seaice_changes = 0
 do iCell = 1, nCellsSolve
    if(xice(iCell) .ge. xice_threshold .or. &
      (landmask(iCell).eq.0 .and. skintemp(iCell).lt.xice_tsk_threshold)) then

       num_seaice_changes = num_seaice_changes + 1
       !... sea-ice points are converted to land points:
       if(landmask(iCell) .eq. 0) tmn(iCell) = 271.4_RKIND
       ivgtyp(iCell) = isice_lu
       isltyp(iCell) = 16
       vegfra(iCell) = 0._RKIND
       xland(iCell)  = 1._RKIND

       !... recalculate the soil temperature and soil moisture:
       do iSoil = 1, nSoilLevels
          mid_point_depth = total_depth/nSoilLevels/2. &
                          + (iSoil-1)*(total_depth/nSoilLevels)
          tslb(iSoil,iCell) = ((total_depth-mid_point_depth) * skintemp(iCell) &
                             +  mid_point_depth * tmn(iCell)) / total_depth
          smois(iSoil,iCell)  = 1._RKIND
          sh2o(iSoil,iCell)   = 0._RKIND
          smcrel(iSoil,iCell) = 0._RKIND
       enddo
       
    elseif(xice(iCell) .lt. xice_threshold) then
       if(xice(iCell) .gt. 0._RKIND) &
          write(0,101) iCell,DBLE(landmask(iCell)),xland(iCell),xice(iCell)
       xice(iCell) = 0._RKIND
    endif

 enddo
 write(0,*)
 write(mess,fmt='(A,i12)') 'number of seaice cells converted to land cells 2 =', &
       num_seaice_changes
 call physics_message(mess)

!finally, update the sea-ice flag:
 do iCell = 1, nCellsSolve
    seaice(iCell) = 0._RKIND
    if(xice(iCell) > 0._RKIND) then
       seaice(iCell) = 1._RKIND
    endif
 enddo
 101 format(i9,5(1x,e15.8))

! write(0,*) '--- end physics_init_seaice:'

 end subroutine physics_init_seaice

!==================================================================================================
 end module mpas_atmphys_initialize_real
!==================================================================================================

