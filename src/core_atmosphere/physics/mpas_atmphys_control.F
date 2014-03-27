! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_control
 use mpas_kind_types
 use mpas_grid_types

 use mpas_atmphys_utilities

 implicit none
 private
 public:: physics_namelist_check, &
          physics_idealized_init, &
          physics_registry_init

 logical,public:: moist_physics


!>\brief MPAS control and initialization routines.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!> mpas_atmphys_control contains subroutines called during the initialization.
!>
!> subroutines called in mpas_atmphys_control:
!> -------------------------------------------
!> physics_namelist_check: checks that physics namelist parameters are defined correctly.
!> physics_registry_init : initializes thicknesses of soil layers for NOAH land-surface scheme.
!> physics_idealized_init: initializes physics variables needed to run idealized cases.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed the namelist option config_eddy_scheme and associated sourcecode.
!>    * removed the namelist option config_conv_shallow_scheme and associated sourcecode.
!>    * removed controls to the updated Kain-Fritsch convection scheme.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-05-29.
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-22.
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-05-15.


 contains


!==================================================================================================
 subroutine physics_namelist_check(mesh,configs)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: config_microp_scheme,    &
                                  config_conv_deep_scheme, &
                                  config_lsm_scheme,       &
                                  config_pbl_scheme,       &
                                  config_gwdo_scheme,      &
                                  config_radt_cld_scheme,  &
                                  config_radt_lw_scheme,   &
                                  config_radt_sw_scheme,   &
                                  config_sfclayer_scheme

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme'   ,config_microp_scheme   )
 call mpas_pool_get_config(configs,'config_conv_deep_scheme',config_conv_deep_scheme)
 call mpas_pool_get_config(configs,'config_lsm_scheme'      ,config_lsm_scheme      )
 call mpas_pool_get_config(configs,'config_pbl_scheme'      ,config_pbl_scheme      )
 call mpas_pool_get_config(configs,'config_gwdo_scheme'     ,config_gwdo_scheme     )
 call mpas_pool_get_config(configs,'config_radt_cld_scheme' ,config_radt_cld_scheme )
 call mpas_pool_get_config(configs,'config_radt_lw_scheme'  ,config_radt_lw_scheme  )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme'  ,config_radt_sw_scheme  )
 call mpas_pool_get_config(configs,'config_sfclayer_scheme' ,config_sfclayer_scheme )

 write(0,*)
 write(0,*) '--- enter subroutine physics_namelist_check:'
 write(0,*) '    config_microp_scheme       = ', trim(config_microp_scheme)
 write(0,*) '    config_conv_deep_scheme    = ', trim(config_conv_deep_scheme)
 write(0,*) '    config_lsm_scheme          = ', trim(config_lsm_scheme)
 write(0,*) '    config_pbl_scheme          = ', trim(config_pbl_scheme)
 write(0,*) '    config_gwdo_scheme         = ', trim(config_gwdo_scheme)
 write(0,*) '    config_radt_cld_scheme     = ', trim(config_radt_cld_scheme)
 write(0,*) '    config_radt_lw_scheme      = ', trim(config_radt_lw_scheme)
 write(0,*) '    config_radt_sw_scheme      = ', trim(config_radt_sw_scheme)
 write(0,*) '    config_sfclayer_scheme     = ', trim(config_sfclayer_scheme)

!cloud microphysics scheme:
 if(.not. (config_microp_scheme .eq. 'off'      .or. &
           config_microp_scheme .eq. 'kessler'  .or. &
           config_microp_scheme .eq. 'wsm6'     )) then
          
    write(mpas_err_message,'(A,A10)') 'illegal value for config_microp_scheme:', &
          trim(config_microp_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!deep convection scheme:
 if(.not. (config_conv_deep_scheme .eq. 'off'                  .or. &
           config_conv_deep_scheme .eq. 'kain_fritsch'         .or. &
           config_conv_deep_scheme .eq. 'tiedtke'              )) then

    write(mpas_err_message,'(A,A10)') 'illegal value for config_deep_conv_scheme: ', &
          trim(config_conv_deep_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!pbl scheme:
 if(.not. (config_pbl_scheme .eq. 'off' .or. &
           config_pbl_scheme .eq. 'ysu')) then

    write(mpas_err_message,'(A,A10)') 'illegal value for pbl_scheme: ', &
          trim(config_pbl_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!gravity wave drag over orography scheme:
 if(.not. (config_gwdo_scheme .eq. 'off' .or. &
           config_gwdo_scheme .eq. 'ysu_gwdo')) then

    write(mpas_err_message,'(A,A10)') 'illegal value for gwdo_scheme: ', &
          trim(config_gwdo_scheme)
    call physics_error_fatal(mpas_err_message)

 elseif(config_gwdo_scheme .eq. 'ysu_gwdo' .and. config_pbl_scheme .ne. 'ysu') then

    write(mpas_err_message,'(A,A10)') 'turn YSU PBL scheme on with config_gwdo = ysu_gwdo:', &
          trim(config_gwdo_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!lw radiation scheme:
 if(.not. (config_radt_lw_scheme .eq. 'off'    .or. &
           config_radt_lw_scheme .eq. 'cam_lw' .or. &
           config_radt_lw_scheme .eq. 'rrtmg_lw')) then
 
    write(mpas_err_message,'(A,A10)') 'illegal value for longwave radiation scheme: ', &
          trim(config_radt_lw_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!sw radiation scheme:
 if(.not. (config_radt_sw_scheme .eq. 'off'    .or. &
           config_radt_sw_scheme .eq. 'cam_sw' .or. &
           config_radt_sw_scheme .eq. 'rrtmg_sw')) then
 
    write(mpas_err_message,'(A,A10)') 'illegal value for shortwave radiation _scheme: ', &
          trim(config_radt_sw_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!cloud fraction for radiation schemes:
 if(.not. (config_radt_cld_scheme .eq. 'off'           .or. &
           config_radt_cld_scheme .eq. 'cld_incidence' .or. &
           config_radt_cld_scheme .eq. 'cld_fraction')) then

    write(mpas_err_message,'(A,A10)') 'illegal value for calculation of cloud fraction: ', &
          trim(config_radt_cld_scheme)
    call physics_error_fatal(mpas_err_message)

 endif
 if((config_radt_lw_scheme.ne.'off' .and. config_radt_cld_scheme.eq.'off') .or. &
    (config_radt_sw_scheme.ne.'off' .and. config_radt_cld_scheme.eq.'off')) then

    write(0,*)
    write(mpas_err_message,'(A,A10)') &
       '    config_radt_cld_scheme is not set for radiation calculation'
    call physics_message(mpas_err_message)
    write(mpas_err_message,'(A,A10)') &
       '    switch calculation of cloud fraction to config_radt_cld_scheme = cld_incidence'
    call physics_message(mpas_err_message)
    config_radt_cld_scheme = "cld_incidence"

 endif

!surface-layer scheme:
 if(.not. (config_sfclayer_scheme .eq. 'off'  .or. &
           config_sfclayer_scheme .eq. 'monin_obukhov')) then
 
    write(mpas_err_message,'(A,A10)') 'illegal value for surface layer scheme: ', &
          trim(config_sfclayer_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!land-surface scheme: note that config_sfclayer_scheme must be defined for the land-surface
!scheme to be called:
 if(config_lsm_scheme .ne. 'off' .and. config_sfclayer_scheme .eq. 'off') then
 
    call physics_error_fatal('land surface scheme: ' // &
                             'set config_sfclayer_scheme different than off')
    
 elseif(.not. (config_lsm_scheme .eq. 'off ' .or. &
               config_lsm_scheme .eq. 'noah')) then
 
    write(mpas_err_message,'(A,A10)') 'illegal value for land surface scheme: ', &
          trim(config_lsm_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!checks if any physics process is called. if not, return:
 moist_physics = .true.
 
 if(config_microp_scheme       .eq. 'off' .and. &
    config_conv_deep_scheme    .eq. 'off' .and. &
    config_lsm_scheme          .eq. 'off' .and. & 
    config_pbl_scheme          .eq. 'off' .and. & 
    config_radt_lw_scheme      .eq. 'off' .and. &
    config_radt_sw_scheme      .eq. 'off' .and. &
    config_sfclayer_scheme     .eq. 'off') moist_physics = .false.

 write(0,*) '--- end subroutine physics_namelist_check:'

 end subroutine physics_namelist_check

!==================================================================================================
 subroutine physics_registry_init(mesh,configs,sfc_input)
!==================================================================================================

!input and inout arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 logical,pointer:: config_do_restart
 character(len=StrKIND),pointer:: config_lsm_scheme
 integer,pointer:: nCells
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:,:),pointer:: dzs

!local variables:
 integer:: iCell

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_do_restart',config_do_restart)
 call mpas_pool_get_config(configs,'config_lsm_scheme',config_lsm_scheme)

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(sfc_input,'landmask',landmask)
 call mpas_pool_get_array(sfc_input,'dzs'     , dzs    )

!initialization of input variables, if needed:

 if(.not. config_do_restart) then

    lsm_select: select case(trim(config_lsm_scheme))

       case("noah")
       !initialize the thickness of the soil layers for the Noah scheme:
          do iCell = 1, nCells
             if(landmask(iCell) == 1) then  
                dzs(1,iCell) = 0.10_RKIND
                dzs(2,iCell) = 0.30_RKIND
                dzs(3,iCell) = 0.60_RKIND
                dzs(4,iCell) = 1.00_RKIND
             endif
          enddo

       case default
    
    end select lsm_select
    
 endif

 end subroutine physics_registry_init

!==================================================================================================
 subroutine physics_idealized_init(mesh,sfc_input)
!==================================================================================================

!input and output arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers: 
 integer,pointer:: nCells,nMonths,nSoilLevels
 integer,dimension(:),pointer:: isltyp,ivgtyp,landmask
 
 real(kind=RKIND),dimension(:),pointer:: shdmin,shdmax,snoalb,snow,snowc,ter,vegfra,xice
 real(kind=RKIND),dimension(:),pointer:: skintemp,sst,tmn
 real(kind=RKIND),dimension(:,:),pointer:: tslb,smcrel,sh2o,smois,dzs
 real(kind=RKIND),dimension(:,:),pointer:: albedo12m,greenfrac

!local variables:
 integer:: iCell,iMonth,iSoil
 
!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nMonths',nMonths)
 call mpas_pool_get_dimension(mesh,'nSoilLevel',nSoilLevels)

 call mpas_pool_get_array(sfc_input,'isltyp'   ,isltyp   )
 call mpas_pool_get_array(sfc_input,'ivgtyp'   ,ivgtyp   )
 call mpas_pool_get_array(sfc_input,'landmask' ,landmask )
 call mpas_pool_get_array(sfc_input,'shdmin'   ,shdmin   )
 call mpas_pool_get_array(sfc_input,'shdmax'   ,shdmax   )
 call mpas_pool_get_array(sfc_input,'snoalb'   ,snoalb   )
 call mpas_pool_get_array(sfc_input,'snow'     ,snow     )
 call mpas_pool_get_array(sfc_input,'snowc'    ,snowc    )
 call mpas_pool_get_array(sfc_input,'ter'      ,ter      )
 call mpas_pool_get_array(sfc_input,'vegfra'   ,vegfra   )
 call mpas_pool_get_array(sfc_input,'xice'     ,xice     )
 call mpas_pool_get_array(sfc_input,'skintemp' ,skintemp )
 call mpas_pool_get_array(sfc_input,'sst'      ,sst      )
 call mpas_pool_get_array(sfc_input,'tmn'      ,tmn      )
 call mpas_pool_get_array(sfc_input,'dzs'      ,dzs      )
 call mpas_pool_get_array(sfc_input,'sh2o'     ,sh2o     )
 call mpas_pool_get_array(sfc_input,'smcrel'   ,smcrel   )
 call mpas_pool_get_array(sfc_input,'smois'    ,smois    )
 call mpas_pool_get_array(sfc_input,'tslb'     ,tslb     )
 call mpas_pool_get_array(sfc_input,'albedo12m',albedo12m)
 call mpas_pool_get_array(sfc_input,'greenfrac',greenfrac)

!initialization of surface input variables that are not needed if we run the current set of
!idealized test cases:

 do iCell = 1, nCells
    !terrain,soil type, and vegetation:
    landmask(iCell) = 0
    ivgtyp(iCell)   = 0
    isltyp(iCell)   = 0
    ter(iCell)      = 0._RKIND
    xice(iCell)     = 0._RKIND
    shdmin(iCell)   = 0._RKIND
    shdmax(iCell)   = 0._RKIND
    vegfra(iCell)   = 0._RKIND

    !snow coverage:    
    snow(iCell)   = 0._RKIND
    snowc(iCell)  = 0._RKIND
    snoalb(iCell) = 0.08_RKIND
    
    !surface and sea-surface temperatures:
    skintemp(iCell) = 288._RKIND
    sst(iCell)      = 288._RKIND
    
    !soil layers:
    tmn(iCell) = 288._RKIND
    do iSoil = 1, nSoilLevels
       dzs(iSoil,iCell)    =   0._RKIND
       sh2o(iSoil,iCell)   =   0._RKIND
       smcrel(iSoil,iCell) =   0._RKIND
       smois(iSoil,iCell)  =   0._RKIND
       tslb(iSoil,iCell)   = 288._RKIND
    enddo

    !monthly climatological surface albedo and greeness fraction:
    do iMonth = 1, nMonths
       albedo12m(iMonth,iCell) = 0.08_RKIND
       greenfrac(iMonth,iCell) = 0._RKIND
    enddo

 enddo   

 end subroutine physics_idealized_init

!==================================================================================================
 end module mpas_atmphys_control
!==================================================================================================


