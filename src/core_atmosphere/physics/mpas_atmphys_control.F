! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=============================================================================================
 module mpas_atmphys_control
 use mpas_configure
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
!>      Laura D. Fowler (birch.ucar.edu) / 2013-05-29.

 contains

!==================================================================================================
 subroutine physics_namelist_check
!==================================================================================================

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
 subroutine physics_registry_init(config_do_restart,mesh,sfc_input)
!==================================================================================================

!input and output arguments:
 logical,intent(in):: config_do_restart
 type(mesh_type),intent(in):: mesh
 type(sfc_input_type),intent(inout):: sfc_input

!local variables:
 integer:: iCell
 
!--------------------------------------------------------------------------------------------------

!initialization of input variables, if needed:

 if(.not. config_do_restart) then

    lsm_select: select case(trim(config_lsm_scheme))

       case("noah")
       !initialize the thickness of the soil layers for the Noah scheme:
          do iCell = 1, mesh % nCells
             if(sfc_input % landmask % array(iCell) == 1) then  
                sfc_input % dzs % array(1,iCell) = 0.10
                sfc_input % dzs % array(2,iCell) = 0.30
                sfc_input % dzs % array(3,iCell) = 0.60
                sfc_input % dzs % array(4,iCell) = 1.00
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
 type(mesh_type),intent(in):: mesh
 type(sfc_input_type),intent(inout):: sfc_input

!local variables:
 integer:: iCell,iMonth,iSoil
 
!--------------------------------------------------------------------------------------------------

!initialization of surface input variables that are not needed if we run the current set of
!idealized test cases:

 do iCell = 1, mesh % nCells
    !terrain,soil type, and vegetation:
    sfc_input % ter      % array(iCell) = 0.
    sfc_input % xice     % array(iCell) = 0.
    sfc_input % landmask % array(iCell) = 0
    sfc_input % ivgtyp   % array(iCell) = 0
    sfc_input % isltyp   % array(iCell) = 0
    sfc_input % shdmin   % array(iCell) = 0.
    sfc_input % shdmax   % array(iCell) = 0.
    sfc_input % vegfra   % array(iCell) = 0.

    !snow coverage:    
    sfc_input % snow     % array(iCell) = 0.
    sfc_input % snowc    % array(iCell) = 0.
    sfc_input % snoalb   % array(iCell) = 0.08 
    
    !surface and sea-surface temperatures:
    sfc_input % skintemp % array(iCell) = 288.0
    sfc_input % sst      % array(iCell) = 288.0
    
    !soil layers:
    sfc_input % tmn % array(iCell) = 288.0
    do iSoil = 1, mesh % nSoilLevels
       sfc_input % tslb % array(iSoil,iCell)   = 288.0
       sfc_input % smcrel % array(iSoil,iCell) =   0.0
       sfc_input % sh2o   % array(iSoil,iCell) =   0.0
       sfc_input % smois  % array(iSoil,iCell) =   0.0
       sfc_input % dzs    % array(iSoil,iCell) =   0.0
    enddo

    !monthly climatological surface albedo and greeness fraction:
    do iMonth = 1, mesh % nMonths
       sfc_input % albedo12m % array(iMonth,iCell) = 0.08
       sfc_input % greenfrac % array(iMonth,iCell) = 0.
    enddo

 enddo   

 end subroutine physics_idealized_init

!==================================================================================================
 end module mpas_atmphys_control
!==================================================================================================


