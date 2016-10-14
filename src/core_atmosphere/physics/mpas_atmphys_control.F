! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_control
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_utilities

 implicit none
 private
 public:: physics_namelist_check, &
          physics_registry_init

 logical,public:: moist_physics


!MPAS control and initialization routines.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines called in mpas_atmphys_control:
! -------------------------------------------
! * physics_namelist_check: checks that physics namelist parameters are defined correctly.
! * physics_registry_init : initializes thicknesses of soil layers for NOAH land-surface scheme.
! * physics_idealized_init: initializes physics variables needed to run idealized cases.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the namelist option config_eddy_scheme and associated sourcecode.
! * removed the namelist option config_conv_shallow_scheme and associated sourcecode.
! * removed controls to the updated Kain-Fritsch convection scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * removed subroutine physics_idealized_init, also available in mpas_init_atm_cases.F in 
!   core_init_atmosphere.
!   Laura D. Fowler (laura@ucar.edu) / 2014-08-11.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * renamed "wsm6" to "mp_wsm6" and "kessler" to "mp_kessler".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-09.
! * renamed "kain_fritsch" to "cu_kain_fritsch" and "tiedtke" to "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * renamed "ysu" to "bl_ysu", "ysu_gwdo" to "bl_gwdo_ysu", and "monin_obukhov" to "sf_monin_obukhov".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25. 
! * added the options mp_thompson_aerosols and mp_thompson.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25.
! * added the options cu_grell_freitas and cu_grell_freitas_momentum.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-31.
! * added the options sf_mynn and bl_mynn and for the MYNN parameterization from WRF version 3.6.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-11.
! * added the options sf_mynn_wrf38 and bl_mynn_wrf38 for the MYNN parameterization from WRF version 3.8.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-14.
! * added the options cu_tiedtke_wrf381 for the Tiedtke parameterization of convection from WRF version 3.8.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-08-18.
! * added the option cu_ntiedtke_wrf381 for the "new" Tiedtke parameterization of convection from WRF
!   version 3.8.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-19.
! * added the option mp_wsm6_wrf381 for the WSM6 parameterization of cloud microphysics from WRF version 3.8.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-16.


 contains


!=================================================================================================================
 subroutine physics_namelist_check(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: config_physics_suite,     &
                                  config_microp_scheme,     &
                                  config_convection_scheme, &
                                  config_lsm_scheme,        &
                                  config_pbl_scheme,        &
                                  config_gwdo_scheme,       &
                                  config_radt_cld_scheme,   &
                                  config_radt_lw_scheme,    &
                                  config_radt_sw_scheme,    &
                                  config_sfclayer_scheme

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_namelist_check:'

 call mpas_pool_get_config(configs,'config_physics_suite'    ,config_physics_suite    )
 call mpas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )
 call mpas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
 call mpas_pool_get_config(configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call mpas_pool_get_config(configs,'config_pbl_scheme'       ,config_pbl_scheme       )
 call mpas_pool_get_config(configs,'config_gwdo_scheme'      ,config_gwdo_scheme      )
 call mpas_pool_get_config(configs,'config_radt_cld_scheme'  ,config_radt_cld_scheme  )
 call mpas_pool_get_config(configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )
 call mpas_pool_get_config(configs,'config_sfclayer_scheme'  ,config_sfclayer_scheme  )

 !
 ! Setup schemes according to the selected physics suite
 !
 if (trim(config_physics_suite) == 'mesoscale_reference') then

    if (trim(config_microp_scheme) == 'suite')     config_microp_scheme     = 'mp_wsm6'
    if (trim(config_convection_scheme) == 'suite') config_convection_scheme = 'cu_tiedtke'
    if (trim(config_pbl_scheme) == 'suite')        config_pbl_scheme        = 'bl_ysu'
    if (trim(config_gwdo_scheme) == 'suite')       config_gwdo_scheme       = 'off'
    if (trim(config_radt_lw_scheme) == 'suite')    config_radt_lw_scheme    = 'rrtmg_lw'
    if (trim(config_radt_sw_scheme) == 'suite')    config_radt_sw_scheme    = 'rrtmg_sw'
    if (trim(config_radt_cld_scheme) == 'suite')   config_radt_cld_scheme   = 'cld_fraction'
    if (trim(config_sfclayer_scheme) == 'suite')   config_sfclayer_scheme   = 'sf_monin_obukhov'
    if (trim(config_lsm_scheme) == 'suite')        config_lsm_scheme        = 'noah'

 else if (trim(config_physics_suite) == 'none') then

    if (trim(config_microp_scheme) == 'suite')     config_microp_scheme     = 'off'
    if (trim(config_convection_scheme) == 'suite') config_convection_scheme = 'off'
    if (trim(config_pbl_scheme) == 'suite')        config_pbl_scheme        = 'off'
    if (trim(config_gwdo_scheme) == 'suite')       config_gwdo_scheme       = 'off'
    if (trim(config_radt_lw_scheme) == 'suite')    config_radt_lw_scheme    = 'off'
    if (trim(config_radt_sw_scheme) == 'suite')    config_radt_sw_scheme    = 'off'
    if (trim(config_radt_cld_scheme) == 'suite')   config_radt_cld_scheme   = 'off'
    if (trim(config_sfclayer_scheme) == 'suite')   config_sfclayer_scheme   = 'off'
    if (trim(config_lsm_scheme) == 'suite')        config_lsm_scheme        = 'off'

 else

    write(mpas_err_message,'(A)') 'Unrecognized choice of physics suite: config_physics_suite = '''// &
                                  trim(config_physics_suite)//''''
    call physics_error_fatal(mpas_err_message)

 end if

!cloud microphysics scheme:
 if(.not. (config_microp_scheme .eq. 'off'                  .or. &
           config_microp_scheme .eq. 'mp_kessler'           .or. &
           config_microp_scheme .eq. 'mp_thompson_aerosols' .or. &
           config_microp_scheme .eq. 'mp_thompson'          .or. &
           config_microp_scheme .eq. 'mp_wsm6_wrf381')) then

    write(mpas_err_message,'(A,A10)') 'illegal value for config_microp_scheme:', &
          trim(config_microp_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!convection scheme:
 if(.not. (config_convection_scheme .eq. 'off'                       .or. &
           config_convection_scheme .eq. 'cu_grell_freitas'          .or. &
           config_convection_scheme .eq. 'cu_grell_freitas_momentum' .or. &
           config_convection_scheme .eq. 'cu_kain_fritsch'           .or. &
           config_convection_scheme .eq. 'cu_tiedtke'                .or. &
           config_convection_scheme .eq. 'cu_tiedtke_wrf381'         .or. &
           config_convection_scheme .eq. 'cu_ntiedtke_wrf381'        )) then

    write(mpas_err_message,'(A,A10)') 'illegal value for config_convection_scheme: ', &
          trim(config_convection_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!pbl scheme:
 if(.not. (config_pbl_scheme .eq. 'off'           .or. &
           config_pbl_scheme .eq. 'bl_mynn'       .or. &
           config_pbl_scheme .eq. 'bl_mynn_wrf38' .or. &
           config_pbl_scheme .eq. 'bl_ysu'        )) then

    write(mpas_err_message,'(A,A10)') 'illegal value for pbl_scheme: ', &
          trim(config_pbl_scheme)
    call physics_error_fatal(mpas_err_message)

 endif

!gravity wave drag over orography scheme:
 if(.not. (config_gwdo_scheme .eq. 'off' .or. &
           config_gwdo_scheme .eq. 'bl_ysu_gwdo')) then

    write(mpas_err_message,'(A,A10)') 'illegal value for gwdo_scheme: ', &
          trim(config_gwdo_scheme)
    call physics_error_fatal(mpas_err_message)

 elseif(config_gwdo_scheme .eq. 'bl_ysu_gwdo' .and. config_pbl_scheme .ne. 'bl_ysu') then

    write(mpas_err_message,'(A,A10)') 'turn YSU PBL scheme on with config_gwdo = bl_ysu_gwdo:', &
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
 if(.not. (config_radt_cld_scheme .eq. 'off'                  .or. &
           config_radt_cld_scheme .eq. 'cld_incidence'        .or. &
           config_radt_cld_scheme .eq. 'cld_fraction'         .or. &
           config_radt_cld_scheme .eq. 'cld_fraction_thompson')) then

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
 if(.not. (config_sfclayer_scheme .eq. 'off'           .or.    &
           config_sfclayer_scheme .eq. 'sf_mynn'       .or. &
           config_sfclayer_scheme .eq. 'sf_mynn_wrf38' .or. &
           config_sfclayer_scheme .eq. 'sf_monin_obukhov')) then
 
    write(mpas_err_message,'(A,A10)') 'illegal value for surface layer scheme: ', &
          trim(config_sfclayer_scheme)
    call physics_error_fatal(mpas_err_message)
 else
    if(config_pbl_scheme == 'bl_mynn') then
       config_sfclayer_scheme = 'sf_mynn'
    elseif(config_pbl_scheme == 'bl_mynn_wrf38') then
       config_sfclayer_scheme = 'sf_mynn_wrf38'
    elseif(config_pbl_scheme == 'bl_ysu') then
       config_sfclayer_scheme = 'sf_monin_obukhov'
    endif
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
 
 if(config_microp_scheme     .eq. 'off' .and. &
    config_convection_scheme .eq. 'off' .and. &
    config_lsm_scheme        .eq. 'off' .and. &
    config_pbl_scheme        .eq. 'off' .and. &
    config_radt_lw_scheme    .eq. 'off' .and. &
    config_radt_sw_scheme    .eq. 'off' .and. &
    config_sfclayer_scheme   .eq. 'off') moist_physics = .false.

 write(0,*) '    config_microp_scheme       = ', trim(config_microp_scheme)
 write(0,*) '    config_convection_scheme   = ', trim(config_convection_scheme)
 write(0,*) '    config_lsm_scheme          = ', trim(config_lsm_scheme)
 write(0,*) '    config_pbl_scheme          = ', trim(config_pbl_scheme)
 write(0,*) '    config_gwdo_scheme         = ', trim(config_gwdo_scheme)
 write(0,*) '    config_radt_cld_scheme     = ', trim(config_radt_cld_scheme)
 write(0,*) '    config_radt_lw_scheme      = ', trim(config_radt_lw_scheme)
 write(0,*) '    config_radt_sw_scheme      = ', trim(config_radt_sw_scheme)
 write(0,*) '    config_sfclayer_scheme     = ', trim(config_sfclayer_scheme)
 write(0,*)

 end subroutine physics_namelist_check

!=================================================================================================================
 subroutine physics_registry_init(mesh,configs,sfc_input)
!=================================================================================================================

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

!-----------------------------------------------------------------------------------------------------------------

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

!=================================================================================================================
 end module mpas_atmphys_control
!=================================================================================================================

