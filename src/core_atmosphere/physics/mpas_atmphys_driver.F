! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=============================================================================================
 module mpas_atmphys_driver
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_driver_cloudiness
 use mpas_atmphys_driver_convection_deep
 use mpas_atmphys_driver_gwdo
 use mpas_atmphys_driver_pbl
 use mpas_atmphys_driver_lsm
 use mpas_atmphys_driver_radiation_sw 
 use mpas_atmphys_driver_radiation_lw
 use mpas_atmphys_driver_sfclayer
 use mpas_atmphys_constants
 use mpas_atmphys_update
 use mpas_atmphys_vars, only: l_camlw,l_conv,l_radtlw,l_radtsw
#ifdef non_hydrostatic_core
 use mpas_atmphys_interface_nhyd
#elif hydrostatic_core
 use mpas_atmphys_interface_hyd
#endif

 implicit none
 private
 public:: physics_driver

 contains

!=============================================================================================
 subroutine physics_driver(domain,itimestep,xtime_s)
!=============================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(domain_type),intent(inout):: domain

!local variables:
 type(block_type),pointer:: block
 integer:: i,j,k

!=============================================================================================

 if(config_conv_deep_scheme    .ne. 'off' .or. &
    config_lsm_scheme          .ne. 'off' .or. &
    config_pbl_scheme          .ne. 'off' .or. & 
    config_radt_lw_scheme      .ne. 'off' .or. &
    config_radt_sw_scheme      .ne. 'off' .or. &
    config_sfclayer_scheme     .ne. 'off') then

 block => domain % blocklist
 do while(associated(block))

    !allocate arrays shared by all physics parameterizations:
    call allocate_forall_physics

    !physics prep step:
#ifdef non_hydrostatic_core
    call MPAS_to_physics(block%mesh,block%state%time_levs(1)%state,block%diag, &
                         block%diag_physics)
#elif hydrostatic_core
    call MPAS_to_physics(block%state%time_levs(1)%state,block%diag)
#endif

    !call to cloud scheme:
    if(l_radtlw .or. l_radtsw) then
       call allocate_cloudiness
       call driver_cloudiness(block%diag_physics)
    endif

    !call to short wave radiation scheme:
    if(l_radtsw) then
       call allocate_radiation_sw(xtime_s)
       call driver_radiation_sw(itimestep,block%mesh,block%state%time_levs(1)%state,  &
                               block%diag_physics,block%sfc_input,block%tend_physics, &
                               xtime_s)
    endif

    !call to long wave radiation scheme:
    if(l_radtlw) then
       call allocate_radiation_lw(xtime_s)
       call driver_radiation_lw(xtime_s,block%mesh,block%state%time_levs(1)%state, &
                              block%diag_physics,block%sfc_input,block%tend_physics)
    endif
!   if(l_camlw .and. config_radt_lw_scheme .eq. 'cam_lw') &
!                             call radiation_camlw_to_MPAS(block%diag_physics)

    !call to accumulate long- and short-wave diagnostics if needed:
     if(config_bucket_update /= 'none' .and. config_bucket_radt .gt. 0._RKIND) &
        call update_radiation_diagnostics(config_bucket_radt,block%mesh,block%diag_physics)

    !deallocate all radiation arrays:
    if(config_radt_sw_scheme.ne.'off' .or. config_radt_lw_scheme.ne.'off') &
       call deallocate_cloudiness
    if(config_radt_sw_scheme.ne.'off') call deallocate_radiation_sw
    if(config_radt_lw_scheme.ne.'off') call deallocate_radiation_lw

    !call to surface-layer scheme:
    if(config_sfclayer_scheme .ne. 'off') then
       call allocate_sfclayer
       call driver_sfclayer(block%mesh,block%diag_physics,block%sfc_input)
       call deallocate_sfclayer
    endif

    !call to land-surface scheme:
    if(config_lsm_scheme .ne. 'off') then
       call allocate_lsm
       call driver_lsm(itimestep,block%mesh,block%diag_physics,block%sfc_input)
       call deallocate_lsm
    endif

    !call to pbl schemes:
    if(config_pbl_scheme .ne. 'off' .and. config_sfclayer_scheme .ne. 'off') then
       call allocate_pbl
       call driver_pbl(block%sfc_input,block%diag_physics,block%tend_physics)
       call deallocate_pbl
    endif

    !call to gravity wave drag over orography scheme:
    if(config_gwdo_scheme .ne. 'off') then
       call allocate_gwdo
       call driver_gwdo(itimestep,block%mesh,block%sfc_input,block%diag_physics, &
                        block%tend_physics)
       call deallocate_gwdo
    endif

    !call to convection scheme:
    call update_convection_step1(block%mesh,block%diag_physics,block%tend_physics)
    if(l_conv) then
       call allocate_convection_deep
       call driver_convection_deep(itimestep,block%mesh,block%sfc_input,block%diag_physics, &
                                   block%tend_physics)
       call deallocate_convection_deep
    endif
    !update diagnostics:
    call update_convection_step2(config_bucket_rainc,block%mesh,block%diag_physics)

    !deallocate arrays shared by all physics parameterizations:
    call deallocate_forall_physics

    block => block % next
 end do 

 endif
 write(0,*)
 write(0,*) '--- end physics_driver:'
 write(0,*)

 end subroutine physics_driver

!=============================================================================================
 end module mpas_atmphys_driver
!=============================================================================================
