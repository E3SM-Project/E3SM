! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_driver_cloudiness
 use mpas_atmphys_driver_convection_deep
 use mpas_atmphys_driver_gwdo
 use mpas_atmphys_driver_lsm
 use mpas_atmphys_driver_pbl
 use mpas_atmphys_driver_radiation_lw
 use mpas_atmphys_driver_radiation_sw 
 use mpas_atmphys_driver_sfclayer
 use mpas_atmphys_constants
 use mpas_atmphys_interface_nhyd
 use mpas_atmphys_update
 use mpas_atmphys_vars, only: l_camlw,l_conv,l_radtlw,l_radtsw

 implicit none
 private
 public:: physics_driver


!>\brief MPAS top physics driver.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!> subroutine physics_driver is the top physics driver from which separate drivers for all physics
!> parameterizations, except cloud microphysics parameterizations are called.
!>
!> subroutines called in mpas_atmphys_driver:
!> ------------------------------------------
!> allocate_forall_physics     : allocate local arrays defining atmospheric soundings (pressure,..)
!> allocate_cloudiness         : allocate all local arrays used in driver_cloudiness. 
!> allocate_convection_deep    : allocate all local arrays used in driver_convection_deep.
!> allocate_gwdo               : allocate all local arrays used in driver_gwdo.
!> allocate_lsm                : allocate all local arrays used in driver_lsm.
!> allocate_pbl                : allocate all local arrays used in driver_pbl.
!> allocate_radiation_lw       : allocate all local arrays used in driver_radiation_lw.
!> allocate_radiation_sw       : allocate all local arrays used in driver_radiation_sw.
!> allocate_sfclayer           : allocate all local arrays used in driver_sfclayer.
!>
!> deallocate_forall_physics   : deallocate local arrays defining atmospheric soundings.
!> deallocate_cloudiness       : dedeallocate all local arrays used in driver_cloudiness.
!> deallocate_convection_deep  : deallocate all local arrays used in driver_convection_deep.
!> deallocate_gwdo             : deallocate all local arrays used in driver_gwdo.
!> deallocate_lsm              : deallocate all local arrays used in driver_lsm.
!> deallocate_pbl              : deallocate all local arrays used in driver_pbl.
!> deallocate_radiation_lw     : deallocate all local arrays used in driver_radiation_lw.
!> deallocate_radiation_sw     : deallocate all local arrays used in driver_radiation_sw.
!> deallocate_sfclayer         : deallocate all local arrays used in driver_sfclayer.
!>
!> MPAS_to_physics             :
!> driver_cloudiness           : driver for parameterization of fractional cloudiness.
!> driver_convection_deep      : driver for parameterization of convection.
!> driver_gwdo                 : driver for parameterization of gravity wave drag over orography.
!> driver_lsm                  : driver for land-surface scheme.
!> driver_pbl                  : driver for planetary boundary layer scheme.
!> driver_radiation_sw         : driver for short wave radiation schemes.
!> driver_radiation_lw         : driver for long wave radiation schemes.
!> driver_sfclayer             : driver for surface layer scheme.
!> update_convection_step1     : updates lifetime of deep convective clouds in Kain-Fritsch scheme.
!> update_convection_step2     : updates accumulated precipitation output from convection schemes.
!> update_radiation_diagnostics: updates accumualted radiation diagnostics from radiation schemes.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed call to calculate atmospheric soundings for the hydrostatic dynamical core.
!>      Laura D. Fowler (2013-05-06).
!>    * removed the namelist option config_eddy_scheme and associated sourcecode.
!>    * removed the namelist option config_conv_shallow_scheme and associated sourcecode.
!>      Laura D. Fowler (birch.ucar.edu) / 2013-05-29.

 contains

!==================================================================================================
 subroutine physics_driver(domain,itimestep,xtime_s)
!==================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(domain_type),intent(inout):: domain

!local variables:
 type(block_type),pointer:: block

!==================================================================================================

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
    call MPAS_to_physics(block%mesh,block%state%time_levs(1)%state,block%diag, &
                         block%diag_physics)

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

!==================================================================================================
 end module mpas_atmphys_driver
!==================================================================================================
