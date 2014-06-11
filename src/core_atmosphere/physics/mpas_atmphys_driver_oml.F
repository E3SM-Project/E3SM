! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_oml
 use mpas_kind_types
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_landuse
 use mpas_atmphys_lsm_noahinit
 use mpas_atmphys_vars

!wrf physics
 use module_sf_oml
 
 implicit none
 private

 integer,private:: i,j,k,n

!>\brief MPAS driver for 1d ocean mixed layer
!>\author Bill Skamarock (send comments to skamaroc@ucar.edu).
!>\date 2014-06-11.
!>

 contains

!==================================================================================================
 subroutine driver_olm(configs,mesh,diag,diag_physics)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics


  logical,pointer:: config_oml1d
  real(kind=RKIND),pointer:: olm_gamma

  real(kind=RKIND),pointer:: t_oml
  real(kind=RKIND),pointer:: t_oml_initial
  real(kind=RKIND),pointer:: t_oml_200m_initial
  real(kind=RKIND),pointer:: h_oml
  real(kind=RKIND),pointer:: h_oml_initial
  real(kind=RKIND),pointer:: hu_oml
  real(kind=RKIND),pointer:: hv_oml

  real(kind=RKIND),pointer:: skintemp
  real(kind=RKIND),pointer:: hfx
  real(kind=RKIND),pointer:: lh
  real(kind=RKIND),pointer:: gws
  real(kind=RKIND),pointer:: glw
  real(kind=RKIND),pointer:: sfc_emiss
  real(kind=RKIND),pointer:: ust
  real(kind=RKIND),pointer:: xland

  real(kind=RKIND),pointer:: uReconstructZonal
  real(kind=RKIND),pointer:: uReconstructmeridional
  real(kind=RKIND),pointer:: latCell

  integer, pointer:: nCells, nCellsSolve

  integer :: iCell
  real(kind=RKIND):: f_coriolis

!---------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine driver_olm:'

! namelist parameters for ocean mixed layer model
 call mpas_pool_get_config(configs,'config_olm1d'       ,config_olm1d )
 call mpas_pool_get_config(configs,'config_olm_gamma'   ,olm_gamma )

! state and initial state for ocean mixed layer model
 call mpas_pool_get_config(diag_physics,'t_oml'              ,t_oml )
 call mpas_pool_get_config(diag_physics,'t_oml_initial'      ,t_oml_initial )
 call mpas_pool_get_config(diag_physics,'t_oml_200m_initial' ,t_oml_200m_initial )
 call mpas_pool_get_config(diag_physics,'h_oml'              ,h_oml )
 call mpas_pool_get_config(diag_physics,'h_oml_initial'      ,h_oml_initial )
 call mpas_pool_get_config(diag_physics,'hu_oml'             ,hu_oml )
 call mpas_pool_get_config(diag_physics,'hv_oml'             ,hv_oml )

! state and diagnostics from mpas
 call mpas_pool_get_config(diag_physics,'skintemp'  ,skintemp )
 call mpas_pool_get_config(diag_physics,'hfx'       ,hfx )
 call mpas_pool_get_config(diag_physics,'lh'        ,lh )
 call mpas_pool_get_config(diag_physics,'gsw'       ,gsw )
 call mpas_pool_get_config(diag_physics,'glw'       ,glw )
 call mpas_pool_get_config(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_config(diag_physics,'ust'       ,ust )
 call mpas_pool_get_config(diag_physics,'xland'     ,xland )

 call mpas_pool_get_config(diag,'uReconstructZonal'      ,uReconstructZonal )
 call mpas_pool_get_config(diag,'uReconstructMeridional' ,uReconstructMeridional )

 call mpas_pool_get_config(mesh,'latCell'      ,latCell )
 call mpas_pool_get_config(mesh,'nCells'       ,nCells )
 call mpas_pool_get_config(mesh,'nCellsSolve'  ,nCellsSolve )

 if (config_olm1d) then
 
   do iCell = 1, nCellsSolve

! if ocean point, call the 1d ocean mixed layer model 
     if( xland(iCell) .gt. 1.5) then
       f_corioils = 2.*omega*cos(latCell(iCell)
       call oml1d( t_oml(iCell), t_oml_initial(iCell), h_oml(iCell), h_oml_initial(iCell),  &
                   hu_oml(iCell), hv_oml(iCell), skintemp(iCell), hfx(iCell), lh(iCell),    &
                   lh(iCell), gsw(iCell), glw(iCell), t_oml_200m_initial(iCell),            &
                   uReconstructZonal(1,iCell), uReconstructMeridional(1,iCell),             &
                   ust(iCell), f_coriolis, sfc_emiss(iCell), stbolt, gravity, dt_pbl,       &
                   olm_gamma                                                               )
     end if
   end do

 end if

 end subroutine driver_olm

!==================================================================================================
 end module mpas_atmphys_driver_olm
!==================================================================================================
