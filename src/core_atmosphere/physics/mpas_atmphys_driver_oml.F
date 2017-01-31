! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_oml
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_constants
 use mpas_atmphys_landuse
 use mpas_atmphys_lsm_noahinit
 use mpas_atmphys_vars
 use mpas_constants

!wrf physics
 use module_sf_oml
 
 implicit none
 private
 public :: driver_oml1d

 integer,private:: i,j,k,n

!MPAS driver for 1d ocean mixed layer
!Bill Skamarock (send comments to skamaroc@ucar.edu).
!2014-06-11.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added update of the sea-surface temperature below call to subroutine oml1d.
!   Laura D. Fowler (laura@ucar.edu) / 2014-07-24.


 contains


!=================================================================================================================
 subroutine driver_oml1d(configs,mesh,diag,diag_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input


  logical,pointer:: config_oml1d
  real(kind=RKIND),pointer:: oml_gamma
  real(kind=RKIND),pointer:: oml_relaxation_time

  real(kind=RKIND), dimension(:), pointer:: t_oml
  real(kind=RKIND), dimension(:), pointer:: t_oml_initial
  real(kind=RKIND), dimension(:), pointer:: t_oml_200m_initial
  real(kind=RKIND), dimension(:), pointer:: h_oml
  real(kind=RKIND), dimension(:), pointer:: h_oml_initial
  real(kind=RKIND), dimension(:), pointer:: hu_oml
  real(kind=RKIND), dimension(:), pointer:: hv_oml

  real(kind=RKIND), dimension(:), pointer:: sst
  real(kind=RKIND), dimension(:), pointer:: skintemp
  real(kind=RKIND), dimension(:), pointer:: hfx
  real(kind=RKIND), dimension(:), pointer:: lh
  real(kind=RKIND), dimension(:), pointer:: gsw
  real(kind=RKIND), dimension(:), pointer:: glw
  real(kind=RKIND), dimension(:), pointer:: sfc_emiss
  real(kind=RKIND), dimension(:), pointer:: ust
  real(kind=RKIND), dimension(:), pointer:: xland

  real(kind=RKIND), dimension(:,:), pointer:: uReconstructZonal
  real(kind=RKIND), dimension(:,:), pointer:: uReconstructmeridional
  real(kind=RKIND), dimension(:)  , pointer:: latCell

  integer, pointer:: nCells, nCellsSolve

  integer :: iCell
  real(kind=RKIND):: f_coriolis

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine driver_oml1d:'

! namelist parameters for ocean mixed layer model
 call mpas_pool_get_config(configs,'config_oml1d'               ,config_oml1d )
 call mpas_pool_get_config(configs,'config_oml_gamma'           ,oml_gamma )
 call mpas_pool_get_config(configs,'config_oml_relaxation_time' ,oml_relaxation_time )

! state and initial state for ocean mixed layer model
 call mpas_pool_get_array(diag_physics,'t_oml'              ,t_oml )
 call mpas_pool_get_array(diag_physics,'t_oml_initial'      ,t_oml_initial )
 call mpas_pool_get_array(diag_physics,'t_oml_200m_initial' ,t_oml_200m_initial )
 call mpas_pool_get_array(diag_physics,'h_oml'              ,h_oml )
 call mpas_pool_get_array(diag_physics,'h_oml_initial'      ,h_oml_initial )
 call mpas_pool_get_array(diag_physics,'hu_oml'             ,hu_oml )
 call mpas_pool_get_array(diag_physics,'hv_oml'             ,hv_oml )

! state and diagnostics from mpas
 call mpas_pool_get_array(diag_physics,'hfx'       ,hfx )
 call mpas_pool_get_array(diag_physics,'lh'        ,lh )
 call mpas_pool_get_array(diag_physics,'gsw'       ,gsw )
 call mpas_pool_get_array(diag_physics,'glw'       ,glw )
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'ust'       ,ust )

 call mpas_pool_get_array(sfc_input,'sst'       ,sst )
 call mpas_pool_get_array(sfc_input,'skintemp'  ,skintemp )
 call mpas_pool_get_array(sfc_input,'xland'     ,xland )

 call mpas_pool_get_array(diag,'uReconstructZonal'      ,uReconstructZonal )
 call mpas_pool_get_array(diag,'uReconstructMeridional' ,uReconstructMeridional )

 call mpas_pool_get_array(mesh,'latCell'      ,latCell )

 call mpas_pool_get_dimension(mesh,'nCells'       ,nCells )
 call mpas_pool_get_dimension(mesh,'nCellsSolve'  ,nCellsSolve )

 if (config_oml1d) then
 
   do iCell = 1, nCellsSolve

! if ocean point, call the 1d ocean mixed layer model 
     if( xland(iCell) .gt. 1.5) then
       f_coriolis = 2.*omega*cos(latCell(iCell))
       call oml1d( t_oml(iCell), t_oml_initial(iCell), h_oml(iCell), h_oml_initial(iCell),  &
                   hu_oml(iCell), hv_oml(iCell), skintemp(iCell), hfx(iCell),               &
                   lh(iCell), gsw(iCell), glw(iCell), t_oml_200m_initial(iCell),            &
                   uReconstructZonal(1,iCell), uReconstructMeridional(1,iCell),             &
                   ust(iCell), f_coriolis, sfc_emiss(iCell), stbolt, gravity, dt_pbl,       &
                   oml_gamma, oml_relaxation_time                                          )

       sst(iCell) = skintemp(iCell)

     end if

   end do

 end if

!write(0,*) '--- end subroutine driver_oml1d.'

 end subroutine driver_oml1d

!=================================================================================================================
 end module mpas_atmphys_driver_oml
!=================================================================================================================
