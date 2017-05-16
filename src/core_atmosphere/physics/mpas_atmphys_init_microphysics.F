! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!=================================================================================================================
 module mpas_atmphys_init_microphysics
 use mpas_dmpar
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_utilities
!use module_mp_thompson, only: is_aerosol_aware,naCCN0,naCCN1,naIN0,naIN1,ntb_arc,ntb_arw,ntb_art,ntb_arr, &
!                              ntb_ark,tnccn_act

 implicit none
 private
 public:: init_thompson_clouddroplets_forMPAS

!MPAS main initialization of the Thompson parameterization of cloud microphysics with nucleation of cloud
!droplets based on distributions of CCNs and INs (aerosol-aware parameterization).
!Laura D. Fowler (send comments to laura@ucar.edu).
!2016-03-28.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added "use mpas_dmpar" at the top of the module.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-04.


 contains


!=================================================================================================================
 subroutine init_thompson_clouddroplets_forMPAS(mesh,sfc_input,diag_physics)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer:: nt_c,mu_c

 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine init_thompson_clouddroplets_forMPAS:')

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(sfc_input,'landmask',landmask)
 call mpas_pool_get_array(diag_physics,'nt_c',nt_c)
 call mpas_pool_get_array(diag_physics,'mu_c',mu_c)

!... initialize the prescribed number of cloud droplets, and mu_c (parameter in the exponential of the generalized
!gamma distribution) as a function of the land-cean mask. as set in the thompson cloud microphysics scheme, nt_c
!is set to 100 per cc (100.E6 m^-3) for maritime cases and 300 per cc (300.E6 m^-3) for continental cases.
 do iCell = 1, nCellsSolve
    if(landmask(iCell) .eq. 1) then
       nt_c(iCell) = 300.e6 
    elseif(landmask(iCell) .eq. 0) then
       nt_c(iCell) = 100.e6
    endif
    mu_c(iCell) = MIN(15., (1000.e6/nt_c(iCell) + 2.))
 enddo

!call mpas_log_write('--- end subroutine init_thompson_clouddroplets_forMPAS.')

 end subroutine init_thompson_clouddroplets_forMPAS 

!=================================================================================================================
 end module mpas_atmphys_init_microphysics
!=================================================================================================================

 
 
