! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_gwdo
 use mpas_configure, only: len_disp => config_len_disp
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_vars

!wrf physics:
 use module_bl_gwdo

 implicit none
 private
 public:: allocate_gwdo,   &
          deallocate_gwdo, &
          driver_gwdo

 integer,private:: i,j,k


!>\brief MPAS driver for parameterization of gravity wave drag over orography.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_driver_gwdo:
!> ----------------------------------------
!> allocate_gwdo  : allocate local arrays for parameterization of gravity wave drag.
!> deallocate_gwdo: deallocate local arrays for parameterization of gravity wave drag.
!> driver_gwdo    : main driver (called from subroutine physics_driver).
!> gwdo_from_MPAS : initialize local arrays.
!> gwdo_to_MPAS   : copy local arrays to MPAS arrays.
!>
!> WRF physics called from driver_gwdo:
!> --------------------------- --------
!>    * module_bl_gwdo       : parameterization of gravity wave drag over orography. 
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutine gwdo. 
!>      Laura D. Fowler (birch.ucar.edu) / 2013-05-29.
!>    * changed the definition of dx_p to the mean distance between cell centers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-08-23.

 contains

!==================================================================================================
 subroutine allocate_gwdo
!==================================================================================================

 if(.not.allocated(dx_p)    ) allocate(dx_p(ims:ime,jms:jme)    )
 if(.not.allocated(var2d_p) ) allocate(var2d_p(ims:ime,jms:jme) )
 if(.not.allocated(con_p)   ) allocate(con_p(ims:ime,jms:jme)   )
 if(.not.allocated(oa1_p)   ) allocate(oa1_p(ims:ime,jms:jme)   )
 if(.not.allocated(oa2_p)   ) allocate(oa2_p(ims:ime,jms:jme)   )
 if(.not.allocated(oa3_p)   ) allocate(oa3_p(ims:ime,jms:jme)   )
 if(.not.allocated(oa4_p)   ) allocate(oa4_p(ims:ime,jms:jme)   )
 if(.not.allocated(ol1_p)   ) allocate(ol1_p(ims:ime,jms:jme)   )
 if(.not.allocated(ol2_p)   ) allocate(ol2_p(ims:ime,jms:jme)   )
 if(.not.allocated(ol3_p)   ) allocate(ol3_p(ims:ime,jms:jme)   )
 if(.not.allocated(ol4_p)   ) allocate(ol4_p(ims:ime,jms:jme)   )
 if(.not.allocated(kpbl_p  )) allocate(kpbl_p(ims:ime,jms:jme)  )
 if(.not.allocated(dusfcg_p)) allocate(dusfcg_p(ims:ime,jms:jme))
 if(.not.allocated(dvsfcg_p)) allocate(dvsfcg_p(ims:ime,jms:jme))
 
 if(.not.allocated(dtaux3d_p)) allocate(dtaux3d_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(dtauy3d_p)) allocate(dtauy3d_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rublten_p)) allocate(rublten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rvblten_p)) allocate(rvblten_p(ims:ime,kms:kme,jms:jme))

 end subroutine allocate_gwdo

!==================================================================================================
 subroutine deallocate_gwdo
!==================================================================================================

 if(allocated(dx_p)    ) deallocate(dx_p    )
 if(allocated(var2d_p) ) deallocate(var2d_p )
 if(allocated(con_p)   ) deallocate(con_p   )
 if(allocated(oa1_p)   ) deallocate(oa1_p   )
 if(allocated(oa2_p)   ) deallocate(oa2_p   )
 if(allocated(oa3_p)   ) deallocate(oa3_p   )
 if(allocated(oa4_p)   ) deallocate(oa4_p   )
 if(allocated(ol1_p)   ) deallocate(ol1_p   )
 if(allocated(ol2_p)   ) deallocate(ol2_p   )
 if(allocated(ol3_p)   ) deallocate(ol3_p   )
 if(allocated(ol4_p)   ) deallocate(ol4_p   )
 if(allocated(kpbl_p  )) deallocate(kpbl_p  )
 if(allocated(dusfcg_p)) deallocate(dusfcg_p)
 if(allocated(dvsfcg_p)) deallocate(dvsfcg_p)
 
 if(allocated(dtaux3d_p)) deallocate(dtaux3d_p)
 if(allocated(dtauy3d_p)) deallocate(dtauy3d_p)
 if(allocated(rublten_p)) deallocate(rublten_p)
 if(allocated(rvblten_p)) deallocate(rvblten_p)

 end subroutine deallocate_gwdo

!==================================================================================================
 subroutine gwdo_from_MPAS(mesh,sfc_input,diag_physics,tend_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 type(sfc_input_type),intent(in)   :: sfc_input
 type(diag_physics_type),intent(in):: diag_physics
 type(tend_physics_type),intent(in):: tend_physics

!--------------------------------------------------------------------------------------------------

 do j = jts,jte
 do i = its,ite
    var2d_p(i,j) = sfc_input % var2d % array(i)
    con_p(i,j)   = sfc_input % con   % array(i)
    oa1_p(i,j)   = sfc_input % oa1   % array(i)
    oa2_p(i,j)   = sfc_input % oa2   % array(i)
    oa3_p(i,j)   = sfc_input % oa3   % array(i)
    oa4_p(i,j)   = sfc_input % oa4   % array(i)
    ol1_p(i,j)   = sfc_input % ol1   % array(i)
    ol2_p(i,j)   = sfc_input % ol2   % array(i)
    ol3_p(i,j)   = sfc_input % ol3   % array(i)
    ol4_p(i,j)   = sfc_input % ol4   % array(i)
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    dx_p(i,j)     = diag_physics % dcEdge_m % array(i)
    kpbl_p(i,j)   = diag_physics % kpbl     % array(i)
    dusfcg_p(i,j) = diag_physics % dusfcg   % array(i)
    dvsfcg_p(i,j) = diag_physics % dvsfcg   % array(i)
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    dtaux3d_p(i,k,j) = diag_physics % dtaux3d % array(k,i)
    dtauy3d_p(i,k,j) = diag_physics % dtauy3d % array(k,i)
    rublten_p(i,k,j) = tend_physics % rublten % array(k,i)
    rvblten_p(i,k,j) = tend_physics % rvblten % array(k,i)
 enddo
 enddo
 enddo

 end subroutine gwdo_from_MPAS
 
!==================================================================================================
 subroutine gwdo_to_MPAS(diag_physics,tend_physics)
!==================================================================================================

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!--------------------------------------------------------------------------------------------------

 do j = jts,jte
 do i = its,ite
    diag_physics % dusfcg % array(i) = dusfcg_p(i,j) 
    diag_physics % dvsfcg % array(i) = dvsfcg_p(i,j)
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    diag_physics % dtaux3d % array(k,i) = dtaux3d_p(i,k,j)
    diag_physics % dtauy3d % array(k,i) = dtauy3d_p(i,k,j)
    diag_physics % rubldiff % array(k,i) = rublten_p(i,k,j)-tend_physics%rublten%array(k,i)
    diag_physics % rvbldiff % array(k,i) = rvblten_p(i,k,j)-tend_physics%rvblten%array(k,i) 

    tend_physics % rublten % array(k,i) = rublten_p(i,k,j)
    tend_physics % rvblten % array(k,i) = rvblten_p(i,k,j)
 enddo
 enddo
 enddo

!write(0,*)
!write(0,*) '--- end subroutine gwdo_to_MPAS:'
!do i = its,ite
!   write(0,101) i,diag_physics%dusfcg%array(i),diag_physics%dvsfcg%array(i)
!enddo
!101 format(i8,2(1x,e15.8))

 end subroutine gwdo_to_MPAS
 
!==================================================================================================
 subroutine driver_gwdo(itimestep,mesh,sfc_input,diag_physics,tend_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 type(sfc_input_type),intent(in):: sfc_input
 integer,intent(in):: itimestep

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!local variables:
 integer:: i,iCell,iEdge
 real(kind=RKIND),dimension(:),allocatable:: dx_max

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine driver_gwdo: dt_pbl=',dt_pbl

!copy MPAS arrays to local arrays:
 call gwdo_from_MPAS(mesh,sfc_input,diag_physics,tend_physics)

 gwdo_select: select case (trim(gwdo_scheme))

    case("ysu_gwdo")
       call gwdo ( &
                  p3d       = pres_hydd_p , p3di      = pres2_hydd_p , pi3d    = pi_p      , &
                  u3d       = u_p         , v3d       = v_p          , t3d     = t_p       , & 
                  qv3d      = qv_p        , z         = z_p          , rublten = rublten_p , &
                  rvblten   = rvblten_p   , dtaux3d   = dtaux3d_p    , dtauy3d = dtauy3d_p , &
                  dusfcg    = dusfcg_p    , dvsfcg    = dvsfcg_p     , kpbl2d  = kpbl_p    , &
                  itimestep = itimestep   , dt        = dt_pbl       , dx      = dx_p      , & 
                  cp        = cp          , g         = g            , rd      = R_d       , & 
                  rv        = R_v         , ep1       = ep_1         , pi      = pii       , & 
                  var2d     = var2d_p     , oc12d     = con_p        , oa2d1   = oa1_p     , & 
                  oa2d2     = oa2_p       , oa2d3     = oa3_p        , oa2d4   = oa4_p     , &
                  ol2d1     = ol1_p       , ol2d2     = ol2_p        , ol2d3   = ol3_p     , & 
                  ol2d4     = ol4_p       ,                                                  &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,    &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,    &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte      &
                 )

     case default

 end select gwdo_select

!copy local arrays to MPAS grid:
 call gwdo_to_MPAS(diag_physics,tend_physics)
 write(0,*) '--- end subroutine driver_gwdo'

 end subroutine driver_gwdo

!==================================================================================================
 end module mpas_atmphys_driver_gwdo
!==================================================================================================
