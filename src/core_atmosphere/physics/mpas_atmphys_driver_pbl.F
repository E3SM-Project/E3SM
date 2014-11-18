! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_pbl
 use mpas_kind_types
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_vars

!wrf physics:
 use module_bl_ysu

 implicit none
 private
 public:: allocate_pbl,   &
          deallocate_pbl, &
          driver_pbl

 integer,private:: i,j,k


!>\brief MPAS driver for parameterization of Planetary Boundary Layer (PBL) processes.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_driver_pbl:
!> ---------------------------------------
!> allocate_pbl  : allocate local arrays for parameterization of PBL processes.
!> deallocate_pbl: deallocate local arrays for parameterization of PBL processes.
!> driver_pbl    : main driver (called from subroutine physics_driver).
!> pbl_from_MPAS : initialize local arrays.
!> pbl_to_MPAS   : copy local arrays to MPAS arrays.
!>
!> WRF physics called from driver_pbl:
!> -----------------------------------
!>    * module_bl_ysu : YSU PBL scheme.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutine ysu.
!>      Laura D. Fowler (birch.ucar.edu) / 2013-05-29.
!>    * in call to subroutine ysu, replaced the variable g (that originally pointed to gravity)
!>      with gravity, for simplicity.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-05-15.


 contains


!==================================================================================================
 subroutine allocate_pbl
!==================================================================================================

!from surface-layer model:
 if(.not.allocated(br_p)      ) allocate(br_p(ims:ime,jms:jme)               )
 if(.not.allocated(ctopo_p)   ) allocate(ctopo_p(ims:ime,jms:jme)            )
 if(.not.allocated(ctopo2_p)  ) allocate(ctopo2_p(ims:ime,jms:jme)           )
 if(.not.allocated(gz1oz0_p)  ) allocate(gz1oz0_p(ims:ime,jms:jme)           )
 if(.not.allocated(hfx_p)     ) allocate(hfx_p(ims:ime,jms:jme)              )
 if(.not.allocated(psih_p)    ) allocate(psih_p(ims:ime,jms:jme)             )
 if(.not.allocated(psim_p)    ) allocate(psim_p(ims:ime,jms:jme)             )
 if(.not.allocated(qfx_p)     ) allocate(qfx_p(ims:ime,jms:jme)              )
 if(.not.allocated(regime_p)  ) allocate(regime_p(ims:ime,jms:jme)           )
 if(.not.allocated(u10_p)     ) allocate(u10_p(ims:ime,jms:jme)              )
 if(.not.allocated(ust_p)     ) allocate(ust_p(ims:ime,jms:jme)              )
 if(.not.allocated(v10_p)     ) allocate(v10_p(ims:ime,jms:jme)              )
 if(.not.allocated(wspd_p)    ) allocate(wspd_p(ims:ime,jms:jme)             )
 if(.not.allocated(xland_p)   ) allocate(xland_p(ims:ime,jms:jme)            )
 if(.not.allocated(znt_p)     ) allocate(znt_p(ims:ime,jms:jme)              )
!for YSU PBL scheme:
 if(.not.allocated(hpbl_p)    ) allocate(hpbl_p(ims:ime,jms:jme)             )
 if(.not.allocated(kpbl_p)    ) allocate(kpbl_p(ims:ime,jms:jme)             )
 if(.not.allocated(exch_p)    ) allocate(exch_p(ims:ime,kms:kme,jms:jme)     )
!tendencies:
 if(.not.allocated(rublten_p) ) allocate(rublten_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(rvblten_p) ) allocate(rvblten_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(rthblten_p)) allocate(rthblten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rqvblten_p)) allocate(rqvblten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rqcblten_p)) allocate(rqcblten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rqiblten_p)) allocate(rqiblten_p(ims:ime,kms:kme,jms:jme) )
!temporary for debugging the YSU PBL scheme:
 if(.not.allocated(kzh_p)     ) allocate(kzh_p(ims:ime,kms:kme,jms:jme)      )
 if(.not.allocated(kzm_p)     ) allocate(kzm_p(ims:ime,kms:kme,jms:jme)      )
 if(.not.allocated(kzq_p)     ) allocate(kzq_p(ims:ime,kms:kme,jms:jme)      )

 end subroutine allocate_pbl

!==================================================================================================
 subroutine deallocate_pbl
!==================================================================================================

!from surface-layer model:
 if(allocated(br_p)      ) deallocate(br_p       )
 if(allocated(ctopo_p)   ) deallocate(ctopo_p    )
 if(allocated(ctopo2_p)  ) deallocate(ctopo2_p   )
 if(allocated(gz1oz0_p)  ) deallocate(gz1oz0_p   )
 if(allocated(hfx_p)     ) deallocate(hfx_p      )
 if(allocated(psih_p)    ) deallocate(psih_p     )
 if(allocated(psim_p)    ) deallocate(psim_p     )
 if(allocated(qfx_p)     ) deallocate(qfx_p      )
 if(allocated(regime_p)  ) deallocate(regime_p   )
 if(allocated(u10_p)     ) deallocate(u10_p      )
 if(allocated(ust_p)     ) deallocate(ust_p      )
 if(allocated(v10_p)     ) deallocate(v10_p      )
 if(allocated(wspd_p)    ) deallocate(wspd_p     )
 if(allocated(xland_p)   ) deallocate(xland_p    )
 if(allocated(znt_p)     ) deallocate(znt_p      )
!for YSU PBL scheme:
 if(allocated(hpbl_p)    ) deallocate(hpbl_p     )
 if(allocated(kpbl_p)    ) deallocate(kpbl_p     )
 if(allocated(exch_p)    ) deallocate(exch_p     )
!tendencies
 if(allocated(rublten_p) ) deallocate(rublten_p  )
 if(allocated(rvblten_p) ) deallocate(rvblten_p  )
 if(allocated(rthblten_p)) deallocate(rthblten_p )
 if(allocated(rqvblten_p)) deallocate(rqvblten_p )
 if(allocated(rqcblten_p)) deallocate(rqcblten_p )
 if(allocated(rqiblten_p)) deallocate(rqiblten_p )
!temporary for debugging the YSU PBL scheme:
 if(allocated(kzh_p)     ) deallocate(kzh_p      )
 if(allocated(kzm_p)     ) deallocate(kzm_p      )
 if(allocated(kzq_p)     ) deallocate(kzq_p      )

 end subroutine deallocate_pbl

!==================================================================================================
 subroutine pbl_from_MPAS(sfc_input,diag_physics)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input

!local pointers:
 real(kind=RKIND),dimension(:),pointer:: xland
 real(kind=RKIND),dimension(:),pointer:: br,gz1oz0,hfx,hpbl,fm,fh,qfx,regime,u10,ust,v10,wspd,znt

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'br'    ,br    )
 call mpas_pool_get_array(diag_physics,'gz1oz0',gz1oz0)
 call mpas_pool_get_array(diag_physics,'hfx'   ,hfx   )
 call mpas_pool_get_array(diag_physics,'hpbl'  ,hpbl  )
 call mpas_pool_get_array(diag_physics,'fm'    ,fm    )
 call mpas_pool_get_array(diag_physics,'fh'    ,fh    )
 call mpas_pool_get_array(diag_physics,'qfx'   ,qfx   )
 call mpas_pool_get_array(diag_physics,'regime',regime)
 call mpas_pool_get_array(diag_physics,'u10'   ,u10   )
 call mpas_pool_get_array(diag_physics,'ust'   ,ust   )
 call mpas_pool_get_array(diag_physics,'v10'   ,v10   )
 call mpas_pool_get_array(diag_physics,'wspd'  ,wspd  )
 call mpas_pool_get_array(diag_physics,'znt'   ,znt   )
 call mpas_pool_get_array(sfc_input   ,'xland' ,xland )

 do j = jts,jte
 do i = its,ite
!from surface-layer model:
    br_p(i,j)     = br(i)
    gz1oz0_p(i,j) = gz1oz0(i)
    hfx_p(i,j)    = hfx(i)
    hpbl_p(i,j)   = hpbl(i)
    psim_p(i,j)   = fm(i)
    psih_p(i,j)   = fh(i)
    qfx_p(i,j)    = qfx(i)
    regime_p(i,j) = regime(i)
    u10_p(i,j)    = u10(i)
    ust_p(i,j)    = ust(i)
    v10_p(i,j)    = v10(i)
    wspd_p(i,j)   = wspd(i)
    znt_p(i,j)    = znt(i)

    xland_p(i,j)  = xland(i)       
!initialization for YSU PBL scheme:
    ctopo_p(i,j)  = 1._RKIND
    ctopo2_p(i,j) = 1._RKIND
    kpbl_p(i,j)   = 1
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    exch_p(i,k,j)     = 0._RKIND
    rublten_p(i,k,j)  = 0._RKIND
    rvblten_p(i,k,j)  = 0._RKIND
    rthblten_p(i,k,j) = 0._RKIND
    rqvblten_p(i,k,j) = 0._RKIND
    rqcblten_p(i,k,j) = 0._RKIND
    rqiblten_p(i,k,j) = 0._RKIND
!temporary for debugging the YSU PBL scheme:
    kzh_p(i,k,j)    = 0._RKIND
    kzm_p(i,k,j)    = 0._RKIND
    kzq_p(i,k,j)    = 0._RKIND
 enddo
 enddo
 enddo

 end subroutine pbl_from_MPAS
 
!==================================================================================================
 subroutine pbl_to_MPAS(diag_physics,tend_physics)
!==================================================================================================

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,dimension(:),pointer:: kpbl
 real(kind=RKIND),dimension(:),pointer  :: hpbl
 real(kind=RKIND),dimension(:,:),pointer:: exch_h,kzh,kzm,kzq
 real(kind=RKIND),dimension(:,:),pointer:: rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'kpbl'  ,kpbl  )
 call mpas_pool_get_array(diag_physics,'hpbl'  ,hpbl  )
 call mpas_pool_get_array(diag_physics,'exch_h',exch_h)
 call mpas_pool_get_array(diag_physics,'kzh'   ,kzh   )
 call mpas_pool_get_array(diag_physics,'kzm'   ,kzm   )
 call mpas_pool_get_array(diag_physics,'kzq'   ,kzq   )

 call mpas_pool_get_array(tend_physics,'rublten',rublten  )
 call mpas_pool_get_array(tend_physics,'rvblten',rvblten  )
 call mpas_pool_get_array(tend_physics,'rthblten',rthblten)
 call mpas_pool_get_array(tend_physics,'rqvblten',rqvblten)
 call mpas_pool_get_array(tend_physics,'rqcblten',rqcblten)
 call mpas_pool_get_array(tend_physics,'rqiblten',rqiblten)

 do j = jts,jte
 do i = its,ite
    hpbl(i) = hpbl_p(i,j)
    kpbl(i) = kpbl_p(i,j)
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    exch_h(k,i) = exch_p(i,k,j)
    rublten(k,i) = rublten_p(i,k,j) 
    rvblten(k,i) = rvblten_p(i,k,j)
    rthblten(k,i) = rthblten_p(i,k,j)
    rqvblten(k,i) = rqvblten_p(i,k,j)
    rqcblten(k,i) = rqcblten_p(i,k,j)
    rqiblten(k,i) = rqiblten_p(i,k,j)
!temporary for debugging the YSU PBL scheme:
    kzh(k,i) = kzh_p(i,k,j)
    kzm(k,i) = kzm_p(i,k,j)
    kzq(k,i) = kzq_p(i,k,j)
 enddo
 enddo
 enddo

 end subroutine pbl_to_MPAS
 
!==================================================================================================
 subroutine driver_pbl(sfc_input,diag_physics,tend_physics)
!==================================================================================================

!input and output arguments:
!---------------------------
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!--------------------------------------------------------------------------------------------------
! write(0,*)
! write(0,*) '--- enter subroutine driver_pbl: dt_pbl=',dt_pbl

!copy MPAS arrays to local arrays:
 call pbl_from_MPAS(sfc_input,diag_physics)

 pbl_select: select case (trim(pbl_scheme))

    case("ysu")
       call ysu ( &
                 p3d      = pres_hyd_p , p3di     = pres2_hyd_p , psfc     = psfc_p     , &
                 th3d     = th_p       , t3d      = t_p         , dz8w     = dz_p       , &
                 pi3d     = pi_p       , u3d      = u_p         , v3d      = v_p        , &
                 qv3d     = qv_p       , qc3d     = qc_p        , qi3d     = qi_p       , &
                 rublten  = rublten_p  , rvblten  = rvblten_p   , rthblten = rthblten_p , &
                 rqvblten = rqvblten_p , rqcblten = rqcblten_p  , rqiblten = rqiblten_p , & 
                 flag_qi  = f_qi       , cp       = cp          , g        = gravity    , &
                 rovcp    = rcp        , rd       = R_d         , rovg     = rdg        , & 
                 ep1      = ep_1       , ep2      = ep_2        , karman   = karman     , &
                 xlv      = xlv        , rv       = R_v         , znt      = znt_p      , &
                 ust      = ust_p      , hpbl     = hpbl_p      , psim     = psim_p     , &
                 psih     = psih_p     , xland    = xland_p     , hfx      = hfx_p      , &
                 qfx      = qfx_p      , wspd     = wspd_p      , br       = br_p       , &
                 dt       = dt_pbl     , kpbl2d   = kpbl_p      , exch_h   = exch_p     , &
                 u10      = u10_p      , v10      = v10_p       , ctopo    = ctopo_p    , & 
                 ctopo2   = ctopo2_p   , regime   = regime_p    , rho      = rho_p      , &
                 kzhout   = kzh_p      , kzmout   = kzm_p       , kzqout   = kzq_p      , &
                 ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde  , &
                 ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme  , &
                 its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte    &
                )

     case default

 end select pbl_select

!copy local arrays to MPAS grid:
 call pbl_to_MPAS(diag_physics,tend_physics)
! write(0,*) '--- end subroutine driver_pbl'

 end subroutine driver_pbl

!==================================================================================================
 end module mpas_atmphys_driver_pbl
!==================================================================================================
