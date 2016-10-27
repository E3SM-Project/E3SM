! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_pbl
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_constants
 use mpas_atmphys_vars

!wrf physics:
 use module_bl_mynn
 use module_bl_ysu

 implicit none
 private
 public:: allocate_pbl,   &
          deallocate_pbl, &
          driver_pbl

!MPAS driver for parameterization of Planetary Boundary Layer (PBL) processes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_pbl:
! ---------------------------------------
! allocate_pbl  : allocate local arrays for parameterization of PBL processes.
! deallocate_pbl: deallocate local arrays for parameterization of PBL processes.
! driver_pbl    : main driver (called from subroutine physics_driver).
! pbl_from_MPAS : initialize local arrays.
! pbl_to_MPAS   : copy local arrays to MPAS arrays.
!
! WRF physics called from driver_pbl:
! -----------------------------------
! * module_bl_ysu : YSU PBL scheme.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutine ysu.
!   Laura D. Fowler (birch.ucar.edu) / 2013-05-29.
! * in call to subroutine ysu, replaced the variable g (that originally pointed to gravity)
!   with gravity, for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed "ysu" with "bl_ysu".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25.
! * added the implementation of the MYNN PBL scheme from WRF 3.6.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
! * corrected the initialization of sh3d for the mynn parameterization.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-13.
! * for the mynn parameterization, change the definition of dx_p to match that used in other physics
!   parameterizations.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-18.
! * updated the call to subroutine ysu in comjunction with updating module_bl_ysu.F from WRF version 3.6.1 to
!   WRF version 3.8.1
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-27.


 contains


!=================================================================================================================
 subroutine allocate_pbl
!=================================================================================================================

 if(.not.allocated(hfx_p)  ) allocate(hfx_p(ims:ime,jms:jme)  )
 if(.not.allocated(qfx_p)  ) allocate(qfx_p(ims:ime,jms:jme)  )
 if(.not.allocated(ust_p)  ) allocate(ust_p(ims:ime,jms:jme)  )
 if(.not.allocated(wspd_p) ) allocate(wspd_p(ims:ime,jms:jme) )
 if(.not.allocated(xland_p)) allocate(xland_p(ims:ime,jms:jme))
 if(.not.allocated(hpbl_p) ) allocate(hpbl_p(ims:ime,jms:jme) )
 if(.not.allocated(kpbl_p) ) allocate(kpbl_p(ims:ime,jms:jme) )
 if(.not.allocated(znt_p)  ) allocate(znt_p(ims:ime,jms:jme)  )
 if(.not.allocated(delta_p)) allocate(delta_p(ims:ime,jms:jme))
 if(.not.allocated(wstar_p)) allocate(wstar_p(ims:ime,jms:jme))
 if(.not.allocated(uoce_p) ) allocate(uoce_p(ims:ime,jms:jme) )
 if(.not.allocated(voce_p) ) allocate(voce_p(ims:ime,jms:jme) )


 !tendencies:
 if(.not.allocated(rublten_p) ) allocate(rublten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rvblten_p) ) allocate(rvblten_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rthblten_p)) allocate(rthblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqvblten_p)) allocate(rqvblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqcblten_p)) allocate(rqcblten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqiblten_p)) allocate(rqiblten_p(ims:ime,kms:kme,jms:jme))

 !exchange coefficients:
 if(.not.allocated(kzh_p)) allocate(kzh_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(kzm_p)) allocate(kzm_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(kzq_p)) allocate(kzq_p(ims:ime,kms:kme,jms:jme))

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       !from surface-layer model:
       if(.not.allocated(br_p)    ) allocate(br_p(ims:ime,jms:jme)          )
       if(.not.allocated(ctopo_p) ) allocate(ctopo_p(ims:ime,jms:jme)       )
       if(.not.allocated(ctopo2_p)) allocate(ctopo2_p(ims:ime,jms:jme)      )
       if(.not.allocated(psih_p)  ) allocate(psih_p(ims:ime,jms:jme)        )
       if(.not.allocated(psim_p)  ) allocate(psim_p(ims:ime,jms:jme)        )
       if(.not.allocated(regime_p)) allocate(regime_p(ims:ime,jms:jme)      )
       if(.not.allocated(u10_p)   ) allocate(u10_p(ims:ime,jms:jme)         )
       if(.not.allocated(v10_p)   ) allocate(v10_p(ims:ime,jms:jme)         )
       if(.not.allocated(exch_p)  ) allocate(exch_p(ims:ime,kms:kme,jms:jme))
       !from radiation schemes:
       if(.not.allocated(rthraten_p)) allocate(rthraten_p(ims:ime,kms:kme,jms:jme))

    case("bl_mynn")
       if(.not.allocated(dx_p)   ) allocate(dx_p(ims:ime,jms:jme)   )
       if(.not.allocated(ch_p)   ) allocate(ch_p(ims:ime,jms:jme)   )
       if(.not.allocated(qcg_p)  ) allocate(qcg_p(ims:ime,jms:jme)  )
       if(.not.allocated(qsfc_p) ) allocate(qsfc_p(ims:ime,jms:jme) )
       if(.not.allocated(rmol_p) ) allocate(rmol_p(ims:ime,jms:jme) )
       if(.not.allocated(tsk_p)  ) allocate(tsk_p(ims:ime,jms:jme)  )
       if(.not.allocated(vdfg_p) ) allocate(vdfg_p(ims:ime,jms:jme) )

       if(.not.allocated(cov_p)   ) allocate(cov_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(qke_p)   ) allocate(qke_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(qsq_p)   ) allocate(qsq_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(tsq_p)   ) allocate(tsq_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(qkeadv_p)) allocate(qkeadv_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(elpbl_p) ) allocate(elpbl_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(tkepbl_p)) allocate(tkepbl_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(sh3d_p)  ) allocate(sh3d_p(ims:ime,kms:kme,jms:jme)  )

       if(.not.allocated(dqke_p)  ) allocate(dqke_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(qbuoy_p) ) allocate(qbuoy_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(qdiss_p) ) allocate(qdiss_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(qshear_p)) allocate(qshear_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qwt_p)   ) allocate(qwt_p(ims:ime,kms:kme,jms:jme)   )

       if(.not.allocated(rniblten_p)) allocate(rniblten_p(ims:ime,kms:kme,jms:jme))

    case default

 end select pbl_select

 end subroutine allocate_pbl

!=================================================================================================================
 subroutine deallocate_pbl
!=================================================================================================================

 if(allocated(hfx_p)  ) deallocate(hfx_p  )
 if(allocated(qfx_p)  ) deallocate(qfx_p  )
 if(allocated(ust_p)  ) deallocate(ust_p  )
 if(allocated(wspd_p) ) deallocate(wspd_p )
 if(allocated(xland_p)) deallocate(xland_p)
 if(allocated(hpbl_p) ) deallocate(hpbl_p )
 if(allocated(kpbl_p) ) deallocate(kpbl_p )
 if(allocated(znt_p)  ) deallocate(znt_p  )
 if(allocated(delta_p)) deallocate(delta_p)
 if(allocated(wstar_p)) deallocate(wstar_p)
 if(allocated(uoce_p) ) deallocate(uoce_p )
 if(allocated(voce_p) ) deallocate(voce_p )

 !tendencies:
 if(allocated(rublten_p) ) deallocate(rublten_p )
 if(allocated(rvblten_p) ) deallocate(rvblten_p )
 if(allocated(rthblten_p)) deallocate(rthblten_p)
 if(allocated(rqvblten_p)) deallocate(rqvblten_p)
 if(allocated(rqcblten_p)) deallocate(rqcblten_p)
 if(allocated(rqiblten_p)) deallocate(rqiblten_p)

 !exchange coefficients:
 if(allocated(kzh_p)) deallocate(kzh_p)
 if(allocated(kzm_p)) deallocate(kzm_p)
 if(allocated(kzq_p)) deallocate(kzq_p)

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       !from surface-layer model:
       if(allocated(br_p)    ) deallocate(br_p    )
       if(allocated(ctopo_p) ) deallocate(ctopo_p )
       if(allocated(ctopo2_p)) deallocate(ctopo2_p)
       if(allocated(psih_p)  ) deallocate(psih_p  )
       if(allocated(psim_p)  ) deallocate(psim_p  )
       if(allocated(regime_p)) deallocate(regime_p)
       if(allocated(u10_p)   ) deallocate(u10_p   )
       if(allocated(v10_p)   ) deallocate(v10_p   )
       if(allocated(exch_p)  ) deallocate(exch_p  )
       !from radiation schemes:
       if(allocated(rthraten_p)) deallocate(rthraten_p)

    case("bl_mynn")
       if(allocated(dx_p)   ) deallocate(dx_p   )
       if(allocated(ch_p)   ) deallocate(ch_p   )
       if(allocated(qcg_p)  ) deallocate(qcg_p  )
       if(allocated(qsfc_p) ) deallocate(qsfc_p )
       if(allocated(rmol_p) ) deallocate(rmol_p )
       if(allocated(tsk_p)  ) deallocate(tsk_p  )
       if(allocated(vdfg_p) ) deallocate(vdfg_p )

       if(allocated(cov_p)   ) deallocate(cov_p   )
       if(allocated(qke_p)   ) deallocate(qke_p   )
       if(allocated(qsq_p)   ) deallocate(qsq_p   )
       if(allocated(tsq_p)   ) deallocate(tsq_p   )
       if(allocated(qkeadv_p)) deallocate(qkeadv_p)
       if(allocated(elpbl_p) ) deallocate(elpbl_p )
       if(allocated(tkepbl_p)) deallocate(tkepbl_p)
       if(allocated(sh3d_p)  ) deallocate(sh3d_p  )
       if(allocated(dqke_p)  ) deallocate(dqke_p  )
       if(allocated(qbuoy_p) ) deallocate(qbuoy_p )
       if(allocated(qdiss_p) ) deallocate(qdiss_p )
       if(allocated(qshear_p)) deallocate(qshear_p)
       if(allocated(qwt_p)   ) deallocate(qwt_p   )

       if(allocated(rniblten_p)) deallocate(rniblten_p)

    case default

 end select pbl_select

 end subroutine deallocate_pbl

!=================================================================================================================
 subroutine pbl_from_MPAS(configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: tend_physics

 integer,intent(in):: its,ite

!local variables:
 integer:: i,k,j

!local pointers:
 real(kind=RKIND),dimension(:),pointer:: hfx,hpbl,qfx,ust,wspd,xland,znt
 real(kind=RKIND),dimension(:),pointer:: delta,wstar

!local pointers for YSU scheme:
 logical,pointer:: config_ysu_pblmix
 real(kind=RKIND),dimension(:),pointer:: br,fh,fm,regime,u10,v10
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw

!local pointers for MYNN scheme:
 real(kind=RKIND),pointer:: len_disp
 real(kind=RKIND),dimension(:),pointer  :: meshDensity
 real(kind=RKIND),dimension(:),pointer  :: ch,qsfc,qcg,rmol,skintemp
 real(kind=RKIND),dimension(:,:),pointer:: cov,qke,qsq,tsq,sh3d,tke_pbl,qke_adv,el_pbl

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'hfx'   ,hfx   )
 call mpas_pool_get_array(diag_physics,'hpbl'  ,hpbl  )
 call mpas_pool_get_array(diag_physics,'qfx'   ,qfx   )
 call mpas_pool_get_array(diag_physics,'ust'   ,ust   )
 call mpas_pool_get_array(diag_physics,'wspd'  ,wspd  )
 call mpas_pool_get_array(diag_physics,'znt'   ,znt   )
 call mpas_pool_get_array(diag_physics,'delta' ,delta )
 call mpas_pool_get_array(diag_physics,'wstar' ,wstar )

 call mpas_pool_get_array(sfc_input   ,'xland' ,xland )

 do j = jts,jte
 do i = its,ite
    !from surface-layer model:
    hfx_p(i,j)    = hfx(i)
    hpbl_p(i,j)   = hpbl(i)
    qfx_p(i,j)    = qfx(i)
    ust_p(i,j)    = ust(i)
    wspd_p(i,j)   = wspd(i)
    xland_p(i,j)  = xland(i)
    kpbl_p(i,j)   = 1
    znt_p(i,j)    = znt(i)
    delta_p(i,j)  = delta(i)
    wstar_p(i,j)  = wstar(i)
    !... ocean currents are set to zero:
    uoce_p(i,j)  = 0._RKIND
    voce_p(i,j)  = 0._RKIND
 enddo
 enddo

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       call mpas_pool_get_config(configs,'config_ysu_pblmix',config_ysu_pblmix)

       call mpas_pool_get_array(diag_physics,'br'    ,br    )
       call mpas_pool_get_array(diag_physics,'fm'    ,fm    )
       call mpas_pool_get_array(diag_physics,'fh'    ,fh    )
       call mpas_pool_get_array(diag_physics,'regime',regime)
       call mpas_pool_get_array(diag_physics,'u10'   ,u10   )
       call mpas_pool_get_array(diag_physics,'v10'   ,v10   )

       call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
       call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)

       ysu_pblmix = 0
       if(config_ysu_pblmix) ysu_pblmix = 1

       do j = jts,jte
       do i = its,ite
          !from surface-layer model:
          br_p(i,j)     = br(i)
          psim_p(i,j)   = fm(i)
          psih_p(i,j)   = fh(i)
          regime_p(i,j) = regime(i)
          u10_p(i,j)    = u10(i)
          v10_p(i,j)    = v10(i)
          !initialization for YSU PBL scheme:
          ctopo_p(i,j)  = 1._RKIND
          ctopo2_p(i,j) = 1._RKIND
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          exch_p(i,k,j)     = 0._RKIND
          rthraten_p(i,k,j) = rthratenlw(k,i) + rthratensw(k,i) 
       enddo
       enddo
       enddo

    case("bl_mynn")
       call mpas_pool_get_config(configs,'config_len_disp',len_disp)
       call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

       call mpas_pool_get_array(sfc_input   ,'skintemp',skintemp)
       call mpas_pool_get_array(diag_physics,'ch'      ,ch      )
       call mpas_pool_get_array(diag_physics,'qcg'     ,qcg     )
       call mpas_pool_get_array(diag_physics,'qsfc'    ,qsfc    )
       call mpas_pool_get_array(diag_physics,'rmol'    ,rmol    )

       call mpas_pool_get_array(diag_physics,'el_pbl' ,el_pbl   )
       call mpas_pool_get_array(diag_physics,'cov'    ,cov      )
       call mpas_pool_get_array(diag_physics,'qke'    ,qke      )
       call mpas_pool_get_array(diag_physics,'qke_adv',qke_adv  )
       call mpas_pool_get_array(diag_physics,'qsq'    ,qsq      )
       call mpas_pool_get_array(diag_physics,'tsq'    ,tsq      )
       call mpas_pool_get_array(diag_physics,'tke_pbl',tke_pbl  )
       call mpas_pool_get_array(diag_physics,'sh3d'   ,sh3d     )

       do j = jts,jte
       do i = its,ite
          dx_p(i,j)    = len_disp / meshDensity(i)**0.25
          ch_p(i,j)    = ch(i)
          qcg_p(i,j)   = qcg(i)
          qsfc_p(i,j)  = qsfc(i)
          rmol_p(i,j)  = rmol(i)
          tsk_p(i,j)   = skintemp(i)
          !... no gravitational settling of fog/cloud droplets (grav_settling = 0):
          vdfg_p(i,j)  = 0._RKIND
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          elpbl_p(i,k,j)    = el_pbl(k,i)
          cov_p(i,k,j)      = cov(k,i)
          qke_p(i,k,j)      = qke(k,i)
          qsq_p(i,k,j)      = qsq(k,i)
          tsq_p(i,k,j)      = tsq(k,i)
          tkepbl_p(i,k,j)   = tke_pbl(k,i)
          qkeadv_p(i,k,j)   = qke_adv(k,i)
          sh3d_p(i,k,j)     = sh3d(k,i)
          rniblten_p(i,k,j) = 0._RKIND

          !... outputs:
          dqke_p(i,k,j)   = 0._RKIND
          qbuoy_p(i,k,j)  = 0._RKIND
          qdiss_p(i,k,j)  = 0._RKIND
          qshear_p(i,k,j) = 0._RKIND
          qwt_p(i,k,j)    = 0._RKIND
       enddo
       enddo
       enddo

    case default

 end select pbl_select

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    rublten_p(i,k,j)  = 0._RKIND
    rvblten_p(i,k,j)  = 0._RKIND
    rthblten_p(i,k,j) = 0._RKIND
    rqvblten_p(i,k,j) = 0._RKIND
    rqcblten_p(i,k,j) = 0._RKIND
    rqiblten_p(i,k,j) = 0._RKIND

    kzh_p(i,k,j)    = 0._RKIND
    kzm_p(i,k,j)    = 0._RKIND
    kzq_p(i,k,j)    = 0._RKIND
 enddo
 enddo
 enddo

 end subroutine pbl_from_MPAS

!=================================================================================================================
 subroutine pbl_to_MPAS(diag_physics,tend_physics,its,ite)
!=================================================================================================================

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

 integer,intent(in):: its,ite

!local variables:
 integer:: i,k,j

!local pointers:
 integer,dimension(:),pointer:: kpbl
 real(kind=RKIND),dimension(:),pointer  :: hpbl
 real(kind=RKIND),dimension(:,:),pointer:: kzh,kzm,kzq
 real(kind=RKIND),dimension(:,:),pointer:: rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten, &
                                           rniblten

!local pointers for YSU scheme:
 real(kind=RKIND),dimension(:,:),pointer:: exch_h

!local pointers for MYNN scheme:
 real(kind=RKIND),dimension(:),pointer  :: delta,wstar
 real(kind=RKIND),dimension(:,:),pointer:: cov,qke,qsq,tsq,sh3d,tke_pbl,qke_adv,el_pbl,dqke,qbuoy, &
                                           qdiss,qshear,qwt

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'kpbl' ,kpbl )
 call mpas_pool_get_array(diag_physics,'hpbl' ,hpbl )
 call mpas_pool_get_array(diag_physics,'kzh'  ,kzh  )
 call mpas_pool_get_array(diag_physics,'kzm'  ,kzm  )
 call mpas_pool_get_array(diag_physics,'kzq'  ,kzq  )
 call mpas_pool_get_array(diag_physics,'delta',delta)
 call mpas_pool_get_array(diag_physics,'wstar',wstar)

 call mpas_pool_get_array(tend_physics,'rublten' ,rublten )
 call mpas_pool_get_array(tend_physics,'rvblten' ,rvblten )
 call mpas_pool_get_array(tend_physics,'rthblten',rthblten)
 call mpas_pool_get_array(tend_physics,'rqvblten',rqvblten)
 call mpas_pool_get_array(tend_physics,'rqcblten',rqcblten)
 call mpas_pool_get_array(tend_physics,'rqiblten',rqiblten)

 do j = jts,jte
 do i = its,ite
    hpbl(i)  = hpbl_p(i,j)
    kpbl(i)  = kpbl_p(i,j)
    delta(i) = delta_p(i,j)
    wstar(i) = wstar_p(i,j)
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    rublten(k,i)  = rublten_p(i,k,j)
    rvblten(k,i)  = rvblten_p(i,k,j)
    rthblten(k,i) = rthblten_p(i,k,j)
    rqvblten(k,i) = rqvblten_p(i,k,j)
    rqcblten(k,i) = rqcblten_p(i,k,j)
    rqiblten(k,i) = rqiblten_p(i,k,j)

    kzh(k,i) = kzh_p(i,k,j)
    kzm(k,i) = kzm_p(i,k,j)
    kzq(k,i) = kzq_p(i,k,j)
 enddo
 enddo
 enddo

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
       call mpas_pool_get_array(diag_physics,'exch_h',exch_h)

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          exch_h(k,i) = exch_p(i,k,j)
       enddo
       enddo
       enddo

    case("bl_mynn")
       call mpas_pool_get_array(diag_physics,'el_pbl'  ,el_pbl  )
       call mpas_pool_get_array(diag_physics,'cov'     ,cov     )
       call mpas_pool_get_array(diag_physics,'qke'     ,qke     )
       call mpas_pool_get_array(diag_physics,'qke_adv' ,qke_adv )
       call mpas_pool_get_array(diag_physics,'qsq'     ,qsq     )
       call mpas_pool_get_array(diag_physics,'tsq'     ,tsq     )
       call mpas_pool_get_array(diag_physics,'tke_pbl' ,tke_pbl )
       call mpas_pool_get_array(diag_physics,'sh3d'    ,sh3d    )
       call mpas_pool_get_array(diag_physics,'dqke'    ,dqke    )
       call mpas_pool_get_array(diag_physics,'qbuoy'   ,qbuoy   )
       call mpas_pool_get_array(diag_physics,'qdiss'   ,qdiss   )
       call mpas_pool_get_array(diag_physics,'qshear'  ,qshear  )
       call mpas_pool_get_array(diag_physics,'qwt'     ,qwt     )
       call mpas_pool_get_array(tend_physics,'rniblten',rniblten)

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          el_pbl(k,i)   = elpbl_p(i,k,j)
          cov(k,i)      = cov_p(i,k,j)
          qke(k,i)      = qke_p(i,k,j)
          qsq(k,i)      = qsq_p(i,k,j)
          tsq(k,i)      = tsq_p(i,k,j)
          sh3d(k,i)     = sh3d_p(i,k,j)
          tke_pbl(k,i)  = tkepbl_p(i,k,j)
          qke_adv(k,i)  = qkeadv_p(i,k,j)
          !... outputs:
          dqke(k,i)     = dqke_p(i,k,j)
          qbuoy(k,i)    = qbuoy_p(i,k,j)
          qdiss(k,i)    = qdiss_p(i,k,j)
          qshear(k,i)   = qshear_p(i,k,j)
          qwt(k,i)      = qwt_p(i,k,j)

          rniblten(k,i) = rniblten_p(i,k,j)
       enddo
       enddo
       enddo

    case default

 end select pbl_select

 end subroutine pbl_to_MPAS
 
!=================================================================================================================
 subroutine driver_pbl(itimestep,configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 logical,pointer:: config_do_restart

!local variables:
 integer:: initflag
 integer:: i,k,j

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine driver_pbl:'

 call mpas_pool_get_config(configs,'config_do_restart',config_do_restart)

!copy MPAS arrays to local arrays:
 call pbl_from_MPAS(configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)

 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0

 pbl_select: select case (trim(pbl_scheme))

    case("bl_ysu")
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
                 wstar    = wstar_p    , delta    = delta_p     , uoce     = uoce_p     , &
                 voce     = voce_p     , rthraten = rthraten_p  , u10      = u10_p      , &
                 v10      = v10_p      , ctopo    = ctopo_p     , ctopo2   = ctopo2_p   , &
                 regime   = regime_p   , rho      = rho_p       , kzhout   = kzh_p      , &
                 kzmout   = kzm_p      , kzqout   = kzq_p       ,                         & 
                 ysu_topdown_pblmix = ysu_pblmix ,                                        &
                 ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde  , &
                 ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme  , &
                 its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte    &
                )

    case("bl_mynn")
       call  mynn_bl_driver ( &
                 p        = pres_hyd_p  , exner    = pi_p        , ps       = psfc_p     , &
                 th       = th_p        , dz       = dz_p        , u        = u_p        , &
                 v        = v_p         , qv       = qv_p        , qc       = qc_p       , &
                 qi       = qi_p        , qni      = ni_p        , rho      = rho_p      , &
                 du       = rublten_p   , dv       = rvblten_p   , dth      = rthblten_p , &
                 dqv      = rqvblten_p  , dqc      = rqcblten_p  , dqi      = rqiblten_p , &
                 dqni     = rniblten_p  , flag_qc  = f_qc        , flag_qnc = f_qnc      , &
                 flag_qi  = f_qi        , flag_qni = f_qni       , kpbl     = kpbl_p     , &
                 pblh     = hpbl_p      , xland    = xland_p     , ts       = tsk_p      , &
                 hfx      = hfx_p       , qfx      = qfx_p       , ch       = ch_p       , &
                 sh3d     = sh3d_p      , tsq      = tsq_p       , qsq      = qsq_p      , &
                 cov      = cov_p       , el_pbl   = elpbl_p     , qsfc     = qsfc_p     , &
                 qcg      = qcg_p       , ust      = ust_p       , rmol     = rmol_p     , &
                 wspd     = wspd_p      , wstar    = wstar_p     , delta    = delta_p    , &
                 delt     = dt_pbl      , k_h      = kzh_p       , k_m      = kzm_p      , &
                 k_q      = kzq_p       , uoce     = uoce_p      , voce     = voce_p     , &
                 qke      = qke_p       , qke_adv  = qkeadv_p    , vdfg     = vdfg_p     , &
                 tke_pbl  = tkepbl_p    , dqke     = dqke_p      , qwt      = qwt_p      , &
                 qshear   = qshear_p    , qbuoy    = qbuoy_p     , qdiss    = qdiss_p    , &
                 initflag          = initflag          ,                                   &
                 grav_settling     = grav_settling     ,                                   &
                 bl_mynn_cloudpdf  = bl_mynn_cloudpdf  ,                                   &
                 bl_mynn_tkeadvect = bl_mynn_tkeadvect ,                                   &
                 bl_mynn_tkebudget = bl_mynn_tkebudget ,                                   &
                 ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde   , &
                 ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme   , &
                 its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     &
                            )

    case default

 end select pbl_select

!copy local arrays to MPAS grid:
 call pbl_to_MPAS(diag_physics,tend_physics,its,ite)

!write(0,*) '--- end subroutine driver_pbl.'

 end subroutine driver_pbl

!=================================================================================================================
 end module mpas_atmphys_driver_pbl
!=================================================================================================================
