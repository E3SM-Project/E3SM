! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_sfclayer
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_constants
 use mpas_atmphys_vars

!wrf physics:
 use module_sf_mynn
 use module_sf_sfclay

 implicit none
 private
 public:: init_sfclayer,       &
          allocate_sfclayer,   &
          deallocate_sfclayer, &
          driver_sfclayer

 integer,parameter,private:: isfflx   = 1        !=1 for surface heat and moisture fluxes.
 integer,parameter,private:: isftcflx = 0        !=0,(Charnock and Carlson-Boland).
 integer,parameter,private:: iz0tlnd  = 0        !=0,(Carlson-Boland).
 integer,parameter,private:: scm_force_flux = 0  !SCM surface forcing by surface fluxes.
                                                 !0=no 1=yes (WRF single column model option only).

!MPAS driver for parameterization of the surface layer.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_sfclayer:
! --------------------------------------------
! allocate_sfclayer    : allocate local arrays for parameterization of surface layer.
! deallocate_sfclayer  : deallocate local arrays for parameterization of surface layer.
! init_sfclayer        : initialization of individual surface layer schemes.
! driver_sfclayer      : main driver (called from subroutine physics_driver).
! sfclayer_from_MPAS   : initialize local arrays.
! sfclayer_to_MPAS     : copy local arrays to MPAS arrays.
!
! WRF physics called from driver_sfclayer:
! ----------------------------------------
! * module_sf_sfclay: Monin-Obukhov surface layer scheme.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to the subroutine sfclay.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * updated the definition of the horizontal resolution to the actual mean distance between cell centers.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-23.
! * in call to subroutine sfclay, replaced the variable g (that originally pointed to gravity)
!   with gravity, for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * in subroutine sfclayer_from_MPAS, added initialization of ustm, cd, cda, ck, and cka. in
!   subroutine sfclayer_to_MPAS, filled diag_physics%ustm with ustm_p after call to subroutine sfclay.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-16. 
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * added initialization of local logical "allowed_to read" in subroutine init_sfclayer. This logical
!   is actually not used in subroutine sfclayinit.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-25. 
! * renamed "monin_obukhov" with "sf_monin_obukhov".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-25.
! * added the implementation of the MYNN surface layer scheme from WRF 3.6.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
! * added the calculation of surface layer variables over seaice cells when config_frac_seaice is set to true.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-03.
! * changed the definition of dx_p to match that used in other physics parameterizations.
!   parameterizations.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-18.

 contains


!=================================================================================================================
 subroutine allocate_sfclayer(config_frac_seaice)
!=================================================================================================================

 logical,intent(in):: config_frac_seaice
!-----------------------------------------------------------------------------------------------------------------

 if(.not.allocated(dx_p)    ) allocate(dx_p(ims:ime,jms:jme)    )
 if(.not.allocated(br_p)    ) allocate(br_p(ims:ime,jms:jme)    )
 if(.not.allocated(cd_p)    ) allocate(cd_p(ims:ime,jms:jme)    )
 if(.not.allocated(cda_p)   ) allocate(cda_p(ims:ime,jms:jme)   )
 if(.not.allocated(chs_p)   ) allocate(chs_p(ims:ime,jms:jme)   )
 if(.not.allocated(chs2_p)  ) allocate(chs2_p(ims:ime,jms:jme)  )
 if(.not.allocated(ck_p)    ) allocate(ck_p(ims:ime,jms:jme)    )
 if(.not.allocated(cka_p)   ) allocate(cka_p(ims:ime,jms:jme)   )
 if(.not.allocated(cpm_p)   ) allocate(cpm_p(ims:ime,jms:jme)   )
 if(.not.allocated(cqs2_p)  ) allocate(cqs2_p(ims:ime,jms:jme)  )
 if(.not.allocated(gz1oz0_p)) allocate(gz1oz0_p(ims:ime,jms:jme))
 if(.not.allocated(flhc_p)  ) allocate(flhc_p(ims:ime,jms:jme)  )
 if(.not.allocated(flqc_p)  ) allocate(flqc_p(ims:ime,jms:jme)  )
 if(.not.allocated(hfx_p)   ) allocate(hfx_p(ims:ime,jms:jme)   )
 if(.not.allocated(hpbl_p)  ) allocate(hpbl_p(ims:ime,jms:jme)  )
 if(.not.allocated(lh_p)    ) allocate(lh_p(ims:ime,jms:jme)    )
 if(.not.allocated(mavail_p)) allocate(mavail_p(ims:ime,jms:jme))
 if(.not.allocated(mol_p)   ) allocate(mol_p(ims:ime,jms:jme)   )
 if(.not.allocated(psih_p)  ) allocate(psih_p(ims:ime,jms:jme)  )
 if(.not.allocated(psim_p)  ) allocate(psim_p(ims:ime,jms:jme)  )
 if(.not.allocated(q2_p)    ) allocate(q2_p(ims:ime,jms:jme)    )
 if(.not.allocated(qfx_p)   ) allocate(qfx_p(ims:ime,jms:jme)   )
 if(.not.allocated(qgh_p)   ) allocate(qgh_p(ims:ime,jms:jme)   )
 if(.not.allocated(qsfc_p)  ) allocate(qsfc_p(ims:ime,jms:jme)  )
 if(.not.allocated(regime_p)) allocate(regime_p(ims:ime,jms:jme))
 if(.not.allocated(rmol_p)  ) allocate(rmol_p(ims:ime,jms:jme)  )
 if(.not.allocated(t2m_p)   ) allocate(t2m_p(ims:ime,jms:jme)   )
 if(.not.allocated(tsk_p)   ) allocate(tsk_p(ims:ime,jms:jme)   )
 if(.not.allocated(th2m_p)  ) allocate(th2m_p(ims:ime,jms:jme)  )
 if(.not.allocated(u10_p)   ) allocate(u10_p(ims:ime,jms:jme)   )
 if(.not.allocated(ust_p)   ) allocate(ust_p(ims:ime,jms:jme)   )
 if(.not.allocated(ustm_p)  ) allocate(ustm_p(ims:ime,jms:jme)  )
 if(.not.allocated(v10_p)   ) allocate(v10_p(ims:ime,jms:jme)   )
 if(.not.allocated(wspd_p)  ) allocate(wspd_p(ims:ime,jms:jme)  )
 if(.not.allocated(xland_p) ) allocate(xland_p(ims:ime,jms:jme) )
 if(.not.allocated(zol_p)   ) allocate(zol_p(ims:ime,jms:jme)   )
 if(.not.allocated(znt_p)   ) allocate(znt_p(ims:ime,jms:jme)   )

 if(config_frac_seaice) then
    if(.not.allocated(sst_p)      ) allocate(sst_p(ims:ime,jms:jme)      )
    if(.not.allocated(xice_p)     ) allocate(xice_p(ims:ime,jms:jme)     )

    if(.not.allocated(br_sea)     ) allocate(br_sea(ims:ime,jms:jme)     )
    if(.not.allocated(chs_sea)    ) allocate(chs_sea(ims:ime,jms:jme)    )
    if(.not.allocated(chs2_sea)   ) allocate(chs2_sea(ims:ime,jms:jme)   )
    if(.not.allocated(cqs2_sea)   ) allocate(cqs2_sea(ims:ime,jms:jme)   )
    if(.not.allocated(cpm_sea)    ) allocate(cpm_sea(ims:ime,jms:jme)    )
    if(.not.allocated(flhc_sea)   ) allocate(flhc_sea(ims:ime,jms:jme)   )
    if(.not.allocated(flqc_sea)   ) allocate(flqc_sea(ims:ime,jms:jme)   )
    if(.not.allocated(gz1oz0_sea) ) allocate(gz1oz0_sea(ims:ime,jms:jme) )
    if(.not.allocated(hfx_sea)    ) allocate(hfx_sea(ims:ime,jms:jme)    )
    if(.not.allocated(qfx_sea)    ) allocate(qfx_sea(ims:ime,jms:jme)    )
    if(.not.allocated(mavail_sea) ) allocate(mavail_sea(ims:ime,jms:jme) )
    if(.not.allocated(mol_sea)    ) allocate(mol_sea(ims:ime,jms:jme)    )
    if(.not.allocated(lh_sea)     ) allocate(lh_sea(ims:ime,jms:jme)     )
    if(.not.allocated(psih_sea)   ) allocate(psih_sea(ims:ime,jms:jme)   )
    if(.not.allocated(psim_sea)   ) allocate(psim_sea(ims:ime,jms:jme)   )
    if(.not.allocated(qgh_sea)    ) allocate(qgh_sea(ims:ime,jms:jme)    )
    if(.not.allocated(qsfc_sea)   ) allocate(qsfc_sea(ims:ime,jms:jme)   )
    if(.not.allocated(regime_sea) ) allocate(regime_sea(ims:ime,jms:jme) )
    if(.not.allocated(rmol_sea)   ) allocate(rmol_sea(ims:ime,jms:jme)   )
    if(.not.allocated(tsk_sea)    ) allocate(tsk_sea(ims:ime,jms:jme)    )
    if(.not.allocated(ust_sea)    ) allocate(ust_sea(ims:ime,jms:jme)    )
    if(.not.allocated(ustm_sea)   ) allocate(ustm_sea(ims:ime,jms:jme)   )
    if(.not.allocated(wspd_sea)   ) allocate(wspd_sea(ims:ime,jms:jme)   )
    if(.not.allocated(xland_sea)  ) allocate(xland_sea(ims:ime,jms:jme)  )
    if(.not.allocated(zol_sea)    ) allocate(zol_sea(ims:ime,jms:jme)    )
    if(.not.allocated(znt_sea)    ) allocate(znt_sea(ims:ime,jms:jme)    )

    if(.not.allocated(cd_sea)     ) allocate(cd_sea(ims:ime,jms:jme)     )
    if(.not.allocated(cda_sea)    ) allocate(cda_sea(ims:ime,jms:jme)    )
    if(.not.allocated(ck_sea)     ) allocate(ck_sea(ims:ime,jms:jme)     )
    if(.not.allocated(cka_sea)    ) allocate(cka_sea(ims:ime,jms:jme)    )
    if(.not.allocated(t2m_sea)    ) allocate(t2m_sea(ims:ime,jms:jme)    )
    if(.not.allocated(th2m_sea)   ) allocate(th2m_sea(ims:ime,jms:jme)   )
    if(.not.allocated(q2_sea)     ) allocate(q2_sea(ims:ime,jms:jme)     )
    if(.not.allocated(u10_sea)    ) allocate(u10_sea(ims:ime,jms:jme)    )
    if(.not.allocated(v10_sea)    ) allocate(v10_sea(ims:ime,jms:jme)    )

    if(.not.allocated(regime_hold)) allocate(regime_hold(ims:ime,jms:jme))
 endif

 sfclayer_select: select case (trim(sfclayer_scheme))

    case("sf_monin_obukhov")
       if(.not.allocated(fh_p)) allocate(fh_p(ims:ime,jms:jme))
       if(.not.allocated(fm_p)) allocate(fm_p(ims:ime,jms:jme))
       if(config_frac_seaice) then
          if(.not.allocated(fh_sea)) allocate(fh_sea(ims:ime,jms:jme))
          if(.not.allocated(fm_sea)) allocate(fm_sea(ims:ime,jms:jme))
       endif

    case("sf_mynn")
       if(.not.allocated(snowh_p)) allocate(snowh_p(ims:ime,jms:jme))
       if(.not.allocated(ch_p)   ) allocate(ch_p(ims:ime,jms:jme)   )
       if(.not.allocated(qcg_p)  ) allocate(qcg_p(ims:ime,jms:jme)  )
       if(config_frac_seaice) then
          if(.not.allocated(ch_sea)) allocate(ch_sea(ims:ime,jms:jme))
       endif

       if(.not.allocated(cov_p)     ) allocate(cov_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(qsq_p)     ) allocate(qsq_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(tsq_p)     ) allocate(tsq_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(sh3d_p)    ) allocate(sh3d_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(elpbl_p)   ) allocate(elpbl_p(ims:ime,kms:kme,jms:jme)   )

    case default

 end select sfclayer_select

 end subroutine allocate_sfclayer

!=================================================================================================================
 subroutine deallocate_sfclayer(config_frac_seaice)
!=================================================================================================================

 logical,intent(in):: config_frac_seaice
!-----------------------------------------------------------------------------------------------------------------

 if(allocated(dx_p)    ) deallocate(dx_p    )
 if(allocated(br_p)    ) deallocate(br_p    )
 if(allocated(cd_p)    ) deallocate(cd_p    )
 if(allocated(cda_p)   ) deallocate(cda_p   )
 if(allocated(chs_p)   ) deallocate(chs_p   )
 if(allocated(chs2_p)  ) deallocate(chs2_p  )
 if(allocated(ck_p)    ) deallocate(ck_p    )
 if(allocated(cka_p)   ) deallocate(cka_p   )
 if(allocated(cpm_p)   ) deallocate(cpm_p   )
 if(allocated(cqs2_p)  ) deallocate(cqs2_p  )
 if(allocated(gz1oz0_p)) deallocate(gz1oz0_p)
 if(allocated(flhc_p)  ) deallocate(flhc_p  )
 if(allocated(flqc_p)  ) deallocate(flqc_p  )
 if(allocated(hfx_p)   ) deallocate(hfx_p   )
 if(allocated(hpbl_p)  ) deallocate(hpbl_p  )
 if(allocated(lh_p)    ) deallocate(lh_p    )
 if(allocated(mavail_p)) deallocate(mavail_p)
 if(allocated(mol_p)   ) deallocate(mol_p   )
 if(allocated(psih_p)  ) deallocate(psih_p  )
 if(allocated(psim_p)  ) deallocate(psim_p  )
 if(allocated(q2_p)    ) deallocate(q2_p    )
 if(allocated(qfx_p)   ) deallocate(qfx_p   )
 if(allocated(qgh_p)   ) deallocate(qgh_p   )
 if(allocated(qsfc_p)  ) deallocate(qsfc_p  )
 if(allocated(regime_p)) deallocate(regime_p)
 if(allocated(rmol_p)  ) deallocate(rmol_p  )
 if(allocated(t2m_p)   ) deallocate(t2m_p   )
 if(allocated(tsk_p)   ) deallocate(tsk_p   )
 if(allocated(th2m_p)  ) deallocate(th2m_p  )
 if(allocated(u10_p)   ) deallocate(u10_p   )
 if(allocated(ust_p)   ) deallocate(ust_p   )
 if(allocated(ustm_p)  ) deallocate(ustm_p  )
 if(allocated(v10_p)   ) deallocate(v10_p   )
 if(allocated(wspd_p)  ) deallocate(wspd_p  )
 if(allocated(xland_p) ) deallocate(xland_p )
 if(allocated(zol_p)   ) deallocate(zol_p   )
 if(allocated(znt_p)   ) deallocate(znt_p   )

 if(config_frac_seaice) then
    if(allocated(sst_p)      ) deallocate(sst_p      )
    if(allocated(xice_p)     ) deallocate(xice_p     )

    if(allocated(br_sea)     ) deallocate(br_sea     )
    if(allocated(flhc_p)     ) deallocate(flhc_sea   )
    if(allocated(flqc_p)     ) deallocate(flqc_sea   )
    if(allocated(gz1oz0_sea) ) deallocate(gz1oz0_sea )
    if(allocated(mol_sea)    ) deallocate(mol_sea    )
    if(allocated(psih_sea)   ) deallocate(psih_sea   )
    if(allocated(psim_sea)   ) deallocate(psim_sea   )
    if(allocated(rmol_sea)   ) deallocate(rmol_sea   )
    if(allocated(ust_sea)    ) deallocate(ust_sea    )
    if(allocated(ustm_sea)   ) deallocate(ustm_sea   )
    if(allocated(wspd_sea)   ) deallocate(wspd_sea   )
    if(allocated(zol_sea)    ) deallocate(zol_sea    )
    if(allocated(cd_sea)     ) deallocate(cd_sea     )
    if(allocated(cda_sea)    ) deallocate(cda_sea    )
    if(allocated(ck_sea)     ) deallocate(ck_sea     )
    if(allocated(cka_sea)    ) deallocate(cka_sea    )
    if(allocated(t2m_sea)    ) deallocate(t2m_sea    )
    if(allocated(th2m_sea)   ) deallocate(th2m_sea   )
    if(allocated(q2_sea)     ) deallocate(q2_sea     )
    if(allocated(u10_sea)    ) deallocate(u10_sea    )
    if(allocated(v10_sea)    ) deallocate(v10_sea    )
    if(allocated(regime_hold)) deallocate(regime_hold)

    if(allocated(mavail_sea) ) deallocate(mavail_sea )
    if(allocated(tsk_sea)    ) deallocate(tsk_sea    )
    if(allocated(xland_sea)  ) deallocate(xland_sea  )
    if(allocated(znt_sea)    ) deallocate(znt_sea    )
 endif

 sfclayer_select: select case (trim(sfclayer_scheme))

    case("sf_monin_obukhov")
       if(allocated(fh_p)) deallocate(fh_p)
       if(allocated(fm_p)) deallocate(fm_p)
       if(config_frac_seaice) then
          if(allocated(fh_sea)) deallocate(fh_sea)
          if(allocated(fm_sea)) deallocate(fm_sea)
       endif

    case("sf_mynn")
       if(allocated(snowh_p)) deallocate(snowh_p)
       if(allocated(ch_p)   ) deallocate(ch_p   )
       if(allocated(qcg_p)  ) deallocate(qcg_p  )
       if(config_frac_seaice) then
          if(allocated(ch_sea)) deallocate(ch_sea)
       endif

       if(allocated(cov_p)     ) deallocate(cov_p     )
       if(allocated(qsq_p)     ) deallocate(qsq_p     )
       if(allocated(tsq_p)     ) deallocate(tsq_p     )
       if(allocated(sh3d_p)    ) deallocate(sh3d_p    )
       if(allocated(elpbl_p)   ) deallocate(elpbl_p   )

    case default

 end select sfclayer_select

 end subroutine deallocate_sfclayer

!=================================================================================================================
 subroutine sfclayer_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics

 integer,intent(in):: its,ite

!local variables:
 integer:: i,j,k

!local pointers:
 logical,pointer:: config_frac_seaice

 real(kind=RKIND),pointer:: len_disp
 real(kind=RKIND),dimension(:),pointer:: meshDensity
 real(kind=RKIND),dimension(:),pointer:: skintemp,sst,xice,xland
 real(kind=RKIND),dimension(:),pointer:: hpbl,mavail
 real(kind=RKIND),dimension(:),pointer:: br,cpm,chs,chs2,cqs2,flhc,flqc,gz1oz0,hfx,qfx,  &
                                         qgh,qsfc,lh,mol,psim,psih,regime,rmol,ust,ustm, &
                                         wspd,znt,zol

!local pointers specific to monin_obukhov:
 real(kind=RKIND),dimension(:),pointer:: fh,fm

!local pointers specific to mynn:
 real(kind=RKIND),dimension(:),pointer:: ch,qcg,snowh
 real(kind=RKIND),dimension(:,:),pointer:: cov,el_pbl,qsq,sh3d,tsq

!-----------------------------------------------------------------------------------------------------------------

!input variables:
 call mpas_pool_get_config(configs,'config_len_disp'   ,len_disp)
 call mpas_pool_get_config(configs,'config_frac_seaice',config_frac_seaice)
 call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

 call mpas_pool_get_array(diag_physics,'hpbl'    ,hpbl    )
 call mpas_pool_get_array(diag_physics,'mavail'  ,mavail  )
 call mpas_pool_get_array(sfc_input   ,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input   ,'xland'   ,xland   )

!inout variables:
 call mpas_pool_get_array(diag_physics,'br'      ,br      )
 call mpas_pool_get_array(diag_physics,'cpm'     ,cpm     )
 call mpas_pool_get_array(diag_physics,'chs'     ,chs     )
 call mpas_pool_get_array(diag_physics,'chs2'    ,chs2    )
 call mpas_pool_get_array(diag_physics,'cqs2'    ,cqs2    )
 call mpas_pool_get_array(diag_physics,'flhc'    ,flhc    )
 call mpas_pool_get_array(diag_physics,'flqc'    ,flqc    )
 call mpas_pool_get_array(diag_physics,'gz1oz0'  ,gz1oz0  )
 call mpas_pool_get_array(diag_physics,'hfx'     ,hfx     )
 call mpas_pool_get_array(diag_physics,'qfx'     ,qfx     )
 call mpas_pool_get_array(diag_physics,'qgh'     ,qgh     )
 call mpas_pool_get_array(diag_physics,'qsfc'    ,qsfc    )
 call mpas_pool_get_array(diag_physics,'lh'      ,lh      )
 call mpas_pool_get_array(diag_physics,'mol'     ,mol     )
 call mpas_pool_get_array(diag_physics,'psih'    ,psih    )
 call mpas_pool_get_array(diag_physics,'psim'    ,psim    )
 call mpas_pool_get_array(diag_physics,'regime'  ,regime  )
 call mpas_pool_get_array(diag_physics,'rmol'    ,rmol    )
 call mpas_pool_get_array(diag_physics,'ust'     ,ust     )
 call mpas_pool_get_array(diag_physics,'ustm'    ,ustm    )
 call mpas_pool_get_array(diag_physics,'wspd'    ,wspd    )
 call mpas_pool_get_array(diag_physics,'znt'     ,znt     )
 call mpas_pool_get_array(diag_physics,'zol'     ,zol     )

 do j = jts,jte
 do i = its,ite
    !input variables:
    dx_p(i,j)     = len_disp / meshDensity(i)**0.25
    hpbl_p(i,j)   = hpbl(i)
    mavail_p(i,j) = mavail(i)
    tsk_p(i,j)    = skintemp(i)
    xland_p(i,j)  = xland(i)

    !inout variables:
    br_p(i,j)     = br(i)
    cpm_p(i,j)    = cpm(i)
    chs_p(i,j)    = chs(i)
    chs2_p(i,j)   = chs2(i)
    cqs2_p(i,j)   = cqs2(i)
    flhc_p(i,j)   = flhc(i)
    flqc_p(i,j)   = flqc(i)
    gz1oz0_p(i,j) = gz1oz0(i)
    hfx_p(i,j)    = hfx(i)
    qfx_p(i,j)    = qfx(i)
    qgh_p(i,j)    = qgh(i)
    qsfc_p(i,j)   = qsfc(i)
    lh_p(i,j)     = lh(i)
    mol_p(i,j)    = mol(i)
    psim_p(i,j)   = psim(i)
    psih_p(i,j)   = psih(i)
    regime_p(i,j) = regime(i)
    rmol_p(i,j)   = rmol(i)
    ust_p(i,j)    = ust(i)
    wspd_p(i,j)   = wspd(i)
    znt_p(i,j)    = znt(i)
    zol_p(i,j)    = zol(i)

    !output variables:
    q2_p(i,j)     = 0._RKIND
    t2m_p(i,j)    = 0._RKIND
    th2m_p(i,j)   = 0._RKIND
    u10_p(i,j)    = 0._RKIND
    v10_p(i,j)    = 0._RKIND

    !output variables (optional):
    cd_p(i,j)     = 0._RKIND
    cda_p(i,j)    = 0._RKIND
    ck_p(i,j)     = 0._RKIND
    cka_p(i,j)    = 0._RKIND
    ustm_p(i,j)   = ustm(i)
 enddo
 enddo

 if(config_frac_seaice) then
    call mpas_pool_get_array(sfc_input,'sst' ,sst)
    call mpas_pool_get_array(sfc_input,'xice',xice)
    do j = jts,jte
    do i = its,ite
       sst_p(i,j)      = sst(i)
       xice_p(i,j)     = xice(i)

       !input variables:
       mavail_sea(i,j) = mavail(i)
       tsk_sea(i,j)    = skintemp(i)
       xland_sea(i,j)  = xland(i)
       !inout variables:
       br_sea(i,j)      = br(i)
       cpm_sea(i,j)     = cpm(i)
       chs_sea(i,j)     = chs(i)
       chs2_sea(i,j)    = chs2(i)
       cqs2_sea(i,j)    = cqs2(i)
       flhc_sea(i,j)    = flhc(i)
       flqc_sea(i,j)    = flqc(i)
       gz1oz0_sea(i,j)  = gz1oz0(i)
       lh_sea(i,j)      = lh(i)
       hfx_sea(i,j)     = hfx(i)
       qfx_sea(i,j)     = qfx(i)
       mol_sea(i,j)     = mol(i)
       psim_sea(i,j)    = psim(i)
       psih_sea(i,j)    = psih(i)
       qgh_sea(i,j)     = qgh(i)
       rmol_sea(i,j)    = rmol(i)
       regime_sea(i,j)  = regime(i)
       ust_sea(i,j)     = ust(i)
       ustm_sea(i,j)    = ustm(i)
       wspd_sea(i,j)    = wspd(i)
       zol_sea(i,j)     = zol(i)
       znt_sea(i,j)     = znt(i)
       regime_hold(i,j) = regime(i)
       !output variables:
       cd_sea(i,j)      = 0._RKIND
       cda_sea(i,j)     = 0._RKIND
       ck_sea(i,j)      = 0._RKIND
       cka_sea(i,j)     = 0._RKIND
       qsfc_sea(i,j)    = 0._RKIND
       q2_sea(i,j)      = 0._RKIND
       t2m_sea(i,j)     = 0._RKIND
       th2m_sea(i,j)    = 0._RKIND
       u10_sea(i,j)     = 0._RKIND
       v10_sea(i,j)     = 0._RKIND

       !overwrite some local variables for sea-ice cells:
       if(xice_p(i,j).ge.xice_threshold .and. xice_p(i,j).le.1._RKIND) then
          xland_sea(i,j)  = 2._RKIND
          mavail_sea(i,j) = 1._RKIND
          znt_sea(i,j)    = 0.0001_RKIND
          tsk_sea(i,j)    = max(sst_p(i,j),271.4_RKIND)
       else
          xland_sea(i,j)  = xland_p(i,j)
          mavail_sea(i,j) = mavail_p(i,j)
          znt_sea(i,j)    = znt_p(i,j)
          tsk_sea(i,j)    = tsk_p(i,j)
       endif
    enddo
    enddo
 endif

 sfclayer_select: select case (trim(sfclayer_scheme))

    case("sf_monin_obukhov")
       call mpas_pool_get_array(diag_physics,'fh',fh)
       call mpas_pool_get_array(diag_physics,'fm',fm)

       do j = jts,jte
       do i = its,ite
          fh_p(i,j) = fh(i)
          fm_p(i,j) = fm(i)
          if(config_frac_seaice) then
             fh_sea(i,j) = fh(i)
             fm_sea(i,j) = fm(i)
          endif
       enddo
       enddo

    case("sf_mynn")
       !input variables:
       call mpas_pool_get_array(diag_physics,'qcg'   ,qcg   )
       call mpas_pool_get_array(sfc_input   ,'snowh' ,snowh )
       call mpas_pool_get_array(diag_physics,'cov'   ,cov   )
       call mpas_pool_get_array(diag_physics,'el_pbl',el_pbl)
       call mpas_pool_get_array(diag_physics,'qsq'   ,qsq   )
       call mpas_pool_get_array(diag_physics,'sh3d'  ,sh3d  )
       call mpas_pool_get_array(diag_physics,'tsq'   ,tsq   )
       !inout variables:
       call mpas_pool_get_array(diag_physics,'ch',ch        )

       do j = jts,jte
       do i = its,ite
          !input variables:
          snowh_p(i,j) = snowh(i)
          qcg_p(i,j)   = qcg(i)
          !inout variables:
          ch_p(i,j)    = ch(i)
          if(config_frac_seaice) then
             ch_sea(i,j) = ch(i)
          endif
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          !input variables:
          cov_p(i,k,j)      = cov(k,i)
          qsq_p(i,k,j)      = qsq(k,i)
          tsq_p(i,k,j)      = tsq(k,i)
          sh3d_p(i,k,j)     = sh3d(k,i)
          elpbl_p(i,k,j)    = el_pbl(k,i)
       enddo
       enddo
       enddo

    case default

 end select sfclayer_select

 end subroutine sfclayer_from_MPAS

!=================================================================================================================
 subroutine sfclayer_to_MPAS(configs,sfc_input,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: sfc_input
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables:
 integer:: i,j

!local pointers:
 logical,pointer:: config_frac_seaice

 real(kind=RKIND),dimension(:),pointer:: br,cpm,chs,chs2,cqs2,flhc,flqc,gz1oz0,hfx,qfx,  &
                                         qgh,qsfc,lh,mol,psim,psih,regime,rmol,ust,wspd, &
                                         znt,zol
 real(kind=RKIND),dimension(:),pointer:: q2,t2m,th2m,u10,v10
 real(kind=RKIND),dimension(:),pointer:: cd,cda,ck,cka,ustm
 real(kind=RKIND),dimension(:),pointer:: xice

!local pointers specific to monin_obukhov:
 real(kind=RKIND),dimension(:),pointer:: fh,fm

!local pointers specific to mynn:
 real(kind=RKIND),dimension(:),pointer:: ch,qcg
 real(kind=RKIND),dimension(:,:),pointer:: cov,el_pbl,qsq,sh3d,tsq

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_frac_seaice',config_frac_seaice)

!inout variables:
 call mpas_pool_get_array(diag_physics,'br'    ,br    )
 call mpas_pool_get_array(diag_physics,'cpm'   ,cpm   )
 call mpas_pool_get_array(diag_physics,'chs'   ,chs   )
 call mpas_pool_get_array(diag_physics,'chs2'  ,chs2  )
 call mpas_pool_get_array(diag_physics,'cqs2'  ,cqs2  )
 call mpas_pool_get_array(diag_physics,'flhc'  ,flhc  )
 call mpas_pool_get_array(diag_physics,'flqc'  ,flqc  )
 call mpas_pool_get_array(diag_physics,'gz1oz0',gz1oz0)
 call mpas_pool_get_array(diag_physics,'hfx'   ,hfx   )
 call mpas_pool_get_array(diag_physics,'qfx'   ,qfx   )
 call mpas_pool_get_array(diag_physics,'qgh'   ,qgh   )
 call mpas_pool_get_array(diag_physics,'qsfc'  ,qsfc  )
 call mpas_pool_get_array(diag_physics,'lh'    ,lh    )
 call mpas_pool_get_array(diag_physics,'mol'   ,mol   )
 call mpas_pool_get_array(diag_physics,'psih'  ,psih  )
 call mpas_pool_get_array(diag_physics,'psim'  ,psim  )
 call mpas_pool_get_array(diag_physics,'regime',regime)
 call mpas_pool_get_array(diag_physics,'rmol'  ,rmol  )
 call mpas_pool_get_array(diag_physics,'ust'   ,ust   )
 call mpas_pool_get_array(diag_physics,'wspd'  ,wspd  )
 call mpas_pool_get_array(diag_physics,'znt'   ,znt   )
 call mpas_pool_get_array(diag_physics,'zol'   ,zol   )

!output variables:
 call mpas_pool_get_array(diag_physics,'q2'    ,q2    )
 call mpas_pool_get_array(diag_physics,'t2m'   ,t2m   )
 call mpas_pool_get_array(diag_physics,'th2m'  ,th2m  )
 call mpas_pool_get_array(diag_physics,'u10'   ,u10   )
 call mpas_pool_get_array(diag_physics,'v10'   ,v10   )

!output variables (optional):
 call mpas_pool_get_array(diag_physics,'cd'    ,cd    )
 call mpas_pool_get_array(diag_physics,'cda'   ,cda   )
 call mpas_pool_get_array(diag_physics,'ck'    ,ck    )
 call mpas_pool_get_array(diag_physics,'cka'   ,cka   )
 call mpas_pool_get_array(diag_physics,'ustm'  ,ustm  )

!output variables (optional):
 call mpas_pool_get_array(diag_physics,'cd'    ,cd    )
 call mpas_pool_get_array(diag_physics,'cda'   ,cda   )
 call mpas_pool_get_array(diag_physics,'ck'    ,ck    )
 call mpas_pool_get_array(diag_physics,'cka'   ,cka   )

 do j = jts,jte
 do i = its,ite
    !inout variables:
    br(i)     = br_p(i,j)
    cpm(i)    = cpm_p(i,j)
    chs(i)    = chs_p(i,j)
    chs2(i)   = chs2_p(i,j)
    cqs2(i)   = cqs2_p(i,j)
    flhc(i)   = flhc_p(i,j)
    flqc(i)   = flqc_p(i,j)
    gz1oz0(i) = gz1oz0_p(i,j)
    hfx(i)    = hfx_p(i,j)
    lh(i)     = lh_p(i,j)
    mol(i)    = mol_p(i,j)
    qfx(i)    = qfx_p(i,j)
    qgh(i)    = qgh_p(i,j)
    qsfc(i)   = qsfc_p(i,j)
    psim(i)   = psim_p(i,j)
    psih(i)   = psih_p(i,j)
    regime(i) = regime_p(i,j)
    rmol(i)   = rmol_p(i,j)
    ust(i)    = ust_p(i,j)
    wspd(i)   = wspd_p(i,j)
    zol(i)    = zol_p(i,j)
    znt(i)    = znt_p(i,j)
    !output variables:
    q2(i)     = q2_p(i,j)
    t2m(i)    = t2m_p(i,j)
    th2m(i)   = th2m_p(i,j)
    u10(i)    = u10_p(i,j)
    v10(i)    = v10_p(i,j)
    !output variables (optional):
    cd(i)     = cd_p(i,j)
    cda(i)    = cda_p(i,j)
    ck(i)     = ck_p(i,j)
    cka(i)    = cka_p(i,j)
    ustm(i)   = ustm_p(i,j)
 enddo
 enddo

 if(config_frac_seaice) then
    call mpas_pool_get_array(sfc_input,'xice',xice)
    do j = jts,jte
    do i = its,ite
       if(xice(i).ge.xice_threshold .and. xice(i).le.1._RKIND) then
          br(i)     = br_p(i,j)*xice(i)     + (1._RKIND-xice(i))*br_sea(i,j)
          flhc(i)   = flhc_p(i,j)*xice(i)   + (1._RKIND-xice(i))*flhc_sea(i,j)
          flqc(i)   = flqc_p(i,j)*xice(i)   + (1._RKIND-xice(i))*flqc_sea(i,j)
          gz1oz0(i) = gz1oz0_p(i,j)*xice(i) + (1._RKIND-xice(i))*gz1oz0_sea(i,j)
          mol(i)    = mol_p(i,j)*xice(i)    + (1._RKIND-xice(i))*mol_sea(i,j)
          psih(i)   = psih_p(i,j)*xice(i)   + (1._RKIND-xice(i))*psih_sea(i,j)
          psim(i)   = psim_p(i,j)*xice(i)   + (1._RKIND-xice(i))*psim_sea(i,j)
          rmol(i)   = rmol_p(i,j)*xice(i)   + (1._RKIND-xice(i))*rmol_sea(i,j)
          ust(i)    = ust_p(i,j)*xice(i)    + (1._RKIND-xice(i))*ust_sea(i,j)
          wspd(i)   = wspd_p(i,j)*xice(i)   + (1._RKIND-xice(i))*wspd_sea(i,j)
          zol(i)    = zol_p(i,j)*xice(i)    + (1._RKIND-xice(i))*zol_sea(i,j)
          if(xice(i) .ge. 0.5_RKIND) regime(i) = regime_hold(i,j)
          !output variables:
          q2(i)     = q2_p(i,j)*xice(i)   + (1._RKIND-xice(i))*q2_sea(i,j)
          t2m(i)    = t2m_p(i,j)*xice(i)  + (1._RKIND-xice(i))*t2m_sea(i,j)
          th2m(i)   = th2m_p(i,j)*xice(i) + (1._RKIND-xice(i))*th2m_sea(i,j)
          u10(i)    = u10_p(i,j)*xice(i)  + (1._RKIND-xice(i))*u10_sea(i,j)
          v10(i)    = v10_p(i,j)*xice(i)  + (1._RKIND-xice(i))*v10_sea(i,j)
          !output variables (optional):
          cd(i)     = cd_p(i,j)*xice(i)   + (1._RKIND-xice(i))*cd_sea(i,j)
          cda(i)    = cda_p(i,j)*xice(i)  + (1._RKIND-xice(i))*cda_sea(i,j)
          ck(i)     = ck_p(i,j)*xice(i)   + (1._RKIND-xice(i))*ck_sea(i,j)
          cka(i)    = cka_p(i,j)*xice(i)  + (1._RKIND-xice(i))*cka_sea(i,j)
          ustm(i)   = ustm_p(i,j)*xice(i)   + (1._RKIND-xice(i))*ustm_sea(i,j)
       endif
    enddo
    enddo
 endif

 sfclayer_select: select case (trim(sfclayer_scheme))

    case("sf_monin_obukhov")
       call mpas_pool_get_array(diag_physics,'fh',fh)
       call mpas_pool_get_array(diag_physics,'fm',fm)

       do j = jts,jte
       do i = its,ite
          fh(i) = fh_p(i,j)
          fm(i) = fm_p(i,j)
       enddo
       enddo
       if(config_frac_seaice) then
          call mpas_pool_get_array(sfc_input,'xice',xice)
          do j = jts,jte
          do i = its,ite
             if(xice(i).ge.xice_threshold .and. xice(i).le.1._RKIND) then
                fh(i) = fh_p(i,j)*xice(i) + (1._RKIND-xice(i))*fh_sea(i,j)
                fm(i) = fm_p(i,j)*xice(i) + (1._RKIND-xice(i))*fm_sea(i,j)
             endif
          enddo
          enddo
       endif

    case("sf_mynn")
       call mpas_pool_get_array(diag_physics,'ch',ch)

       do j = jts,jte
       do i = its,ite
          ch(i) = ch_p(i,j)
       enddo
       enddo
       if(config_frac_seaice) then
          call mpas_pool_get_array(sfc_input,'xice',xice)
          do j = jts,jte
          do i = its,ite
             if(xice(i).ge.xice_threshold .and. xice(i).le.1._RKIND) then
                ch(i) = ch_p(i,j)*xice(i) + (1._RKIND-xice(i))*ch_sea(i,j)
             endif
          enddo
          enddo
       endif

    case default

 end select sfclayer_select

 end subroutine sfclayer_to_MPAS

!=================================================================================================================
 subroutine init_sfclayer
!=================================================================================================================

!local variables:
 logical, parameter:: allowed_to_read = .false. !actually not used in subroutine sfclayinit.

!-----------------------------------------------------------------------------------------------------------------

 sfclayer_select: select case (trim(sfclayer_scheme))

    case("sf_monin_obukhov")
       call sfclayinit(allowed_to_read)

    case("sf_mynn")
       call mynn_sf_init_driver(allowed_to_read)

    case default

 end select sfclayer_select


 end subroutine init_sfclayer

!=================================================================================================================
 subroutine driver_sfclayer(itimestep,configs,mesh,diag_physics,sfc_input,its,ite)
!=================================================================================================================

!input and inout arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: sfc_input

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 logical,pointer:: config_do_restart,config_frac_seaice
 real(kind=RKIND),dimension(:),pointer:: areaCell

!local variables:
 integer:: initflag
 real(kind=RKIND):: dx

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine driver_sfclayer:'

 call mpas_pool_get_config(configs,'config_do_restart' ,config_do_restart )
 call mpas_pool_get_config(configs,'config_frac_seaice',config_frac_seaice)

 call mpas_pool_get_array(mesh,'areaCell',areaCell)

!copy all MPAS arrays to rectanguler grid:
 call sfclayer_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite)

 dx = sqrt(maxval(areaCell))

 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0

 sfclayer_select: select case (trim(sfclayer_scheme))

    case("sf_monin_obukhov")
       call sfclay( &
                   p3d      = pres_hyd_p , psfc     = psfc_p     , t3d      = t_p        , &
                   u3d      = u_p        , v3d      = v_p        , qv3d     = qv_p       , &
                   dz8w     = dz_p       , cp       = cp         , g        = gravity    , &
                   rovcp    = rcp        , R        = R_d        , xlv      = xlv        , & 
                   chs      = chs_p      , chs2     = chs2_p     , cqs2     = cqs2_p     , &
                   cpm      = cpm_p      , znt      = znt_p      , ust      = ust_p      , &
                   pblh     = hpbl_p     , mavail   = mavail_p   , zol      = zol_p      , &
                   mol      = mol_p      , regime   = regime_p   , psim     = psim_p     , &
                   psih     = psih_p     , fm       = fm_p       , fh       = fh_p       , &
                   xland    = xland_p    , hfx      = hfx_p      , qfx      = qfx_p      , &
                   lh       = lh_p       , tsk      = tsk_p      , flhc     = flhc_p     , &
                   flqc     = flqc_p     , qgh      = qgh_p      , qsfc     = qsfc_p     , &
                   rmol     = rmol_p     , u10      = u10_p      , v10      = v10_p      , &
                   th2      = th2m_p     , t2       = t2m_p      , q2       = q2_p       , &
                   gz1oz0   = gz1oz0_p   , wspd     = wspd_p     , br       = br_p       , &
                   isfflx   = isfflx     , dx       = dx         , svp1     = svp1       , &
                   svp2     = svp2       , svp3     = svp3       , svpt0    = svpt0      , &
                   ep1      = ep_1       , ep2      = ep_2       , karman   = karman     , &
                   eomeg    = eomeg      , stbolt   = stbolt     , P1000mb  = P0         , &
                   dxCell   = dx_p       , ustm     = ustm_p     , ck       = ck_p       , & 
                   cka      = cka_p      , cd       = cd_p       , cda      = cda_p      , & 
                   isftcflx = isftcflx   , iz0tlnd  = iz0tlnd    ,                         &
                   scm_force_flux = scm_force_flux               ,                         & 
                   ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
                   ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme , &
                   its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                  )
        
       if(config_frac_seaice) then
          call sfclay( &
                      p3d      = pres_hyd_p , psfc     = psfc_p     , t3d      = t_p        , &
                      u3d      = u_p        , v3d      = v_p        , qv3d     = qv_p       , &
                      dz8w     = dz_p       , cp       = cp         , g        = gravity    , &
                      rovcp    = rcp        , R        = R_d        , xlv      = xlv        , &
                      chs      = chs_sea    , chs2     = chs2_sea   , cqs2     = cqs2_sea   , &
                      cpm      = cpm_sea    , znt      = znt_sea    , ust      = ust_sea    , &
                      pblh     = hpbl_p     , mavail   = mavail_sea , zol      = zol_sea    , &
                      mol      = mol_sea    , regime   = regime_sea , psim     = psim_sea   , &
                      psih     = psih_sea   , fm       = fm_sea     , fh       = fh_sea     , &
                      xland    = xland_sea  , hfx      = hfx_sea    , qfx      = qfx_sea    , &
                      lh       = lh_sea     , tsk      = tsk_sea    , flhc     = flhc_sea   , &
                      flqc     = flqc_sea   , qgh      = qgh_sea    , qsfc     = qsfc_sea   , &
                      rmol     = rmol_sea   , u10      = u10_sea    , v10      = v10_sea    , &
                      th2      = th2m_sea    , t2      = t2m_sea    , q2       = q2_sea     , &
                      gz1oz0   = gz1oz0_sea , wspd     = wspd_sea   , br       = br_sea     , &
                      isfflx   = isfflx     , dx       = dx         , svp1     = svp1       , &
                      svp2     = svp2       , svp3     = svp3       , svpt0    = svpt0      , &
                      ep1      = ep_1       , ep2      = ep_2       , karman   = karman     , &
                      eomeg    = eomeg      , stbolt   = stbolt     , P1000mb  = P0         , &
                      dxCell   = dx_p       , ustm     = ustm_sea   , ck       = ck_sea     , &
                      cka      = cka_sea    , cd       = cd_sea     , cda      = cda_sea    , &
                      isftcflx = isftcflx   , iz0tlnd  = iz0tlnd    ,                         &
                      scm_force_flux = scm_force_flux               ,                         &
                      ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
                      ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme , &
                      its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                     )
       endif

    case("sf_mynn")
       call sfclay_mynn( &
                   p3d      = pres_hyd_p , pi3d     = pi_p       , psfcpa    = psfc_p     , &
                   th3d     = th_p       , t3d      = t_p        , u3d       = u_p        , &
                   v3d      = v_p        , qv3d     = qv_p       , qc3d      = qc_p       , &
                   rho3d    = rho_p      , dz8w     = dz_p       , cp        = cp         , &
                   g        = gravity    , rovcp    = rcp        , R         = R_d        , &
                   xlv      = xlv        , chs      = chs_p      , chs2      = chs2_p     , &
                   cqs2     = cqs2_p     , cpm      = cpm_p      , znt       = znt_p      , &
                   ust      = ust_p      , pblh     = hpbl_p     , mavail    = mavail_p   , &
                   zol      = zol_p      , mol      = mol_p      , regime    = regime_p   , &
                   psim     = psim_p     , psih     = psih_p     , xland     = xland_p    , &
                   hfx      = hfx_p      , qfx      = qfx_p      , lh        = lh_p       , &
                   tsk      = tsk_p      , flhc     = flhc_p     , flqc      = flqc_p     , &
                   qgh      = qgh_p      , qsfc     = qsfc_p     , rmol      = rmol_p     , &
                   u10      = u10_p      , v10      = v10_p      , th2       = th2m_p     , &
                   t2       = t2m_p      , q2       = q2_p       , snowh     = snowh_p    , &
                   gz1oz0   = gz1oz0_p   , wspd     = wspd_p     , br        = br_p       , &
                   isfflx   = isfflx     , dx       = dx         , svp1      = svp1       , &
                   svp2     = svp2       , svp3     = svp3       , svpt0     = svpt0      , &
                   ep1      = ep_1       , ep2      = ep_2       , karman    = karman     , &
                   dxCell   = dx_p       , ustm     = ustm_p     , ck        = ck_p       , &
                   cka      = cka_p      , cd       = cd_p       , cda       = cda_p      , &
                   isftcflx = isftcflx   , iz0tlnd  = iz0tlnd    , itimestep = initflag   , &
                   ch       = ch_p       , cov      = cov_p      , tsq       = tsq_p      , &
                   qsq      = qsq_p      , sh3d     = sh3d_p     , el_pbl    = elpbl_p    , &
                   qcg      = qcg_p      , bl_mynn_cloudpdf = bl_mynn_cloudpdf            , &
                   ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde  , &
                   ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme  , &
                   its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte    &
                       )

       if(config_frac_seaice) then
          call sfclay_mynn( &
                      p3d      = pres_hyd_p , pi3d     = pi_p       , psfcpa    = psfc_p     , &
                      th3d     = th_p       , t3d      = t_p        , u3d       = u_p        , &
                      v3d      = v_p        , qv3d     = qv_p       , qc3d      = qc_p       , &
                      rho3d    = rho_p      , dz8w     = dz_p       , cp        = cp         , &
                      g        = gravity    , rovcp    = rcp        , R         = R_d        , &
                      xlv      = xlv        , chs      = chs_sea    , chs2      = chs2_sea   , &
                      cqs2     = cqs2_sea   , cpm      = cpm_sea    , znt       = znt_sea    , &
                      ust      = ust_sea    , pblh     = hpbl_p     , mavail    = mavail_sea , &
                      zol      = zol_sea    , mol      = mol_sea    , regime    = regime_sea , &
                      psim     = psim_sea   , psih     = psih_sea   , xland     = xland_sea  , &
                      hfx      = hfx_sea    , qfx      = qfx_sea    , lh        = lh_sea     , &
                      tsk      = tsk_sea    , flhc     = flhc_sea   , flqc      = flqc_sea   , &
                      qgh      = qgh_sea    , qsfc     = qsfc_sea   , rmol      = rmol_sea   , &
                      u10      = u10_sea    , v10      = v10_sea    , th2       = th2m_sea   , &
                      t2       = t2m_sea    , q2       = q2_sea     , snowh     = snowh_p    , &
                      gz1oz0   = gz1oz0_sea , wspd     = wspd_sea   , br        = br_sea     , &
                      isfflx   = isfflx     , dx       = dx         , svp1      = svp1       , &
                      svp2     = svp2       , svp3     = svp3       , svpt0     = svpt0      , &
                      ep1      = ep_1       , ep2      = ep_2       , karman    = karman     , &
                      dxCell   = dx_p       , ustm     = ustm_sea   , ck        = ck_sea     , &
                      cka      = cka_sea    , cd       = cd_sea     , cda       = cda_sea    , &
                      isftcflx = isftcflx   , iz0tlnd  = iz0tlnd    , itimestep = initflag   , &
                      ch       = ch_sea     , cov      = cov_p      , tsq       = tsq_p      , &
                      qsq      = qsq_p      , sh3d     = sh3d_p     , el_pbl    = elpbl_p    , &
                      qcg      = qcg_p      , bl_mynn_cloudpdf = bl_mynn_cloudpdf            , &
                      ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde  , &
                      ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme  , &
                      its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte    &
                          )
       endif

    case default

 end select sfclayer_select

!copy local arrays to MPAS grid:
 call sfclayer_to_MPAS(configs,sfc_input,diag_physics,its,ite)

!write(0,*) '--- end subroutine driver_sfclayer.'

 end subroutine driver_sfclayer

!=================================================================================================================
 end module mpas_atmphys_driver_sfclayer
!=================================================================================================================
