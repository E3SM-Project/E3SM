! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_lsm
 use mpas_kind_types
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_landuse
 use mpas_atmphys_lsm_noahinit
 use mpas_atmphys_vars

!wrf physics
 use module_sf_noahdrv
 
 implicit none
 private
 public:: init_lsm,       &
          allocate_lsm,   &
          deallocate_lsm, &
          driver_lsm
 
 logical,parameter:: rdmaxalb = .false. !use snow albedo from geogrid;false use table values
 logical,parameter:: myj      = .false. !true if using Mellor-Yamada PBL scheme.
 logical,parameter:: frpcpn   = .false.
 logical,parameter:: rdlai2d  = .false.

!urban physics: MPAS does not plan to run the urban physics option.
 integer,parameter:: sf_urban_physics = 0 !activate urban canopy model (=0: no urban canopy) 

 integer,private:: i,j,k,n

!>\brief MPAS driver for parameterization of land surface processes.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_driver_lsm:
!> ---------------------------------------
!> allocate_lsm  : allocate local arrays for land surface parameterization.
!> deallocate_lsm: deallocate local arrays for land surface parameterization.
!> driver_lsm    : main driver (called from subroutine physics_driver).
!> lsm_from_MPAS : initialize local arrays.
!> lsm_to_MPAS   : copy local arrays to MPAS arrays.
!>
!> WRF physics called from driver_lsm:
!> ------------------------ ----------
!>    * module_sf_noahdrv : NOAH 4-layers land surface scheme.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutine lsm.
!>      Laura D. Fowler (birch.ucar.edu) / 2013-05-29.
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
!>    * modified sourcecode to use pools.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
!>    * moved the definition of isurban to landuse_init_forMPAS in mpas_atmphys_landuse.F.
!>      isurban is now defined as a function of the input landuse data file.
!>      Dominikus Heinzeller (IMK) / 2014-07-24.

!>
!> DOCUMENTATION:
!> ./physics_wrf/module_sf_noahdrv.F: main driver for the "NOAH" land-surface parameterization.
!> In the argument list,I added "OPTIONAL" to the declaration of a few arrays to avoid compiling
!> with the "urban physics" option. These arrays are:
!> .. num_roof_layers; num_wall_layers; num_road_layers;num_urban_layers.
!> .. ust_urb2d;frc_urb2d;utype_urb2d.
!> Laura D. Fowler (01-18-2011).


 contains


!==================================================================================================
 subroutine allocate_lsm
!==================================================================================================

!arrays for soil layer properties:
 if(.not.allocated(dzs_p)    ) allocate(dzs_p(1:num_soils)                    )
 if(.not.allocated(smcrel_p) ) allocate(smcrel_p(ims:ime,1:num_soils,jms:jme) )
 if(.not.allocated(sh2o_p)   ) allocate(sh2o_p(ims:ime,1:num_soils,jms:jme)   )
 if(.not.allocated(smois_p)  ) allocate(smois_p(ims:ime,1:num_soils,jms:jme)  )
 if(.not.allocated(tslb_p)   ) allocate(tslb_p(ims:ime,1:num_soils,jms:jme)   )

!other arrays:
 if(.not.allocated(acsnom_p)     ) allocate(acsnom_p(ims:ime,jms:jme)       )
 if(.not.allocated(acsnow_p)     ) allocate(acsnow_p(ims:ime,jms:jme)       )
 if(.not.allocated(canwat_p)     ) allocate(canwat_p(ims:ime,jms:jme)       )
 if(.not.allocated(chs_p)        ) allocate(chs_p(ims:ime,jms:jme)          )
 if(.not.allocated(chs2_p)       ) allocate(chs2_p(ims:ime,jms:jme)         )
 if(.not.allocated(chklowq_p)    ) allocate(chklowq_p(ims:ime,jms:jme)      )
 if(.not.allocated(cpm_p)        ) allocate(cpm_p(ims:ime,jms:jme)          )
 if(.not.allocated(cqs2_p)       ) allocate(cqs2_p(ims:ime,jms:jme)         )
 if(.not.allocated(isltyp_p)     ) allocate(isltyp_p(ims:ime,jms:jme)       )
 if(.not.allocated(ivgtyp_p)     ) allocate(ivgtyp_p(ims:ime,jms:jme)       )
 if(.not.allocated(glw_p)        ) allocate(glw_p(ims:ime,jms:jme)          )
 if(.not.allocated(grdflx_p)     ) allocate(grdflx_p(ims:ime,jms:jme)       )
 if(.not.allocated(gsw_p)        ) allocate(gsw_p(ims:ime,jms:jme)          )
 if(.not.allocated(hfx_p)        ) allocate(hfx_p(ims:ime,jms:jme)          )
 if(.not.allocated(lai_p)        ) allocate(lai_p(ims:ime,jms:jme)          )
 if(.not.allocated(lh_p)         ) allocate(lh_p(ims:ime,jms:jme)           ) 
 if(.not.allocated(noahres_p)    ) allocate(noahres_p(ims:ime,jms:jme)      )
 if(.not.allocated(potevp_p)     ) allocate(potevp_p(ims:ime,jms:jme)       )
 if(.not.allocated(qfx_p)        ) allocate(qfx_p(ims:ime,jms:jme)          )
 if(.not.allocated(qgh_p)        ) allocate(qgh_p(ims:ime,jms:jme)          )
 if(.not.allocated(qsfc_p)       ) allocate(qsfc_p(ims:ime,jms:jme)         )
 if(.not.allocated(qz0_p)        ) allocate(qz0_p(ims:ime,jms:jme)          )
 if(.not.allocated(rainbl_p)     ) allocate(rainbl_p(ims:ime,jms:jme)       )
 if(.not.allocated(br_p)         ) allocate(br_p(ims:ime,jms:jme)          )
 if(.not.allocated(sfc_albbck_p) ) allocate(sfc_albbck_p(ims:ime,jms:jme)   )
 if(.not.allocated(sfc_albedo_p) ) allocate(sfc_albedo_p(ims:ime,jms:jme)   )
 if(.not.allocated(sfc_emibck_p) ) allocate(sfc_emibck_p(ims:ime,jms:jme)   )
 if(.not.allocated(sfc_emiss_p)  ) allocate(sfc_emiss_p(ims:ime,jms:jme)    )
 if(.not.allocated(sfcrunoff_p)  ) allocate(sfcrunoff_p(ims:ime,jms:jme)    )
 if(.not.allocated(shdmin_p)     ) allocate(shdmin_p(ims:ime,jms:jme)       )
 if(.not.allocated(shdmax_p)     ) allocate(shdmax_p(ims:ime,jms:jme)       )
 if(.not.allocated(smstav_p)     ) allocate(smstav_p(ims:ime,jms:jme)       )
 if(.not.allocated(smstot_p)     ) allocate(smstot_p(ims:ime,jms:jme)       )
 if(.not.allocated(snoalb_p)     ) allocate(snoalb_p(ims:ime,jms:jme)       )
 if(.not.allocated(snotime_p)    ) allocate(snotime_p(ims:ime,jms:jme)      )
 if(.not.allocated(snopcx_p)     ) allocate(snopcx_p(ims:ime,jms:jme)       )
 if(.not.allocated(snow_p)       ) allocate(snow_p(ims:ime,jms:jme)         )
 if(.not.allocated(snowc_p)      ) allocate(snowc_p(ims:ime,jms:jme)        )
 if(.not.allocated(snowh_p)      ) allocate(snowh_p(ims:ime,jms:jme)        )
 if(.not.allocated(sr_p)         ) allocate(sr_p(ims:ime,jms:jme)           )
 if(.not.allocated(swdown_p)     ) allocate(swdown_p(ims:ime,jms:jme)       )
 if(.not.allocated(tmn_p)        ) allocate(tmn_p(ims:ime,jms:jme)          )
 if(.not.allocated(tsk_p)        ) allocate(tsk_p(ims:ime,jms:jme)          )
 if(.not.allocated(udrunoff_p)   ) allocate(udrunoff_p(ims:ime,jms:jme)     )
 if(.not.allocated(vegfra_p)     ) allocate(vegfra_p(ims:ime,jms:jme)       )
 if(.not.allocated(xice_p)       ) allocate(xice_p(ims:ime,jms:jme)         )
 if(.not.allocated(xland_p)      ) allocate(xland_p(ims:ime,jms:jme)        )
 if(.not.allocated(z0_p)         ) allocate(z0_p(ims:ime,jms:jme)           )
 if(.not.allocated(znt_p)        ) allocate(znt_p(ims:ime,jms:jme)          )

 end subroutine allocate_lsm

!==================================================================================================
 subroutine deallocate_lsm
!==================================================================================================

!arrays for soil layer properties:
 if(allocated(dzs_p)    ) deallocate(dzs_p    )
 if(allocated(smcrel_p) ) deallocate(smcrel_p )
 if(allocated(sh2o_p)   ) deallocate(sh2o_p   )
 if(allocated(smois_p)  ) deallocate(smois_p  )
 if(allocated(tslb_p)   ) deallocate(tslb_p   )

!other arrays:
 if(allocated(acsnom_p)     ) deallocate(acsnom_p     )
 if(allocated(acsnow_p)     ) deallocate(acsnow_p     )
 if(allocated(canwat_p)     ) deallocate(canwat_p     )
 if(allocated(chs_p)        ) deallocate(chs_p        )
 if(allocated(chs2_p)       ) deallocate(chs2_p       )
 if(allocated(chklowq_p)    ) deallocate(chklowq_p    )
 if(allocated(cpm_p)        ) deallocate(cpm_p        )
 if(allocated(cqs2_p)       ) deallocate(cqs2_p       )
 if(allocated(glw_p)        ) deallocate(glw_p        )
 if(allocated(grdflx_p)     ) deallocate(grdflx_p     )
 if(allocated(gsw_p)        ) deallocate(gsw_p        )
 if(allocated(hfx_p)        ) deallocate(hfx_p        )
 if(allocated(isltyp_p)     ) deallocate(isltyp_p     )
 if(allocated(ivgtyp_p)     ) deallocate(ivgtyp_p     )
 if(allocated(lai_p)        ) deallocate(lai_p        )
 if(allocated(lh_p)         ) deallocate(lh_p         )
 if(allocated(noahres_p)    ) deallocate(noahres_p    )
 if(allocated(potevp_p)     ) deallocate(potevp_p     )
 if(allocated(qfx_p)        ) deallocate(qfx_p        )
 if(allocated(qgh_p)        ) deallocate(qgh_p        )
 if(allocated(qsfc_p)       ) deallocate(qsfc_p       )
 if(allocated(qz0_p)        ) deallocate(qz0_p        )
 if(allocated(rainbl_p)     ) deallocate(rainbl_p     )
 if(allocated(br_p)         ) deallocate(br_p         )
 if(allocated(sfc_albbck_p) ) deallocate(sfc_albbck_p )
 if(allocated(sfc_albedo_p) ) deallocate(sfc_albedo_p )
 if(allocated(sfc_emibck_p) ) deallocate(sfc_emibck_p )
 if(allocated(sfc_emiss_p)  ) deallocate(sfc_emiss_p  )
 if(allocated(sfcrunoff_p)  ) deallocate(sfcrunoff_p  )
 if(allocated(shdmin_p)     ) deallocate(shdmin_p     )
 if(allocated(shdmax_p)     ) deallocate(shdmax_p     )
 if(allocated(smstav_p)     ) deallocate(smstav_p     )
 if(allocated(smstot_p)     ) deallocate(smstot_p     )
 if(allocated(snoalb_p)     ) deallocate(snoalb_p     )
 if(allocated(snotime_p)    ) deallocate(snotime_p    )
 if(allocated(snopcx_p)     ) deallocate(snopcx_p     )
 if(allocated(snow_p)       ) deallocate(snow_p       )
 if(allocated(snowc_p)      ) deallocate(snowc_p      )
 if(allocated(snowh_p)      ) deallocate(snowh_p      )
 if(allocated(sr_p)         ) deallocate(sr_p         )
 if(allocated(swdown_p)     ) deallocate(swdown_p     )
 if(allocated(tmn_p)        ) deallocate(tmn_p        )
 if(allocated(tsk_p)        ) deallocate(tsk_p        )
 if(allocated(udrunoff_p)   ) deallocate(udrunoff_p   )
 if(allocated(vegfra_p)     ) deallocate(vegfra_p     )
 if(allocated(xice_p)       ) deallocate(xice_p       )
 if(allocated(xland_p)      ) deallocate(xland_p      )
 if(allocated(z0_p)         ) deallocate(z0_p         )
 if(allocated(znt_p)        ) deallocate(znt_p        )

 end subroutine deallocate_lsm

!==================================================================================================
 subroutine lsm_from_MPAS(mesh,diag_physics,sfc_input)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 integer,dimension(:),pointer:: isltyp,ivgtyp

 real(kind=RKIND),dimension(:),pointer  :: acsnom,acsnow,canwat,chs,chs2,chklowq,cpm,cqs2,glw, &
                                           grdflx,gsw,hfx,lai,lh,noahres,potevp,qfx,qgh,qsfc,  &
                                           qz0,br,sfc_albedo,sfc_emibck,sfc_emiss,sfcrunoff,   &
                                           smstav,smstot,snotime,snopcx,sr,udrunoff,z0,znt
 real(kind=RKIND),dimension(:),pointer  :: shdmin,shdmax,snoalb,sfc_albbck,snow,snowc,snowh,tmn, &
                                           skintemp,vegfra,xice,xland
 real(kind=RKIND),dimension(:),pointer  :: raincv,rainncv
 real(kind=RKIND),dimension(:,:),pointer:: sh2o,smcrel,smois,tslb,dzs

!local variables and arrays:
 logical:: do_fill
 integer:: ip,iEdg

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'acsnom'    ,acsnom    )
 call mpas_pool_get_array(diag_physics,'acsnow'    ,acsnow    )
 call mpas_pool_get_array(diag_physics,'canwat'    ,canwat    )
 call mpas_pool_get_array(diag_physics,'chs'       ,chs       )
 call mpas_pool_get_array(diag_physics,'chs2'      ,chs2      )
 call mpas_pool_get_array(diag_physics,'chklowq'   ,chklowq   )
 call mpas_pool_get_array(diag_physics,'cpm'       ,cpm       )
 call mpas_pool_get_array(diag_physics,'cqs2'      ,cqs2      )
 call mpas_pool_get_array(diag_physics,'glw'       ,glw       )
 call mpas_pool_get_array(diag_physics,'grdflx'    ,grdflx    )
 call mpas_pool_get_array(diag_physics,'gsw'       ,gsw       )
 call mpas_pool_get_array(diag_physics,'hfx'       ,hfx       )
 call mpas_pool_get_array(diag_physics,'lai'       ,lai       )
 call mpas_pool_get_array(diag_physics,'lh'        ,lh        )
 call mpas_pool_get_array(diag_physics,'noahres'   ,noahres   )
 call mpas_pool_get_array(diag_physics,'potevp'    ,potevp    )
 call mpas_pool_get_array(diag_physics,'qfx'       ,qfx       )
 call mpas_pool_get_array(diag_physics,'qgh'       ,qgh       )
 call mpas_pool_get_array(diag_physics,'qsfc'      ,qsfc      )
 call mpas_pool_get_array(diag_physics,'qz0'       ,qz0       )
 call mpas_pool_get_array(diag_physics,'br'        ,br        )
 call mpas_pool_get_array(diag_physics,'raincv'    ,raincv    )
 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'sfc_emibck',sfc_emibck)
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'sfcrunoff' ,sfcrunoff )
 call mpas_pool_get_array(diag_physics,'smstav'    ,smstav    )
 call mpas_pool_get_array(diag_physics,'smstot'    ,smstot    )
 call mpas_pool_get_array(diag_physics,'snotime'   ,snotime   )
 call mpas_pool_get_array(diag_physics,'snopcx'    ,snopcx    )
 call mpas_pool_get_array(diag_physics,'sr'        ,sr        ) 
 call mpas_pool_get_array(diag_physics,'udrunoff'  ,udrunoff  )
 call mpas_pool_get_array(diag_physics,'z0'        ,z0        )
 call mpas_pool_get_array(diag_physics,'znt'       ,znt       )

 call mpas_pool_get_array(sfc_input,'isltyp'       ,isltyp    )
 call mpas_pool_get_array(sfc_input,'ivgtyp'       ,ivgtyp    )
 call mpas_pool_get_array(sfc_input,'shdmin'       ,shdmin    )
 call mpas_pool_get_array(sfc_input,'shdmax'       ,shdmax    )
 call mpas_pool_get_array(sfc_input,'snoalb'       ,snoalb    )
 call mpas_pool_get_array(sfc_input,'sfc_albbck'   ,sfc_albbck)
 call mpas_pool_get_array(sfc_input,'snow'         ,snow      )
 call mpas_pool_get_array(sfc_input,'snowc'        ,snowc     )
 call mpas_pool_get_array(sfc_input,'snowh'        ,snowh     )
 call mpas_pool_get_array(sfc_input,'tmn'          ,tmn       )
 call mpas_pool_get_array(sfc_input,'skintemp'     ,skintemp  )
 call mpas_pool_get_array(sfc_input,'vegfra'       ,vegfra    )
 call mpas_pool_get_array(sfc_input,'xice'         ,xice      )
 call mpas_pool_get_array(sfc_input,'xland'        ,xland     )
 call mpas_pool_get_array(sfc_input,'dzs'          ,dzs       )
 call mpas_pool_get_array(sfc_input,'sh2o'         ,sh2o      )
 call mpas_pool_get_array(sfc_input,'smcrel'       ,smcrel    )
 call mpas_pool_get_array(sfc_input,'smois'        ,smois     )
 call mpas_pool_get_array(sfc_input,'tslb'         ,tslb      )

 do n = 1,num_soils
    dzs_p(n) = maxval(dzs(n,:))
 enddo

 do j = jts,jte
 do n = 1,num_soils
 do i = its,ite
    sh2o_p(i,n,j)   = sh2o(n,i)
    smcrel_p(i,n,j) = smcrel(n,i)
    smois_p(i,n,j)  = smois(n,i)
    tslb_p(i,n,j)   = tslb(n,i)
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    acsnom_p(i,j)     = acsnom(i)
    acsnow_p(i,j)     = acsnow(i)
    canwat_p(i,j)     = canwat(i)
    chs_p(i,j)        = chs(i)
    chs2_p(i,j)       = chs2(i)
    chklowq_p(i,j)    = chklowq(i)
    cpm_p(i,j)        = cpm(i)
    cqs2_p(i,j)       = cqs2(i)
    glw_p(i,j)        = glw(i)
    grdflx_p(i,j)     = grdflx(i)
    gsw_p(i,j)        = gsw(i)
    hfx_p(i,j)        = hfx(i)
    lai_p(i,j)        = lai(i)
    lh_p(i,j)         = lh(i)
    noahres_p(i,j)    = noahres(i)
    potevp_p(i,j)     = potevp(i)
    qfx_p(i,j)        = qfx(i)
    qgh_p(i,j)        = qgh(i)
    qsfc_p(i,j)       = qsfc(i)
    qz0_p(i,j)        = qz0(i)
    br_p(i,j)         = br(i)
    sfc_albedo_p(i,j) = sfc_albedo(i)
    sfc_emibck_p(i,j) = sfc_emibck(i)
    sfc_emiss_p(i,j)  = sfc_emiss(i)
    sfcrunoff_p(i,j)  = sfcrunoff(i)
    smstav_p(i,j)     = smstav(i)
    smstot_p(i,j)     = smstot(i)
    snotime_p(i,j)    = snotime(i)
    snopcx_p(i,j)     = snopcx(i)
    sr_p(i,j)         = sr(i)
    udrunoff_p(i,j)   = udrunoff(i)
    z0_p(i,j)         = z0(i)
    znt_p(i,j)        = znt(i)

    isltyp_p(i,j)     = isltyp(i)
    ivgtyp_p(i,j)     = ivgtyp(i)
    shdmin_p(i,j)     = shdmin(i)
    shdmax_p(i,j)     = shdmax(i)
    snoalb_p(i,j)     = snoalb(i)
    sfc_albbck_p(i,j) = sfc_albbck(i)
    snow_p(i,j)       = snow(i)
    snowc_p(i,j)      = snowc(i)
    snowh_p(i,j)      = snowh(i)
    tmn_p(i,j)        = tmn(i)
    tsk_p(i,j)        = skintemp(i)
    vegfra_p(i,j)     = vegfra(i)
    xice_p(i,j)       = xice(i)
    xland_p(i,j)      = xland(i)

 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    rainbl_p(i,j) = raincv(i) + rainncv(i)
    swdown_p(i,j) = gsw(i) / (1._RKIND - sfc_albedo(i))
 enddo
 enddo

 end subroutine lsm_from_MPAS
 
!==================================================================================================
 subroutine lsm_to_MPAS(mesh,diag_physics,sfc_input)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 integer,dimension(:),pointer:: isltyp,ivgtyp

 real(kind=RKIND),dimension(:),pointer  :: acsnom,acsnow,canwat,chs,chs2,chklowq,cpm,cqs2,glw, &
                                           grdflx,gsw,hfx,lai,lh,noahres,potevp,qfx,qgh,qsfc,  &
                                           qz0,br,sfc_albedo,sfc_emibck,sfc_emiss,sfcrunoff,   &
                                           smstav,smstot,snotime,snopcx,sr,udrunoff,z0,znt
 real(kind=RKIND),dimension(:),pointer  :: shdmin,shdmax,snoalb,sfc_albbck,snow,snowc,snowh,tmn, &
                                           skintemp,vegfra,xice,xland
 real(kind=RKIND),dimension(:),pointer  :: raincv,rainncv
 real(kind=RKIND),dimension(:,:),pointer:: sh2o,smcrel,smois,tslb

!local variables and arrays:
 integer:: ip,iEdg

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'acsnom'    ,acsnom    )
 call mpas_pool_get_array(diag_physics,'acsnow'    ,acsnow    )
 call mpas_pool_get_array(diag_physics,'canwat'    ,canwat    )
 call mpas_pool_get_array(diag_physics,'chs'       ,chs       )
 call mpas_pool_get_array(diag_physics,'chs2'      ,chs2      )
 call mpas_pool_get_array(diag_physics,'chklowq'   ,chklowq   )
 call mpas_pool_get_array(diag_physics,'cpm'       ,cpm       )
 call mpas_pool_get_array(diag_physics,'cqs2'      ,cqs2      )
 call mpas_pool_get_array(diag_physics,'glw'       ,glw       )
 call mpas_pool_get_array(diag_physics,'grdflx'    ,grdflx    )
 call mpas_pool_get_array(diag_physics,'gsw'       ,gsw       )
 call mpas_pool_get_array(diag_physics,'hfx'       ,hfx       )
 call mpas_pool_get_array(diag_physics,'lai'       ,lai       )
 call mpas_pool_get_array(diag_physics,'lh'        ,lh        )
 call mpas_pool_get_array(diag_physics,'noahres'   ,noahres   )
 call mpas_pool_get_array(diag_physics,'potevp'    ,potevp    )
 call mpas_pool_get_array(diag_physics,'qfx'       ,qfx       )
 call mpas_pool_get_array(diag_physics,'qgh'       ,qgh       )
 call mpas_pool_get_array(diag_physics,'qsfc'      ,qsfc      )
 call mpas_pool_get_array(diag_physics,'qz0'       ,qz0       )
 call mpas_pool_get_array(diag_physics,'br'        ,br        )
 call mpas_pool_get_array(diag_physics,'raincv'    ,raincv    )
 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'sfc_emibck',sfc_emibck)
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'sfcrunoff' ,sfcrunoff )
 call mpas_pool_get_array(diag_physics,'smstav'    ,smstav    )
 call mpas_pool_get_array(diag_physics,'smstot'    ,smstot    )
 call mpas_pool_get_array(diag_physics,'snotime'   ,snotime   )
 call mpas_pool_get_array(diag_physics,'snopcx'    ,snopcx    )
 call mpas_pool_get_array(diag_physics,'sr'        ,sr        ) 
 call mpas_pool_get_array(diag_physics,'udrunoff'  ,udrunoff  )
 call mpas_pool_get_array(diag_physics,'z0'        ,z0        )
 call mpas_pool_get_array(diag_physics,'znt'       ,znt       )

 call mpas_pool_get_array(sfc_input,'isltyp'       ,isltyp    )
 call mpas_pool_get_array(sfc_input,'ivgtyp'       ,ivgtyp    )
 call mpas_pool_get_array(sfc_input,'shdmin'       ,shdmin    )
 call mpas_pool_get_array(sfc_input,'shdmax'       ,shdmax    )
 call mpas_pool_get_array(sfc_input,'snoalb'       ,snoalb    )
 call mpas_pool_get_array(sfc_input,'sfc_albbck'   ,sfc_albbck)
 call mpas_pool_get_array(sfc_input,'snow'         ,snow      )
 call mpas_pool_get_array(sfc_input,'snowc'        ,snowc     )
 call mpas_pool_get_array(sfc_input,'snowh'        ,snowh     )
 call mpas_pool_get_array(sfc_input,'tmn'          ,tmn       )
 call mpas_pool_get_array(sfc_input,'skintemp'     ,skintemp  )
 call mpas_pool_get_array(sfc_input,'vegfra'       ,vegfra    )
 call mpas_pool_get_array(sfc_input,'xice'         ,xice      )
 call mpas_pool_get_array(sfc_input,'xland'        ,xland     )
 call mpas_pool_get_array(sfc_input,'sh2o'         ,sh2o      )
 call mpas_pool_get_array(sfc_input,'smcrel'       ,smcrel    )
 call mpas_pool_get_array(sfc_input,'smois'        ,smois     )
 call mpas_pool_get_array(sfc_input,'tslb'         ,tslb      )

 do j = jts,jte
 do n = 1,num_soils
 do i = its,ite
    sh2o(n,i)   = sh2o_p(i,n,j)
    smcrel(n,i) = smcrel_p(i,n,j)
    smois(n,i)  = smois_p(i,n,j)
    tslb(n,i)   = tslb_p(i,n,j)
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    acsnom(i)     = acsnom_p(i,j)
    acsnow(i)     = acsnow_p(i,j)
    canwat(i)     = canwat_p(i,j)
    chs(i)        = chs_p(i,j)
    chs2(i)       = chs2_p(i,j)
    chklowq(i)    = chklowq_p(i,j)
    cpm(i)        = cpm_p(i,j)
    cqs2(i)       = cqs2_p(i,j)
    glw(i)        = glw_p(i,j)
    grdflx(i)     = grdflx_p(i,j)
    gsw(i)        = gsw_p(i,j)
    hfx(i)        = hfx_p(i,j)
    lai(i)        = lai_p(i,j)
    lh(i)         = lh_p(i,j)
    noahres(i)    = noahres_p(i,j)
    potevp(i)     = potevp_p(i,j)
    qfx(i)        = qfx_p(i,j)
    qgh(i)        = qgh_p(i,j)
    qsfc(i)       = qsfc_p(i,j)
    qz0(i)        = qz0_p(i,j)
    br(i)         = br_p(i,j)
    sfc_albedo(i) = sfc_albedo_p(i,j)
    sfc_emibck(i) = sfc_emibck_p(i,j)
    sfc_emiss(i)  = sfc_emiss_p(i,j)
    sfcrunoff(i)  = sfcrunoff_p(i,j)
    smstav(i)     = smstav_p(i,j)
    smstot(i)     = smstot_p(i,j)
    snotime(i)    = snotime_p(i,j)
    snopcx(i)     = snopcx_p(i,j)
    sr(i)         = sr_p(i,j)
    udrunoff(i)   = udrunoff_p(i,j)
    z0(i)         = z0_p(i,j)
    znt(i)        = znt_p(i,j)

    isltyp(i)     = isltyp_p(i,j)
    ivgtyp(i)     = ivgtyp_p(i,j)
    shdmin(i)     = shdmin_p(i,j)
    shdmax(i)     = shdmax_p(i,j)
    snoalb(i)     = snoalb_p(i,j)
    sfc_albbck(i) = sfc_albbck_p(i,j)
    snow(i)       = snow_p(i,j)
    snowc(i)      = snowc_p(i,j)
    snowh(i)      = snowh_p(i,j)
    skintemp(i)   = tsk_p(i,j)
    tmn(i)        = tmn_p(i,j)
    vegfra(i)     = vegfra_p(i,j)
    xice(i)       = xice_p(i,j)
    xland(i)      = xland_p(i,j)
 enddo
 enddo

 end subroutine lsm_to_MPAS
 
!==================================================================================================
 subroutine init_lsm(dminfo,mesh,configs,diag_physics,sfc_input)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type):: mesh
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input

!--------------------------------------------------------------------------------------------------
! write(0,*)
! write(0,*) '--- enter land surface model initialization:'

 lsm_select: select case (trim(lsm_scheme))

    case ("noah")
!       write(0,*) '    enter subroutine noah_init_forMPAS:'
       call noah_init_forMPAS(dminfo,mesh,configs,diag_physics,sfc_input)
!       write(0,*) '    end subroutine noah_init_forMPAS'
    
    case default
 
 end select lsm_select

! write(0,*) '--- end land surface model initialization'

 end subroutine init_lsm

!==================================================================================================
 subroutine driver_lsm(itimestep,configs,mesh,diag_physics,sfc_input)
!==================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 logical,pointer:: config_sfc_albedo
 character(len=StrKIND),pointer:: mminlu

!---------------------------------------------------------------------------------------------
! write(0,*)
! write(0,*) '--- enter subroutine driver_lsm:'
! write(0,*) '--- isice   = ',isice
! write(0,*) '--- iswater = ', iswater
! write(0,*) '--- isurban = ', isurban

 call mpas_pool_get_config(configs,'config_sfc_albedo' ,config_sfc_albedo )
 call mpas_pool_get_array(sfc_input,'mminlu',mminlu)

!formats:
 101 format(2i6,8(1x,e15.8))
 102 format(3i6,8(1x,e15.8))

!copy MPAS arrays to local arrays:
 call lsm_from_MPAS(mesh,diag_physics,sfc_input)
! write(0,*) '--- end lsm_from_MPAS'

!call to land-surface scheme:
 lsm_select: select case (trim(lsm_scheme))

    case("noah")
       call lsm( &
                dz8w      = dz_p      , p8w3d     = pres2_hyd_p  , t3d       = t_p          , &  
                qv3d      = qv_p      , xland     = xland_p      , xice      = xice_p       , &
                ivgtyp    = ivgtyp_p  , isltyp    = isltyp_p     , tmn       = tmn_p        , &
                vegfra    = vegfra_p  , shdmin    = shdmin_p     , shdmax    = shdmax_p     , &
                snoalb    = snoalb_p  , glw       = glw_p        , gsw       = gsw_p        , &
                swdown    = swdown_p  , rainbl    = rainbl_p     , embck     = sfc_emibck_p , & 
                sr        = sr_p      , qgh       = qgh_p        , cpm       = cpm_p        , &
                qz0       = qz0_p     , tsk       = tsk_p        , hfx       = hfx_p        , & 
                qfx       = qfx_p     , lh        = lh_p         , grdflx    = grdflx_p     , &
                qsfc      = qsfc_p    , cqs2      = cqs2_p       , chs       = chs_p        , &
                chs2      = chs2_p    , snow      = snow_p       , snowc     = snowc_p      , &
                snowh     = snowh_p   , canwat    = canwat_p     , smstav    = smstav_p     , &
                smstot    = smstot_p  , sfcrunoff = sfcrunoff_p  , udrunoff  = udrunoff_p   , &               
                acsnom    = acsnom_p  , acsnow    = acsnow_p     , snotime   = snotime_p    , &
                snopcx    = snopcx_p  , emiss     = sfc_emiss_p  , rib       = br_p         , &
                potevp    = potevp_p  , albedo    = sfc_albedo_p , albbck    = sfc_albbck_p , &
                z0        = z0_p      , znt       = znt_p        , lai       = lai_p        , &
                noahres   = noahres_p , chklowq   = chklowq_p    , sh2o      = sh2o_p       , &
                smois     = smois_p   , tslb      = tslb_p       , smcrel    = smcrel_p     , &
                dzs       = dzs_p     , isurban   = isurban      , isice     = isice        , &                
                rovcp     = rcp       , dt        = dt_pbl       , myj       = myj          , &
                itimestep = itimestep , frpcpn    = frpcpn       , rdlai2d   = rdlai2d      , &
                xice_threshold   = xice_threshold     ,                                       &
                usemonalb        = config_sfc_albedo  ,                                       &
                mminlu           = mminlu ,                                                   &
                num_soil_layers  = num_soils          ,                                       &         
                num_roof_layers  = num_soils    ,                                             &
                num_wall_layers  = num_soils    ,                                             &
                num_road_layers  = num_soils    ,                                             &
                num_urban_layers = num_soils    ,                                             &
                sf_urban_physics = sf_urban_physics   ,                                       &
                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,       &
                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,       &
                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte         &
               )

    case default

 end select lsm_select

!copy local arrays to MPAS grid:
 call lsm_to_MPAS(mesh,diag_physics,sfc_input)
! write(0,*) '--- end subroutine driver_lsm'

 end subroutine driver_lsm

!==================================================================================================
 end module mpas_atmphys_driver_lsm
!==================================================================================================
