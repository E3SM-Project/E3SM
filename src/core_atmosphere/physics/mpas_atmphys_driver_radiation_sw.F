! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_radiation_sw
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_timer

 use mpas_atmphys_constants
 use mpas_atmphys_manager, only: gmt,curr_julday,julday,year
 use mpas_atmphys_camrad_init
 use mpas_atmphys_rrtmg_swinit
 use mpas_atmphys_vars
 
!wrf physics:
 use module_ra_cam
 use module_ra_rrtmg_sw

 implicit none
 private
 public:: allocate_radiation_sw,   &
          deallocate_radiation_sw, &
          driver_radiation_sw,     &
          init_radiation_sw,       &
          radconst

 integer,private:: i,j,k,kk,n

!>\brief MPAS driver for parameterization of shortwave radiation codes.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_driver_radiation_sw:
!> ------------------------------------------------
!> allocate_radiation_sw     : allocate local arrays for parameterization of sw radiation codes.
!> deallocate_radiation_sw   : deallocate local arrays for parameterization of sw radiation codes.
!> init_radiation_sw         : initialization of individual sw radiation codes.
!> driver_radiation_sw       : main driver (called from subroutine physics_driver).
!> radiation_sw_from_MPAS    : initialize local arrays.
!> radiation_sw_to_MPAS      : copy local arrays to MPAS arrays.
!> radconst                  : calculate solar declination,...
!>
!> WRF physics called from driver_radiation_sw:
!> --------------------------------------------
!>    * module_ra_cam        : CAM short wave radiation code.
!>    * module_ra_rrtmg_sw   : RRTMG short wave radiation code.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutines
!>      rrtmg_sw and camrad.
!>      Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
!>    * added structure diag in the call to subroutine init_radiation_sw and call to subroutine
!>      camradinit for initialization of variable mxaerl.
!>      Laura D. Fowler (laura@ucar.edu) / 2013-07-01.
!>    * modified the call to subroutine rrtmg_swrad to include the option of using the same ozone
!>      climatology as the one used in the CAM radiation codes.
!>      Laura D. Fowler (laura@ucar.edu) / 2013-07-17.
!>    * in call to subroutine rrtmg_swrad, replaced the variable g (that originally pointed to
!>      gravity) with gravity, for simplicity.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
!>    * modified sourcecode to use pools.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
!>    * cleaned-up the call to rrtmg_swrad after cleaning up subroutine rrtmg_swrad in module_ra_rrtmg_sw.F.
!>      Laura D. Fowler (laura@ucar.edu) / 2016-07-05.
!>    * added the cloud radii for cloud water, cloud ice, and snow calculated in the Thompson cloud microphysics
!>      scheme in the call to subroutine rrtmg_swrad.
!>      Laura D. Fowler (laura@ucar.edu) / 2016-07-07.
!>    * removed qr_p, and qg_p in the call to rrtmg_swrad since not used in the calculation of the cloud optical
!>      properties.
!>      Laura D. Fowler (laura@ucar.edu) / 2016-07-08.
!>    * in the call to rrtmg_swrad, substituted the variables qv_p, qc_p, qi_p, and qs_p with qvrad_p, qcrad_p,
!>      qirad_p, and qsrad_p initialized in subroutine cloudiness_from_MPAS.
!>      Laura D. Fowler (laura@ucar.edu) / 2016-07-09.


 contains


!=================================================================================================================
 subroutine allocate_radiation_sw(xtime_s)
!=================================================================================================================

!input arguments:
 real(kind=RKIND),intent(in):: xtime_s

!-----------------------------------------------------------------------------------------------------------------

 if(.not.allocated(f_ice)        ) allocate(f_ice(ims:ime,kms:kme,jms:jme)        )
 if(.not.allocated(f_rain)       ) allocate(f_rain(ims:ime,kms:kme,jms:jme)       )

 if(.not.allocated(xlat_p)       ) allocate(xlat_p(ims:ime,jms:jme)               )
 if(.not.allocated(xlon_p)       ) allocate(xlon_p(ims:ime,jms:jme)               )

 if(.not.allocated(sfc_albedo_p) ) allocate(sfc_albedo_p(ims:ime,jms:jme)         )
 if(.not.allocated(snow_p)       ) allocate(snow_p(ims:ime,jms:jme)               )
 if(.not.allocated(tsk_p)        ) allocate(tsk_p(ims:ime,jms:jme)                )
 if(.not.allocated(xice_p)       ) allocate(xice_p(ims:ime,jms:jme)               )
 if(.not.allocated(xland_p)      ) allocate(xland_p(ims:ime,jms:jme)              )

 if(.not.allocated(coszr_p)      ) allocate(coszr_p(ims:ime,jms:jme)              )
 if(.not.allocated(gsw_p)        ) allocate(gsw_p(ims:ime,jms:jme)                )
 if(.not.allocated(swcf_p)       ) allocate(swcf_p(ims:ime,jms:jme)               )
 if(.not.allocated(swdnb_p)      ) allocate(swdnb_p(ims:ime,jms:jme)              )
 if(.not.allocated(swdnbc_p)     ) allocate(swdnbc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swdnt_p)      ) allocate(swdnt_p(ims:ime,jms:jme)              )
 if(.not.allocated(swdntc_p)     ) allocate(swdntc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swupb_p)      ) allocate(swupb_p(ims:ime,jms:jme)              )
 if(.not.allocated(swupbc_p)     ) allocate(swupbc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swupt_p)      ) allocate(swupt_p(ims:ime,jms:jme)              )
 if(.not.allocated(swuptc_p)     ) allocate(swuptc_p(ims:ime,jms:jme)             )
 
 if(.not.allocated(rthratensw_p) ) allocate(rthratensw_p(ims:ime,kms:kme,jms:jme) )

 radiation_sw_select: select case (trim(radt_sw_scheme))

    case("rrtmg_sw")
       if(.not.allocated(recloud_p)    ) allocate(recloud_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(reice_p)      ) allocate(reice_p(ims:ime,kms:kme,jms:jme)         )
       if(.not.allocated(resnow_p)     ) allocate(resnow_p(ims:ime,kms:kme,jms:jme)        )

       if(.not.allocated(alswvisdir_p) ) allocate(alswvisdir_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswvisdif_p) ) allocate(alswvisdif_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswnirdir_p) ) allocate(alswnirdir_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswnirdif_p) ) allocate(alswnirdif_p(ims:ime,jms:jme)            )
       if(.not.allocated(swvisdir_p)   ) allocate(swvisdir_p(ims:ime,jms:jme)              )
       if(.not.allocated(swvisdif_p)   ) allocate(swvisdif_p(ims:ime,jms:jme)              )
       if(.not.allocated(swnirdir_p)   ) allocate(swnirdir_p(ims:ime,jms:jme)              )
       if(.not.allocated(swnirdif_p)   ) allocate(swnirdif_p(ims:ime,jms:jme)              )

       if(.not.allocated(swdnflx_p)    ) allocate(swdnflx_p(ims:ime,kms:kme+1,jms:jme)     )
       if(.not.allocated(swdnflxc_p)   ) allocate(swdnflxc_p(ims:ime,kms:kme+1,jms:jme)    )
       if(.not.allocated(swupflx_p)    ) allocate(swupflx_p(ims:ime,kms:kme+1,jms:jme)     )
       if(.not.allocated(swupflxc_p)   ) allocate(swupflxc_p(ims:ime,kms:kme+1,jms:jme)    )

       if(.not.allocated(pin_p)        ) allocate(pin_p(num_oznlevels)                     )
       if(.not.allocated(o3clim_p)     ) allocate(o3clim_p(ims:ime,1:num_oznlevels,jms:jme))

    case("cam_sw")
       if(.not.allocated(glw_p)        ) allocate(glw_p(ims:ime,jms:jme)                )
       if(.not.allocated(lwcf_p)       ) allocate(lwcf_p(ims:ime,jms:jme)               )
       if(.not.allocated(lwdnb_p)      ) allocate(lwdnb_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwdnbc_p)     ) allocate(lwdnbc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwdnt_p)      ) allocate(lwdnt_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwdntc_p)     ) allocate(lwdntc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwupb_p)      ) allocate(lwupb_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwupbc_p)     ) allocate(lwupbc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwupt_p)      ) allocate(lwupt_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwuptc_p)     ) allocate(lwuptc_p(ims:ime,jms:jme)             )
       if(.not.allocated(olrtoa_p)     ) allocate(olrtoa_p(ims:ime,jms:jme)             )
       if(.not.allocated(sfc_emiss_p)  ) allocate(sfc_emiss_p(ims:ime,jms:jme)          )
       if(.not.allocated(rthratenlw_p) ) allocate(rthratenlw_p(ims:ime,kms:kme,jms:jme) )

       if(.not.allocated(cemiss_p)     ) allocate(cemiss_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(taucldc_p)    ) allocate(taucldc_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(taucldi_p)    ) allocate(taucldi_p(ims:ime,kms:kme,jms:jme)    )

       if(.not.allocated(pin_p)        ) allocate(pin_p(num_oznlevels)                  )
       if(.not.allocated(ozmixm_p) ) &
          allocate(ozmixm_p(ims:ime,1:num_oznlevels,jms:jme,num_months) )
       
       if(.not.allocated(m_hybi_p)     ) allocate(m_hybi_p(num_aerlevels)               )
       if(.not.allocated(m_psn_p)      ) allocate(m_psn_p(ims:ime,jms:jme)              )
       if(.not.allocated(m_psp_p)      ) allocate(m_psp_p(ims:ime,jms:jme)              )
       if(.not.allocated(aerosolcn_p)  ) &
          allocate(aerosolcn_p(ims:ime,1:num_aerlevels,jms:jme,num_aerosols) )
       if(.not.allocated(aerosolcp_p)  ) &
          allocate(aerosolcp_p(ims:ime,1:num_aerlevels,jms:jme,num_aerosols) )

       !allocate these arrays on the first time step, only:
       if(xtime_s .lt. 1.e-12) then

          if(.not.allocated(emstot_p) ) allocate(emstot_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(abstot_p) ) &
             allocate(abstot_p(ims:ime,kms:kme,cam_abs_dim2,jms:jme) )
          if(.not.allocated(absnxt_p) ) &
             allocate(absnxt_p(ims:ime,kms:kme,cam_abs_dim1,jms:jme) )

       endif

    case default

 end select radiation_sw_select

 end subroutine allocate_radiation_sw

!=================================================================================================================
 subroutine deallocate_radiation_sw
!=================================================================================================================

 if(allocated(f_ice)        ) deallocate(f_ice        )
 if(allocated(f_rain)       ) deallocate(f_rain       )
 if(allocated(xlat_p)       ) deallocate(xlat_p       )
 if(allocated(xlon_p)       ) deallocate(xlon_p       )
 if(allocated(sfc_albedo_p) ) deallocate(sfc_albedo_p )
 if(allocated(snow_p)       ) deallocate(snow_p       )
 if(allocated(tsk_p)        ) deallocate(tsk_p        )
 if(allocated(xice_p)       ) deallocate(xice_p       )
 if(allocated(xland_p)      ) deallocate(xland_p      )
 if(allocated(coszr_p)      ) deallocate(coszr_p      )
 if(allocated(gsw_p)        ) deallocate(gsw_p        )
 if(allocated(swcf_p)       ) deallocate(swcf_p       )
 if(allocated(swdnb_p)      ) deallocate(swdnb_p      )
 if(allocated(swdnbc_p)     ) deallocate(swdnbc_p     )
 if(allocated(swdnt_p)      ) deallocate(swdnt_p      )
 if(allocated(swdntc_p)     ) deallocate(swdntc_p     )
 if(allocated(swupb_p)      ) deallocate(swupb_p      )
 if(allocated(swupbc_p)     ) deallocate(swupbc_p     )
 if(allocated(swupt_p)      ) deallocate(swupt_p      )
 if(allocated(swuptc_p)     ) deallocate(swuptc_p     )
 
 if(allocated(rthratensw_p) ) deallocate(rthratensw_p )

 radiation_sw_select: select case (trim(radt_sw_scheme))

    case("rrtmg_sw")
       if(allocated(recloud_p)    ) deallocate(recloud_p    )
       if(allocated(reice_p)      ) deallocate(reice_p      )
       if(allocated(resnow_p)     ) deallocate(resnow_p     )

       if(allocated(alswvisdir_p) ) deallocate(alswvisdir_p )
       if(allocated(alswvisdif_p) ) deallocate(alswvisdif_p )
       if(allocated(alswnirdir_p) ) deallocate(alswnirdir_p )
       if(allocated(alswnirdif_p) ) deallocate(alswnirdif_p )

       if(allocated(swdnflx_p)    ) deallocate(swdnflx_p    )
       if(allocated(swdnflxc_p)   ) deallocate(swdnflxc_p   )
       if(allocated(swupflx_p)    ) deallocate(swupflx_p    )
       if(allocated(swupflxc_p)   ) deallocate(swupflxc_p   )

       if(allocated(pin_p)        ) deallocate(pin_p        )
       if(allocated(o3clim_p)     ) deallocate(o3clim_p     )

    case("cam_sw")
       if(allocated(pin_p)        ) deallocate(pin_p        )
       if(allocated(m_hybi_p)     ) deallocate(m_hybi_p     )

       if(allocated(xlat_p)       ) deallocate(xlat_p       )
       if(allocated(xlon_p)       ) deallocate(xlon_p       )

       if(allocated(glw_p)        ) deallocate(glw_p        )
       if(allocated(lwcf_p)       ) deallocate(lwcf_p       )
       if(allocated(lwdnb_p)      ) deallocate(lwdnb_p      )
       if(allocated(lwdnbc_p)     ) deallocate(lwdnbc_p     )
       if(allocated(lwdnt_p)      ) deallocate(lwdnt_p      )
       if(allocated(lwdntc_p)     ) deallocate(lwdntc_p     )
       if(allocated(lwupb_p)      ) deallocate(lwupb_p      )
       if(allocated(lwupbc_p)     ) deallocate(lwupbc_p     )
       if(allocated(lwupt_p)      ) deallocate(lwupt_p      )
       if(allocated(lwuptc_p)     ) deallocate(lwuptc_p     )
       if(allocated(olrtoa_p)     ) deallocate(olrtoa_p     )
       if(allocated(sfc_emiss_p)  ) deallocate(sfc_emiss_p  )
       if(allocated(rthratenlw_p) ) deallocate(rthratenlw_p )

       if(allocated(cemiss_p)     ) deallocate(cemiss_p     )
       if(allocated(ozmixm_p)     ) deallocate(ozmixm_p     )
       if(allocated(taucldc_p)    ) deallocate(taucldc_p    )
       if(allocated(taucldi_p)    ) deallocate(taucldi_p    )

       if(allocated(m_psn_p)      ) deallocate(m_psn_p      )
       if(allocated(m_psp_p)      ) deallocate(m_psp_p      )
       if(allocated(aerosolcn_p)  ) deallocate(aerosolcn_p  )
       if(allocated(aerosolcp_p)  ) deallocate(aerosolcp_p  )

    case default

 end select radiation_sw_select

 end subroutine deallocate_radiation_sw

!=================================================================================================================
 subroutine radiation_sw_from_MPAS(configs,mesh,state,time_lev,diag_physics,atm_input, &
                                   sfc_input,xtime_s)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: atm_input
 type(mpas_pool_type),intent(in):: sfc_input

 integer,intent(in):: time_lev
 real(kind=RKIND),intent(in):: xtime_s

!local pointers:
 logical,pointer:: config_o3climatology

 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: skintemp,snow,xice,xland
 real(kind=RKIND),dimension(:),pointer    :: m_ps,pin
 real(kind=RKIND),dimension(:),pointer    :: sfc_albedo,sfc_emiss
 real(kind=RKIND),dimension(:,:),pointer  :: cldfrac,m_hybi,o3clim
 real(kind=RKIND),dimension(:,:),pointer  :: re_cloud,re_ice,re_snow
 real(kind=RKIND),dimension(:,:,:),pointer:: aerosols,ozmixm

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_o3climatology',config_o3climatology)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'m_hybi' ,m_hybi )

 call mpas_pool_get_array(state,'m_ps'    ,m_ps    ,time_lev)
 call mpas_pool_get_array(state,'aerosols',aerosols,time_lev)

 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input,'xice'    ,xice    )
 call mpas_pool_get_array(sfc_input,'xland'   ,xland   )

 call mpas_pool_get_array(atm_input,'pin'     ,pin     )
 call mpas_pool_get_array(atm_input,'ozmixm'  ,ozmixm  )

 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'cldfrac'   ,cldfrac   )
 call mpas_pool_get_array(diag_physics,'o3clim'    ,o3clim    )

 do j = jts,jte
 do i = its,ite
    xlat_p(i,j)       = latCell(i) / degrad
    xlon_p(i,j)       = lonCell(i) / degrad

    sfc_albedo_p(i,j) = sfc_albedo(i)
    snow_p(i,j)       = snow(i)
    tsk_p(i,j)        = skintemp(i)
    xice_p(i,j)       = xice(i)
    xland_p(i,j)      = xland(i)
 enddo
 enddo
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac_p(i,k,j) = cldfrac(k,i)
 enddo
 enddo
 enddo

!initialization:
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    f_ice(i,k,j)  = 0.0_RKIND
    f_rain(i,k,j) = 0.0_RKIND
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    coszr_p(i,j)    = 0.0_RKIND
    gsw_p(i,j)      = 0.0_RKIND
    swcf_p(i,j)     = 0.0_RKIND
    swdnb_p(i,j)    = 0.0_RKIND
    swdnbc_p(i,j)   = 0.0_RKIND
    swdnt_p(i,j)    = 0.0_RKIND
    swdntc_p(i,j)   = 0.0_RKIND
    swupb_p(i,j)    = 0.0_RKIND
    swupbc_p(i,j)   = 0.0_RKIND
    swupt_p(i,j)    = 0.0_RKIND
    swuptc_p(i,j)   = 0.0_RKIND
 enddo

 do k = kts,kte
 do i = its,ite
    rthratensw_p(i,k,j) = 0.0_RKIND
 enddo
 enddo
 enddo

 radiation_sw_select: select case (trim(radt_sw_scheme))

    case("rrtmg_sw")
       call mpas_pool_get_array(diag_physics,'re_cloud',re_cloud)
       call mpas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
       call mpas_pool_get_array(diag_physics,'re_snow' ,re_snow )

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          recloud_p(i,k,j) = re_cloud(k,i)
          reice_p(i,k,j)   = re_ice(k,i)
          resnow_p(i,k,j)  = re_snow(k,i)
       enddo
       enddo
       enddo

       do j = jts,jte
       do k = kts,kte+2
       do i = its,ite
          swdnflx_p(i,k,j)  = 0.0_RKIND
          swdnflxc_p(i,k,j) = 0.0_RKIND
          swupflx_p(i,k,j)  = 0.0_RKIND
          swupflxc_p(i,k,j) = 0.0_RKIND
       enddo
       enddo
       enddo

       !ozone volum mixing ratio:
       if(config_o3climatology) then
          do k = 1, num_oznLevels
             pin_p(k) = pin(k)
          enddo
          do j = jts,jte
          do k = 1, num_oznLevels
             do i = its,ite
                o3clim_p(i,k,j) = o3clim(k,i)
             enddo
          enddo
          enddo
       else
          do k = 1, num_oznLevels
             pin_p(k) = 0.0_RKIND
          enddo
          do j = jts,jte
          do k = 1, num_oznLevels
             do i = its,ite
                o3clim_p(i,k,j) = 0.0_RKIND
             enddo
          enddo
          enddo
       endif

    case("cam_sw")
       do j = jts,jte
       do i = its,ite
          sfc_emiss_p(i,j)  = sfc_emiss(i)

          olrtoa_p(i,j)     = 0.0_RKIND
          glw_p(i,j)        = 0.0_RKIND
          lwcf_p(i,j)       = 0.0_RKIND
          lwdnb_p(i,j)      = 0.0_RKIND
          lwdnbc_p(i,j)     = 0.0_RKIND
          lwdnt_p(i,j)      = 0.0_RKIND
          lwdntc_p(i,j)     = 0.0_RKIND
          lwupb_p(i,j)      = 0.0_RKIND
          lwupbc_p(i,j)     = 0.0_RKIND
          lwupt_p(i,j)      = 0.0_RKIND
          lwuptc_p(i,j)     = 0.0_RKIND
       enddo
       do k = kts,kte
       do i = its,ite
          rthratenlw_p(i,k,j) = 0.0_RKIND
          cemiss_p(i,k,j)     = 0.0_RKIND
          taucldc_p(i,k,j)    = 0.0_RKIND
          taucldi_p(i,k,j)    = 0.0_RKIND
       enddo
       enddo
       enddo
       !infrared absorption:
       if(xtime_s .lt. 1.e-12) then
          do j = jts,jte
          do n = 1,cam_abs_dim1
          do k = kts,kte
          do i = its,ite
             absnxt_p(i,k,n,j) = 0.0_RKIND
          enddo
          enddo
          enddo
          do n = 1,cam_abs_dim2
          do k = kts,kte+1
          do i = its,ite
             abstot_p(i,k,n,j) = 0.0_RKIND
          enddo
          enddo
          enddo
          do k = kts,kte+1
          do i = its,ite
             emstot_p(i,k,j) = 0.0_RKIND
          enddo
          enddo
          enddo
       endif
       !ozone mixing ratio:
       do k = 1, num_oznlevels
          pin_p(k) = pin(k)
       enddo
       do n = 1, num_months
          do j = jts,jte
          do k = 1, num_oznlevels
          do i = its,ite
             ozmixm_p(i,k,j,n) = ozmixm(n,k,i)
          enddo
          enddo
          enddo
       enddo
       !aerosol mixing ratio:
       do k = 1, num_aerlevels
          m_hybi_p(k) = m_hybi(k,1)
       enddo
       do i = its,ite
       do j = jts,jte
          m_psp_p(i,j) = m_ps(i)
          m_psn_p(i,j) = m_ps(i)
       enddo
       enddo
       do n = 1,num_aerosols
       do j = jts,jte
       do k = 1, num_aerlevels
       do i = its,ite
          aerosolcp_p(i,k,j,n) = aerosols(n,k,i)
          aerosolcn_p(i,k,j,n) = aerosols(n,k,i)
       enddo
       enddo
       enddo
       enddo

    case default

 end select radiation_sw_select

!format:
 101 format(i3,2i6,12(1x,e15.8))
 102 format(i6,12(1x,e15.8))

 end subroutine radiation_sw_from_MPAS

!=================================================================================================================
 subroutine radiation_sw_to_MPAS(diag_physics,tend_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 real(kind=RKIND),dimension(:),pointer  :: coszr,gsw,swcf,swdnb,swdnbc,swdnt,swdntc, &
                                           swupb,swupbc,swupt,swuptc
!real(kind=RKIND),dimension(:,:),pointer:: swdnflx,swdnflxc,swupflx,swupflxc
 real(kind=RKIND),dimension(:,:),pointer:: rthratensw

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'coszr'     ,coszr     )
 call mpas_pool_get_array(diag_physics,'gsw'       ,gsw       )
 call mpas_pool_get_array(diag_physics,'swcf'      ,swcf      )
 call mpas_pool_get_array(diag_physics,'swdnb'     ,swdnb     )
 call mpas_pool_get_array(diag_physics,'swdnbc'    ,swdnbc    )
 call mpas_pool_get_array(diag_physics,'swdnt'     ,swdnt     )
 call mpas_pool_get_array(diag_physics,'swdntc'    ,swdntc    )
 call mpas_pool_get_array(diag_physics,'swupb'     ,swupb     )
 call mpas_pool_get_array(diag_physics,'swupbc'    , swupbc   )
 call mpas_pool_get_array(diag_physics,'swupt'     ,swupt     )
 call mpas_pool_get_array(diag_physics,'swuptc'    ,swuptc    )
!call mpas_pool_get_array(diag_physics,'swdnflx'   ,swdnflx   )
!call mpas_pool_get_array(diag_physics,'swdnflxc'  ,swdnflxc  )
!call mpas_pool_get_array(diag_physics,'swupflx'   ,swupflx   )
!call mpas_pool_get_array(diag_physics,'swupflxc'  ,swupflxc  )
 call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)

 do j = jts,jte

 do i = its,ite
    coszr(i) = coszr_p(i,j)
    gsw(i)    = gsw_p(i,j)
    swcf(i)   = swcf_p(i,j)
    swdnb(i)  = swdnb_p(i,j)
    swdnbc(i) = swdnbc_p(i,j)
    swdnt(i)  = swdnt_p(i,j)
    swdntc(i) = swdntc_p(i,j)
    swupb(i)  = swupb_p(i,j)
    swupbc(i) = swupbc_p(i,j)
    swupt(i)  = swupt_p(i,j)
    swuptc(i) = swuptc_p(i,j)
 enddo
!not needed:
!do k = kts,kte+2
!do i = its,ite
!   swdnflx(k,i)  = swdnflx_p(i,k,j)
!   swdnflxc(k,i) = swdnflxc_p(i,k,j)
!   swupflx(k,i)  = swupflx_p(i,k,j)
!   swupflxc(k,i) = swupflxc_p(i,k,j)
!enddo
!enddo

 do k = kts,kte
 do i = its,ite
    rthratensw(k,i) = rthratensw_p(i,k,j)
 enddo
 enddo

 enddo

 end subroutine radiation_sw_to_MPAS

!=================================================================================================================
 subroutine init_radiation_sw(dminfo,mesh,atm_input,diag,state,time_lev)
!=================================================================================================================

!input arguments:
 type(dm_info), intent(in):: dminfo
 type(mpas_pool_type),intent(in),optional:: mesh
 type(mpas_pool_type),intent(in),optional:: diag

 integer,intent(in),optional:: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout),optional:: atm_input
 type(mpas_pool_type),intent(inout),optional:: state

!-----------------------------------------------------------------------------------------------------------------
! write(0,*)
! write(0,*) '--- enter radiation_sw initialization:'

!call to shortwave radiation scheme:
 radiation_sw_select: select case (trim(radt_sw_scheme))

    case ("rrtmg_sw")
!       write(0,*) '    enter subroutine rrtmg_swinit:'
       call rrtmg_initsw_forMPAS(dminfo)
!       write(0,*) '    end subroutine rrtmg_swinit'

    case("cam_sw")
!       write(0,*) '    enter subroutine camradinit:'
       call camradinit(dminfo,mesh,atm_input,diag,state,time_lev)
!       write(0,*) '    end subroutine camradinit'

    case default

 end select radiation_sw_select

! write(0,*) '--- end radiation_sw initialization'

 end subroutine init_radiation_sw

!=================================================================================================================
 subroutine driver_radiation_sw(itimestep,configs,mesh,state,time_lev,diag_physics,atm_input, &
                                sfc_input,tend_physics,xtime_s)
!=================================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: time_lev
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: atm_input
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: tend_physics


!local pointers:
 logical,pointer:: config_o3climatology

!local variables:
 integer:: o3input
 real(kind=RKIND):: radt,xtime_m

!-----------------------------------------------------------------------------------------------------------------
! write(0,100) itimestep

 call mpas_pool_get_config(configs,'config_o3climatology',config_o3climatology)

!formats:
 100 format(/,' --- enter subroutine driver_radiation_sw: ',i6)
 101 format('     ITIMESTEP   =', i6   ,/, &
            '     YEAR        =', i6   ,/, &
            '     JULDAY      =', i6   ,/, &
            '     GMT         =', f16.9,/, &
            '     XTIME_M     =', f16.9,/, &
            '     CURR_JULDAY =', f16.9,/, &
            '     SOLCON      =', f16.9,/, &
            '     DECLIN      =', f16.9)

 xtime_m = xtime_s/60.

!copy MPAS arrays to local arrays:
 call radiation_sw_from_MPAS(configs,mesh,state,time_lev,diag_physics,atm_input,sfc_input,xtime_s)

!... calculates solar declination:
!call radconst(declin,solcon,julday,degrad,dpd)
 call radconst(declin,solcon,curr_julday,degrad,dpd)
! write(0,101) itimestep,year,julday,gmt,xtime_m,curr_julday,solcon,declin

!... convert the radiation time_step to minutes:
 radt = dt_radtsw/60.
    
!call to shortwave radiation scheme:
 radiation_sw_select: select case (trim(radt_sw_scheme))

    case ("rrtmg_sw")

!      write(0,*) '--- enter subroutine rrtmg_swrad:'
       o3input = 0
       if(config_o3climatology) o3input = 2

       call rrtmg_swrad( &
              p3d        = pres_hyd_p   , p8w        = pres2_hyd_p   , pi3d     = pi_p     , &
              t3d        = t_p          , t8w        = t2_p          , dz8w     = dz_p     , &
!             qv3d       = qv_p         , qc3d       = qc_p          , qi3d     = qi_p     , &
!             qs3d       = qs_p         , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              qv3d       = qvrad_p      , qc3d       = qcrad_p       , qi3d     = qirad_p  , &
              qs3d       = qsrad_p      , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              albedo     = sfc_albedo_p , xland      = xland_p       , xice     = xice_p   , &
              snow       = snow_p       , coszr      = coszr_p       , xtime    = xtime_m  , &
              gmt        = gmt          , julday     = julday        , radt     = radt     , &
              degrad     = degrad       , declin     = declin        , solcon   = solcon   , &
              xlat       = xlat_p       , xlong      = xlon_p        , icloud   = icloud   , &
              o3input    = o3input      , noznlevels = num_oznlevels , pin      = pin_p    , &
              o3clim     = o3clim_p     , gsw        = gsw_p         , swcf     = swcf_p   , &
              rthratensw = rthratensw_p , has_reqc   = has_reqc      , has_reqi = has_reqi , &
              has_reqs   = has_reqs     , re_cloud   = recloud_p     , re_ice   = reice_p  , &
              re_snow    = resnow_p     , swupt      = swupt_p       , swuptc   = swuptc_p , &
              swdnt      = swdnt_p      , swdntc     = swdntc_p      , swupb    = swupb_p  , &
              swupbc     = swupbc_p     , swdnb      = swdnb_p       , swdnbc   = swdnbc_p , &
              ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,        &
              ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,        &
              its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                       )

    case ("cam_sw")

!      write(0,*) '--- enter subroutine camrad_sw:'
       call camrad( dolw = .false. , dosw = .true. ,                                         &
                p_phy         = pres_hyd_p    , p8w           = pres2_hyd_p   ,              &
                pi_phy        = pi_p          , t_phy         = t_p           ,              &
                z             = zmid_p        , dz8w          = dz_p          ,              &            
                rthratenlw    = rthratenlw_p  , rthratensw    = rthratensw_p  ,              &
                swupt         = swupt_p       , swuptc        = swuptc_p      ,              &
                swdnt         = swdnt_p       , swdntc        = swdntc_p      ,              &
                lwupt         = lwupt_p       , lwuptc        = lwuptc_p      ,              &
                lwdnt         = lwdnt_p       , lwdntc        = lwdntc_p      ,              &
                swupb         = swupb_p       , swupbc        = swupbc_p      ,              &
                swdnb         = swdnb_p       , swdnbc        = swdnbc_p      ,              &
                lwupb         = lwupb_p       , lwupbc        = lwupbc_p      ,              &
                lwdnb         = lwdnb_p       , lwdnbc        = lwdnbc_p      ,              &
                swcf          = swcf_p        , lwcf          = lwcf_p        ,              &
                gsw           = gsw_p         , glw           = glw_p         ,              &
                olr           = olrtoa_p      , cemiss        = cemiss_p      ,              &
                taucldc       = taucldc_p     , taucldi       = taucldi_p     ,              & 
                coszr         = coszr_p       , albedo        = sfc_albedo_p  ,              & 
                emiss         = sfc_emiss_p   , tsk           = tsk_p         ,              & 
                xlat          = xlat_p        , xlong         = xlon_p        ,              &
                rho_phy       = rho_p         , qv3d          = qv_p          ,              & 
                qc3d          = qc_p          , qr3d          = qr_p          ,              &
                qi3d          = qi_p          , qs3d          = qs_p          ,              &
                qg3d          = qg_p          , f_qv          = f_qv          ,              &
                f_qc          = f_qc          , f_qr          = f_qr          ,              &
                f_qi          = f_qi          , f_qs          = f_qs          ,              &
                f_qg          = f_qg          , f_ice_phy     = f_ice         ,              &
                f_rain_phy    = f_rain        , cldfra        = cldfrac_p     ,              &
                xland         = xland_p       , xice          = xice_p        ,              &
                num_months    = num_months    , levsiz        = num_oznlevels ,              & 
                pin0          = pin_p         , ozmixm        = ozmixm_p      ,              &
                paerlev       = num_aerlevels , naer_c        = num_aerosols  ,              &
                m_psp         = m_psp_p       , m_psn         = m_psn_p       ,              &
                aerosolcp     = aerosolcp_p   , aerosolcn     = aerosolcn_p   ,              &
                m_hybi0       = m_hybi_p      , snow          = snow_p        ,              &
                cam_abs_dim1  = cam_abs_dim1  , cam_abs_dim2  = cam_abs_dim2  ,              &
                gmt           = gmt           , yr            = year          ,              &
                julday        = julday        , julian        = curr_julday   ,              &
                dt            = dt_dyn        , xtime         = xtime_m       ,              &
                declin        = declin        , solcon        = solcon        ,              &
                radt          = radt          , degrad        = degrad        ,              &
                n_cldadv      = 3             , abstot_3d     = abstot_p      ,              &
                absnxt_3d     = absnxt_p      , emstot_3d     = emstot_p      ,              &
                doabsems      = doabsems      ,                                              &
                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,      &
                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,      &
                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte        &
                  )

    case default

 end select radiation_sw_select

!copy local arrays to MPAS grid:
 call radiation_sw_to_MPAS(diag_physics,tend_physics)

! write(0,*) '--- end subroutine driver_radiation_sw'

!formats:
 200 format(i3,i6,8(1x,e15.8))

 end subroutine driver_radiation_sw

!=================================================================================================================
 subroutine radconst(declin,solcon,julian,degrad,dpd)
!=================================================================================================================

!input arguments:
!integer,intent(in):: julian
 real(kind=RKIND),intent(in):: julian
 real(kind=RKIND),intent(in):: degrad,dpd

!output arguments:
 real(kind=RKIND),intent(out):: declin,solcon

!local variables:
 real(kind=RKIND):: obecl,sinob,sxlong,arg,decdeg,djul,rjul,eccfac

!-----------------------------------------------------------------------------------------------------------------

   declin=0.
   solcon=0.

!obecl : obliquity = 23.5 degree.
        
   obecl=23.5*degrad
   sinob=sin(obecl)
        
!calculate longitude of the sun from vernal equinox:        

   if(julian.ge.80.)sxlong=dpd*(julian-80.)
   if(julian.lt.80.)sxlong=dpd*(julian+285.)
   sxlong=sxlong*degrad
   arg=sinob*sin(sxlong)
   declin=asin(arg)
   decdeg=declin/degrad

!solar constant eccentricity factor (paltridge and platt 1976)

   djul=julian*360./365.
   rjul=djul*degrad
   eccfac=1.000110+0.034221*cos(rjul)+0.001280*sin(rjul)+0.000719*  &
          cos(2*rjul)+0.000077*sin(2*rjul)
   solcon=solcon_0*eccfac

 end subroutine radconst

!=================================================================================================================
 end module mpas_atmphys_driver_radiation_sw
!=================================================================================================================
