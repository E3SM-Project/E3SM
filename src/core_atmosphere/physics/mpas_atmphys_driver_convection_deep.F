! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_convection_deep
 use mpas_kind_types
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use mpas_atmphys_vars

!wrf physics:
 use module_cu_kfeta
 use module_cu_tiedtke

 implicit none
 private
 public:: allocate_convection_deep,   &
          deallocate_convection_deep, &
          init_convection_deep,       &
          driver_convection_deep,     &
          update_convection_step1,    &
          update_convection_step2

 integer, private:: i,k,j


!>\brief MPAS driver for parameterization of convection.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_driver_convection_deep:
!> ---------------------------------------------------
!> allocate_convection_deep  : allocate local arrays for parameterization of convection.
!> deallocate_convection_deep: deallocate local arrays for parameterization of convection.
!> init_convection_deep      : initialization of individual convection scheme.
!> driver_convection_deep    : main driver (called from subroutine physics_driver).
!> convection_from_MPAS      : initialize local arrays.
!> convection_to_MPAS        : copy local arrays to MPAS arrays.
!> update_convection_step1   : updates lifetime of deep convective clouds in Kain-Fritsch scheme.
!> update_convection_step2   : updates accumulated precipitation output from convection schemes.
!>
!> WRF physics called from driver_convection_deep:
!> -----------------------------------------------
!>    * module_cu_kfeta      : Kain-Fritsch convection scheme.
!>    * module_cu_tiedtke    : Tiedtke convection scheme.

!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed the pre-processor option "do_hydrostatic_pressure" before call to the subroutines 
!>      kf_eta_cps and tiedtke.
!>    * removed call to the updated Kain-Fritsch convection scheme.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-05-29.
!>    * added the mean distance between cell centers in the call to the Kain-Fritsch convection
!>      parameterization of convection.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-08-22. 
!>    * in call to subroutine kf_eta_cps, replaced the variable g (that originally pointed to
!>      gravity) with gravity, for simplicity.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-03-21.
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-22.
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-05-15.


 contains


!==================================================================================================
 subroutine allocate_convection_deep
!==================================================================================================

 if(.not.allocated(cu_act_flag)) allocate(cu_act_flag(ims:ime,jms:jme)       )
 if(.not.allocated(rthcuten_p) ) allocate(rthcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqvcuten_p) ) allocate(rqvcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqccuten_p) ) allocate(rqccuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqicuten_p) ) allocate(rqicuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(pratec_p)   ) allocate(pratec_p(ims:ime,jms:jme)          )
 if(.not.allocated(raincv_p)   ) allocate(raincv_p(ims:ime,jms:jme)          )

 do i = its,ite
 do j = jts,jte
    pratec_p(i,j) = 0._RKIND
    raincv_p(i,j) = 0._RKIND
 enddo
 enddo

 do i = its,ite
 do k = kts,kte
 do j = jts,jte
    rthcuten_p(i,k,j) = 0._RKIND
    rqvcuten_p(i,k,j) = 0._RKIND
    rqccuten_p(i,k,j) = 0._RKIND
    rqicuten_p(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       if(.not.allocated(dx_p)       ) allocate(dx_p(ims:ime,jms:jme)              )
       if(.not.allocated(area_p)     ) allocate(area_p(ims:ime,jms:jme)            )
       if(.not.allocated(nca_p)      ) allocate(nca_p(ims:ime,jms:jme)             )
       if(.not.allocated(cubot_p)    ) allocate(cubot_p(ims:ime,jms:jme)           )
       if(.not.allocated(cutop_p)    ) allocate(cutop_p(ims:ime,jms:jme)           ) 
       if(.not.allocated(w0avg_p)    ) allocate(w0avg_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(rqrcuten_p) ) allocate(rqrcuten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rqscuten_p) ) allocate(rqscuten_p(ims:ime,kms:kme,jms:jme))

       do i = its,ite
       do j = jts,jte
          cubot_p(i,j) = DBLE(kte+1)
          cutop_p(i,j) = DBLE(kts)
       enddo
       enddo

       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rqrcuten_p(i,k,j) = 0._RKIND
          rqscuten_p(i,k,j) = 0._RKIND
       enddo
       enddo
       enddo

    case ("tiedtke")
       if(.not.allocated(qfx_p)        ) allocate(qfx_p(ims:ime,jms:jme)                )
       if(.not.allocated(xland_p)      ) allocate(xland_p(ims:ime,jms:jme)              )
       if(.not.allocated(rqvdynten_p)  ) allocate(rqvdynten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rqvdynblten_p)) allocate(rqvdynblten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rucuten_p)    ) allocate(rucuten_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(rvcuten_p)    ) allocate(rvcuten_p(ims:ime,kms:kme,jms:jme)    )

       do i = its,ite
       do j = jts,jte
          qfx_p(i,j)   = 0._RKIND
          xland_p(i,j) = 0._RKIND
       enddo
       enddo

       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rqvdynten_p(i,k,j)   = 0._RKIND
          rqvdynblten_p(i,k,j) = 0._RKIND
          rucuten_p(i,k,j)     = 0._RKIND
          rvcuten_p(i,k,j)     = 0._RKIND
       enddo
       enddo
       enddo

    case default

 end select convection_select

 end subroutine allocate_convection_deep

!==================================================================================================
 subroutine deallocate_convection_deep
!==================================================================================================

 if(allocated(cu_act_flag)) deallocate(cu_act_flag)
 if(allocated(rthcuten_p) ) deallocate(rthcuten_p )
 if(allocated(rqvcuten_p) ) deallocate(rqvcuten_p )
 if(allocated(rqccuten_p) ) deallocate(rqccuten_p )
 if(allocated(rqicuten_p) ) deallocate(rqicuten_p )
 if(allocated(pratec_p)   ) deallocate(pratec_p   )
 if(allocated(raincv_p)   ) deallocate(raincv_p   )

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       if(allocated(dx_p)         ) deallocate(dx_p         )
       if(allocated(area_p)       ) deallocate(area_p       )
       if(allocated(nca_p)        ) deallocate(nca_p        )
       if(allocated(cubot_p)      ) deallocate(cubot_p      )
       if(allocated(cutop_p)      ) deallocate(cutop_p      ) 
       if(allocated(w0avg_p)      ) deallocate(w0avg_p      )
       if(allocated(rqrcuten_p)   ) deallocate(rqrcuten_p   )
       if(allocated(rqscuten_p)   ) deallocate(rqscuten_p   )

    case ("tiedtke")
       if(allocated(qfx_p)        ) deallocate(qfx_p        )
       if(allocated(xland_p)      ) deallocate(xland_p      )
       if(allocated(rqvdynten_p)  ) deallocate(rqvdynten_p  )
       if(allocated(rqvdynblten_p)) deallocate(rqvdynblten_p)
       if(allocated(rucuten_p)    ) deallocate(rucuten_p    )
       if(allocated(rvcuten_p)    ) deallocate(rvcuten_p    )

    case default

 end select convection_select

 end subroutine deallocate_convection_deep

!==================================================================================================
 subroutine init_convection_deep(mesh,configs,diag_physics)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 logical,pointer:: config_do_restart
 integer,pointer:: nCells
 real(kind=RKIND),dimension(:),pointer:: nca

!local variables:
 integer:: iCell

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter convection_deep initialization:'

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_config(configs,'config_do_restart',config_do_restart)
 call mpas_pool_get_array(diag_physics,'nca',nca)

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       write(0,*) '    enter kain-fritsch initialization:'
       if(.not. config_do_restart) then
          do iCell = 1, nCells
             nca(iCell) = -100._RKIND
          enddo
       endif
       call kf_lutab(svp1,svp2,svp3,svpt0)
       write(0,*) '    end kain-kritsch initialization'

    case ("tiedtke")
       write(0,*) '    enter tiedtke initialization:'
!      write(mpas_err_message,'(A,A10)') &
!        'Tiedtke is being tested. Do not use right now. Thanks '
!      call physics_error_fatal(mpas_err_message)

    case default

 end select convection_select

 write(0,*) '--- end convection_deep initialization'

 end subroutine init_convection_deep

!==================================================================================================
 subroutine driver_convection_deep(itimestep,mesh,sfc_input,diag_physics,tend_physics)
!==================================================================================================

!input and output arguments:
!---------------------------
 integer,intent(in):: itimestep
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local variables and arrays:
!---------------------------
 logical:: log_convection
 integer:: icount
 real(kind=RKIND):: dx

!variables specific to Kain_Fritsch parameterization:
 logical:: warm_rain,adapt_step_flag
 integer:: ktau
 real(kind=RKIND):: curr_secs
 real(kind=RKIND):: cudt
 real(kind=RKIND):: cudtacttime

!temp:
 real(kind=RKIND):: max_rthcuten
 real(kind=RKIND):: min_rthcuten

!==================================================================================================
 write(0,*)
 write(0,*) '--- enter convection_driver: dt_cu=',dt_cu

!initialize instantaneous precipitation, and copy convective tendencies from the dynamics to
!the physics grid:
 call convection_from_MPAS(dt_dyn,mesh,sfc_input,diag_physics,tend_physics)

!... convert the convection time-step to minutes:
 cudt = dt_cu/60.

!... call to convection schemes:
 curr_secs   = -1
 cudtacttime = -1
 adapt_step_flag = .false.
 do j = jts, jte
 do i = its, ite
    cu_act_flag(i,j) = .false.
 enddo
 enddo

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       if(itimestep == 1) then
          ktau = itimestep
       else
          ktau = itimestep + 1
       endif

       call  kf_eta_cps ( &
             pcps            = pres_hyd_p      , t         = t_p        ,            &
             dt              = dt_dyn          , ktau      = ktau       ,            &
             dxCell          = dx_p            , areaCell  = area_p     ,            &
             cudt            = cudt            , curr_secs = curr_secs  ,            & 
             adapt_step_flag = adapt_step_flag , rho       = rho_p      ,            &
             raincv          = raincv_p        , pratec    = pratec_p   ,            &
             nca             = nca_p           , u         = u_p        ,            &
             v               = v_p             , th        = th_p       ,            &
             pi              = pi_p            , w         = w_p        ,            &
             dz8w            = dz_p            , w0avg     = w0avg_p    ,            &
             xlv0            = xlv0            , xlv1      = xlv1       ,            &
             xls0            = xls0            , xls1      = xls1       ,            &
             cp              = cp              , r         = r_d        ,            &
             g               = gravity         , ep1       = ep_1       ,            &
             ep2             = ep_2            , svp1      = svp1       ,            &
             svp2            = svp2            , svp3      = svp3       ,            &
             svpt0           = svpt0           , stepcu    = n_cu       ,            &
             cu_act_flag     = cu_act_flag     , warm_rain = warm_rain  ,            &
             cutop           = cutop_p         , cubot     = cubot_p    ,            &
             qv              = qv_p            , f_qv      = f_qv       ,            &
             f_qc            = f_qc            , f_qr      = f_qr       ,            &
             f_qi            = f_qi            , f_qs      = f_qs       ,            &
             rthcuten        = rthcuten_p      , rqvcuten  = rqvcuten_p ,            &
             rqccuten        = rqccuten_p      , rqrcuten  = rqrcuten_p ,            &
             rqicuten        = rqicuten_p      , rqscuten  = rqscuten_p ,            &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                        )

    case("tiedtke")
       write(0,*) '--- enter subroutine cu_tiedtke:'
       call cu_tiedtke ( &
             pcps        = pres_hyd_p  , p8w             = pres2_hyd_p     ,         &
             znu         = znu_hyd_p   , t3d             = t_p             ,         &              
             dt          = dt_dyn      , itimestep       = itimestep       ,         &
             stepcu      = n_cu        , raincv          = raincv_p        ,         &
             pratec      = pratec_p    , qfx             = qfx_p           ,         &
             u3d         = u_p         , v3d             = v_p             ,         & 
             w           = w_p         , qv3d            = qv_p            ,         &
             qc3d        = qc_p        , qi3d            = qi_p            ,         &
             pi3d        = pi_p        , rho3d           = rho_p           ,         &
             qvften      = rqvdynten_p , qvpblten        = rqvdynblten_p   ,         &
             dz8w        = dz_p        , xland           = xland_p         ,         &
             cu_act_flag = cu_act_flag , cudt            = dt_cu           ,         &
             f_qv        = f_qv        , f_qc            = f_qc            ,         &
             f_qr        = f_qr        , f_qi            = f_qi            ,         &
             f_qs        = f_qs        , rthcuten        = rthcuten_p      ,         &
             rqvcuten    = rqvcuten_p  , rqccuten        = rqccuten_p      ,         &
             rqicuten    = rqicuten_p  , rucuten         = rucuten_p       ,         &
             rvcuten     = rvcuten_p   ,                                             &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                    )

    case default

 end select convection_select

!copy instantaneous and accumulated precipitation, convective tendencies, and "other" arrays
!specific to convection parameterization back to the dynamics grid:
 call convection_to_MPAS(diag_physics,tend_physics)

 write(0,*) '--- end subroutine convection_driver'

 end subroutine driver_convection_deep

!==================================================================================================
 subroutine convection_from_MPAS(dt_dyn,mesh,sfc_input,diag_physics,tend_physics)
!==================================================================================================
!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in)   :: sfc_input
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: tend_physics
 real(kind=RKIND),intent(in):: dt_dyn

!local pointers:
 real(kind=RKIND),dimension(:),pointer  :: areaCell,dcEdge_m
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,raincv
 real(kind=RKIND),dimension(:),pointer  :: qfx,xland
 real(kind=RKIND),dimension(:,:),pointer:: w0avg
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten 
 real(kind=RKIND),dimension(:,:),pointer:: rqvblten,rqvdynten,rucuten,rvcuten 

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine convection_from_MPAS:'

 call mpas_pool_get_array(mesh,'areaCell',areaCell)

 call mpas_pool_get_array(sfc_input,'xland',xland)

 call mpas_pool_get_array(diag_physics,'dcEdge_m' ,dcEdge_m )
 call mpas_pool_get_array(diag_physics,'nca'      ,nca      )
 call mpas_pool_get_array(diag_physics,'cubot'    ,cubot    )
 call mpas_pool_get_array(diag_physics,'cutop'    ,cutop    )
 call mpas_pool_get_array(diag_physics,'cuprec'   ,cuprec   )
 call mpas_pool_get_array(diag_physics,'raincv'   ,raincv   )
 call mpas_pool_get_array(diag_physics,'w0avg'    ,w0avg    )
 call mpas_pool_get_array(diag_physics,'qfx'      ,qfx      )

 call mpas_pool_get_array(tend_physics,'rthcuten' ,rthcuten )
 call mpas_pool_get_array(tend_physics,'rqvcuten' ,rqvcuten )
 call mpas_pool_get_array(tend_physics,'rqccuten' ,rqccuten )
 call mpas_pool_get_array(tend_physics,'rqicuten' ,rqicuten )
 call mpas_pool_get_array(tend_physics,'rqrcuten' ,rqrcuten )
 call mpas_pool_get_array(tend_physics,'rqscuten' ,rqscuten )
 call mpas_pool_get_array(tend_physics,'rucuten'  ,rucuten  )
 call mpas_pool_get_array(tend_physics,'rvcuten'  ,rvcuten  )
 call mpas_pool_get_array(tend_physics,'rqvblten' ,rqvblten )
 call mpas_pool_get_array(tend_physics,'rqvdynten',rqvdynten)

 do j = jts,jte
 do i = its,ite
    raincv_p(i,j) = raincv(i)
    pratec_p(i,j) = cuprec(i)
    do k = kts,kte
       rthcuten_p(i,k,j) = rthcuten(k,i)
       rqvcuten_p(i,k,j) = rqvcuten(k,i)
       rqccuten_p(i,k,j) = rqccuten(k,i)
       rqicuten_p(i,k,j) = rqicuten(k,i)
    enddo
 enddo
 enddo

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")

       do j = jts,jte
       do i = its,ite
          !area of grid-cell:       
          area_p(i,j)  = areaCell(i)
          dx_p(i,j)    = dcEdge_m(i)
          cubot_p(i,j) = cubot(i)
          cutop_p(i,j) = cutop(i)

          do k = kts,kte
             rqrcuten_p(i,k,j) = rqrcuten(k,i)
             rqscuten_p(i,k,j) = rqscuten(k,i)
          enddo

          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          nca_p(i,j) = nca(i)

          do k = kts,kte
             w0avg_p(i,k,j) = w0avg(k,i)
          enddo
       enddo
       enddo

    case ("tiedtke")
       do j = jts,jte
       do i = its,ite
          xland_p(i,j) = xland(i)
          qfx_p(i,j)   = qfx(i)
       enddo

       do k = kts,kte
       do i = its,ite
          rqvdynblten_p(i,k,j) = rqvblten(k,i) 
          rqvdynten_p(i,k,j)   = rqvdynten(k,i)
          rucuten_p(i,k,j)     = rucuten(k,i)
          rvcuten_p(i,k,j)     = rvcuten(k,i)
       enddo
       enddo
       enddo
 
    case default

 end select convection_select
 
 end subroutine convection_from_MPAS

!==================================================================================================
 subroutine convection_to_MPAS(diag_physics,tend_physics)
!==================================================================================================
!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,raincv
 real(kind=RKIND),dimension(:,:),pointer:: w0avg
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten
 real(kind=RKIND),dimension(:,:),pointer:: rucuten,rvcuten

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag_physics,'nca'     ,nca     )
 call mpas_pool_get_array(diag_physics,'cubot'   ,cubot   )
 call mpas_pool_get_array(diag_physics,'cutop'   ,cutop   )
 call mpas_pool_get_array(diag_physics,'cuprec'  ,cuprec  )
 call mpas_pool_get_array(diag_physics,'raincv'  ,raincv  )
 call mpas_pool_get_array(diag_physics,'w0avg'   ,w0avg   )
 call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
 call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
 call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)
 call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
 call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)
 call mpas_pool_get_array(tend_physics,'rucuten' ,rucuten )
 call mpas_pool_get_array(tend_physics,'rvcuten' ,rvcuten )

 do j = jts,jte
 do i = its,ite
    raincv(i) = raincv_p(i,j)
    cuprec(i) = pratec_p(i,j)
    do k = kts, kte
       rthcuten(k,i) = rthcuten_p(i,k,j)
       rqvcuten(k,i) = rqvcuten_p(i,k,j)
       rqccuten(k,i) = rqccuten_p(i,k,j)
       rqicuten(k,i) = rqicuten_p(i,k,j)
    enddo
 enddo
 enddo

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       do j = jts,jte
       do i = its,ite
          cubot(i) = cubot_p(i,j)
          cutop(i) = cutop_p(i,j)
          nca(i)   = nca_p(i,j)
          do k = kts, kte
             w0avg(k,i)    = w0avg_p(i,k,j)
             rqrcuten(k,i) = rqrcuten_p(i,k,j)
             rqscuten(k,i) = rqscuten_p(i,k,j)
          enddo                          
       enddo
       enddo

    case ("tiedtke")
       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          rucuten(k,i) = rucuten_p(i,k,j)
          rvcuten(k,i) = rvcuten_p(i,k,j)
       enddo
       enddo
       enddo

    case default

 end select convection_select

 end subroutine convection_to_MPAS

!==================================================================================================
 subroutine update_convection_step1(mesh,diag_physics,tend_physics)
!==================================================================================================
!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,raincv
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten

!local variables and arrays:
 integer:: iCell,k

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(diag_physics,'nca'   ,nca   )
 call mpas_pool_get_array(diag_physics,'cubot' ,cubot )
 call mpas_pool_get_array(diag_physics,'cutop' ,cutop )
 call mpas_pool_get_array(diag_physics,'cuprec',cuprec)
 call mpas_pool_get_array(diag_physics,'raincv',raincv)

 call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
 call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
 call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)
 call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
 call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")

       do iCell = 1, nCellsSolve
          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          if(nca(iCell) .gt. 0.) then
             nca(iCell) = nca(iCell) - dt_dyn
             
             if(nca(iCell) .lt. 0.5*dt_dyn) then
                do k = 1, nVertLevels
                   rthcuten(k,iCell) = 0._RKIND
                   rqvcuten(k,iCell) = 0._RKIND
                   rqccuten(k,iCell) = 0._RKIND
                   rqrcuten(k,iCell) = 0._RKIND
                   rqicuten(k,iCell) = 0._RKIND
                   rqscuten(k,iCell) = 0._RKIND
                enddo
                raincv(iCell) = 0._RKIND
                cuprec(iCell) = 0._RKIND
                cubot(iCell)  = kte+1
                cutop(iCell)  = kts
             endif
          endif
       enddo

    case default

 end select convection_select

 end subroutine update_convection_step1

!==================================================================================================
 subroutine update_convection_step2(configs,mesh,diag_physics)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: i_rainc

 real(kind=RKIND),pointer:: bucket_rainc
 real(kind=RKIND),dimension(:),pointer:: cuprec,rainc

!local variables and arrays:
 integer:: iCell

!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_bucket_rainc',bucket_rainc)

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(diag_physics,'i_rainc',i_rainc)
 call mpas_pool_get_array(diag_physics,'cuprec' ,cuprec )
 call mpas_pool_get_array(diag_physics,'rainc'  ,rainc  )

!update the accumulated precipitation at the end of each dynamic time-step:
 do iCell = 1, nCellsSolve
    rainc(iCell) = rainc(iCell) + cuprec(iCell) * dt_dyn

    if(l_acrain .and. bucket_rainc.gt.0._RKIND .and. &
       rainc(iCell).gt.bucket_rainc) then
       i_rainc(iCell) = i_rainc(iCell) + 1
       rainc(iCell)   = rainc(iCell) - bucket_rainc
    endif

 enddo

 end subroutine update_convection_step2

!==================================================================================================
 end module mpas_atmphys_driver_convection_deep
!==================================================================================================
