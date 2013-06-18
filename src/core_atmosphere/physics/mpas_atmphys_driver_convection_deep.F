! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver_convection_deep
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
!>      Laura D. Fowler (birch.ucar.edu) / 2013-05-29.

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
 subroutine init_convection_deep(config_do_restart,mesh,diag_physics)
!==================================================================================================

!input arguments:
 logical,intent(in):: config_do_restart
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: iCell

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter convection_deep initialization:'

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       write(0,*) '    enter kain-fritsch initialization:'
       if(.not. config_do_restart) then
          do iCell = 1, mesh % nCells
             diag_physics % nca % array(iCell) = -100._RKIND
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
 type(mesh_type),intent(in):: mesh
 type(sfc_input_type),intent(in):: sfc_input
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

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
             pcps      = pres_hyd_p , t               = t_p             ,            &
             dt        = dt_dyn     , ktau            = ktau            ,            &
             areaCell  = area_p     , cudt            = cudt            ,            &
             curr_secs = curr_secs  , adapt_step_flag = adapt_step_flag ,            &
             rho       = rho_p      , raincv          = raincv_p        ,            &
             pratec    = pratec_p   , nca             = nca_p           ,            &
             u         = u_p        , v               = v_p             ,            &
             th        = th_p       , pi              = pi_p            ,            &
             w         = w_p        , dz8w            = dz_p            ,            &
             w0avg     = w0avg_p    , xlv0            = xlv0            ,            &
             xlv1      = xlv1       , xls0            = xls0            ,            &
             xls1      = xls1       , cp              = cp              ,            &
             r         = r_d        , g               = g               ,            &
             ep1       = ep_1       , ep2             = ep_2            ,            &
             svp1      = svp1       , svp2            = svp2            ,            &
             svp3      = svp3       , svpt0           = svpt0           ,            &
             stepcu    = n_cu       , cu_act_flag     = cu_act_flag     ,            &
             warm_rain = warm_rain  , cutop           = cutop_p         ,            &
             cubot     = cubot_p    , qv              = qv_p            ,            &
             f_qv      = f_qv       , f_qc            = f_qc            ,            &
             f_qr      = f_qr       , f_qi            = f_qi            ,            &
             f_qs      = f_qs       , rthcuten        = rthcuten_p      ,            &
             rqvcuten  = rqvcuten_p , rqccuten        = rqccuten_p      ,            &
             rqrcuten  = rqrcuten_p , rqicuten        = rqicuten_p      ,            &
             rqscuten  = rqscuten_p ,                                                &
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
 type(mesh_type),intent(in):: mesh
 type(sfc_input_type),intent(in)   :: sfc_input
 type(diag_physics_type),intent(in):: diag_physics
 type(tend_physics_type),intent(in):: tend_physics
 real(kind=RKIND),intent(in):: dt_dyn

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine convection_from_MPAS:'

 do j = jts,jte
 do i = its,ite
    raincv_p(i,j) = diag_physics % raincv % array(i)
    pratec_p(i,j) = diag_physics % cuprec % array(i)
    do k = kts,kte
       rthcuten_p(i,k,j) = tend_physics % rthcuten % array(k,i)
       rqvcuten_p(i,k,j) = tend_physics % rqvcuten % array(k,i)
       rqccuten_p(i,k,j) = tend_physics % rqccuten % array(k,i)
       rqicuten_p(i,k,j) = tend_physics % rqicuten % array(k,i)
    enddo
 enddo
 enddo

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")

       do j = jts,jte
       do i = its,ite
          !area of grid-cell:       
          area_p(i,j)  = mesh % areaCell % array(i)
          cubot_p(i,j) = diag_physics % cubot % array(i)
          cutop_p(i,j) = diag_physics % cutop % array(i)

          do k = kts,kte
             rqrcuten_p(i,k,j) = tend_physics % rqrcuten % array(k,i)
             rqscuten_p(i,k,j) = tend_physics % rqscuten % array(k,i)
          enddo

          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          nca_p(i,j) = diag_physics % nca % array(i)

          do k = kts,kte
             w0avg_p(i,k,j) = diag_physics % w0avg % array(k,i)
          enddo
       enddo
       enddo

    case ("tiedtke")
       do j = jts,jte
       do i = its,ite
          xland_p(i,j) = sfc_input % xland % array(i)
          qfx_p(i,j)   = diag_physics % qfx % array(i)
       enddo

       do k = kts,kte
       do i = its,ite
          rqvdynblten_p(i,k,j) = tend_physics % rqvblten  % array(k,i) 
          rqvdynten_p(i,k,j)   = tend_physics % rqvdynten % array(k,i)
          rucuten_p(i,k,j)     = tend_physics % rucuten % array(k,i)
          rvcuten_p(i,k,j)     = tend_physics % rvcuten % array(k,i)
       enddo
       enddo
       enddo
!      write(0,*) '--- max rqvdynblten = ',maxval(rqvdynblten_p(its:ite,kts:kte,jts:jte))
!      write(0,*) '--- min rqvdynblten = ',minval(rqvdynblten_p(its:ite,kts:kte,jts:jte))
!      write(0,*) '--- max rqvdynten   = ',maxval(rqvdynten_p(its:ite,kts:kte,jts:jte))
!      write(0,*) '--- min rqvdynten   = ',minval(rqvdynten_p(its:ite,kts:kte,jts:jte))
       
    case default

 end select convection_select
 
 end subroutine convection_from_MPAS

!==================================================================================================
 subroutine convection_to_MPAS(diag_physics,tend_physics)
!==================================================================================================
!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!--------------------------------------------------------------------------------------------------

 do j = jts,jte
 do i = its,ite
    diag_physics % raincv % array(i) = raincv_p(i,j)
    diag_physics % cuprec % array(i) = pratec_p(i,j)
    do k = kts, kte
       tend_physics % rthcuten % array(k,i) = rthcuten_p(i,k,j)
       tend_physics % rqvcuten % array(k,i) = rqvcuten_p(i,k,j)
       tend_physics % rqccuten % array(k,i) = rqccuten_p(i,k,j)
       tend_physics % rqicuten % array(k,i) = rqicuten_p(i,k,j)
    enddo
 enddo
 enddo

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       do j = jts,jte
       do i = its,ite
          diag_physics % cubot % array(i) = cubot_p(i,j)
          diag_physics % cutop % array(i) = cutop_p(i,j)
          diag_physics % nca   % array(i) = nca_p(i,j)
          do k = kts, kte
             diag_physics % w0avg % array(k,i) = w0avg_p(i,k,j)
             tend_physics % rqrcuten % array(k,i) = rqrcuten_p(i,k,j)
             tend_physics % rqscuten % array(k,i) = rqscuten_p(i,k,j)
          enddo                          
       enddo
       enddo

    case ("tiedtke")
       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          tend_physics % rucuten % array(k,i) = rucuten_p(i,k,j)
          tend_physics % rvcuten % array(k,i) = rvcuten_p(i,k,j)
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
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!local variables:
 integer:: iCell,k

!--------------------------------------------------------------------------------------------------
 
 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")

       do iCell = 1, mesh%nCellsSolve
          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          if(diag_physics % nca % array(iCell) .gt. 0.) then
             diag_physics % nca % array(iCell) = diag_physics % nca % array(iCell) - dt_dyn
             
             if(diag_physics % nca % array(iCell) .lt. 0.5*dt_dyn) then
                do k = 1, mesh%nVertLevels
                   tend_physics % rthcuten % array(k,iCell) = 0._RKIND
                   tend_physics % rqvcuten % array(k,iCell) = 0._RKIND
                   tend_physics % rqccuten % array(k,iCell) = 0._RKIND
                   tend_physics % rqrcuten % array(k,iCell) = 0._RKIND
                   tend_physics % rqicuten % array(k,iCell) = 0._RKIND
                   tend_physics % rqscuten % array(k,iCell) = 0._RKIND
                enddo
                diag_physics % raincv % array(iCell) = 0._RKIND
                diag_physics % cuprec % array(iCell) = 0._RKIND
                diag_physics % cubot  % array(iCell) = kte+1
                diag_physics  % cutop  % array(iCell) = kts
             endif
          endif
       enddo

    case default

 end select convection_select

 end subroutine update_convection_step1

!==================================================================================================
 subroutine update_convection_step2(bucket_rainc,mesh,diag_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 real(kind=RKIND),intent(in):: bucket_rainc

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: iCell

!--------------------------------------------------------------------------------------------------
 
!update the accumulated precipitation at the end of each dynamic time-step:
 do iCell = 1, mesh % nCellsSolve
    diag_physics % rainc % array(iCell) = diag_physics % rainc % array(iCell) &
                                    + diag_physics % cuprec % array(iCell) * dt_dyn

    if(l_acrain .and. bucket_rainc.gt.0._RKIND .and. &
       diag_physics%rainc%array(iCell).gt.bucket_rainc) then
       diag_physics % i_rainc % array(iCell) = diag_physics % i_rainc % array(iCell) + 1
       diag_physics % rainc % array(iCell) = diag_physics % rainc % array(iCell) &
                                           - bucket_rainc
    endif

 enddo

 end subroutine update_convection_step2

!==================================================================================================
 end module mpas_atmphys_driver_convection_deep
!==================================================================================================
