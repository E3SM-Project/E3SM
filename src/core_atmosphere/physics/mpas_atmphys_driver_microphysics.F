! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_microphysics
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 
 use mpas_atmphys_constants
 use mpas_atmphys_init_microphysics
 use mpas_atmphys_interface
 use mpas_atmphys_vars

!wrf physics:
 use module_mp_kessler
 use module_mp_thompson
 use module_mp_wsm6
 use module_mp_wsm6_wrf381

 implicit none
 private
 public:: allocate_microphysics,   &
          deallocate_microphysics, &
          microphysics_driver  ,   &
          microphysics_init


!MPAS driver for parameterization of cloud microphysics processes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_microphysics:
! ------------------------------------------------
! allocate_microphysics  : allocate local arrays for parameterization of cloud microphysics.
! deallocate_microphysics: deallocate local arrays for parameterization of cloud microphysics.
! microphysics_init      : initialization of individual cloud microphysics schemes.
! microphysics_driver    : main driver (called from mpas_atm_time_integration).
! precip_from_MPAS       : initialize timestep local arrays for precipitation.
! precip_to_MPAS         : copy local arrays to MPAS arrays.
! compute_radar_reflectivity: compute radar reflectivities.
!
! WRF physics called from microphysics_driver:
! --------------------------------------------
!    * module_mp_kessler : Kessler cloud microphysics.
!    * module_mp_wsm6    : WSM6 cloud microphysics.
!
! comments:
! unlike all the other physics parameterizations, parameterizations of cloud microphysics schemes
! are called at the bottom of subroutine atm_srk3 in module atm_time_integration.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed call to the Thompson cloud microphysics scheme until scheme is updated to that in WRF revision 3.5.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added subroutine compute_relhum to calculate the relative humidity using the functions rslf and rsif from 
!   the Thompson cloud microphysics scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-12. 
! * removed the argument tend from the call to microphysics_from_MPAS (not needed).
!   Laura D. Fowler (laura@ucar.edu) / 2013-11-07.
! * in call to subroutine wsm6, replaced the variable g (that originally pointed to gravity) with gravity, 
!   for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * moved the variable relhum from the diag_physics to the diag pool. Changed the argument
!   list for the subroutine compute_relhum accordingly.
!   Laura D. Fowler (laura@ucar.edu) / 2015-04-22.
! * added parameterization of the Thompson cloud microphysics from WRF version 3.8.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-28.
! * in subroutine compute_relhum, multiply relhum by 100. so that it has the same unit as in the initial
!   conditions.
!   Laura D. Fowler (laura@ucar.edu) / 2016-06-20.
! * added parameterization of the WSM6 cloud microphysics from WRF version 3.8.1. To initialize WSM6 as in its
!   original version, set the hail_option to 0.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-19.

!--- initialization option for WSM6 from WRF version 3.8.1. this option could also be set as a namelist parameter.
 integer,parameter:: hail_opt = 0


 contains


!=================================================================================================================
 subroutine allocate_microphysics
!=================================================================================================================

!sounding variables:
 if(.not.allocated(rho_p)    ) allocate(rho_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(th_p)     ) allocate(th_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pi_p)     ) allocate(pi_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pres_p)   ) allocate(pres_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(z_p)      ) allocate(z_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(dz_p)     ) allocate(dz_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(w_p)      ) allocate(w_p(ims:ime,kms:kme,jms:jme)   )

!mass mixing ratios:
 if(.not.allocated(qv_p)     ) allocate(qv_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qc_p)     ) allocate(qc_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qr_p)     ) allocate(qr_p(ims:ime,kms:kme,jms:jme))

 !surface precipitation:
 if(.not.allocated(rainnc_p) ) allocate(rainnc_p(ims:ime,jms:jme) )
 if(.not.allocated(rainncv_p)) allocate(rainncv_p(ims:ime,jms:jme))

 microp_select: select case(microp_scheme)

    case ("mp_thompson","mp_thompson_aerosols","mp_wsm6","mp_wsm6_wrf381")
       !mass mixing ratios:
       if(.not.allocated(qi_p)        ) allocate(qi_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qs_p)        ) allocate(qs_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qg_p)        ) allocate(qg_p(ims:ime,kms:kme,jms:jme))

       !surface precipitation:
       if(.not.allocated(sr_p)        ) allocate(sr_p(ims:ime,jms:jme)        )
       if(.not.allocated(snownc_p)    ) allocate(snownc_p(ims:ime,jms:jme)    )
       if(.not.allocated(snowncv_p)   ) allocate(snowncv_p(ims:ime,jms:jme)   )
       if(.not.allocated(graupelnc_p) ) allocate(graupelnc_p(ims:ime,jms:jme) )
       if(.not.allocated(graupelncv_p)) allocate(graupelncv_p(ims:ime,jms:jme))

       !cloud water,cloud ice,and snow effective radii:
       if(.not.allocated(recloud_p) ) allocate(recloud_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(reice_p)   ) allocate(reice_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(resnow_p)  ) allocate(resnow_p(ims:ime,kms:kme,jms:jme)  )

    microp2_select: select case(microp_scheme)

       case("mp_thompson","mp_thompson_aerosols")
          !number concentrations:
          if(.not.allocated(ntc_p)) allocate(ntc_p(ims:ime,jms:jme))
          if(.not.allocated(muc_p)) allocate(muc_p(ims:ime,jms:jme))
          if(.not.allocated(ni_p) ) allocate(ni_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(nr_p) ) allocate(nr_p(ims:ime,kms:kme,jms:jme))

          if(.not.allocated(rainprod_p)) allocate(rainprod_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(evapprod_p)) allocate(evapprod_p(ims:ime,kms:kme,jms:jme))

       microp3_select: select case(microp_scheme)

          case("mp_thompson_aerosols")
             !number concentrations:
             if(.not.allocated(nc_p)) allocate(nc_p(ims:ime,kms:kme,jms:jme))
             !aerosol concentrations:
             if(.not.allocated(nifa_p)  ) allocate(nifa_p(ims:ime,kms:kme,jms:jme))
             if(.not.allocated(nwfa_p)  ) allocate(nwfa_p(ims:ime,kms:kme,jms:jme))
             if(.not.allocated(nwfa2d_p)) allocate(nwfa2d_p(ims:ime,jms:jme)      )

          case default

       end select microp3_select

       case default

    end select microp2_select

    case default

 end select microp_select

 end subroutine allocate_microphysics

!=================================================================================================================
 subroutine deallocate_microphysics
!=================================================================================================================

!sounding variables:
 if(allocated(rho_p)     ) deallocate(rho_p     )
 if(allocated(th_p)      ) deallocate(th_p      )
 if(allocated(pi_p)      ) deallocate(pi_p      )
 if(allocated(pres_p)    ) deallocate(pres_p    )
 if(allocated(z_p)       ) deallocate(z_p       )
 if(allocated(dz_p)      ) deallocate(dz_p      )
 if(allocated(w_p)       ) deallocate(w_p       )

!mass mixing ratios:
 if(allocated(qv_p)      ) deallocate(qv_p      )
 if(allocated(qc_p)      ) deallocate(qc_p      )
 if(allocated(qr_p)      ) deallocate(qr_p      )

 !surface precipitation:
 if(allocated(rainnc_p)  ) deallocate(rainnc_p  )
 if(allocated(rainncv_p) ) deallocate(rainncv_p )

 microp_select: select case(microp_scheme)

    case ("mp_thompson","mp_thompson_aerosols","mp_wsm6","mp_wsm6_wrf381")
       !mass mixing ratios:
       if(allocated(qi_p)         ) deallocate(qi_p         )
       if(allocated(qs_p)         ) deallocate(qs_p         )
       if(allocated(qg_p)         ) deallocate(qg_p         )

       !surface precipitation:
       if(allocated(sr_p)         ) deallocate(sr_p         )
       if(allocated(snownc_p)     ) deallocate(snownc_p     )
       if(allocated(snowncv_p)    ) deallocate(snowncv_p    )
       if(allocated(graupelnc_p)  ) deallocate(graupelnc_p  )
       if(allocated(graupelncv_p) ) deallocate(graupelncv_p )

       !cloud water,cloud ice,and snow effective radii:
       if(.not.allocated(recloud_p) ) allocate(recloud_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(reice_p)   ) allocate(reice_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(resnow_p)  ) allocate(resnow_p(ims:ime,kms:kme,jms:jme)  )

    microp2_select: select case(microp_scheme)

       case("mp_thompson","mp_thompson_aerosols")
          !number concentrations:
          if(allocated(ntc_p)) deallocate(ntc_p)
          if(allocated(muc_p)) deallocate(muc_p)
          if(allocated(ni_p) ) deallocate(ni_p )
          if(allocated(nr_p) ) deallocate(nr_p )

          if(allocated(rainprod_p)) deallocate(rainprod_p)
          if(allocated(evapprod_p)) deallocate(evapprod_p)

       microp3_select: select case(microp_scheme)

          case("mp_thompson_aerosols")
             !number concentrations:
             if(allocated(nc_p)) deallocate(nc_p)
             !aerosol concentrations:
             if(allocated(nifa_p)  ) deallocate(nifa_p  )
             if(allocated(nwfa_p)  ) deallocate(nwfa_p  )
             if(allocated(nwfa2d_p)) deallocate(nwfa2d_p)

          case default

       end select microp3_select

       case default

    end select microp2_select

    case default

 end select microp_select

 end subroutine deallocate_microphysics

!=================================================================================================================
 subroutine microphysics_init(dminfo,mesh,state,time_lev,sfc_input,diag_physics)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: state

!-----------------------------------------------------------------------------------------------------------------
! write(0,*)
! write(0,*) '--- enter cloud microphysics initialization:'

  microp_select: select case(microp_scheme)

     case("mp_thompson")
!       write(0,*) '--- enter subroutine thompson_init:'
        call thompson_init
        call init_thompson_clouddroplets_forMPAS(mesh,sfc_input,diag_physics)
!       write(0,*) '--- end subroutine thompson'

     case("mp_thompson_aerosols")
        call thompson_init
        call init_thompson_clouddroplets_forMPAS(mesh,sfc_input,diag_physics)
        call init_thompson_aerosols_forMPAS(dminfo,mesh,state,time_lev,diag_physics)

     case("mp_wsm6")
        call wsm6init(rho_a,rho_r,rho_s,cliq,cpv,.false.)

     case("mp_wsm6_wrf381")
        write(0,*) '--- enter subroutine wsm6init_wrf381:' 
        call wsm6init_wrf381(rho_a,rho_r,rho_s,cliq,cpv,hail_opt,.false.)
        write(0,*) '--- end subroutine wsm6init_wrf381:' 

     case default
    
  end select microp_select

 !write(0,*) '--- end cloud microphysics initialization'

 end subroutine microphysics_init

!=================================================================================================================
 subroutine microphysics_driver(configs,mesh,state,time_lev,diag,diag_physics,tend,itimestep,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: time_lev
 integer,intent(in):: itimestep
 integer,intent(in):: its,ite
    
!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend

!local variables and arrays:
 logical:: log_microphysics
 integer:: i,icell,icount,istep,j,k,kk

!-----------------------------------------------------------------------------------------------------------------
! write(0,*)
! write(0,*) '--- enter subroutine microphysics_driver: itimestep=', itimestep
! write(0,*) '    dt_microp=',dt_microp
! write(0,*) '    n_microp =',n_microp
 
!initialization:
! write(0,*) '    ims= ',ims,' ime=',ime
! write(0,*) '    jms= ',jms,' jme=',jme
! write(0,*) '    kms= ',kms,' kme=',kme
! write(0,*)
! write(0,*) '    ids= ',ids,' ide=',ide
! write(0,*) '    jds= ',jds,' jde=',jde
! write(0,*) '    kds= ',kds,' kde=',kde
! write(0,*)
! write(0,*) '    its= ',its,' ite=',ite
! write(0,*) '    jts= ',jts,' jte=',jte
! write(0,*) '    kts= ',kts,' kte=',kte

!... allocation of microphysics arrays:
!$OMP MASTER
 call allocate_microphysics
!$OMP END MASTER
!$OMP BARRIER

!... initialization of precipitation related arrays:
 call precip_from_MPAS(diag_physics,mesh,its,ite)

!... initialization of soundings for non-hydrostatic dynamical cores.
 call microphysics_from_MPAS(mesh,state,time_lev,diag,diag_physics,its,ite)

!-----------------------------------------------------------------------------------------------------------------

!... call to different cloud microphysics schemes:

!-----------------------------------------------------------------------------------------------------------------

 microp_select: select case(microp_scheme)
    
    case ("mp_kessler")
       call kessler( &
                  t        = th_p      , qv    = qv_p  , qc     = qc_p     ,                &
                  qr       = qr_p      , rho   = rho_p , pii    = pi_p     ,                &
                  dt_in    = dt_microp , z     = z_p   , xlv    = xlv      ,                &
                  cp       = cp        , ep2   = ep_2  , svp1   = svp1     ,                &
                  svp2     = svp2      , svp3  = svp3  , svpt0  = svpt0    ,                &
                  rhowater = rho_w     , dz8w  = dz_p  , rainnc = rainnc_p ,                &
                  rainncv  = rainncv_p ,                                                    &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde   , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme   , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     &
                   )

       case ("mp_thompson")
          istep = 1
          do while (istep .le. n_microp)
          call mp_gt_driver( &
                  th        = th_p        , qv         = qv_p         , qc         = qc_p         , &
                  qr        = qr_p        , qi         = qi_p         , qs         = qs_p         , &
                  qg        = qg_p        , ni         = ni_p         , nr         = nr_p         , &
                  pii       = pi_p        , p          = pres_p       , dz         = dz_p         , & 
                  w         = w_p         , dt_in      = dt_microp    , itimestep  = itimestep    , &
                  rainnc    = rainnc_p    , rainncv    = rainncv_p    , snownc     = snownc_p     , &
                  snowncv   = snowncv_p   , graupelnc  = graupelnc_p  , graupelncv = graupelncv_p , &
                  sr        = sr_p        , rainprod   = rainprod_p   , evapprod   = evapprod_p   , &
                  re_cloud  = recloud_p   , re_ice     = reice_p      , re_snow    = resnow_p     , &
                  has_reqc  = has_reqc    , has_reqi   = has_reqi     , has_reqs   = has_reqs     , &
                  ntc       = ntc_p       , muc        = muc_p        ,                             &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde           , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme           , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte             &
                           )

          istep = istep + 1
          enddo

       case ("mp_thompson_aerosols")
          istep = 1
          do while (istep .le. n_microp)
          call mp_gt_driver( &
                  th        = th_p        , qv         = qv_p         , qc         = qc_p         , &
                  qr        = qr_p        , qi         = qi_p         , qs         = qs_p         , &
                  qg        = qg_p        , ni         = ni_p         , nr         = nr_p         , &
                  pii       = pi_p        , p          = pres_p       , dz         = dz_p         , & 
                  w         = w_p         , dt_in      = dt_microp    , itimestep  = itimestep    , &
                  rainnc    = rainnc_p    , rainncv    = rainncv_p    , snownc     = snownc_p     , &
                  snowncv   = snowncv_p   , graupelnc  = graupelnc_p  , graupelncv = graupelncv_p , &
                  sr        = sr_p        , rainprod   = rainprod_p   , evapprod   = evapprod_p   , &
                  re_cloud  = recloud_p   , re_ice     = reice_p      , re_snow    = resnow_p     , &
                  has_reqc  = has_reqc    , has_reqi   = has_reqi     , has_reqs   = has_reqs     , &
                  nc        = nc_p        , nifa       = nifa_p       , nwfa       = nwfa_p       , & 
                  nwfa2d    = nwfa2d_p    , ntc        = ntc_p        , muc        = muc_p        , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde           , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme           , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte             &
                           )
          write(0,*) '--- end mp_thompson_aerosols:',istep
          istep = istep + 1
          enddo

    case ("mp_wsm6")
      call wsm6( &
                  th      = th_p        , q          = qv_p         , qc      = qc_p      , &
                  qr      = qr_p        , qi         = qi_p         , qs      = qs_p      , &
                  qg      = qg_p        , den        = rho_p        , pii     = pi_p      , &
                  p       = pres_p      , delz       = dz_p         , delt    = dt_microp , &
                  g       = gravity     , cpd        = cp           , cpv     = cpv       , &
                  rd      = R_d         , rv         = R_v          , t0c     = svpt0     , &
                  ep1     = ep_1        , ep2        = ep_2         , qmin    = epsilon   , &
                  xls     = xls         , xlv0       = xlv          , xlf0    = xlf       , &
                  den0    = rho_a       , denr       = rho_w        , cliq    = cliq      , &
                  cice    = cice        , psat       = psat         , rain    = rainnc_p  , &
                  rainncv = rainncv_p   , snow       = snownc_p     , snowncv = snowncv_p , &
                  graupel = graupelnc_p , graupelncv = graupelncv_p , sr      = sr_p      , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde   , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme   , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     &
                   )

    case ("mp_wsm6_wrf381")
      write(0,*) '--- enter subroutine mp_wsm6_wrf381:'
      call wsm6_wrf381( &
                  th       = th_p        , q          = qv_p         , qc       = qc_p      , &
                  qr       = qr_p        , qi         = qi_p         , qs       = qs_p      , &
                  qg       = qg_p        , den        = rho_p        , pii      = pi_p      , &
                  p        = pres_p      , delz       = dz_p         , delt     = dt_microp , &
                  g        = gravity     , cpd        = cp           , cpv      = cpv       , &
                  rd       = R_d         , rv         = R_v          , t0c      = svpt0     , &
                  ep1      = ep_1        , ep2        = ep_2         , qmin     = epsilon   , &
                  xls      = xls         , xlv0       = xlv          , xlf0     = xlf       , &
                  den0     = rho_a       , denr       = rho_w        , cliq     = cliq      , &
                  cice     = cice        , psat       = psat         , rain     = rainnc_p  , &
                  rainncv  = rainncv_p   , snow       = snownc_p     , snowncv  = snowncv_p , &
                  graupel  = graupelnc_p , graupelncv = graupelncv_p , sr       = sr_p      , &
                  re_cloud = recloud_p   , re_ice     = reice_p      , re_snow  = resnow_p  , &
                  has_reqc = has_reqc    , has_reqi   = has_reqi     , has_reqs = has_reqs  , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde     , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme     , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte       &
                      )
      write(0,*) '--- end subroutine mp_wsm6_wrf381.'

    case default
       
 end select microp_select

!... calculate the 10cm radar reflectivity and relative humidity, if needed:
 if (l_diags) then
 
    ! Ensure that we only call compute_radar_reflectivity() if we are using an MPS that
    !    supports the computation of simulated radar reflectivity
    if(trim(microp_scheme) == "mp_wsm6"        .or. &
       trim(microp_scheme) == "mp_wsm6_wrf381" .or. &
       trim(microp_scheme) == "mp_thompson"    .or. &
       trim(microp_scheme) == "mp_thompson_aerosols") then
       call compute_radar_reflectivity(diag_physics,its,ite)
    else
       write(0,*) '*** NOTICE: NOT computing simulated radar reflectivity'
       write(0,*) '            since WSM6 or Thompson microphysics scheme was not selected'
    end if

    ! calculate the relative humidity over water if the temperature is strictly greater than 0.C,
    ! over ice otherwise.
    call compute_relhum(mesh,diag,its,ite)

 end if
 write(0,*) '--- l_diags =', l_diags
!call compute_relhum(diag_physics)

!... copy updated precipitation from the wrf-physics grid back to the geodesic-dynamics grid:
 call precip_to_MPAS(configs,mesh,diag_physics,its,ite)

!... copy updated cloud microphysics variables from the wrf-physics grid back to the geodesic-
!    dynamics grid:
 call microphysics_to_MPAS(mesh,state,time_lev,diag,diag_physics,tend,itimestep,its,ite)

!... deallocation of all microphysics arrays:
!$OMP BARRIER
!$OMP MASTER
 call deallocate_microphysics
!$OMP END MASTER

!formats:
 200 format(i3,i6,10(1x,e15.8))
 201 format(3i6,10(1x,e15.8))
 203 format('microphysics begins:',3i6,2(1x,f6.1))
 204 format('microphysics precip:',3i6,8(1x,e15.8))

 end subroutine microphysics_driver

!=================================================================================================================
 subroutine precip_from_MPAS(diag_physics,mesh,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type) ,intent(in):: mesh
 integer,intent(in):: its,ite

!output variables:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer,pointer:: nCellsSolve
 real,dimension(:),pointer:: graupelncv,rainncv,snowncv,sr 

!local variables and arrays:
 integer:: i,iCell,j

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(diag_physics,'graupelncv',graupelncv)
 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call mpas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
 call mpas_pool_get_array(diag_physics,'sr'        ,sr        )

!variables common to all cloud microphysics schemes:
 do j = jts, jte
 do i = its, ite
    rainncv_p(i,j) = 0._RKIND
    rainnc_p(i,j)  = 0._RKIND
 enddo
 enddo

 do iCell = its,ite
    rainncv(iCell) = 0._RKIND
 enddo

!variables specific to different cloud microphysics schemes:

 microp_select_init: select case(microp_scheme)

    case ("mp_thompson","mp_thompson_aerosols","mp_wsm6","mp_wsm6_wrf381")
       do j = jts, jte
       do i = its, ite
          snowncv_p(i,j)    = 0._RKIND
          graupelncv_p(i,j) = 0._RKIND
          snownc_p(i,j)     = 0._RKIND
          graupelnc_p(i,j)  = 0._RKIND
          sr_p(i,j)         = 0._RKIND
       enddo
       enddo

       do iCell = its,ite
          snowncv(iCell)    = 0._RKIND
          graupelncv(iCell) = 0._RKIND
          sr(iCell)         = 0._RKIND
       enddo 

    case default

 end select microp_select_init

 end subroutine precip_from_MPAS

!=================================================================================================================
 subroutine precip_to_MPAS(configs,mesh,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: i_rainnc

 real(kind=RKIND),pointer:: config_bucket_rainnc
 real(kind=RKIND),dimension(:),pointer:: precipw
 real(kind=RKIND),dimension(:),pointer:: graupelnc,rainnc,snownc
 real(kind=RKIND),dimension(:),pointer:: graupelncv,rainncv,snowncv,sr

!local variables and arrays:
 integer:: i,j,k
 real(kind=RKIND):: rho_a

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_bucket_rainnc',config_bucket_rainnc)

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(diag_physics,'i_rainnc'  ,i_rainnc  )
 call mpas_pool_get_array(diag_physics,'precipw'   ,precipw   )
 call mpas_pool_get_array(diag_physics,'graupelnc' ,graupelnc )
 call mpas_pool_get_array(diag_physics,'graupelncv',graupelncv)
 call mpas_pool_get_array(diag_physics,'rainnc'    ,rainnc    )
 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call mpas_pool_get_array(diag_physics,'snownc'    ,snownc    )
 call mpas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
 call mpas_pool_get_array(diag_physics,'sr'        ,sr        )

 do i = its,ite
    precipw(i) = 0._RKIND    
 enddo

!variables common to all cloud microphysics schemes:

 do j = jts,jte
 do i = its,ite

    !precipitable water:
    do k = kts,kte
       rho_a = rho_p(i,k,j) / (1._RKIND + qv_p(i,k,j))
       precipw(i) = precipw(i) + qv_p(i,k,j) * rho_a * dz_p(i,k,j)
    enddo

    !time-step precipitation:
    rainncv(i) = rainnc_p(i,j)
    
    !accumulated precipitation:
    rainnc(i) = rainnc(i) + rainncv(i)

   if(l_acrain .and. config_bucket_rainnc.gt.0._RKIND .and. &
      rainnc(i).gt.config_bucket_rainnc) then
      i_rainnc(i) = i_rainnc(i) + 1
      rainnc(i)   = rainnc(i) - config_bucket_rainnc
   endif
 
 enddo
 enddo

!variables specific to different cloud microphysics schemes:

 microp_select_init: select case(microp_scheme)

    case ("mp_thompson","mp_thompson_aerosols","mp_wsm6","mp_wsm6_wrf381")
       do j = jts,jte
       do i = its,ite
          !time-step precipitation:
          snowncv(i)    = snownc_p(i,j)
          graupelncv(i) = graupelnc_p(i,j)
          sr(i) = (snownc_p(i,j) + graupelnc_p(i,j)) / (rainnc_p(i,j)+1.e-12)

          !accumulated precipitation:
          snownc(i)    = snownc(i) + snowncv(i)
          graupelnc(i) = graupelnc(i) + graupelncv(i)
       enddo
       enddo

    case default

 end select microp_select_init

 end subroutine precip_to_MPAS

!=================================================================================================================
 subroutine compute_radar_reflectivity(diag_physics,its,ite)
!=================================================================================================================

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

 integer,intent(in):: its,ite

!local pointers:
 real(kind=RKIND),dimension(:),pointer:: refl10cm_max, refl10cm_1km, refl10cm_1km_max

!local variables and arrays:
 integer:: i,j,k,kp
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qc1d,qr1d,qs1d,qg1d,t1d,p1d,nr1d,dBZ1d,zp
 real(kind=RKIND):: w1,w2

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine COMPUTE_RADAR_REFLECTIVITY:'

 call mpas_pool_get_array(diag_physics,'refl10cm_max',refl10cm_max)
 call mpas_pool_get_array(diag_physics,'refl10cm_1km',refl10cm_1km)
 call mpas_pool_get_array(diag_physics,'refl10cm_1km_max',refl10cm_1km_max)

 microp_select: select case(microp_scheme)

    case ("mp_kessler")
       call physics_error_fatal('--- calculation of radar reflectivity is not available' // &
                                 'with kessler cloud microphysics')

    case ("mp_wsm6","mp_wsm6_wrf381")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp))    allocate(zp(kts:kte))

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j)+0.5*dz_p(i,1,j) ! height AGL
          enddo

          call refl10cm_wsm6(qv1d,qr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp))    deallocate(zp)

    case ("mp_thompson","mp_thompson_aerosols")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qc1d) ) allocate(qc1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(nr1d) ) allocate(nr1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp))    allocate(zp(kts:kte))

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             qc1d(k)  = qc_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             nr1d(k)  = nr_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j)+0.5*dz_p(i,1,j) ! height AGL
          enddo

          call calc_refl10cm(qv1d,qc1d,qr1d,nr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qc1d) ) deallocate(qc1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(nr1d) ) deallocate(nr1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp))    deallocate(zp)

    case default

 end select microp_select
!write(0,*) '--- end subroutine COMPUTE_RADAR_REFLECTIVITY'

 end subroutine compute_radar_reflectivity

!=================================================================================================================
 subroutine compute_relhum(mesh,diag,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 real(kind=RKIND),dimension(:,:),pointer:: relhum

!local variables and arrays:
 integer:: i,j,k

 real(kind=RKIND):: tempc, rh
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qvs1d,t1d,p1d

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(diag,'relhum',relhum)

 relhum(1:nVertLevels,its:ite) = 0._RKIND

 if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
 if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
 if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
 if(.not.allocated(qvs1d)) allocate(qvs1d(kts:kte))

 do j = jts,jte
 do i = its,ite

    do k = kts,kte
       p1d(k)  = pres_p(i,k,j)
       t1d(k)  = th_p(i,k,j) * pi_p(i,k,j)
       tempc   = t1d(k) - 273.16_RKIND
       qvs1d(k) = rslf(p1d(k),t1d(k))
       if(tempc .le. 0._RKIND) qvs1d(k) = rsif(p1d(k),t1d(k))
       qv1d(k) = qv_p(i,k,j)
       relhum(k,i) = qv1d(k) / qvs1d(k)
       relhum(k,i) = relhum(k,i) * 100._RKIND
    enddo
!   write(0,*)

 enddo
 enddo
 101 format(2(1x,i8),5(1x,e15.8))

 if(allocated(p1d)  ) deallocate(p1d  )
 if(allocated(t1d)  ) deallocate(t1d  )
 if(allocated(qv1d) ) deallocate(qv1d )
 if(allocated(qvs1d)) deallocate(qvs1d)

 end subroutine compute_relhum

!=================================================================================================================
 end module mpas_atmphys_driver_microphysics
!=================================================================================================================
