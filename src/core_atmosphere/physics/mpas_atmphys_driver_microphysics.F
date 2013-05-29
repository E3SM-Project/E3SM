! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=============================================================================================
 module mpas_atmphys_driver_microphysics
 use mpas_configure
 use mpas_grid_types
 
 use mpas_atmphys_constants
 use mpas_atmphys_vars
 use mpas_atmphys_interface_nhyd

!wrf physics:
 use module_mp_kessler
 use module_mp_wsm6

 implicit none
 private
 public:: allocate_microphysics,   &
          deallocate_microphysics, &
          microphysics_driver  ,   &
          microphysics_init


!>\brief MPAS driver for parameterization of cloud microphysics processes.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_driver_microphysics:
!> ------------------------------------------------
!> allocate_microphysics  : allocate local arrays for parameterization of cloud microphysics.
!> deallocate_microphysics: deallocate local arrays for parameterization of cloud microphysics.
!> microphysics_init      : initialization of individual cloud microphysics schemes.
!> microphysics_driver    : main driver (called from mpas_atm_time_integration).
!> precip_from_MPAS       : initialize timestep local arrays for precipitation.
!> precip_to_MPAS         : copy local arrays to MPAS arrays.
!> compute_radar_reflectivity: compute radar reflectivities.
!>
!> WRF physics called from microphysics_driver:
!> --------------------------------------------
!>    * module_mp_kessler : Kessler cloud microphysics.
!>    * module_mp_wsm6    : WSM6 cloud microphysics.
!>
!> comments:
!> unlike all the other physics parameterizations, parameterizations of cloud microphysics schemes
!> are called at the bottom of subroutine atm_srk3 in module atm_time_integration.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed call to the Thompson cloud microphysics scheme until the scheme is updated to that
!>      in WRF revision 3.5.
!>      Laura D. Fowler (birch.ucar.edu) / 2013-05-29.


 contains

!==================================================================================================
 subroutine allocate_microphysics
!==================================================================================================

!sounding variables:
 if(.not.allocated(rho_p)     ) allocate(rho_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(th_p)      ) allocate(th_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(pi_p)      ) allocate(pi_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(pres_p)    ) allocate(pres_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(z_p)       ) allocate(z_p(ims:ime,kms:kme,jms:jme)     )
 if(.not.allocated(dz_p)      ) allocate(dz_p(ims:ime,kms:kme,jms:jme)    )

!mass mixing ratios:
 if(.not.allocated(qv_p)      ) allocate(qv_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(qc_p)      ) allocate(qc_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(qr_p)      ) allocate(qr_p(ims:ime,kms:kme,jms:jme) )

 !surface precipitation:
 if(.not.allocated(rainnc_p)  ) allocate(rainnc_p(ims:ime,jms:jme)     )
 if(.not.allocated(rainncv_p) ) allocate(rainncv_p(ims:ime,jms:jme)    )

 microp_select: select case(microp_scheme)

    case ("wsm6")
       !mass mixing ratios:
       if(.not.allocated(qi_p)         ) allocate(qi_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(qs_p)         ) allocate(qs_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(qg_p)         ) allocate(qg_p(ims:ime,kms:kme,jms:jme)  )

       !surface precipitation:
       if(.not.allocated(sr_p)         ) allocate(sr_p(ims:ime,jms:jme)          )
       if(.not.allocated(snownc_p)     ) allocate(snownc_p(ims:ime,jms:jme)      )
       if(.not.allocated(snowncv_p)    ) allocate(snowncv_p(ims:ime,jms:jme)     )
       if(.not.allocated(graupelnc_p)  ) allocate(graupelnc_p(ims:ime,jms:jme)   )
       if(.not.allocated(graupelncv_p) ) allocate(graupelncv_p(ims:ime,jms:jme)  )

    case default

 end select microp_select

 end subroutine allocate_microphysics

!==================================================================================================
 subroutine deallocate_microphysics
!==================================================================================================

!sounding variables:
 if(allocated(rho_p)     ) deallocate(rho_p     )
 if(allocated(th_p)      ) deallocate(th_p      )
 if(allocated(pi_p)      ) deallocate(pi_p      )
 if(allocated(pres_p)    ) deallocate(pres_p    )
 if(allocated(z_p)       ) deallocate(z_p       )
 if(allocated(dz_p)      ) deallocate(dz_p      )

!mass mixing ratios:
 if(allocated(qv_p)      ) deallocate(qv_p      )
 if(allocated(qc_p)      ) deallocate(qc_p      )
 if(allocated(qr_p)      ) deallocate(qr_p      )

 !surface precipitation:
 if(allocated(rainnc_p)  ) deallocate(rainnc_p  )
 if(allocated(rainncv_p) ) deallocate(rainncv_p )

 microp_select: select case(microp_scheme)

    case ("wsm6")
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

    case default

 end select microp_select

 end subroutine deallocate_microphysics

!==================================================================================================
 subroutine microphysics_init
!==================================================================================================
 write(0,*)
 write(0,*) '--- enter cloud microphysics initialization:'

 microp_select: select case(microp_scheme)

    case("wsm6")
       write(0,*) '--- enter subroutine wsm6init:'
       call wsm6init(rho_a,rho_r,rho_s,cliq,cpv,.false.)
       write(0,*) '--- end subroutine wsm6init'

    case default

 end select microp_select

 write(0,*) '--- end cloud microphysics initialization'

 end subroutine microphysics_init

!==================================================================================================
 subroutine microphysics_driver(state,diag,diag_physics,tend,mesh,itimestep)
!==================================================================================================

!input arguments:
!----------------
 type(mesh_type),intent(in):: mesh
 integer,intent(in):: itimestep
    
!inout arguments:
!----------------
 type(state_type),intent(inout):: state
 type(diag_type),intent(inout):: diag
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_type),intent(inout):: tend

!local variables and arrays:
!---------------------------
 logical:: log_microphysics
 integer:: i,icell,icount,istep,j,k,kk

!==================================================================================================
 write(0,*)
 write(0,*) '--- enter subroutine microphysics_driver: itimestep=', itimestep
 write(0,*) '    dt_microp=',dt_microp
 write(0,*) '    n_microp =',n_microp
 
!initialization:
 write(0,*) '    ims= ',ims,' ime=',ime
 write(0,*) '    jms= ',jms,' jme=',jme
 write(0,*) '    kms= ',kms,' kme=',kme
 write(0,*)
 write(0,*) '    ids= ',ids,' ide=',ide
 write(0,*) '    jds= ',jds,' jde=',jde
 write(0,*) '    kds= ',kds,' kde=',kde
 write(0,*)
 write(0,*) '    its= ',its,' ite=',ite
 write(0,*) '    jts= ',jts,' jte=',jte
 write(0,*) '    kts= ',kts,' kte=',kte

!... allocation of microphysics arrays:
 call allocate_microphysics

!... initialization of precipitation related arrays:
 call precip_from_MPAS(diag_physics,mesh)

!... initialization of soundings for non-hydrostatic dynamical cores.
 call microphysics_from_MPAS(mesh,state,tend,diag)

!--------------------------------------------------------------------------------------------------

!... call to different cloud microphysics schemes:

!--------------------------------------------------------------------------------------------------

 istep = 1

 do while (istep .le. n_microp)

    microp_select: select case(microp_scheme)

       case ("kessler")

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
          write(0,*) '--- end kessler:',istep

       case ("wsm6")
          call wsm6( &
                  th      = th_p        , q          = qv_p         , qc      = qc_p      , &
                  qr      = qr_p        , qi         = qi_p         , qs      = qs_p      , &
                  qg      = qg_p        , den        = rho_p        , pii     = pi_p      , &
                  p       = pres_p      , delz       = dz_p         , delt    = dt_microp , &
                  g       = g           , cpd        = cp           , cpv     = cpv       , &
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
          write(0,*) '--- end wsm6:',istep

       case default
       
    end select microp_select
 
    istep = istep + 1
    
 end do

!... calculate the 10cm radar reflectivity, if needed:
 if (l_diags) then
 
    ! Ensure that we only call compute_radar_reflectivity() if we are using an MPS that
    !    supports the computation of simulated radar reflectivity
    if (trim(microp_scheme) == "wsm6") then
       call compute_radar_reflectivity(diag_physics)
    else
       write(0,*) '*** NOTICE: NOT computing simulated radar reflectivity'
       write(0,*) '            since WSM6 microphysics scheme was not selected'
    end if
 end if

!... copy updated precipitation from the wrf-physics grid back to the geodesic-dynamics grid:
 call precip_to_MPAS(config_bucket_rainnc,diag_physics)

!... copy updated cloud microphysics variables from the wrf-physics grid back to the geodesic-
!    dynamics grid:
 call microphysics_to_MPAS(mesh,state,diag,tend,itimestep)

!... deallocation of all microphysics arrays:
 call deallocate_microphysics

!formats:
 200 format(i3,i6,10(1x,e15.8))
 201 format(3i6,10(1x,e15.8))
 203 format('microphysics begins:',3i6,2(1x,f6.1))
 204 format('microphysics precip:',3i6,8(1x,e15.8))

 end subroutine microphysics_driver

!==================================================================================================
 subroutine precip_from_MPAS(diag_physics,mesh)
!==================================================================================================

!input variables:
 type(mesh_type) ,intent(in):: mesh

!output variables:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: i,iCell,j

!--------------------------------------------------------------------------------------------------

!variables common to all cloud microphysics schemes:
 do j = jts, jte
 do i = its, ite
    rainncv_p(i,j) = 0._RKIND
    rainnc_p(i,j)  = 0._RKIND
 enddo
 enddo

 do iCell = 1, mesh % nCellsSolve
    diag_physics % rainncv % array(iCell) = 0._RKIND
 enddo

!variables specific to different cloud microphysics schemes:

 microp_select_init: select case(microp_scheme)

    case ("wsm6")
 
       do j = jts, jte
       do i = its, ite
          snowncv_p(i,j)    = 0._RKIND
          graupelncv_p(i,j) = 0._RKIND
          snownc_p(i,j)     = 0._RKIND
          graupelnc_p(i,j)  = 0._RKIND
          sr_p(i,j)         = 0._RKIND
       enddo
       enddo

       do iCell = 1, mesh % nCellsSolve
          diag_physics % snowncv % array(iCell)    = 0._RKIND
          diag_physics % graupelncv % array(iCell) = 0._RKIND
          diag_physics % sr % array(iCell)         = 0._RKIND
       enddo 

    case default

 end select microp_select_init

 end subroutine precip_from_MPAS

!==================================================================================================
 subroutine precip_to_MPAS(bucket_rainnc,diag_physics)
!==================================================================================================

!output variables:
 real(kind=RKIND),intent(in):: bucket_rainnc
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: i,j,k
 real(kind=RKIND):: rho_a

!--------------------------------------------------------------------------------------------------

!variables common to all cloud microphysics schemes:

 do j = jts,jte
 do i = its,ite

    !precipitable water:
    diag_physics % precipw % array(i) = 0._RKIND
    do k = kts,kte
       rho_a = rho_p(i,k,j) / (1._RKIND + qv_p(i,k,j))
       diag_physics % precipw % array(i) = &
            diag_physics % precipw % array(i) + qv_p(i,k,j) * rho_a * dz_p(i,k,j)
    enddo

    !time-step precipitation:
    diag_physics % rainncv % array(i) = rainnc_p(i,j)
    
    !accumulated precipitation:
    diag_physics % rainnc % array(i) = diag_physics % rainnc % array(i) &
                                     + diag_physics % rainncv % array(i)

   if(l_acrain .and. bucket_rainnc.gt.0._RKIND .and. &
      diag_physics%rainnc%array(i).gt.bucket_rainnc) then
      diag_physics % i_rainnc % array(i) = diag_physics % i_rainnc % array(i) + 1
      diag_physics % rainnc % array(i) = diag_physics % rainnc % array(i) - bucket_rainnc
   endif
 
 enddo
 enddo

!variables specific to different cloud microphysics schemes:

 microp_select_init: select case(microp_scheme)

    case ("wsm6")

       do j = jts,jte
       do i = its,ite
       
          !time-step precipitation:
          diag_physics % snowncv % array(i)    = snownc_p(i,j)
          diag_physics % graupelncv % array(i) = graupelnc_p(i,j)
          diag_physics % sr % array(i) = &
                         (snownc_p(i,j) + graupelnc_p(i,j)) / (rainnc_p(i,j)+1.e-12)

          !accumulated precipitation:
          diag_physics % snownc % array(i)    = diag_physics % snownc  % array(i)    &
                                              + diag_physics % snowncv % array(i)
          diag_physics % graupelnc % array(i) = diag_physics % graupelnc  % array(i) &
                                              + diag_physics % graupelncv % array(i)

       enddo
       enddo

    case default

 end select microp_select_init

 end subroutine precip_to_MPAS

!==================================================================================================
 subroutine compute_radar_reflectivity(diag_physics)
!==================================================================================================

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: i,j,k
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qr1d,qs1d,qg1d,t1d,p1d,dBZ1d

!--------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine COMPUTE_RADAR_REFLECTIVITY:'

 microp_select: select case(microp_scheme)

    case ("kessler")
       call physics_error_fatal('--- calculation of radar reflectivity is not available' // &
                                 'with kessler cloud microphysics')

    case ("wsm6")

       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))

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
          enddo

          call refl10cm_wsm6(qv1d,qr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)

          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
!            write(0,201) i,k,dBZ1d(k)
          enddo
          diag_physics % refl10cm_max % array(i) = maxval(dBZ1d(:))
!         if(diag_physics % refl10cm_max % array(i) .gt. 0.) &
!            write(0,201) j,i,diag_physics % refl10cm_max % array(i)
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)

    case default

 end select microp_select
!write(0,*) '--- end subroutine COMPUTE_RADAR_REFLECTIVITY'

 201 format(2i6,e15.8)

 end subroutine compute_radar_reflectivity

!==================================================================================================
 end module mpas_atmphys_driver_microphysics
!==================================================================================================
