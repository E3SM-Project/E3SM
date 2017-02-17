! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_interface
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_constants
 use mpas_atmphys_vars

 implicit none
 private
 public:: allocate_forall_physics,   &
          deallocate_forall_physics, &
          MPAS_to_physics,           &
          microphysics_from_MPAS,    &
          microphysics_to_MPAS

!Interface for conversion between variables used in the MPAS dynamical core and variables needed in the
!physics parameterizations.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_interface_nhyd:
! -------------------------------------------
! allocate_forall_physics  : allocation of all meteorological variables needed in the physics.
! deallocate_forall_physics: deallocation of all meteorological variables needed in the physics.
! MPAS_to_physics          : conversion of input "MPAS" variables to "physics" variables.
! microphysics_from_MPAS   : initialize local arrays needed in cloud microphysics schemes.
! microphysics_to_MPAS     : copy local arrays needed in cloud microphysics schemesto MPAS arrays.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * in subroutine MPAS_to_physics, moved the calculation of the local arrays qv_p,qc_p, and qr_p above the 
!   calculation of th_p so that we do not need to use the pointer qv.
! * in subroutine microphysics_from_MPAS, moved the calculation of the local arrays qv_p,qc_p, and qr_p above 
!   the calculation of th_p so that we do not need to use the pointer qv.
! * in subroutine microphysics_to_MPAS, since microphysics schemes update the temperature and water vapor
!   mixing ratio, we first update the total pressure and exner functions. Then, we update the modified
!   potential temperature and calculate the diabatic tendency due to cloud microphysics processes.
!   Laura D. Fowler (laura@ucar.edu) / 2013-11-07.
! * replaced the variable g (that originally pointed to gravity) with gravity, for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * added calculation of the surface pressure tendency. Moved the calculation of znu_p below the calculation 
!   of the surface pressure to avoid dividing by zero if the surface pressure is not output in the init file.
!   Laura D. Fowler (birch.mmm.ucar.ecu) / 2014-06-23. 
! * renamed module mpas_atmphys_interface_nhyd to mpas_atmphys_interface.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-19.
! * in subroutine microphysics_to_MPAS, reverted the calculation of cloud microphysics tendency rt_diabatic_tend, 
!   and update of the state variables to git hash identifier 7a66f273e182f4. This change reflects the fact that
!   we want to compute rt_diabatic_tend at constant volume.
!   Laura D. Fowler (laura@ucar.edu) / 2014-01-015. 
! * added the initialization of local variables needed in the parameterization of the MYNN surface layer scheme
!   and planetary boundary layer scheme from WRF 3.6.1.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-11.
! * corrected the calculation of the surface pressure, mainly extrapolation of the air density to the surface.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-25.


 contains


!=================================================================================================================
 subroutine allocate_forall_physics(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 if(.not.allocated(psfc_p) ) allocate(psfc_p(ims:ime,jms:jme)         )
 if(.not.allocated(ptop_p) ) allocate(ptop_p(ims:ime,jms:jme)         )

 if(.not.allocated(u_p)    ) allocate(u_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(v_p)    ) allocate(v_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(fzm_p)  ) allocate(fzm_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(fzp_p)  ) allocate(fzp_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(zz_p)   ) allocate(zz_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(pres_p) ) allocate(pres_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(pi_p)   ) allocate(pi_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(z_p)    ) allocate(z_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(zmid_p) ) allocate(zmid_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(dz_p)   ) allocate(dz_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(t_p)    ) allocate(t_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(th_p)   ) allocate(th_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(al_p)   ) allocate(al_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(rho_p)  ) allocate(rho_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(rh_p)   ) allocate(rh_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(znu_p)  ) allocate(znu_p(ims:ime,kms:kme,jms:jme)  )

 if(.not.allocated(w_p)    ) allocate(w_p(ims:ime,kms:kme,jms:jme)    )
 if(.not.allocated(pres2_p)) allocate(pres2_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(t2_p)   ) allocate(t2_p(ims:ime,kms:kme,jms:jme)   )
 
 if(.not.allocated(qv_p)   ) allocate(qv_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qc_p)   ) allocate(qc_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qr_p)   ) allocate(qr_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qi_p)   ) allocate(qi_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qs_p)   ) allocate(qs_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(qg_p)   ) allocate(qg_p(ims:ime,kms:kme,jms:jme)   )

 pbl_select: select case (trim(pbl_scheme))
    case("bl_mynn")
       if(.not.allocated(ni_p)) allocate(ni_p(ims:ime,kms:kme,jms:jme))

    case default

 end select pbl_select

!... arrays used for calculating the hydrostatic pressure and exner function:
 if(.not.allocated(psfc_hyd_p)  ) allocate(psfc_hyd_p(ims:ime,jms:jme)          )
 if(.not.allocated(psfc_hydd_p) ) allocate(psfc_hydd_p(ims:ime,jms:jme)         )
 if(.not.allocated(pres_hyd_p)  ) allocate(pres_hyd_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pres_hydd_p) ) allocate(pres_hydd_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(pres2_hyd_p) ) allocate(pres2_hyd_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(pres2_hydd_p)) allocate(pres2_hydd_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(znu_hyd_p)   ) allocate(znu_hyd_p(ims:ime,kms:kme,jms:jme)   )
 
 end subroutine allocate_forall_physics

!=================================================================================================================
 subroutine deallocate_forall_physics(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 if(allocated(psfc_p)  ) deallocate(psfc_p  )
 if(allocated(ptop_p)  ) deallocate(ptop_p  )

 if(allocated(u_p)     ) deallocate(u_p     )
 if(allocated(v_p)     ) deallocate(v_p     )
 if(allocated(fzm_p)   ) deallocate(fzm_p   )
 if(allocated(fzp_p)   ) deallocate(fzp_p   )
 if(allocated(zz_p)    ) deallocate(zz_p    )
 if(allocated(pres_p)  ) deallocate(pres_p  )
 if(allocated(pi_p)    ) deallocate(pi_p    )
 if(allocated(z_p)     ) deallocate(z_p     )
 if(allocated(zmid_p)  ) deallocate(zmid_p  )
 if(allocated(dz_p)    ) deallocate(dz_p    )
 if(allocated(t_p)     ) deallocate(t_p     )
 if(allocated(th_p)    ) deallocate(th_p    )
 if(allocated(al_p)    ) deallocate(al_p    )
 if(allocated(rho_p)   ) deallocate(rho_p   ) 
 if(allocated(rh_p)    ) deallocate(rh_p    ) 
 if(allocated(znu_p)   ) deallocate(znu_p   )

 if(allocated(w_p)     ) deallocate(w_p     )
 if(allocated(pres2_p) ) deallocate(pres2_p )
 if(allocated(t2_p)    ) deallocate(t2_p    )

 if(allocated(qv_p)    ) deallocate(qv_p    )
 if(allocated(qc_p)    ) deallocate(qc_p    )
 if(allocated(qr_p)    ) deallocate(qr_p    )
 if(allocated(qi_p)    ) deallocate(qi_p    )
 if(allocated(qs_p)    ) deallocate(qs_p    )
 if(allocated(qg_p)    ) deallocate(qg_p    )

 pbl_select: select case (trim(pbl_scheme))
    case("bl_mynn")
       if(allocated(ni_p)) deallocate(ni_p)

    case default

 end select pbl_select

 if(allocated(psfc_hyd_p)  ) deallocate(psfc_hyd_p  )
 if(allocated(psfc_hydd_p) ) deallocate(psfc_hydd_p )
 if(allocated(pres_hyd_p)  ) deallocate(pres_hyd_p  )
 if(allocated(pres_hydd_p) ) deallocate(pres_hydd_p )
 if(allocated(pres2_hyd_p) ) deallocate(pres2_hyd_p )
 if(allocated(pres2_hydd_p)) deallocate(pres2_hydd_p)
 if(allocated(znu_hyd_p)   ) deallocate(znu_hyd_p   )
 
 end subroutine deallocate_forall_physics

!=================================================================================================================
 subroutine MPAS_to_physics(configs,mesh,state,time_lev,diag,diag_physics,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag

 integer,intent(in):: its,ite
 integer,intent(in):: time_lev

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: pbl_scheme

 integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
 integer,pointer:: index_ni

 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: fzm,fzp,rdzw
 real(kind=RKIND),dimension(:),pointer    :: surface_pressure,plrad
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid
 real(kind=RKIND),dimension(:,:),pointer  :: zz,exner,pressure_b,rtheta_p,rtheta_b
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p,u,v,w
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qr,qi,qs,qg
 real(kind=RKIND),dimension(:,:),pointer  :: ni
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!local variables:
 integer:: i,k,j

 real(kind=RKIND):: z0,z1,z2,w1,w2
 real(kind=RKIND):: rho_a,rho1,rho2,tem1,tem2

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine MPAS_to_phys:'
!write(0,*) 'ims=',ims,' ime=',ime
!write(0,*) 'jms=',jms,' jme=',jme
!write(0,*) 'kms=',kms,' kme=',kme
!write(0,*)
!write(0,*) 'its=',its,' ite=',ite
!write(0,*) 'jts=',jts,' jte=',jte
!write(0,*) 'kts=',kts,' kte=',kte

!initialization:
 call mpas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'fzm'    ,fzm    )
 call mpas_pool_get_array(mesh,'fzp'    ,fzp    )
 call mpas_pool_get_array(mesh,'rdzw'   ,rdzw   )
 call mpas_pool_get_array(mesh,'zgrid'  ,zgrid  )
 call mpas_pool_get_array(mesh,'zz'     ,zz     )

 call mpas_pool_get_array(diag,'surface_pressure'      ,surface_pressure)
 call mpas_pool_get_array(diag,'exner'                 ,exner           )
 call mpas_pool_get_array(diag,'pressure_base'         ,pressure_b      )
 call mpas_pool_get_array(diag,'pressure_p'            ,pressure_p      )
 call mpas_pool_get_array(diag,'rtheta_base'           ,rtheta_b        )
 call mpas_pool_get_array(diag,'rtheta_p'              ,rtheta_p        )
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,u               )
 call mpas_pool_get_array(diag,'uReconstructMeridional',v               )

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'w'      ,w      ,time_lev)

 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_dimension(state,'index_qc',index_qc)
 call mpas_pool_get_dimension(state,'index_qr',index_qr)
 call mpas_pool_get_dimension(state,'index_qi',index_qi)
 call mpas_pool_get_dimension(state,'index_qs',index_qs)
 call mpas_pool_get_dimension(state,'index_qg',index_qg)

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)
 qc => scalars(index_qc,:,:)
 qr => scalars(index_qr,:,:)
 qi => scalars(index_qi,:,:)
 qs => scalars(index_qs,:,:)
 qg => scalars(index_qg,:,:)

 call mpas_pool_get_array(diag_physics,'plrad',plrad)

 do j = jts, jte
 do k = kts, kte
 do i = its, ite

    !water vapor and moist arrays:
    qv_p(i,k,j) = max(0.,qv(k,i))
    qc_p(i,k,j) = max(0.,qc(k,i))
    qr_p(i,k,j) = max(0.,qr(k,i))
    qi_p(i,k,j) = max(0.,qi(k,i))
    qs_p(i,k,j) = max(0.,qs(k,i))
    qg_p(i,k,j) = max(0.,qg(k,i))

    !arrays located at theta points:
    u_p(i,k,j) = u(k,i)
    v_p(i,k,j) = v(k,i)

    zz_p(i,k,j)  = zz(k,i)
    rho_p(i,k,j) = zz(k,i) * rho_zz(k,i)
    rho_p(i,k,j) = rho_p(i,k,j)*(1._RKIND + qv_p(i,k,j))
    th_p(i,k,j)  = theta_m(k,i) / (1._RKIND + R_v/R_d * qv_p(i,k,j))
    t_p(i,k,j)   = th_p(i,k,j)*exner(k,i)

    pi_p(i,k,j)   = exner(k,i)
    pres_p(i,k,j) = pressure_p(k,i) + pressure_b(k,i)

    zmid_p(i,k,j) = 0.5*(zgrid(k+1,i)+zgrid(k,i))
    dz_p(i,k,j)   = zgrid(k+1,i)-zgrid(k,i)

 enddo
 enddo
 enddo

 pbl_select: select case (trim(pbl_scheme))
    case("bl_mynn")
       call mpas_pool_get_dimension(state,'index_ni',index_ni)
       ni => scalars(index_ni,:,:)

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          ni_p(i,k,j) = max(0.,ni(k,i))
       enddo
       enddo
       enddo

    case default

 end select pbl_select

!calculation of the surface pressure using hydrostatic assumption down to the surface::
 do j = jts,jte
 do i = its,ite
    tem1 = zgrid(2,i)-zgrid(1,i)
    tem2 = zgrid(3,i)-zgrid(2,i)
    rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv_p(i,1,j))
    rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv_p(i,2,j))
!   surface_pressure(i) = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
!                   * (rho1 + 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i) = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
                    * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i) = surface_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
 enddo
 
 do k = kts,kte
 do i = its,ite
    znu_p(i,k,j)  = pres_p(i,k,j) / surface_pressure(i)
 enddo
 enddo
 enddo

!arrays located at w points:
 do j = jts, jte
 do k = kts,kte+1
 do i = its,ite
    w_p(i,k,j) = w(k,i)
    z_p(i,k,j) = zgrid(k,i)
 enddo
 enddo
 enddo

!check that the pressure in the layer above the surface is greater than that in the layer
!above it:
 201 format(3i8,10(1x,e15.8))
 do j = jts,jte
 do i = its,ite
    if(pres_p(i,1,j) .le. pres_p(i,2,j)) then
       write(0,*)
       write(0,*) '--- subroutine MPAS_to_phys - pressure(1) < pressure(2):'
       write(0,*) 'i      =', i
       write(0,*) 'latCell=', latCell(i)/degrad
       write(0,*) 'lonCell=', lonCell(i)/degrad
       do k = kts,kte
          write(0,201) j,i,k,dz_p(i,k,j),pressure_b(k,i),pressure_p(k,i),pres_p(i,k,j), &
             rho_p(i,k,j),th_p(i,k,j),t_p(i,k,j),qv_p(i,k,j)
       enddo
!      call mpas_dmpar_global_abort('ERROR: pressure increasing with height')
    endif
 enddo
 enddo

!interpolation of pressure and temperature from theta points to w points:
 do j = jts,jte
 do k = kts+1,kte
 do i = its,ite
    tem1 = 1./(zgrid(k+1,i)-zgrid(k-1,i))
    fzm_p(i,k,j) = (zgrid(k,i)-zgrid(k-1,i)) * tem1
    fzp_p(i,k,j) = (zgrid(k+1,i)-zgrid(k,i)) * tem1
    t2_p(i,k,j)    = fzm_p(i,k,j)*t_p(i,k,j) + fzp_p(i,k,j)*t_p(i,k-1,j)
    pres2_p(i,k,j) = fzm_p(i,k,j)*pres_p(i,k,j) + fzp_p(i,k,j)*pres_p(i,k-1,j)
 enddo
 enddo
 enddo

!interpolation of pressure and temperature from theta points to the top-of-the-model: follows
!the calculation of the top-of-the-model pressure and temperature in WRF (subroutine phy_prep
!in ./dyn_em/module_big_step_utilities.F):
 k = kte+1
 do j = jts,jte
 do i = its,ite
    z0 = zgrid(k,i)
    z1 = 0.5*(zgrid(k,i)+zgrid(k-1,i)) 
    z2 = 0.5*(zgrid(k-1,i)+zgrid(k-2,i))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    t2_p(i,k,j) = w1*t_p(i,k-1,j) + w2*t_p(i,k-2,j)
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pres2_p(i,k,j) = exp(w1*log(pres_p(i,k-1,j))+w2*log(pres_p(i,k-2,j)))
 enddo
 enddo

!ldf (2012-06-22): recalculates the pressure at the surface as an extrapolation of the
!pressures in the 2 layers above the surface, as was originally done:
 k = kts
 do j = jts,jte
 do i = its,ite
    z0 = zgrid(k,i)
    z1 = 0.5*(zgrid(k,i)+zgrid(k+1,i)) 
    z2 = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    t2_p(i,k,j)    = w1*t_p(i,k,j)+w2*t_p(i,k+1,j)
    pres2_p(i,k,j) = w1*pres_p(i,k,j)+w2*pres_p(i,k+1,j)
!   psfc_p(i,j) = pres2_p(i,k,j)
    psfc_p(i,j) = surface_pressure(i)
 enddo
 enddo

!calculation of the hydrostatic pressure:
 do j = jts,jte
 do i = its,ite
    !pressure at w-points:
    k = kte+1
    pres2_hyd_p(i,k,j)  = pres2_p(i,k,j)
    pres2_hydd_p(i,k,j) = pres2_p(i,k,j)
    do k = kte,1,-1
       rho_a = rho_p(i,k,j) / (1.+qv_p(i,k,j))
       pres2_hyd_p(i,k,j)  = pres2_hyd_p(i,k+1,j)  + gravity*rho_p(i,k,j)*dz_p(i,k,j)
       pres2_hydd_p(i,k,j) = pres2_hydd_p(i,k+1,j) + gravity*rho_a*dz_p(i,k,j)
    enddo
    !pressure at theta-points:
    do k = kte,1,-1
       pres_hyd_p(i,k,j)  = 0.5*(pres2_hyd_p(i,k+1,j)+pres2_hyd_p(i,k,j))
       pres_hydd_p(i,k,j) = 0.5*(pres2_hydd_p(i,k+1,j)+pres2_hydd_p(i,k,j))
    enddo
    !surface pressure:
    psfc_hyd_p(i,j) = pres2_hyd_p(i,1,j)
    psfc_hydd_p(i,j) = pres2_hydd_p(i,1,j)
    !znu:
    do k = kte,1,-1
       znu_hyd_p(i,k,j) = pres_hyd_p(i,k,j) / psfc_hyd_p(i,j) 
    enddo
 enddo
 enddo

!save the model-top pressure:
 do j = jts,jte
 do i = its,ite
    plrad(i) = pres2_p(i,kte+1,j) 
 enddo
 enddo

 end subroutine MPAS_to_physics

!=================================================================================================================
 subroutine microphysics_from_MPAS(configs,mesh,state,time_lev,diag,diag_physics,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics

 integer,intent(in):: its,ite
 integer:: time_lev

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
 integer,pointer:: index_ni,index_nr
 real(kind=RKIND),dimension(:),pointer    :: nt_c,mu_c
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid,w
 real(kind=RKIND),dimension(:,:),pointer  :: zz,exner,pressure_b
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qr,qi,qs,qg
 real(kind=RKIND),dimension(:,:),pointer  :: ni,nr
 real(kind=RKIND),dimension(:,:),pointer  :: rainprod,evapprod
 real(kind=RKIND),dimension(:,:),pointer  :: re_cloud,re_ice,re_snow
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!local variables:
 integer:: i,k,j

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

 call mpas_pool_get_array(mesh,'zgrid',zgrid)
 call mpas_pool_get_array(mesh,'zz'   ,zz   )

 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)

 call mpas_pool_get_array(diag_physics,'nt_c'    ,nt_c    )
 call mpas_pool_get_array(diag_physics,'mu_c'    ,mu_c    )
 call mpas_pool_get_array(diag_physics,'rainprod',rainprod)
 call mpas_pool_get_array(diag_physics,'evapprod',evapprod)
 call mpas_pool_get_array(diag_physics,'re_cloud',re_cloud)
 call mpas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
 call mpas_pool_get_array(diag_physics,'re_snow' ,re_snow )

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'w'      ,w      ,time_lev)

 call mpas_pool_get_dimension(state,'index_qv'  ,index_qv  )
 call mpas_pool_get_dimension(state,'index_qc'  ,index_qc  )
 call mpas_pool_get_dimension(state,'index_qr'  ,index_qr  )
 call mpas_pool_get_dimension(state,'index_qi'  ,index_qi  )
 call mpas_pool_get_dimension(state,'index_qs'  ,index_qs  )
 call mpas_pool_get_dimension(state,'index_qg'  ,index_qg  )
 call mpas_pool_get_dimension(state,'index_ni'  ,index_ni  )
 call mpas_pool_get_dimension(state,'index_nr'  ,index_nr  )

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv   => scalars(index_qv,:,:)
 qc   => scalars(index_qc,:,:)
 qr   => scalars(index_qr,:,:)
 qi   => scalars(index_qi,:,:)
 qs   => scalars(index_qs,:,:)
 qg   => scalars(index_qg,:,:)
 ni   => scalars(index_ni,:,:)
 nr   => scalars(index_nr,:,:)

!initialize variables needed in the cloud microphysics schemes:
 do j = jts, jte
 do k = kts, kte
 do i = its, ite
    qv_p(i,k,j) = qv(k,i)
    qc_p(i,k,j) = qc(k,i)
    qr_p(i,k,j) = qr(k,i)

    rho_p(i,k,j)  = zz(k,i) * rho_zz(k,i)
    th_p(i,k,j)   = theta_m(k,i) / (1._RKIND + R_v/R_d * max(0._RKIND,qv_p(i,k,j)))

    pi_p(i,k,j)   = exner(k,i)
    pres_p(i,k,j) = pressure_b(k,i) + pressure_p(k,i)

    z_p(i,k,j)    = zgrid(k,i)
    dz_p(i,k,j)   = zgrid(k+1,i) - zgrid(k,i)
    w_p(i,k,j)    = w(k,i)
 enddo
 enddo
 enddo

!additional initialization as function of cloud microphysics scheme:
 microp_select_init: select case(microp_scheme)

    case ("mp_thompson","mp_wsm6")
       do j = jts, jte
       do k = kts, kte
       do i = its, ite
          qi_p(i,k,j) = qi(k,i)
          qs_p(i,k,j) = qs(k,i)
          qg_p(i,k,j) = qg(k,i)
       enddo
       enddo
       enddo

    microp2_select: select case(microp_scheme)

       case("mp_thompson")
          do j = jts,jte
          do i = its,ite
             muc_p(i,j) = mu_c(i)
             ntc_p(i,j) = nt_c(i)
          enddo
          enddo
          do j = jts, jte
          do k = kts, kte
          do i = its, ite
             ni_p(i,k,j) = ni(k,i)
             nr_p(i,k,j) = nr(k,i)
             rainprod_p(i,k,j) = rainprod(k,i)
             evapprod_p(i,k,j) = evapprod(k,i)
             recloud_p(i,k,j)  = re_cloud(k,i)
             reice_p(i,k,j)    = re_ice(k,i)
             resnow_p(i,k,j)   = re_snow(k,i)
          enddo
          enddo
          enddo

       case default

    end select microp2_select

    case default

 end select microp_select_init

 end subroutine microphysics_from_MPAS

!=================================================================================================================
 subroutine microphysics_to_MPAS(configs,mesh,state,time_lev,diag,diag_physics,tend,itimestep,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: itimestep,time_lev
 integer,intent(in):: its,ite

!output variables:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: tend
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
 integer,pointer:: index_ni,index_nr
 real(kind=RKIND),dimension(:),pointer    :: surface_pressure,tend_sfc_pressure
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid
 real(kind=RKIND),dimension(:,:),pointer  :: zz,exner,exner_b,pressure_b,rtheta_p,rtheta_b
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p
 real(kind=RKIND),dimension(:,:),pointer  :: rt_diabatic_tend
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qr,qi,qs,qg
 real(kind=RKIND),dimension(:,:),pointer  :: ni,nr
 real(kind=RKIND),dimension(:,:),pointer  :: rainprod,evapprod
 real(kind=RKIND),dimension(:,:),pointer  :: re_cloud,re_ice,re_snow
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!local variables:
 integer:: icount
 integer:: i,k,j
 real(kind=RKIND):: rho1,rho2,tem1,tem2

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

 call mpas_pool_get_array(mesh,'zz'   ,zz   )
 call mpas_pool_get_array(mesh,'zgrid',zgrid)

 call mpas_pool_get_array(diag,'exner'           ,exner           )
 call mpas_pool_get_array(diag,'exner_base'      ,exner_b         )
 call mpas_pool_get_array(diag,'pressure_base'   ,pressure_b      )
 call mpas_pool_get_array(diag,'pressure_p'      ,pressure_p      )
 call mpas_pool_get_array(diag,'rtheta_base'     ,rtheta_b        )
 call mpas_pool_get_array(diag,'rtheta_p'        ,rtheta_p        )
 call mpas_pool_get_array(diag,'surface_pressure',surface_pressure)

 call mpas_pool_get_array(diag_physics,'rainprod',rainprod)
 call mpas_pool_get_array(diag_physics,'evapprod',evapprod)
 call mpas_pool_get_array(diag_physics,'re_cloud',re_cloud)
 call mpas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
 call mpas_pool_get_array(diag_physics,'re_snow' ,re_snow )

 call mpas_pool_get_array(tend,'tend_sfc_pressure',tend_sfc_pressure)

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)

 call mpas_pool_get_dimension(state,'index_qv'  ,index_qv  )
 call mpas_pool_get_dimension(state,'index_qc'  ,index_qc  )
 call mpas_pool_get_dimension(state,'index_qr'  ,index_qr  )
 call mpas_pool_get_dimension(state,'index_qi'  ,index_qi  )
 call mpas_pool_get_dimension(state,'index_qs'  ,index_qs  )
 call mpas_pool_get_dimension(state,'index_qg'  ,index_qg  )
 call mpas_pool_get_dimension(state,'index_ni'  ,index_ni  )
 call mpas_pool_get_dimension(state,'index_nr'  ,index_nr  )

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv   => scalars(index_qv,:,:)
 qc   => scalars(index_qc,:,:)
 qr   => scalars(index_qr,:,:)
 qi   => scalars(index_qi,:,:)
 qs   => scalars(index_qs,:,:)
 qg   => scalars(index_qg,:,:)
 ni   => scalars(index_ni,:,:)
 nr   => scalars(index_nr,:,:)

 call mpas_pool_get_array(tend,'rt_diabatic_tend',rt_diabatic_tend)

!update variables needed in the dynamical core:
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    qv(k,i) = qv_p(i,k,j)
    qc(k,i) = qc_p(i,k,j)
    qr(k,i) = qr_p(i,k,j)

    !potential temperature and diabatic forcing:
    rt_diabatic_tend(k,i) = theta_m(k,i)
    theta_m(k,i) = th_p(i,k,j) * (1. + R_v/R_d * qv_p(i,k,j))
    rt_diabatic_tend(k,i) = (theta_m(k,i) - rt_diabatic_tend(k,i)) / dt_dyn

    !density-weighted perturbation potential temperature:
    rtheta_p(k,i) = rho_zz(k,i) * theta_m(k,i) - rtheta_b(k,i)

    !exner function:
    exner(k,i) = (zz(k,i)*(R_d/P0)*(rtheta_p(k,i)+rtheta_b(k,i)))**rcv

    !pertubation pressure:
    pressure_p(k,i) = zz(k,i)*R_d*(exner(k,i)*rtheta_p(k,i) &
                    + (exner(k,i)-exner_b(k,i))*rtheta_b(k,i))

 enddo
 enddo
 enddo

!updates the surface pressure and calculates the surface pressure tendency:
 do j = jts,jte
 do i = its,ite
    tem1 = zgrid(2,i)-zgrid(1,i)
    tem2 = zgrid(3,i)-zgrid(2,i)
    rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv_p(i,1,j))
    rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv_p(i,2,j))

    tend_sfc_pressure(i) = surface_pressure(i)
!   surface_pressure(i)  = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
!                        * (rho1 + 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i)  = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
                         * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i)  = surface_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
    tend_sfc_pressure(i)   = (surface_pressure(i)-tend_sfc_pressure(i)) / dt_dyn
 enddo
 enddo

!variables specific to different cloud microphysics schemes:
 microp_select_init: select case(microp_scheme)

    case ("mp_thompson","mp_wsm6")
       do j = jts, jte
       do k = kts, kte
       do i = its, ite
          qi(k,i) = qi_p(i,k,j)
          qs(k,i) = qs_p(i,k,j)
          qg(k,i) = qg_p(i,k,j)
       enddo
       enddo
       enddo

    microp2_select: select case(microp_scheme)

       case("mp_thompson")
          do j = jts, jte
          do k = kts, kte
          do i = its, ite
             ni(k,i) = ni_p(i,k,j)
             nr(k,i) = nr_p(i,k,j)
             rainprod(k,i) = rainprod_p(i,k,j)
             evapprod(k,i) = evapprod_p(i,k,j)
             re_cloud(k,i) = recloud_p(i,k,j)
             re_ice(k,i)   = reice_p(i,k,j)
             re_snow(k,i)  = resnow_p(i,k,j)
          enddo
          enddo
          enddo

       case default

    end select microp2_select

    case default

 end select microp_select_init

 end subroutine microphysics_to_MPAS

!=================================================================================================================
 end module mpas_atmphys_interface
!=================================================================================================================
