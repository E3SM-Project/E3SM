! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_interface_nhyd
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_vars

 implicit none
 private
 public:: allocate_forall_physics,   &
          deallocate_forall_physics, &
          MPAS_to_physics,           &
          microphysics_from_MPAS,    &
          microphysics_to_MPAS

 integer:: i,j,k


!>\brief interface for conversion between variables used in the MPAS dynamical core and variables
!> needed in the physics parameterizations.  
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_interface_nhyd:
!> -------------------------------------------
!> allocate_forall_physics  : allocation of all meteorological variables needed in the physics.
!> deallocate_forall_physics: deallocation of all meteorological variables needed in the physics.
!> MPAS_to_physics          : conversion of input "MPAS" variables to "physics" variables.
!> microphysics_from_MPAS   : initialize local arrays needed in cloud microphysics schemes.
!> microphysics_to_MPAS     : copy local arrays needed in cloud microphysics schemesto MPAS arrays.


 contains

!==================================================================================================
 subroutine allocate_forall_physics
!==================================================================================================

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

!... arrays used for calculating the hydrostatic pressure and exner function:
 if(.not.allocated(psfc_hyd_p)  ) allocate(psfc_hyd_p(ims:ime,jms:jme)          )
 if(.not.allocated(psfc_hydd_p) ) allocate(psfc_hydd_p(ims:ime,jms:jme)         )
 if(.not.allocated(pres_hyd_p)  ) allocate(pres_hyd_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pres_hydd_p) ) allocate(pres_hydd_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(pres2_hyd_p) ) allocate(pres2_hyd_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(pres2_hydd_p)) allocate(pres2_hydd_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(znu_hyd_p)   ) allocate(znu_hyd_p(ims:ime,kms:kme,jms:jme)   )
 
 end subroutine allocate_forall_physics

!==================================================================================================
 subroutine deallocate_forall_physics
!==================================================================================================

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

 if(allocated(psfc_hyd_p)  ) deallocate(psfc_hyd_p  )
 if(allocated(psfc_hydd_p) ) deallocate(psfc_hydd_p )
 if(allocated(pres_hyd_p)  ) deallocate(pres_hyd_p  )
 if(allocated(pres_hydd_p) ) deallocate(pres_hydd_p )
 if(allocated(pres2_hyd_p) ) deallocate(pres2_hyd_p )
 if(allocated(pres2_hydd_p)) deallocate(pres2_hydd_p)
 if(allocated(znu_hyd_p)   ) deallocate(znu_hyd_p   )
 
 end subroutine deallocate_forall_physics

!==================================================================================================
 subroutine MPAS_to_physics(mesh,state,diag,diag_physics)
!==================================================================================================

!input variables:
 type(mesh_type) ,intent(in):: mesh
 type(state_type),intent(in):: state
 type(diag_type) ,intent(in):: diag

!inout variables:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: i,k,j
 real(kind=RKIND):: z0,z1,z2,w1,w2

 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer  :: fzm,fzp,rdzw
 real(kind=RKIND),dimension(:),pointer  :: sfc_pressure
 real(kind=RKIND),dimension(:,:),pointer:: zgrid
 real(kind=RKIND),dimension(:,:),pointer:: zz,exner,pressure_b,rtheta_p,rtheta_b
 real(kind=RKIND),dimension(:,:),pointer:: rho_zz,theta_m,qv,pressure_p,u,v,w
 real(kind=RKIND),dimension(:,:),pointer:: qvs,rh

 real(kind=RKIND):: rho_a,rho1,rho2,tem1,tem2

!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine MPAS_to_phys:'
 write(0,*) 'ims=',ims,' ime=',ime
 write(0,*) 'jms=',jms,' jme=',jme
 write(0,*) 'kms=',kms,' kme=',kme
 write(0,*)
 write(0,*) 'its=',its,' ite=',ite
 write(0,*) 'jts=',jts,' jte=',jte
 write(0,*) 'kts=',kts,' kte=',kte

!initialization:
 latCell   => mesh % latCell % array
 lonCell   => mesh % lonCell % array

 fzm          => mesh % fzm % array
 fzp          => mesh % fzp % array
 rdzw         => mesh % rdzw % array
 zgrid        => mesh % zgrid % array
 zz           => mesh % zz % array
 sfc_pressure => diag % surface_pressure % array
 exner        => diag % exner % array
 pressure_b   => diag % pressure_base % array
 pressure_p   => diag % pressure_p % array
 rtheta_p     => diag % rtheta_p % array
 rtheta_b     => diag % rtheta_base % array

 rho_zz  => state % rho_zz % array
 theta_m => state % theta_m % array
 qv      => state % scalars % array(state%index_qv,:,:)

 w => state % w % array
 u => diag  % uReconstructZonal % array
 v => diag  % uReconstructMeridional % array

!calculation of the surface pressure using hydrostatic assumption down to the surface::
 do i = its,ite
    tem1 = zgrid(2,i)-zgrid(1,i)
    tem2 = zgrid(3,i)-zgrid(2,i)
    rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv(1,i))
    rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv(2,i))
    sfc_pressure(i) = 0.5*g*(zgrid(2,i)-zgrid(1,i)) &
                    * (rho1 + 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    sfc_pressure(i) = sfc_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
 enddo

!arrays located at theta points:
!do j = jts, jte
!do i = its, ite
!   psfc_p(i,j) = diag % surface_pressure % array(i)
!enddo
!enddo

 do j = jts, jte
 do k = kts, kte
 do i = its, ite

    !moist arrays:
    qv_p(i,k,j) = max(0.,state % scalars % array(state%index_qv,k,i))
    qc_p(i,k,j) = max(0.,state % scalars % array(state%index_qc,k,i))
    qr_p(i,k,j) = max(0.,state % scalars % array(state%index_qr,k,i))
    qi_p(i,k,j) = max(0.,state % scalars % array(state%index_qi,k,i))
    qs_p(i,k,j) = max(0.,state % scalars % array(state%index_qs,k,i))
    qg_p(i,k,j) = max(0.,state % scalars % array(state%index_qg,k,i))
    
    !arrays located at theta points:
    u_p(i,k,j) = u(k,i)
    v_p(i,k,j) = v(k,i)

    zz_p(i,k,j)  = zz(k,i)
    rho_p(i,k,j) = zz(k,i) * rho_zz(k,i)
    rho_p(i,k,j) = rho_p(i,k,j)*(1.+qv_p(i,k,j))
    th_p(i,k,j)  = theta_m(k,i) / (1. + R_v/R_d * qv(k,i))
    t_p(i,k,j)   = theta_m(k,i) * exner(k,i) / (1. + R_v/R_d * qv(k,i))

    pi_p(i,k,j)   = exner(k,i)
    pres_p(i,k,j) = pressure_p(k,i) + pressure_b(k,i)
    znu_p(i,k,j)  = pres_p(i,k,j) / sfc_pressure(i)

    zmid_p(i,k,j) = 0.5*(zgrid(k+1,i)+zgrid(k,i))
    dz_p(i,k,j)   = zgrid(k+1,i)-zgrid(k,i)

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
!      stop
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
    psfc_p(i,j) = pres2_p(i,k,j)
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
       pres2_hyd_p(i,k,j)  = pres2_hyd_p(i,k+1,j)  + g*rho_p(i,k,j)*dz_p(i,k,j)
       pres2_hydd_p(i,k,j) = pres2_hydd_p(i,k+1,j) + g*rho_a*dz_p(i,k,j)
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
    diag_physics % plrad % array(i) = pres2_p(i,kte+1,j) 
 enddo
 enddo

!formats: 
 201 format(3i8,10(1x,e15.8))
 202 format(2i6,10(1x,e15.8))
 203 format(i6,10(1x,e15.8))

 end subroutine MPAS_to_physics

!==================================================================================================
 subroutine microphysics_from_MPAS(mesh,state,tend,diag)
!==================================================================================================

!input variables:
 type(state_type),intent(in):: state
 type(diag_type) ,intent(in):: diag
 type(mesh_type) ,intent(in):: mesh

 type(tend_type),intent(inout):: tend

!local variables:
 integer:: i,k,j

 real(kind=RKIND):: min_theta,min_thp,min_tp
 real(kind=RKIND):: max_theta,max_thp,max_tp
 real(kind=RKIND):: min_qv,min_qc,min_qr,min_qi,min_qs,min_qg,min_qnr,min_qni
 real(kind=RKIND):: max_qv,max_qc,max_qr,max_qi,max_qs,max_qg,max_qnr,max_qni

 real(kind=RKIND),dimension(:,:),pointer:: zgrid
 real(kind=RKIND),dimension(:,:),pointer:: zz,exner,pressure_b,rtheta_p,rtheta_b
 real(kind=RKIND),dimension(:,:),pointer:: rho_zz,theta_m,qv,rh,pressure_p
 real(kind=RKIND),dimension(:,:),pointer:: rt_diabatic_tend

!--------------------------------------------------------------------------------------------------

!initialization:
 write(0,*)
 write(0,*) '--- enter subroutine microphysics_from_MPAS:'

 zgrid => mesh % zgrid % array
 zz    => mesh % zz % array

 exner      => diag % exner % array
 pressure_b => diag % pressure_base % array
 pressure_p => diag % pressure_p % array
 rtheta_p   => diag % rtheta_p % array
 rtheta_b   => diag % rtheta_base % array
 
 rho_zz     => state % rho_zz % array
 theta_m    => state % theta_m % array
 rh         => diag % rh % array

 rt_diabatic_tend => tend % rt_diabatic_tend % array

 qv => state % scalars % array(state%index_qv,:,:)

 do k = kts,kte
 do i = its,ite
    rt_diabatic_tend(k,i) = 0.
 enddo
 enddo

!copy sounding variables from the geodesic grid to the wrf-physics grid:
 do j = jts, jte
 do k = kts, kte
 do i = its, ite
    rho_p(i,k,j)  = zz(k,i) * rho_zz(k,i)
    th_p(i,k,j)   = theta_m(k,i) / (1. + R_v/R_d * max(0.,qv(k,i)))

    pi_p(i,k,j)   = exner(k,i)
    pres_p(i,k,j) = pressure_b(k,i) + pressure_p(k,i)

    z_p(i,k,j)    = zgrid(k,i)
    dz_p(i,k,j)   = zgrid(k+1,i) - zgrid(k,i)

!   qv_p(i,k,j) = max(0.,state % scalars % array(state%index_qv,k,i))
!   qc_p(i,k,j) = max(0.,state % scalars % array(state%index_qc,k,i))
!   qr_p(i,k,j) = max(0.,state % scalars % array(state%index_qr,k,i))
    qv_p(i,k,j) = state % scalars % array(state%index_qv,k,i)
    qc_p(i,k,j) = state % scalars % array(state%index_qc,k,i)
    qr_p(i,k,j) = state % scalars % array(state%index_qr,k,i)
 enddo
 enddo
 enddo

!additional initialization as function of cloud microphysics scheme:
 microp_select_init: select case(microp_scheme)

    case ("wsm6")
       do j = jts, jte
       do k = kts, kte
       do i = its, ite
!         qi_p(i,k,j) = max(0.,vars % scalars % array(vars%index_qi,k,i))
!         qs_p(i,k,j) = max(0.,vars % scalars % array(vars%index_qs,k,i))
!         qg_p(i,k,j) = max(0.,vars % scalars % array(vars%index_qg,k,i))
          qi_p(i,k,j) = state % scalars % array(state%index_qi,k,i)
          qs_p(i,k,j) = state % scalars % array(state%index_qs,k,i)
          qg_p(i,k,j) = state % scalars % array(state%index_qg,k,i)
       enddo
       enddo
       enddo

    case default

 end select microp_select_init

!formats:
 201 format(3i6,10(1x,e15.8))

 end subroutine microphysics_from_MPAS

!==================================================================================================
 subroutine microphysics_to_MPAS(mesh,state,diag,tend,itimestep)
!==================================================================================================

!input variables:
 integer,intent(in):: itimestep
 type(mesh_type),intent(in):: mesh

!output variables:
 type(state_type),intent(inout):: state
 type(diag_type),intent(inout):: diag
 type(tend_type),intent(inout):: tend

 real(kind=RKIND):: min_theta,min_thp,min_tp
 real(kind=RKIND):: max_theta,max_thp,max_tp
 real(kind=RKIND):: min_qv,min_qc,min_qr,min_qi,min_qs,min_qg,min_qnr,min_qni
 real(kind=RKIND):: max_qv,max_qc,max_qr,max_qi,max_qs,max_qg,max_qnr,max_qni

!local variables:
 integer:: icount
 integer:: i,k,j
 real(kind=RKIND),dimension(:,:),pointer:: zz,exner,exner_b,pressure_b,rtheta_p,rtheta_b
 real(kind=RKIND),dimension(:,:),pointer:: rho_zz,theta_m,pressure_p
 real(kind=RKIND),dimension(:,:),pointer:: rt_diabatic_tend

!ldf(2011-11-12): surface pressure.
 real(kind=RKIND):: rho1,rho2,tem1,tem2
 real(kind=RKIND),dimension(:),pointer:: rdzw
 real(kind=RKIND),dimension(:),pointer:: sfc_pressure
 real(kind=RKIND),dimension(:,:),pointer:: zgrid
!ldf end.

!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine microphysics_to_MPAS:'

!initialization:
 zz         => mesh % zz % array
 zgrid      => mesh % zgrid % array
 exner      => diag % exner % array
 exner_b    => diag % exner_base % array
 pressure_b => diag % pressure_base % array
 pressure_p => diag % pressure_p % array
 rtheta_p   => diag % rtheta_p % array
 rtheta_b   => diag % rtheta_base % array

 rho_zz     => state % rho_zz % array
 theta_m    => state % theta_m % array

 rt_diabatic_tend => tend % rt_diabatic_tend % array

!ldf (2011-11-12): update surface pressure.
 rdzw => mesh % rdzw % array
 sfc_pressure => diag % surface_pressure % array
!ldf end.

!variables common to all cloud microphysics schemes:

 do j = jts, jte
 do k = kts, kte
 do i = its, ite

    !potential temperature and diabatic forcing:
    rt_diabatic_tend(k,i) = theta_m(k,i)
    theta_m(k,i) = th_p(i,k,j) * (1. + R_v/R_d * qv_p(i,k,j))
    rt_diabatic_tend(k,i) = (theta_m(k,i) - rt_diabatic_tend(k,i)) / dt_dyn
!   rt_diabatic_tend(k,i) = 0.

    !density-weigthed perturbation potential temperature:
    rtheta_p(k,i) = rho_zz(k,i) * theta_m(k,i) - rtheta_b(k,i)

    !exner function:
    exner(k,i) = (zz(k,i)*(R_d/P0)*(rtheta_p(k,i)+rtheta_b(k,i)))**rcv

    !pertubation pressure:
    pressure_p(k,i) = zz(k,i)*R_d*(exner(k,i)*rtheta_p(k,i) &
                  + (exner(k,i)-exner_b(k,i))*rtheta_b(k,i))

    !mass mixing ratios:
    state % scalars % array(state%index_qv,k,i) = qv_p(i,k,j)
    state % scalars % array(state%index_qc,k,i) = qc_p(i,k,j)
    state % scalars % array(state%index_qr,k,i) = qr_p(i,k,j)

 enddo
 enddo
 enddo

!updates the surface pressure.
 do j = jts,jte
 do i = its,ite
    tem1 = zgrid(2,i)-zgrid(1,i)
    tem2 = zgrid(3,i)-zgrid(2,i)
    rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv_p(i,1,j))
    rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv_P(i,2,j))
    sfc_pressure(i) = 0.5*g*(zgrid(2,i)-zgrid(1,i)) &
                    * (rho1 + 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    sfc_pressure(i) = sfc_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
 enddo
 enddo
!ldf end.

!variables specific to different cloud microphysics schemes:

 microp_select_init: select case(microp_scheme)

    case ("wsm6")

       do j = jts, jte
       do k = kts, kte
       do i = its, ite

          !mass mixing ratios:
          state % scalars % array(state%index_qi,k,i) = qi_p(i,k,j)
          state % scalars % array(state%index_qs,k,i) = qs_p(i,k,j)
          state % scalars % array(state%index_qg,k,i) = qg_p(i,k,j)

       enddo
       enddo
       enddo

    case default

 end select microp_select_init

!formats:
 201 format(3i6,10(1x,e15.8))

 end subroutine microphysics_to_MPAS

!==================================================================================================
 end module mpas_atmphys_interface_nhyd
!==================================================================================================
