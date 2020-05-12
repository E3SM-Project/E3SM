module sat_mod
  implicit none

contains

  ! Saturation vapor pressure and mixing ratio.
  ! Based on Flatau et.al, (JAM, 1992:1507) - valid for T > -80C
  ! sat. vapor over ice below -80C - used Murphy and Koop (2005)
  ! For water below -80C simply assumed esw/esi = 2.
  ! des/dT below -80C computed as a finite difference of es

  real(crm_rknd) function esatw_crm(t)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd) t  ! temperature (K)
    real(crm_rknd) a0,a1,a2,a3,a4,a5,a6,a7,a8
    data  a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
    6.105851, 0.4440316, 0.1430341e-1, &
    0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
    0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
    !       6.11239921, 0.443987641, 0.142986287e-1, &
    !       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
    !       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
    real(crm_rknd) dt
    dt = t-273.16
    if(dt.gt.-80.) then
      esatw_crm = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else
      esatw_crm = 2.*0.01*exp(9.550426 - 5723.265/t + 3.53068*Log(t) - 0.00728332*t)
    end if
  end



  real(crm_rknd) function qsatw_crm(t,p)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd) t  ! temperature (K)
    real(crm_rknd) p  ! pressure    (mb)
    real(crm_rknd) esat_crm
    esat_crm = esatw_crm(t)
    qsatw_crm = 0.622 * esat_crm/max(esat_crm,p-esat_crm)
  end


  real(crm_rknd) function dtesatw_crm(t)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd) t  ! temperature (K)
    real(crm_rknd) a0,a1,a2,a3,a4,a5,a6,a7,a8
    data  a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
    0.443956472, 0.285976452e-1, 0.794747212e-3, &
    0.121167162e-4, 0.103167413e-6, 0.385208005e-9, &
    -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
    real(crm_rknd) dt
    dt = t-273.16
    if(dt.gt.-80.) then
      dtesatw_crm = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else
      dtesatw_crm = esatw_crm(t+1)-esatw_crm(t)
    end if

  end


  real(crm_rknd) function dtqsatw_crm(t,p)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd) t  ! temperature (K)
    real(crm_rknd) p  ! pressure    (mb)
    dtqsatw_crm = 0.622*dtesatw_crm(t)/p
  end


  real(crm_rknd) function esati_crm(t)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd) t  ! temperature (K)
    real(crm_rknd) a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
    6.11147274, 0.503160820, 0.188439774e-1, &
    0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
    0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
    real(crm_rknd) dt
    dt = t-273.16
    if(dt.gt.-80.) then
      esati_crm = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else
      esati_crm = 0.01*exp(9.550426 - 5723.265/t + 3.53068*Log(t) - 0.00728332*t)
    end if
  end



  real(crm_rknd) function qsati_crm(t,p)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd) t  ! temperature (K)
    real(crm_rknd) p  ! pressure    (mb)
    real(crm_rknd) esat_crm
    esat_crm=esati_crm(t)
    qsati_crm=0.622 * esat_crm/max(esat_crm,p-esat_crm)
  end


  real(crm_rknd) function dtesati_crm(t)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd) t  ! temperature (K)
    real(crm_rknd) a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 / &
    0.503223089, 0.377174432e-1,0.126710138e-2, &
    0.249065913e-4, 0.312668753e-6, 0.255653718e-8, &
    0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
    real(crm_rknd) dt
    dt = t-273.16
    if(dt.gt.-80.) then
      dtesati_crm = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else
      dtesati_crm = esati_crm(t+1.)-esati_crm(t)
    end if
  end


  real(crm_rknd) function dtqsati_crm(t,p)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd) t  ! temperature (K)
    real(crm_rknd) p  ! pressure    (mb)
    dtqsati_crm=0.622*dtesati_crm(t)/p
  end

end module sat_mod
