module micro_p3_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from micro_p3 fortran to scream c++.
!

interface

  subroutine find_lookuptable_indices_1a_f(dumi,dumjj,dumii,dumzz,dum1,dum4,dum5,dum6,      &
       qitot,nitot,qirim,rhop) bind(C)
    use iso_c_binding

    ! arguments:
    integer(kind=c_int), intent(out) :: dumi,dumjj,dumii,dumzz
    real(kind=c_real),   intent(out) :: dum1,dum4,dum5,dum6
    real(kind=c_real),   value, intent(in)  :: qitot,nitot,qirim,rhop
  end subroutine find_lookuptable_indices_1a_f

  subroutine find_lookuptable_indices_1b_f(dumj,dum3,qr,nr) bind(C)
    use iso_c_binding

    integer(kind=c_int), intent(out) :: dumj
    real(kind=c_real),   intent(out) :: dum3
    real(kind=c_real),   value, intent(in) :: qr, nr
  end subroutine find_lookuptable_indices_1b_f

  subroutine access_lookup_table_f(dumjj,dumii,dumi,index,dum1,dum4,dum5,proc) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: dumjj, dumii, dumi, index
    real(kind=c_real),   value, intent(in) :: dum1, dum4, dum5
    real(kind=c_real),   intent(out) :: proc
  end subroutine access_lookup_table_f

  subroutine access_lookup_table_coll_f(dumjj,dumii,dumj,dumi,index,dum1,dum3,dum4,dum5,proc) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: dumjj,dumii,dumj,dumi,index
    real(kind=c_real),   value, intent(in) :: dum1,dum3,dum4,dum5
    real(kind=c_real),   intent(out) :: proc
  end subroutine access_lookup_table_coll_f

  subroutine back_to_cell_average_f(lcldm,rcldm,icldm, qcacc,qrevp,qcaut,&
    ncacc,ncslf,ncautc,nrslf,nrevp,ncautr,qisub,nrshdr,qcheti,&
    qrcol,qcshd,qimlt,qccol,qrheti,nimlt,nccol,ncshdc,ncheti,nrcol,nislf,&
    qidep,nrheti,nisub,qinuc,ninuc,qiberg) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: lcldm, rcldm, icldm
    real(kind=c_real), intent(inout) :: qcacc, qrevp, qcaut, ncacc, ncslf, ncautc,  &
                                        nrslf, nrevp, ncautr, qisub,  &
                                        nrshdr, qcheti, qrcol, qcshd, qimlt, qccol, &
                                        qrheti, nimlt, nccol, ncshdc, ncheti, nrcol,&
                                        nislf, qidep, nrheti, nisub, qinuc, ninuc,  &
                                        qiberg
  end subroutine back_to_cell_average_f

  subroutine prevent_ice_overdepletion_f(pres,t,qv,xxls,odt,    &
     qidep,qisub) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: pres, t, qv, xxls, odt
    real(kind=c_real), intent(inout) :: qidep, qisub

  end subroutine prevent_ice_overdepletion_f

  subroutine cloud_water_conservation_f(qc,dt,qcaut,qcacc,qccol,qcheti,qcshd,     &
    qiberg,qisub,qidep) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: qc, dt
    real(kind=c_real), intent(inout) :: qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep
  end subroutine cloud_water_conservation_f

  subroutine rain_water_conservation_f(qr,qcaut,qcacc,qimlt,qcshd,dt,    &
    qrevp,qrcol,qrheti) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: qr, qcaut, qcacc, qimlt, qcshd, dt
    real(kind=c_real), intent(inout) :: qrevp, qrcol, qrheti
  end subroutine rain_water_conservation_f

  subroutine ice_water_conservation_f(qitot,qidep,qinuc,qiberg,qrcol,qccol,qrheti,qcheti,dt,    &
    qisub,qimlt) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: qitot, qidep, qinuc, qrcol, qccol, qrheti, qcheti, qiberg, dt
    real(kind=c_real), intent(inout) :: qisub, qimlt

  end subroutine ice_water_conservation_f

  subroutine get_cloud_dsd2_f(qc,nc,mu_c,rho,nu,lamc,cdist,cdist1,lcldm) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)         :: qc,rho,lcldm
    real(kind=c_real), intent(inout)             :: nc
    real(kind=c_real), intent(out)               :: mu_c,nu,lamc,cdist,cdist1
  end subroutine get_cloud_dsd2_f

  subroutine get_rain_dsd2_f(qr,nr,mu_r,lamr,cdistr,logn0r,rcldm) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: qr,rcldm
    real(kind=c_real), intent(inout)     :: nr
    real(kind=c_real), intent(out)       :: lamr,mu_r,cdistr,logn0r
  end subroutine get_rain_dsd2_f

  subroutine calc_rime_density_f(t,rhofaci,f1pr02,acn,lamc,mu_c,qc_incld,qccol,vtrmi1,rhorime_c) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: t, rhofaci, f1pr02, acn, lamc, mu_c, qc_incld, qccol
    real(kind=c_real), intent(out) :: vtrmi1, rhorime_c
  end subroutine calc_rime_density_f

  subroutine cldliq_immersion_freezing_f(t,lamc,mu_c,cdist1,qc_incld,qc_relvar,qcheti,ncheti) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: t, lamc, mu_c, cdist1, qc_incld, qc_relvar
    real(kind=c_real), intent(out) :: qcheti, ncheti
  end subroutine cldliq_immersion_freezing_f

  subroutine rain_immersion_freezing_f(t,lamr,mu_r,cdistr,qr_incld,qrheti,nrheti) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: t, lamr, mu_r, cdistr, qr_incld
    real(kind=c_real), intent(out) :: qrheti, nrheti
  end subroutine rain_immersion_freezing_f

  subroutine droplet_self_collection_f(rho,inv_rho,qc_incld,mu_c,nu,ncautc,ncslf) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho, inv_rho, qc_incld, mu_c, nu, ncautc
    real(kind=c_real), intent(out) :: ncslf
  end subroutine droplet_self_collection_f

  subroutine cloud_rain_accretion_f(rho,inv_rho,qc_incld,nc_incld,qr_incld,qc_relvar,qcacc,ncacc) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho, inv_rho, qc_incld, nc_incld, qr_incld, qc_relvar
    real(kind=c_real), intent(out) :: qcacc, ncacc
  end subroutine cloud_rain_accretion_f

  subroutine cloud_water_autoconversion_f(rho, qc_incld, nc_incld, qc_relvar, qcaut, ncautc, ncautr) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho, qc_incld, nc_incld, qc_relvar
    real(kind=c_real), intent(inout) :: qcaut, ncautc, ncautr
  end subroutine cloud_water_autoconversion_f

  subroutine rain_self_collection_f(rho, qr_incld, nr_incld, nrslf) bind(C)
    use iso_c_binding

    !arguments;
    real(kind=c_real), value, intent(in) :: rho, qr_incld, nr_incld
    real(kind=c_real), intent(out) :: nrslf

  end subroutine rain_self_collection_f

  subroutine ice_melting_f(rho,t,pres,rhofaci,f1pr05,f1pr14,xxlv,xlf,dv,sc,mu,kap,qv,qitot_incld,nitot_incld,qimlt,nimlt) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho,t,pres,rhofaci,f1pr05,f1pr14,xxlv,xlf,dv,sc,mu,kap,qv,qitot_incld,nitot_incld
    real(kind=c_real), intent(out) :: qimlt,nimlt
  end subroutine ice_melting_f

  subroutine impose_max_total_ni_f(nitot_local, max_total_Ni, inv_rho_local) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), intent(inout) :: nitot_local
    real(kind=c_real), value, intent(in) :: max_total_Ni, inv_rho_local
  end subroutine impose_max_total_ni_f

  subroutine calc_first_order_upwind_step_f(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dzq, num_arrays, fluxes, vs, qnx) bind(C)
    use iso_c_binding

    !arguments:
    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, kbot, k_qxtop, num_arrays
    real(kind=c_real), value, intent(in) :: dt_sub
    real(kind=c_real), dimension(kts:kte), intent(in) :: rho, inv_rho, inv_dzq
    type(c_ptr), intent(in), dimension(num_arrays) :: fluxes, vs, qnx
  end subroutine calc_first_order_upwind_step_f

  subroutine generalized_sedimentation_f(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dzq, inv_rho, rho, num_arrays, vs, fluxes, qnx) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, k_qxtop, kbot, num_arrays
    integer(kind=c_int), intent(inout) :: k_qxbot
    real(kind=c_real), value, intent(in) :: Co_max
    real(kind=c_real), intent(inout) :: dt_left, prt_accum
    real(kind=c_real), dimension(kts:kte), intent(in) :: inv_dzq, inv_rho, rho

    type(c_ptr), intent(in), dimension(num_arrays) :: vs, fluxes, qnx
  end subroutine generalized_sedimentation_f

  subroutine cloud_sedimentation_f(kts,kte,ktop,kbot,kdir,   &
       qc_incld,rho,inv_rho,lcldm,acn,inv_dzq,&
       dt,odt,log_predictNc, &
       qc, nc, nc_incld,mu_c,lamc,prt_liq,qc_tend,nc_tend) bind(C)

    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: qc_incld
    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: lcldm
    real(kind=c_real), intent(in), dimension(kts:kte) :: acn
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dzq

    real(kind=c_real),    value, intent(in) :: dt
    real(kind=c_real),    value, intent(in) :: odt
    logical(kind=c_bool), value, intent(in) :: log_predictNc

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: mu_c
    real(kind=c_real), intent(inout), dimension(kts:kte) :: lamc
    real(kind=c_real), intent(inout) :: prt_liq
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc_tend
  end subroutine cloud_sedimentation_f

  subroutine ice_sedimentation_f(kts,kte,ktop,kbot,kdir,    &
       rho,inv_rho,rhofaci,icldm,inv_dzq,dt,odt,  &
       qitot,qitot_incld,nitot,qirim,qirim_incld,birim,birim_incld,nitot_incld,prt_sol,qi_tend,ni_tend) bind(C)

    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: rhofaci
    real(kind=c_real), intent(in), dimension(kts:kte) :: icldm
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dzq
    real(kind=c_real), value, intent(in) :: dt, odt

    real(kind=c_real), intent(inout), dimension(kts:kte), target :: qitot
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qitot_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: nitot
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nitot_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: qirim
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qirim_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: birim
    real(kind=c_real), intent(inout), dimension(kts:kte) :: birim_incld

    real(kind=c_real), intent(inout) :: prt_sol
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qi_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: ni_tend
  end subroutine ice_sedimentation_f

  subroutine rain_sedimentation_f(kts,kte,ktop,kbot,kdir,   &
       qr_incld,rho,inv_rho,rhofacr,rcldm,inv_dzq,dt,odt,  &
       qr,nr,nr_incld,mu_r,lamr,prt_liq,rflx,qr_tend,nr_tend) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: qr_incld

    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: rhofacr
    real(kind=c_real), intent(in), dimension(kts:kte) :: rcldm
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dzq
    real(kind=c_real), value, intent(in) :: dt, odt

    real(kind=c_real), intent(inout), target, dimension(kts:kte) :: qr
    real(kind=c_real), intent(inout), target, dimension(kts:kte) :: nr
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: mu_r
    real(kind=c_real), intent(inout), dimension(kts:kte) :: lamr
    real(kind=c_real), intent(inout) :: prt_liq
    real(kind=c_real), intent(inout), dimension(kts:kte+1) :: rflx
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qr_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr_tend

  end subroutine rain_sedimentation_f

  subroutine calc_bulk_rho_rime_f(qi_tot, qi_rim, bi_rim, rho_rime) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real),   value, intent(in)  :: qi_tot
    real(kind=c_real),   intent(inout) :: qi_rim, bi_rim
    real(kind=c_real),   intent(out) :: rho_rime
  end subroutine calc_bulk_rho_rime_f

  subroutine homogeneous_freezing_f(kts,kte,ktop,kbot,kdir,t,exner,xlf,    &
   qc,nc,qr,nr,qitot,nitot,qirim,birim,th) bind(C)
    use iso_c_binding

    ! arguments:
    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir
    real(kind=c_real), intent(in), dimension(kts:kte) :: t
    real(kind=c_real), intent(in), dimension(kts:kte) :: exner
    real(kind=c_real), intent(in), dimension(kts:kte) :: xlf

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qr
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qitot
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nitot
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qirim
    real(kind=c_real), intent(inout), dimension(kts:kte) :: birim
    real(kind=c_real), intent(inout), dimension(kts:kte) :: th
  end subroutine homogeneous_freezing_f

  subroutine compute_rain_fall_velocity_f(qr_incld, rcldm, rhofacr, nr, nr_incld, mu_r, lamr, V_qr, V_nr) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: qr_incld, rcldm, rhofacr
    real(kind=c_real), intent(inout) :: nr, nr_incld
    real(kind=c_real), intent(out) :: mu_r, lamr, V_qr, V_nr
  end subroutine compute_rain_fall_velocity_f

  subroutine get_time_space_phys_variables_f(t, pres, rho, xxlv, xxls, qvs, qvi, mu, dv, sc, dqsdt, dqsidt, &
       ab, abi, kap, eii) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: t, pres, rho, xxlv, xxls, qvs, qvi
    real(kind=c_real), intent(out) :: mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii
  end subroutine get_time_space_phys_variables_f

subroutine  update_prognostic_ice_f(qcheti,qccol,qcshd,nccol,ncheti,ncshdc,qrcol,nrcol,qrheti,nrheti,nrshdr, &
       qimlt,nimlt,qisub,qidep,qinuc,ninuc,nislf,nisub,qiberg,exner,xxls,xlf,log_predictNc,log_wetgrowth, &
       dt,nmltratio,rhorime_c,th,qv,qitot,nitot,qirim,birim,qc,nc,qr,nr) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: qcheti, qccol, qcshd, nccol, ncheti, ncshdc, qrcol, nrcol, &
         qrheti, nrheti, nrshdr, qimlt, nimlt, qisub, qidep, qinuc, ninuc, nislf, nisub, qiberg, exner, &
         xlf, xxls, dt, nmltratio, rhorime_c

    logical(kind=c_bool), value, intent(in) :: log_predictNc
    logical(kind=c_bool), value, intent(in) :: log_wetgrowth

    real(kind=c_real), intent(inout) :: th, qv, qc, nc, qr, nr, qitot, nitot, qirim, birim

  end subroutine update_prognostic_ice_f

  subroutine ice_cldliq_collection_f(rho, t, rhofaci, f1pr04, qitot_incld, qc_incld, nitot_incld, &
                                     nc_incld, qccol, nccol, qcshd, ncshdc) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, t, rhofaci, f1pr04
    real(kind=c_real), value, intent(in) :: qitot_incld, qc_incld, nitot_incld, nc_incld
    real(kind=c_real), intent(out) :: qccol, nccol, qcshd, ncshdc
  end subroutine ice_cldliq_collection_f

  subroutine ice_rain_collection_f(rho, t, rhofaci, logn0r, f1pr07, f1pr08, &
                                   qitot_incld, nitot_incld, qr_incld, qrcol, nrcol) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, t, rhofaci, logn0r, f1pr07, f1pr08
    real(kind=c_real), value, intent(in) :: qitot_incld, nitot_incld, qr_incld
    real(kind=c_real), intent(out) :: qrcol, nrcol
  end subroutine ice_rain_collection_f

  subroutine ice_self_collection_f(rho, rhofaci, f1pr03, eii, qirim_incld, qitot_incld, nitot_incld, nislf) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, rhofaci, f1pr03, eii, qirim_incld, qitot_incld, nitot_incld
    real(kind=c_real), intent(out) :: nislf
  end subroutine ice_self_collection_f

  subroutine evaporate_sublimate_precip_f(qr_incld, qc_incld, nr_incld, qitot_incld,  lcldm, rcldm, qvs, ab, &
       epsr, qv, qrevp, nrevp) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: qr_incld, qc_incld, nr_incld, qitot_incld,  lcldm, rcldm, qvs, ab, &
         epsr, qv
    real(kind=c_real), intent(out) :: qrevp, nrevp

  end subroutine evaporate_sublimate_precip_f

  subroutine update_prognostic_liquid_f(qcacc, ncacc, qcaut,ncautc, ncautr, ncslf, &
       qrevp, nrevp, nrslf, log_predictNc, inv_rho, exner, xxlv, dt, th, qv, qc, nc, qr, nr) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: qcacc, ncacc, qcaut, ncautc, ncautr, ncslf, &
         qrevp, nrevp, nrslf

    logical(kind=c_bool), value, intent(in) :: log_predictNc

    real(kind=c_real), value, intent(in) :: inv_rho, exner, xxlv, dt

    real(kind=c_real), intent(inout) :: th, qv, qc, nc, qr, nr

  end subroutine update_prognostic_liquid_f

  subroutine ice_deposition_sublimation_f(qitot_incld, nitot_incld, t, qvs, qvi, epsi, abi, qv, &
       qidep, qisub, nisub, qiberg) bind(C)

    use iso_c_binding
    !arguments

    real(kind=c_real), value, intent(in) :: qitot_incld, nitot_incld, t, qvs, qvi, epsi, abi, qv
    real(kind=c_real), intent(out) :: qidep, qisub, nisub, qiberg

  end subroutine ice_deposition_sublimation_f

  subroutine ice_relaxation_timescale_f(rho, temp, rhofaci, f1pr05, f1pr14,   &
                                        dv, mu, sc, qitot_incld, nitot_incld, &
                                        epsi, epsi_tot) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: rho, temp, rhofaci, f1pr05, f1pr14, &
                                            dv, mu, sc, qitot_incld, nitot_incld
    real(kind=c_real), intent(out) :: epsi
    real(kind=c_real), intent(inout) :: epsi_tot
  end subroutine ice_relaxation_timescale_f

  subroutine calc_liq_relaxation_timescale_f(rho, f1r, f2r, dv, mu, sc, mu_r, &
                                             lamr, cdistr, cdist, qr_incld,   &
                                             qc_incld, epsr, epsc) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: rho,f1r,f2r,dv,mu,sc,mu_r,lamr, &
                                            cdistr,cdist,qr_incld,qc_incld
    real(kind=c_real), intent(out) :: epsr
    real(kind=c_real), intent(out) :: epsc

  end subroutine calc_liq_relaxation_timescale_f

  subroutine ice_nucleation_f(temp, inv_rho, nitot, naai, supi, odt, &
                              log_predictNc, qinuc, ninuc) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: temp, inv_rho, nitot, naai, supi, odt
    logical(kind=c_bool), value, intent(in) :: log_predictNc
    real(kind=c_real), intent(inout) :: qinuc, ninuc
 end subroutine ice_nucleation_f

 subroutine ice_cldliq_wet_growth_f(rho, temp, pres, rhofaci, f1pr05, f1pr14, xxlv, xlf, dv, kap, &
                                    mu, sc, qv, qc_incld,qitot_incld, nitot_incld, qr_incld,     &
                                    log_wetgrowth, qrcol, qccol, qwgrth, nrshdr, qcshd) bind(C)
   use iso_c_binding

   ! argmens
   real(kind=c_real), value, intent(in) :: rho, temp ,pres, rhofaci, f1pr05,f1pr14, xxlv, xlf, dv, &
                                           kap, mu, sc, qv, qc_incld, qitot_incld, nitot_incld,qr_incld
   logical(kind=c_bool), intent(inout) :: log_wetgrowth
   real(kind=c_real), intent(inout) :: qrcol, qccol, qwgrth, nrshdr, qcshd
 end subroutine ice_cldliq_wet_growth_f

 subroutine get_latent_heat_f(its, ite, kts, kte, v, s, f) bind(C)
   use iso_c_binding

   ! arguments
   integer(kind=c_int), value, intent(in) :: its, ite, kts, kte
   real(kind=c_real), dimension(its:ite, kts:kte), intent(out) :: v, s, f
 end subroutine get_latent_heat_f

 real(kind=c_real) function subgrid_variance_scaling_f(relvar,expon) bind(C)
   use iso_c_binding

   ! arguments
   real(kind=c_real), value, intent(in) :: relvar,expon
   ! return
   !real(kind=c_real) :: res

 end function subgrid_variance_scaling_f

 subroutine check_values_f(qv, temp, kts, kte, timestepcount, force_abort, source_ind, col_loc) bind(C)
   use iso_c_binding

   integer(kind=c_int), value, intent(in) :: kts, kte, timestepcount, source_ind
   logical(kind=c_bool), value, intent(in) :: force_abort

   ! arguments
   real(kind=c_real), intent(in) :: qv(kts:kte)
   real(kind=c_real), intent(in) :: temp(kts:kte)
   real(kind=c_real), intent(in) :: col_loc(3)

 end subroutine check_values_f

 subroutine calculate_incloud_mixingratios_f(qc, qr, qitot, qirim, nc, nr, nitot, birim,  &
                                             inv_lcldm, inv_icldm, inv_rcldm, &
                                             qc_incld, qr_incld, qitot_incld, qirim_incld, &
                                             nc_incld, nr_incld, nitot_incld, birim_incld) bind(C)
   use iso_c_binding

   ! argumens
   real(kind=c_real), value, intent(in) :: qc, qr, qitot, qirim, nc, nr, nitot, birim, inv_lcldm, inv_icldm, inv_rcldm
   real(kind=c_real), intent(out) :: qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld
 end subroutine calculate_incloud_mixingratios_f

 subroutine p3_main_pre_main_loop_f(kts, kte, kbot, ktop, kdir, log_predictNc, dt, &
       pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, xxlv, xxls, xlf, &
       t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn, qv, th, qc, nc, qr, nr, &
       qitot, nitot, qirim, birim, qc_incld, qr_incld, qitot_incld, qirim_incld, &
       nc_incld, nr_incld, nitot_incld, birim_incld, log_nucleationPossible, log_hydrometeorsPresent) bind(C)

   use iso_c_binding

   ! arguments
   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   logical(kind=c_bool), value, intent(in) :: log_predictNc
   real(kind=c_real), value, intent(in) :: dt

   real(kind=c_real), intent(in), dimension(kts:kte) :: pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, xxlv, xxls, xlf

   real(kind=c_real), intent(inout), dimension(kts:kte) :: t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, &
        acn, qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, qc_incld, qr_incld, qitot_incld, &
        qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld

   logical(kind=c_bool), intent(out) :: log_nucleationPossible, log_hydrometeorsPresent
 end subroutine p3_main_pre_main_loop_f

 subroutine p3_main_main_loop_f(kts, kte, kbot, ktop, kdir, log_predictNc, dt, odt, &
       pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, naai, qc_relvar, icldm, lcldm, rcldm,&
       t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn, qv, th, qc, nc, qr, nr, qitot, nitot, &
       qirim, birim, xxlv, xxls, xlf, qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, &
       nitot_incld, birim_incld, mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, cmeiout, prain, &
       nevapr, prer_evap, vap_liq_exchange, vap_ice_exchange, liq_ice_exchange, pratot, &
       prctot, log_hydrometeorsPresent) bind(C)

   use iso_c_binding

   !arguments
   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   logical(kind=c_bool), value, intent(in) :: log_predictNc
   real(kind=c_real), value, intent(in) :: dt, odt

   real(kind=c_real), intent(in), dimension(kts:kte) :: pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, &
        inv_rcldm, naai, qc_relvar, icldm, lcldm, rcldm

   real(kind=c_real), intent(inout), dimension(kts:kte) :: t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn, &
        qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, xlf, qc_incld, qr_incld, &
        qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld, mu_c, nu, lamc, cdist, cdist1, &
        cdistr, mu_r, lamr, logn0r, cmeiout, prain, nevapr, prer_evap, vap_liq_exchange, &
        vap_ice_exchange, liq_ice_exchange, pratot, prctot

   logical(kind=c_bool), intent(out) :: log_hydrometeorsPresent

 end subroutine p3_main_main_loop_f

 subroutine p3_main_post_main_loop_f(kts, kte, kbot, ktop, kdir, &
      exner, lcldm, rcldm, &
      rho, inv_rho, rhofaci, qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, &
      mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange, &
      ze_rain, ze_ice, diag_vmi, diag_effi, diag_di, diag_rhoi, diag_ze, diag_effc) bind(C)

   use iso_c_binding

   ! args

   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   real(kind=c_real), intent(in), dimension(kts:kte) :: exner, lcldm, rcldm
   real(kind=c_real), intent(inout), dimension(kts:kte) :: rho, inv_rho, rhofaci, &
        qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, &
        mu_c, nu, lamc, mu_r, &
        lamr, vap_liq_exchange, &
        ze_rain, ze_ice, diag_vmi, diag_effi, diag_di, diag_rhoi, diag_ze, diag_effc

 end subroutine p3_main_post_main_loop_f

  !
  ! These are some routine math operations that are not BFB between
  ! fortran and C++ on all platforms, so fortran will need to use
  ! the C++ versions in order to stay BFB.
  !

  function cxx_pow(base, exp) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base
    real(kind=c_real), value, intent(in)  :: exp

    ! return
    real(kind=c_real)               :: cxx_pow
  end function cxx_pow

  function cxx_sqrt(base) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base

    ! return
    real(kind=c_real)               :: cxx_sqrt
  end function cxx_sqrt

  function cxx_cbrt(base) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base

    ! return
    real(kind=c_real)               :: cxx_cbrt
  end function cxx_cbrt

  function cxx_gamma(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_gamma
  end function cxx_gamma

  function cxx_log(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_log
  end function cxx_log

  function cxx_log10(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_log10
  end function cxx_log10

  function cxx_exp(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_exp
  end function cxx_exp

end interface

end module micro_p3_iso_f
