/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImpl.hpp"
#include "compose_test.hpp"
#include "profiling.hpp"
#include "mpi/Comm.hpp"

namespace Homme {
using CTI = ComposeTransportImpl;

static void fill_ics (const ComposeTransportImpl& cti, const int n0_qdp, const int np1 = -1) {
  const auto pll = CTI::cmvdc(cti.m_geometry.m_sphere_latlon);
  const auto qdp = CTI::cmvdc(cti.m_tracers.qdp);
  const auto dp3d = Kokkos::create_mirror_view(cti.m_state.m_dp3d);
  if (np1 >= 0) Kokkos::deep_copy(dp3d, cti.m_state.m_dp3d);
  const auto f = [&] (int ie, int lev, int i, int j) {
    Real lat = pll(ie,i,j,0), lon = pll(ie,i,j,1);
    compose::test::offset_latlon(cti.num_phys_lev, lev, lat, lon);
    const int p = lev / cti.packn, s = lev % cti.packn;
    for (int q = 0; q < cti.m_data.qsize; ++q)
      compose::test::InitialCondition::init(compose::test::get_ic(cti.m_data.qsize, lev, q),
                                            1, &lat, &lon,
                                            &qdp(ie,n0_qdp,q,i,j,p)[s]);
    if (np1 >= 0) dp3d(ie,np1,i,j,p)[s] = 1;
  };
  cti.loop_host_ie_plev_ij(f);
  Kokkos::deep_copy(cti.m_tracers.qdp, qdp);
  if (np1 >= 0) Kokkos::deep_copy(cti.m_state.m_dp3d, dp3d);
}

static void cp_v_to_vstar (const ComposeTransportImpl& cti, const int np1) {
  const auto vstar = cti.m_derived.m_vstar;
  const auto v = cti.m_state.m_v;
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    CTI::idx_ie_packlev_ij(idx, ie, lev, i, j);
    for (int d = 0; d < 2; ++d)
      vstar(ie,d,i,j,lev) = v(ie,np1,d,i,j,lev);
  };
  cti.launch_ie_packlev_ij(f);
}

static void fill_v (const ComposeTransportImpl& cti, const Real t, const int np1,
                    const bool bfb) {
  const int packn = cti.packn;
  const int num_phys_lev = cti.num_phys_lev;
  const compose::test::NonDivergentWindField wf;
  if (bfb) { // On host b/c trig isn't BFB between host and device.
    const auto pll = CTI::cmvdc(cti.m_geometry.m_sphere_latlon);
    const auto v = CTI::cmvdc(cti.m_state.m_v);
    const auto f = [&] (const int ie, const int lev, const int i, const int j) {
      Real latlon[] = {pll(ie,i,j,0), pll(ie,i,j,1)};
      compose::test::offset_latlon(num_phys_lev, lev, latlon[0], latlon[1]);
      Real uv[2];
      wf.eval(t, latlon, uv);
      const int p = lev / packn, s = lev % packn;
      for (int d = 0; d < 2; ++d) v(ie,np1,d,i,j,p)[s] = uv[d];  
    };
    cti.loop_host_ie_plev_ij(f);
    Kokkos::deep_copy(cti.m_state.m_v, v);
  } else {
    const auto pll = cti.m_geometry.m_sphere_latlon;
    const auto v = cti.m_state.m_v;
    const auto f = KOKKOS_LAMBDA (const int idx) {
      int ie, lev, i, j;
      CTI::idx_ie_physlev_ij(idx, ie, lev, i, j);
      Real latlon[] = {pll(ie,i,j,0), pll(ie,i,j,1)};
      compose::test::offset_latlon(num_phys_lev, lev, latlon[0], latlon[1]);
      Real uv[2];
      wf.eval(t, latlon, uv);
      const int p = lev / packn, s = lev % packn;
      for (int d = 0; d < 2; ++d) v(ie,np1,d,i,j,p)[s] = uv[d];  
    };
    cti.launch_ie_physlev_ij(f);
  }
}

static void finish (const ComposeTransportImpl& cti, const Comm& comm,
                    const int n0_qdp, const int np1, std::vector<Real>& eval) {
  const int nelemd = cti.m_data.nelemd, qsize = cti.m_data.qsize, nlev = cti.num_phys_lev,
    nother_qdp = (n0_qdp + 1) % 2;
  fill_ics(cti, nother_qdp);
  const auto qdp = CTI::cmvdc(cti.m_tracers.qdp);
  const auto dp3d = CTI::cmvdc(cti.m_state.m_dp3d);
  const auto spheremp = CTI::cmvdc(cti.m_geometry.m_spheremp);
  Kokkos::View<Real***, HostMemSpace>
    l2_num("l2_num", nlev, qsize, nelemd), l2_den("l2_den", nlev, qsize, nelemd);
  Kokkos::View<Real**, HostMemSpace>
    mass0("mass0", qsize, nelemd), massf("massf", qsize, nelemd);
  for (int ie = 0; ie < nelemd; ++ie) {
    for (int q = 0; q < qsize; ++q) {
      Real m0 = 0, mf = 0;
      for (int k = 0; k < nlev; ++k) {
        const int p = k / cti.packn, s = k % cti.packn;
        Real num = 0, den = 0;
        for (int i = 0; i < cti.np; ++i)
          for (int j = 0; j < cti.np; ++j) {
            num += spheremp(ie,i,j)*square(qdp(ie,n0_qdp,q,i,j,p)[s]/dp3d(ie,np1,i,j,p)[s] -
                                           qdp(ie,nother_qdp,q,i,j,p)[s]);
            den += spheremp(ie,i,j)*square(qdp(ie,nother_qdp,q,i,j,p)[s]);
            m0 += spheremp(ie,i,j)*qdp(ie,nother_qdp,q,i,j,p)[s] /* times rho = 1 */;
            mf += spheremp(ie,i,j)*qdp(ie,n0_qdp,q,i,j,p)[s];
          }
        l2_num(k,q,ie) = num;
        l2_den(k,q,ie) = den;
      }
      mass0(q,ie) = m0;
      massf(q,ie) = mf;
    }
  }
  const auto mpi_comm = comm.mpi_comm();
  const auto am_root = comm.root();
  const auto fcomm = MPI_Comm_c2f(mpi_comm);
  eval.resize((nlev+1)*qsize);
  int cnt = 0;
  {
    Kokkos::View<Real**, HostMemSpace>
      l2_num_red("l2_num_red", nlev, qsize), l2_den_red("l2_den_red", nlev, qsize);
    const auto nr = nlev*qsize;
    compose_repro_sum(l2_num.data(), l2_num_red.data(), nelemd, nr, fcomm);
    compose_repro_sum(l2_den.data(), l2_den_red.data(), nelemd, nr, fcomm);
    if (am_root)
      for (int q = 0; q < qsize; ++q) {
        printf("COMPOSE (hxx)>");
        for (int k = 0; k < nlev; ++k) {
          const auto err = std::sqrt(l2_num_red(k,q)/l2_den_red(k,q));
          eval[cnt++] = err;
          printf("%23.16e", err);
        }
        printf("\n");
      }
  }
  {
    std::vector<Real> mass0_red(qsize), massf_red(qsize);
    compose_repro_sum(mass0.data(), mass0_red.data(), nelemd, qsize, fcomm);
    compose_repro_sum(massf.data(), massf_red.data(), nelemd, qsize, fcomm);
    if (am_root)
      for (int q = 0; q < qsize; ++q) {
        const auto err = (massf_red[q] - mass0_red[q])/mass0_red[q];
        eval[cnt++] = err;
        printf("COMPOSE (hxx)> mass0 %8.2e mass re %9.2e\n", mass0_red[q], err);
      }
  }
}

void ComposeTransportImpl::test_2d (const bool bfb, const int nstep, std::vector<Real>& eval) {
  SimulationParams& params = Context::singleton().get<SimulationParams>();
  params.qsplit = 1;

  TimeLevel& tl = Context::singleton().get<TimeLevel>();
  tl.nstep = 0;
  tl.update_tracers_levels(params.qsplit);
  const Real twelve_days = 3600 * 24 * 12, dt = twelve_days/nstep;

  fill_ics(*this, tl.n0_qdp, tl.np1);

  GPTLstart("compose_stt_step");
  for (int i = 0; i < nstep; ++i) {
    const auto tprev = dt*i;
    const auto t = dt*(i+1);
    if (i == 0) {
      fill_v(*this, tprev, tl.np1, bfb);
      Kokkos::fence();
    }
    cp_v_to_vstar(*this, tl.np1);
    Kokkos::fence();
    fill_v(*this, t, tl.np1, bfb);
    Kokkos::fence();
    run(tl, dt);
    Kokkos::fence();
    tl.nstep += params.qsplit;
    tl.update_tracers_levels(params.qsplit);
  }
  GPTLstop("compose_stt_step");

  finish(*this, Context::singleton().get<Comm>(), tl.n0_qdp, tl.np1, eval);
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
