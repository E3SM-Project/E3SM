/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Elements.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/TestUtils.hpp"
#include "HybridVCoord.hpp"

#include <limits>
#include <random>
#include <assert.h>

namespace Homme {

void Elements::init(const int num_elems, const bool consthv) {
  m_num_elems = num_elems;

  buffers.init(num_elems);

  m_fcor = ExecViewManaged<Real * [NP][NP]>("FCOR", m_num_elems);
  m_mp = ExecViewManaged<Real * [NP][NP]>("MP", m_num_elems);
  m_spheremp = ExecViewManaged<Real * [NP][NP]>("SPHEREMP", m_num_elems);
  m_rspheremp = ExecViewManaged<Real * [NP][NP]>("RSPHEREMP", m_num_elems);
  m_metinv = ExecViewManaged<Real * [2][2][NP][NP]>("METINV", m_num_elems);
  m_metdet = ExecViewManaged<Real * [NP][NP]>("METDET", m_num_elems);

  if(!consthv){
    m_tensorvisc = ExecViewManaged<Real * [2][2][NP][NP]>("TENSORVISC", m_num_elems);
    m_vec_sph2cart = ExecViewManaged<Real * [2][3][NP][NP]>("VEC_SPH2CART", m_num_elems);
  }

  m_phis = ExecViewManaged<Real * [NP][NP]>("PHIS", m_num_elems);

  //matrix D and its derivatives 
  m_d =
      ExecViewManaged<Real * [2][2][NP][NP]>("matrix D", m_num_elems);
  m_dinv = ExecViewManaged<Real * [2][2][NP][NP]>(
      "DInv - inverse of matrix D", m_num_elems);

  m_omega_p =
      ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Omega P", m_num_elems);
  m_phi = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("PHI", m_num_elems);

  m_derived_vn0 = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>(
      "Derived Lateral Velocities", m_num_elems);

  m_v = ExecViewManaged<Scalar * [NUM_TIME_LEVELS][2][NP][NP][NUM_LEV]>(
      "Horizontal Velocity", m_num_elems);
  m_t = ExecViewManaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]>(
      "Temperature", m_num_elems);
  m_dp3d = ExecViewManaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]>(
      "DP3D", m_num_elems);

  m_ps_v = ExecViewManaged<Real * [NUM_TIME_LEVELS][NP][NP]>("PS_V", m_num_elems);

  m_fm = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("F_Momentum", m_num_elems);
  m_ft = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("F_Temperature", m_num_elems);

  m_eta_dot_dpdn =
      ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("eta_dot_dpdn", m_num_elems);
  m_derived_dpdiss_ave = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
      "mean dp used to compute psdiss_tens", m_num_elems);

  m_derived_dp = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
    "derived_dp", m_num_elems);
  m_derived_divdp = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
    "derived_divdp", m_num_elems);
  m_derived_divdp_proj = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
    "derived_divdp_proj", m_num_elems);
  m_derived_dpdiss_biharmonic = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
    "derived_dpdiss_biharmonic", m_num_elems);
  m_derived_dpdiss_ave = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
    "derived_dpdiss_ave", m_num_elems);
}

void Elements::init_2d(const int ie, CF90Ptr &D, CF90Ptr &Dinv, CF90Ptr &fcor,
                       CF90Ptr &mp, CF90Ptr &spheremp, CF90Ptr &rspheremp,
                       CF90Ptr &metdet, CF90Ptr &metinv, CF90Ptr &phis,
                       CF90Ptr &tensorvisc, 
                       CF90Ptr &vec_sph2cart,
                       const bool consthv) {

  using ScalarView   = ExecViewUnmanaged<Real [NP][NP]>;
  using TensorView   = ExecViewUnmanaged<Real [2][2][NP][NP]>;
  using Tensor23View = ExecViewUnmanaged<Real [2][3][NP][NP]>;

  using ScalarViewF90   = HostViewUnmanaged<const Real [NP][NP]>;
  using TensorViewF90   = HostViewUnmanaged<const Real [2][2][NP][NP]>;
  using Tensor23ViewF90 = HostViewUnmanaged<const Real [2][3][NP][NP]>;

  ScalarView::HostMirror h_fcor      = Kokkos::create_mirror_view(Homme::subview(m_fcor,ie));
  ScalarView::HostMirror h_metdet    = Kokkos::create_mirror_view(Homme::subview(m_metdet,ie));
  ScalarView::HostMirror h_mp        = Kokkos::create_mirror_view(Homme::subview(m_mp,ie));
  ScalarView::HostMirror h_spheremp  = Kokkos::create_mirror_view(Homme::subview(m_spheremp,ie));
  ScalarView::HostMirror h_rspheremp = Kokkos::create_mirror_view(Homme::subview(m_rspheremp,ie));
  ScalarView::HostMirror h_phis      = Kokkos::create_mirror_view(Homme::subview(m_phis,ie));
  TensorView::HostMirror h_metinv    = Kokkos::create_mirror_view(Homme::subview(m_metinv,ie));
  TensorView::HostMirror h_d         = Kokkos::create_mirror_view(Homme::subview(m_d,ie));
  TensorView::HostMirror h_dinv      = Kokkos::create_mirror_view(Homme::subview(m_dinv,ie));

  TensorView::HostMirror h_tensorvisc;
  Tensor23View::HostMirror h_vec_sph2cart;
  if( !consthv ){
    h_tensorvisc   = Kokkos::create_mirror_view(Homme::subview(m_tensorvisc,ie));
    h_vec_sph2cart = Kokkos::create_mirror_view(Homme::subview(m_vec_sph2cart,ie));
  }

  ScalarViewF90 h_fcor_f90         (fcor);
  ScalarViewF90 h_metdet_f90       (metdet);
  ScalarViewF90 h_mp_f90           (mp);
  ScalarViewF90 h_spheremp_f90     (spheremp);
  ScalarViewF90 h_rspheremp_f90    (rspheremp);
  ScalarViewF90 h_phis_f90         (phis);
  TensorViewF90 h_metinv_f90       (metinv);
  TensorViewF90 h_d_f90            (D);
  TensorViewF90 h_dinv_f90         (Dinv);
  TensorViewF90 h_tensorvisc_f90   (tensorvisc);
  Tensor23ViewF90 h_vec_sph2cart_f90 (vec_sph2cart);
  
  // 2d scalars
  for (int igp = 0; igp < NP; ++igp) {
    for (int jgp = 0; jgp < NP; ++jgp) {
      h_fcor      (igp, jgp) = h_fcor_f90      (igp,jgp);
      h_mp        (igp, jgp) = h_mp_f90        (igp, jgp);
      h_spheremp  (igp, jgp) = h_spheremp_f90  (igp,jgp);
      h_rspheremp (igp, jgp) = h_rspheremp_f90 (igp,jgp);
      h_metdet    (igp, jgp) = h_metdet_f90    (igp,jgp);
      h_phis      (igp, jgp) = h_phis_f90      (igp,jgp);
    }
  }

  // 2d tensors
  for (int idim = 0; idim < 2; ++idim) {
    for (int jdim = 0; jdim < 2; ++jdim) {
      for (int igp = 0; igp < NP; ++igp) {
        for (int jgp = 0; jgp < NP; ++jgp) {
          h_d      (idim,jdim,igp,jgp) = h_d_f90      (idim,jdim,igp,jgp);
          h_dinv   (idim,jdim,igp,jgp) = h_dinv_f90   (idim,jdim,igp,jgp);
          h_metinv (idim,jdim,igp,jgp) = h_metinv_f90 (idim,jdim,igp,jgp);
        }
      }
    }
  }
  
  if(!consthv) {
    for (int idim = 0; idim < 2; ++idim) {
      for (int jdim = 0; jdim < 2; ++jdim) {
        for (int igp = 0; igp < NP; ++igp) {
          for (int jgp = 0; jgp < NP; ++jgp) {
            h_tensorvisc   (idim,jdim,igp,jgp) = h_tensorvisc_f90   (idim,jdim,igp,jgp);
          }
        }
      }
    }
    for (int idim = 0; idim < 2; ++idim) {
      for (int jdim = 0; jdim < 3; ++jdim) {
        for (int igp = 0; igp < NP; ++igp) {
          for (int jgp = 0; jgp < NP; ++jgp) {
            h_vec_sph2cart (idim,jdim,igp,jgp) = h_vec_sph2cart_f90 (idim,jdim,igp,jgp);
          }
        }
      }
    }
  }//end if consthv

  Kokkos::deep_copy(Homme::subview(m_fcor,ie), h_fcor);
  Kokkos::deep_copy(Homme::subview(m_metinv,ie), h_metinv);
  Kokkos::deep_copy(Homme::subview(m_metdet,ie), h_metdet);
  Kokkos::deep_copy(Homme::subview(m_mp,ie), h_mp);
  Kokkos::deep_copy(Homme::subview(m_spheremp,ie), h_spheremp);
  Kokkos::deep_copy(Homme::subview(m_rspheremp,ie), h_rspheremp);
  Kokkos::deep_copy(Homme::subview(m_phis,ie), h_phis);
  Kokkos::deep_copy(Homme::subview(m_d,ie), h_d);
  Kokkos::deep_copy(Homme::subview(m_dinv,ie), h_dinv);
  if( !consthv ) {
    Kokkos::deep_copy(Homme::subview(m_tensorvisc,ie), h_tensorvisc);
    Kokkos::deep_copy(Homme::subview(m_vec_sph2cart,ie), h_vec_sph2cart);
  }
}

//test for tensor hv is needed
void Elements::random_init(int num_elems, Real max_pressure) {
  std::random_device rd;
  HybridVCoord hv;
  hv.random_init(rd());
  random_init(num_elems,max_pressure,hv);
}

void Elements::random_init(int num_elems, Real max_pressure, const HybridVCoord& hvcoord) {
  // arbitrary minimum value to generate and minimum determinant allowed
  constexpr const Real min_value = 0.015625;
  // 1 is for const hv
  init(num_elems, 1);
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);

  genRandArray(m_fcor, engine, random_dist);
  genRandArray(m_mp, engine, random_dist);
  genRandArray(m_spheremp, engine, random_dist);
  genRandArray(m_rspheremp, engine, random_dist);
  genRandArray(m_metdet, engine, random_dist);
  genRandArray(m_metinv, engine, random_dist);
  genRandArray(m_phis, engine, random_dist);

  genRandArray(m_omega_p, engine, random_dist);
  genRandArray(m_phi, engine, random_dist);
  genRandArray(m_derived_vn0, engine, random_dist);

  genRandArray(m_v, engine, random_dist);
  genRandArray(m_t, engine, random_dist);

  // Generate ps_v so that it is >> ps0.
  // Note: make sure you init hvcoord before calling this method!
  genRandArray(m_ps_v, engine, std::uniform_real_distribution<Real>(100*hvcoord.ps0,1000*hvcoord.ps0));

  // This ensures the pressure in a single column is monotonically increasing
  // and has fixed upper and lower values
  const auto make_pressure_partition = [=](
      HostViewUnmanaged<Scalar[NUM_LEV]> pt_pressure) {
    // Put in monotonic order
    std::sort(
        reinterpret_cast<Real *>(pt_pressure.data()),
        reinterpret_cast<Real *>(pt_pressure.data() + pt_pressure.size()));
    // Ensure none of the values are repeated
    for (int level = NUM_PHYSICAL_LEV - 1; level > 0; --level) {
      const int prev_ilev = (level - 1) / VECTOR_SIZE;
      const int prev_vlev = (level - 1) % VECTOR_SIZE;
      const int cur_ilev = level / VECTOR_SIZE;
      const int cur_vlev = level % VECTOR_SIZE;
      // Need to try again if these are the same or if the thickness is too
      // small
      if (pt_pressure(cur_ilev)[cur_vlev] <=
          pt_pressure(prev_ilev)[prev_vlev] +
              min_value * std::numeric_limits<Real>::epsilon()) {
        return false;
      }
    }
    // We know the minimum thickness of a layer is min_value * epsilon
    // (due to floating point), so set the bottom layer thickness to that,
    // and subtract that from the top layer
    // This ensures that the total sum is max_pressure
    pt_pressure(0)[0] = min_value * std::numeric_limits<Real>::epsilon();
    const int top_ilev = (NUM_PHYSICAL_LEV - 1) / VECTOR_SIZE;
    const int top_vlev = (NUM_PHYSICAL_LEV - 1) % VECTOR_SIZE;
    // Note that this may not actually change the top level pressure
    // This is okay, because we only need to approximately sum to max_pressure
    pt_pressure(top_ilev)[top_vlev] = max_pressure - pt_pressure(0)[0];
    for (int e_vlev = top_vlev + 1; e_vlev < VECTOR_SIZE; ++e_vlev) {
      pt_pressure(top_ilev)[e_vlev] = std::numeric_limits<Real>::quiet_NaN();
    }
    // Now compute the interval thicknesses
    for (int level = NUM_PHYSICAL_LEV - 1; level > 0; --level) {
      const int prev_ilev = (level - 1) / VECTOR_SIZE;
      const int prev_vlev = (level - 1) % VECTOR_SIZE;
      const int cur_ilev = level / VECTOR_SIZE;
      const int cur_vlev = level % VECTOR_SIZE;
      pt_pressure(cur_ilev)[cur_vlev] -= pt_pressure(prev_ilev)[prev_vlev];
    }
    return true;
  };

  std::uniform_real_distribution<Real> pressure_pdf(min_value, max_pressure);

  // Lambdas used to constrain the metric tensor and its inverse
  const auto compute_det = [](HostViewUnmanaged<Real[2][2]> mtx) {
    return mtx(0, 0) * mtx(1, 1) - mtx(0, 1) * mtx(1, 0);
  };

  const auto constrain_det = [=](HostViewUnmanaged<Real[2][2]> mtx) {
    Real determinant = compute_det(mtx);
    // We want to ensure both the metric tensor and its inverse have reasonable
    // determinants
    if (determinant > min_value && determinant < 1.0 / min_value) {
      return true;
    } else {
      return false;
    }
  };

  // 2d tensors
  // Generating lots of matrices with reasonable determinants can be difficult
  // So instead of generating them all at once and verifying they're correct,
  // generate them one at a time, verifying them individually
  HostViewManaged<Real[2][2]> h_matrix("single host metric matrix");

  ExecViewManaged<Real *[2][2][NP][NP]>::HostMirror h_d =
      Kokkos::create_mirror_view(m_d);
  ExecViewManaged<Real *[2][2][NP][NP]>::HostMirror h_dinv =
      Kokkos::create_mirror_view(m_dinv);
  ExecViewManaged<Real *[NUM_TIME_LEVELS][NP][NP]>::HostMirror h_ps_v=
      Kokkos::create_mirror_view(m_ps_v);
  Kokkos::deep_copy(h_ps_v,m_ps_v);

  Real dp3d_min = std::numeric_limits<Real>::max();
  for (int ie = 0; ie < m_num_elems; ++ie) {
    // Because this constraint is difficult to satisfy for all of the tensors,
    // incrementally generate the view
    for (int igp = 0; igp < NP; ++igp) {
      for (int jgp = 0; jgp < NP; ++jgp) {
        for (int tl = 0; tl < NUM_TIME_LEVELS; ++tl) {
          ExecViewUnmanaged<Scalar[NUM_LEV]> pt_dp3d =
              Homme::subview(m_dp3d, ie, tl, igp, jgp);
          genRandArray(pt_dp3d, engine, pressure_pdf, make_pressure_partition);
          auto h_dp3d = Kokkos::create_mirror_view(pt_dp3d);
          Kokkos::deep_copy(h_dp3d,pt_dp3d);
          for (int ilev=0; ilev<NUM_LEV; ++ilev) {
            for (int iv=0; iv<VECTOR_SIZE; ++iv) {
              dp3d_min = std::min(dp3d_min,h_dp3d(ilev)[iv]);
            }
          }
        }
        genRandArray(h_matrix, engine, random_dist, constrain_det);
        for (int i = 0; i < 2; ++i) {
          for (int j = 0; j < 2; ++j) {
            h_d(ie, i, j, igp, jgp) = h_matrix(i, j);
          }
        }
        const Real determinant = compute_det(h_matrix);
        h_dinv(ie, 0, 0, igp, jgp) = h_matrix(1, 1) / determinant;
        h_dinv(ie, 1, 0, igp, jgp) = -h_matrix(1, 0) / determinant;
        h_dinv(ie, 0, 1, igp, jgp) = -h_matrix(0, 1) / determinant;
        h_dinv(ie, 1, 1, igp, jgp) = h_matrix(0, 0) / determinant;
      }
    }
  }

  // Generate eta_dot_dpdn so that it is << dp3d
  genRandArray(m_eta_dot_dpdn, engine, std::uniform_real_distribution<Real>(0.01*dp3d_min,0.1*dp3d_min));

  Kokkos::deep_copy(m_d, h_d);
  Kokkos::deep_copy(m_dinv, h_dinv);
  return;
}

void Elements::pull_from_f90_pointers(
    CF90Ptr &state_v, CF90Ptr &state_t, CF90Ptr &state_dp3d,
    CF90Ptr &derived_phi, CF90Ptr &derived_omega_p,
    CF90Ptr &derived_v, CF90Ptr &derived_eta_dot_dpdn) {
  pull_3d(derived_phi, derived_omega_p, derived_v);
  pull_4d(state_v, state_t, state_dp3d);
  pull_eta_dot(derived_eta_dot_dpdn);
}

void Elements::pull_3d(CF90Ptr &derived_phi, CF90Ptr &derived_omega_p, CF90Ptr &derived_v) {
  HostViewUnmanaged<const Real *[NUM_PHYSICAL_LEV]   [NP][NP]> derived_phi_f90(derived_phi,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_PHYSICAL_LEV]   [NP][NP]> derived_omega_p_f90(derived_omega_p,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_PHYSICAL_LEV][2][NP][NP]> derived_v_f90(derived_v,m_num_elems);

  sync_to_device(derived_phi_f90,     m_phi);
  sync_to_device(derived_omega_p_f90, m_omega_p);
  sync_to_device(derived_v_f90,       m_derived_vn0);
}

void Elements::pull_4d(CF90Ptr &state_v, CF90Ptr &state_t, CF90Ptr &state_dp3d) {
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV]   [NP][NP]> state_t_f90    (state_t,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV]   [NP][NP]> state_dp3d_f90 (state_dp3d,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]> state_v_f90    (state_v,m_num_elems);

  sync_to_device(state_t_f90,    m_t);
  sync_to_device(state_dp3d_f90, m_dp3d);
  sync_to_device(state_v_f90,    m_v);
}

void Elements::pull_eta_dot(CF90Ptr &derived_eta_dot_dpdn) {
  HostViewUnmanaged<const Real *[NUM_INTERFACE_LEV][NP][NP]> eta_dot_dpdn_f90(derived_eta_dot_dpdn,m_num_elems);
  sync_to_device_i2p(eta_dot_dpdn_f90,m_eta_dot_dpdn);
}

void Elements::push_to_f90_pointers(F90Ptr &state_v, F90Ptr &state_t,
                                    F90Ptr &state_dp3d, F90Ptr &derived_phi,
                                    F90Ptr &derived_omega_p, F90Ptr &derived_v,
                                    F90Ptr &derived_eta_dot_dpdn) const {
  push_3d(derived_phi, derived_omega_p, derived_v);
  push_4d(state_v, state_t, state_dp3d);
  push_eta_dot(derived_eta_dot_dpdn);
}

void Elements::push_3d(F90Ptr &derived_phi, F90Ptr &derived_omega_p, F90Ptr &derived_v) const {
  HostViewUnmanaged<Real *[NUM_PHYSICAL_LEV]   [NP][NP]> derived_phi_f90(derived_phi,m_num_elems);
  HostViewUnmanaged<Real *[NUM_PHYSICAL_LEV]   [NP][NP]> derived_omega_p_f90(derived_omega_p,m_num_elems);
  HostViewUnmanaged<Real *[NUM_PHYSICAL_LEV][2][NP][NP]> derived_v_f90(derived_v,m_num_elems);

  sync_to_host(m_phi,         derived_phi_f90);
  sync_to_host(m_omega_p,     derived_omega_p_f90);
  sync_to_host(m_derived_vn0, derived_v_f90);
}

void Elements::push_4d(F90Ptr &state_v, F90Ptr &state_t, F90Ptr &state_dp3d) const {
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV]   [NP][NP]> state_t_f90    (state_t,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV]   [NP][NP]> state_dp3d_f90 (state_dp3d,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]> state_v_f90    (state_v,m_num_elems);

  sync_to_host(m_t,    state_t_f90);
  sync_to_host(m_dp3d, state_dp3d_f90);
  sync_to_host(m_v,    state_v_f90);
}

void Elements::push_eta_dot(F90Ptr &derived_eta_dot_dpdn) const {
  HostViewUnmanaged<Real *[NUM_INTERFACE_LEV][NP][NP]> eta_dot_dpdn_f90(derived_eta_dot_dpdn,m_num_elems);
  sync_to_host_p2i(m_eta_dot_dpdn,eta_dot_dpdn_f90);
}

void Elements::d(Real *d_ptr, int ie) const {
  ExecViewUnmanaged<Real[2][2][NP][NP]> d_device = Homme::subview(m_d, ie);
  decltype(d_device)::HostMirror d_host = Kokkos::create_mirror_view(d_device);
  HostViewUnmanaged<Real[2][2][NP][NP]> d_wrapper(d_ptr);
  Kokkos::deep_copy(d_host, d_device);
  Kokkos::deep_copy(d_wrapper,d_host);
}

void Elements::dinv(Real *dinv_ptr, int ie) const {
  ExecViewUnmanaged<Real[2][2][NP][NP]> dinv_device = Homme::subview(m_dinv,ie);
  decltype(dinv_device)::HostMirror dinv_host = Kokkos::create_mirror_view(dinv_device);
  HostViewUnmanaged<Real[2][2][NP][NP]> dinv_wrapper(dinv_ptr);
  Kokkos::deep_copy(dinv_host, dinv_device);
  Kokkos::deep_copy(dinv_wrapper,dinv_host);
}

void Elements::BufferViews::init(const int num_elems) {
  pressure =
      ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Pressure buffer", num_elems);
  pressure_grad = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>(
      "Gradient of pressure", num_elems);
  temperature_virt = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
      "Virtual Temperature", num_elems);
  temperature_grad = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>(
      "Gradient of temperature", num_elems);
  omega_p = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
      "Omega_P = omega/pressure = (Dp/Dt)/pressure", num_elems);
  vdp = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("(u,v)*dp", num_elems);
  div_vdp = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
      "Divergence of dp3d * (u,v)", num_elems);
  ephi = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
      "Kinetic Energy + Geopotential Energy", num_elems);
  energy_grad = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>(
      "Gradient of ephi", num_elems);
  vorticity =
      ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Vorticity", num_elems);

  ttens  = ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]>("Temporary for temperature",num_elems);
  dptens = ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]>("Temporary for dp3d",num_elems);
  vtens  = ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]>("Temporary for velocity",num_elems);

  vstar = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("buffer for (flux v)/dp",
       num_elems);
  dpdissk = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
      "dpdissk", num_elems);

  preq_buf = ExecViewManaged<Real * [NP][NP]>("Preq Buffer", num_elems);

  sdot_sum = ExecViewManaged<Real * [NP][NP]>("Sdot sum buffer", num_elems);

  div_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Divergence Buffer",
                                                           num_elems);
  grad_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Gradient Buffer",
                                                            num_elems);
  curl_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Vorticity Buffer",
                                                            num_elems);

  sphere_vector_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("laplacian vector Buffer", num_elems);

  divergence_temp = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Divergence temporary",
                                                            num_elems);
  vorticity_temp = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Vorticity temporary",
                                                            num_elems);
  lapl_buf_1 = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Scalar laplacian Buffer", num_elems);
  lapl_buf_2 = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Scalar laplacian Buffer", num_elems);
  lapl_buf_3 = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Scalar laplacian Buffer", num_elems);
  v_vadv_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("v_vadv buffer",
                                                              num_elems);
  t_vadv_buf = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("t_vadv buffer",
                                                           num_elems);
  eta_dot_dpdn_buf = ExecViewManaged<Scalar * [NP][NP][NUM_LEV_P]>("eta_dot_dpdpn buffer",
                                                                   num_elems);

  kernel_start_times = ExecViewManaged<clock_t *>("Start Times", num_elems);
  kernel_end_times = ExecViewManaged<clock_t *>("End Times", num_elems);
}

} // namespace Homme
