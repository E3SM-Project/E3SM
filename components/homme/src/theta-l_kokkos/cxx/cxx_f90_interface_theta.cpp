/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "CaarFunctor.hpp"
#include "Context.hpp"
#include "Diagnostics.hpp"
#include "DirkFunctor.hpp"
#include "Elements.hpp"
#include "ErrorDefs.hpp"
#include "EulerStepFunctor.hpp"
#include "ForcingFunctor.hpp"
#include "FunctorsBuffersManager.hpp"
#include "HommexxEnums.hpp"
#include "HybridVCoord.hpp"
#include "HyperviscosityFunctor.hpp"
#include "ReferenceElement.hpp"
#include "SimulationParams.hpp"
#include "SphereOperators.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "VerticalRemapManager.hpp"
#include "mpi/BoundaryExchange.hpp"
#include "mpi/MpiBuffersManager.hpp"

#include "utilities/SyncUtils.hpp"

#include "profiling.hpp"

namespace Homme
{

extern "C"
{

void init_simulation_params_c (const int& remap_alg, const int& limiter_option, const int& rsplit, const int& qsplit,
                               const int& time_step_type, const int& qsize, const int& state_frequency,
                               const Real& nu, const Real& nu_p, const Real& nu_q, const Real& nu_s, const Real& nu_div, const Real& nu_top,
                               const int& hypervis_order, const int& hypervis_subcycle, const double& hypervis_scaling, const double& dcmip16_mu,
                               const int& ftype, const int& theta_adv_form, const bool& prescribed_wind, const bool& moisture, const bool& disable_diagnostics,
                               const bool& use_cpstar, const bool& use_semi_lagrangian_transport, const bool& theta_hydrostatic_mode, const char** test_case)
{
  // Check that the simulation options are supported. This helps us in the future, since we
  // are currently 'assuming' some option have/not have certain values. As we support for more
  // options in the C++ build, we will remove some checks
  Errors::check_option("init_simulation_params_c","vert_remap_q_alg",remap_alg,{1,3});
  Errors::check_option("init_simulation_params_c","prescribed_wind",prescribed_wind,{false});
  Errors::check_option("init_simulation_params_c","hypervis_order",hypervis_order,{2});
  Errors::check_option("init_simulation_params_c","use_semi_lagrangian_transport",use_semi_lagrangian_transport,{false});
  Errors::check_option("init_simulation_params_c","time_step_type",time_step_type,{1,4,5,6,7,9,10});
  Errors::check_option("init_simulation_params_c","qsize",qsize,0,Errors::ComparisonOp::GE);
  Errors::check_option("init_simulation_params_c","qsize",qsize,QSIZE_D,Errors::ComparisonOp::LE);
  Errors::check_option("init_simulation_params_c","limiter_option",limiter_option,{8,9});
  Errors::check_option("init_simulation_params_c","ftype",ftype, {-1, 0, 2});
  Errors::check_option("init_simulation_params_c","nu_p",nu_p,0.0,Errors::ComparisonOp::GT);
  Errors::check_option("init_simulation_params_c","nu",nu,0.0,Errors::ComparisonOp::GT);
  Errors::check_option("init_simulation_params_c","nu_div",nu_div,0.0,Errors::ComparisonOp::GT);
  Errors::check_option("init_simulation_params_c","theta_advection_form",theta_adv_form,{0,1});

  // Get the simulation params struct
  SimulationParams& params = Context::singleton().create<SimulationParams>();

  if (remap_alg==1) {
    params.remap_alg = RemapAlg::PPM_MIRRORED;
  } else if (remap_alg == 2) {
    params.remap_alg = RemapAlg::PPM_FIXED_PARABOLA;
  } else if (remap_alg == 3) {
    params.remap_alg = RemapAlg::PPM_FIXED_MEANS;
  }

  if (theta_adv_form==0) {
    params.theta_adv_form = AdvectionForm::Conservative;
  } else {
    params.theta_adv_form = AdvectionForm::NonConservative;
  }

  params.limiter_option                = limiter_option;
  params.rsplit                        = rsplit;
  params.qsplit                        = qsplit;
  params.prescribed_wind               = prescribed_wind;
  params.state_frequency               = state_frequency;
  params.qsize                         = qsize;
  params.nu                            = nu;
  params.nu_p                          = nu_p;
  params.nu_q                          = nu_q;
  params.nu_s                          = nu_s;
  params.nu_div                        = nu_div;
  params.nu_top                        = nu_top;
  params.hypervis_order                = hypervis_order;
  params.hypervis_subcycle             = hypervis_subcycle;
  params.hypervis_scaling              = hypervis_scaling;
  params.disable_diagnostics           = disable_diagnostics;
  params.moisture                      = (moisture ? MoistDry::MOIST : MoistDry::DRY);
  params.use_cpstar                    = use_cpstar;
  params.use_semi_lagrangian_transport = use_semi_lagrangian_transport;
  params.theta_hydrostatic_mode        = theta_hydrostatic_mode;
  params.dcmip16_mu                    = dcmip16_mu;
  if (time_step_type==0) {
    params.time_step_type = TimeStepType::LF;
  } else if (time_step_type==1) {
    params.time_step_type = TimeStepType::RK2;
  } else if (time_step_type==4) {
    params.time_step_type = TimeStepType::IMEX_KG254_EX;
  } else if (time_step_type==5) {
    params.time_step_type = TimeStepType::ULLRICH_RK35;
  } else if (time_step_type==6) {
    params.time_step_type = TimeStepType::IMEX_KG243;
  } else if (time_step_type==7) {
    params.time_step_type = TimeStepType::IMEX_KG254;
  } else if (time_step_type==9) {
    params.time_step_type = TimeStepType::IMEX_KG355;
  } else if (time_step_type==10) {
    params.time_step_type = TimeStepType::IMEX_KG255;
  }

  //set nu_ratios values
  if (params.nu != params.nu_div) {
    Real ratio = params.nu_div / params.nu;
    if (params.hypervis_scaling != 0.0) {
      params.nu_ratio1 = ratio * ratio;
      params.nu_ratio2 = 1.0;
    }else{
      params.nu_ratio1 = ratio;
      params.nu_ratio2 = ratio;
    }
  }else{
    params.nu_ratio1 = 1.0;
    params.nu_ratio2 = 1.0;
  }

  if (ftype == -1) {
    params.ftype = ForcingAlg::FORCING_OFF;
  } else if (ftype == 0) {
    params.ftype = ForcingAlg::FORCING_DEBUG;
  } else if (ftype == 2) {
    params.ftype = ForcingAlg::FORCING_2;
  }

  // TODO Parse a fortran string and set this properly. For now, our code does
  // not depend on this except to throw an error in apply_test_forcing.
  std::string test_name(*test_case);
  if (test_name=="jw_baroclinic") {
    params.test_case = TestCase::JW_BAROCLINIC;
  } else {
    Errors::runtime_abort("Error! Unknown test case '" + test_name + "'.\n");
  }

  // Now this structure can be used safely
  params.params_set = true;

}

void init_hvcoord_c (const Real& ps0, CRCPtr& hybrid_am_ptr, CRCPtr& hybrid_ai_ptr,
                                      CRCPtr& hybrid_bm_ptr, CRCPtr& hybrid_bi_ptr)
{
  HybridVCoord& hvcoord = Context::singleton().create<HybridVCoord>();
  hvcoord.init(ps0,hybrid_am_ptr,hybrid_ai_ptr,hybrid_bm_ptr,hybrid_bi_ptr);
}

void cxx_push_results_to_f90(F90Ptr &elem_state_v_ptr,         F90Ptr &elem_state_w_i_ptr,
                             F90Ptr &elem_state_vtheta_dp_ptr, F90Ptr &elem_state_phinh_i_ptr,
                             F90Ptr &elem_state_dp3d_ptr,      F90Ptr &elem_state_ps_v_ptr,
                             F90Ptr &elem_state_Qdp_ptr,       F90Ptr &elem_Q_ptr,
                             F90Ptr &elem_derived_omega_p_ptr) {
  ElementsState &state = Context::singleton().get<ElementsState>();
  const int num_elems = state.num_elems();

  state.push_to_f90_pointers(elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr,
                             elem_state_phinh_i_ptr, elem_state_dp3d_ptr);

  Tracers &tracers = Context::singleton().get<Tracers>();
  tracers.push_qdp(elem_state_Qdp_ptr);

  // F90 ptrs to arrays (np,np,num_time_levels,nelemd) can be stuffed directly
  // in an unmanaged view
  // with scalar Real*[NUM_TIME_LEVELS][NP][NP] (with runtime dimension nelemd)
  HostViewUnmanaged<Real * [NUM_TIME_LEVELS][NP][NP]> ps_v_f90(
      elem_state_ps_v_ptr, num_elems);

  auto ps_v_host = Kokkos::create_mirror_view(state.m_ps_v);

  Kokkos::deep_copy(ps_v_host, state.m_ps_v);
  Kokkos::deep_copy(ps_v_f90, ps_v_host);

  ElementsDerivedState &derived = Context::singleton().get<ElementsDerivedState>();
  sync_to_host(derived.m_omega_p,
               HostViewUnmanaged<Real * [NUM_PHYSICAL_LEV][NP][NP]>(
                   elem_derived_omega_p_ptr, num_elems));
  sync_to_host(tracers.Q,
               HostViewUnmanaged<Real * [QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]>(
                   elem_Q_ptr, num_elems));

  Diagnostics& diags = Context::singleton().get<Diagnostics>();
  diags.sync_diags_to_host();
}

void push_forcing_to_c (F90Ptr elem_derived_FM,
                        F90Ptr elem_derived_FVTheta,
                        F90Ptr elem_derived_FT,
                        F90Ptr elem_derived_FPHI,
                        F90Ptr elem_derived_FQ) {
  ElementsForcing &forcing = Context::singleton().get<ElementsForcing>();
  const int num_elems = forcing.num_elems();

  HostViewUnmanaged<Real *[NUM_PHYSICAL_LEV][3][NP][NP]> fm_f90(
      elem_derived_FM, num_elems);
  sync_to_device<3>(fm_f90, forcing.m_fm);

  HostViewUnmanaged<Real * [NUM_PHYSICAL_LEV][NP][NP]> fvtheta_f90(
      elem_derived_FVTheta, num_elems);
  sync_to_device(fvtheta_f90, forcing.m_fvtheta);

  HostViewUnmanaged<Real * [NUM_PHYSICAL_LEV][NP][NP]> ft_f90(
      elem_derived_FT, num_elems);
  sync_to_device(ft_f90, forcing.m_ft);

  HostViewUnmanaged<Real * [NUM_INTERFACE_LEV][NP][NP]> fphi_f90(
      elem_derived_FPHI, num_elems);
  sync_to_device(fphi_f90, forcing.m_fphi);

  const SimulationParams &params = Context::singleton().get<SimulationParams>();
  Tracers &tracers = Context::singleton().get<Tracers>();
  if (params.ftype == ForcingAlg::FORCING_DEBUG) {
    if (tracers.fq.data() == nullptr) {
      tracers.fq = decltype(tracers.fq)("fq", num_elems);
    }
    HostViewUnmanaged<Real * [QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]> fq_f90(
        elem_derived_FQ, num_elems);
    sync_to_device(fq_f90, tracers.fq);
  }
}

void init_reference_element_c (CF90Ptr& deriv, CF90Ptr& mass)
{
  ReferenceElement& ref_FE = Context::singleton().create<ReferenceElement> ();
  ref_FE.init(deriv,mass);
}

void init_time_level_c (const int& nm1, const int& n0, const int& np1,
                        const int& nstep, const int& nstep0)
{
  TimeLevel& tl = Context::singleton().create<TimeLevel>();
  tl.nm1    = nm1-1;
  tl.n0     = n0-1;
  tl.np1    = np1-1;
  tl.nstep  = nstep;
  tl.nstep0 = nstep0;
}

void init_elements_c (const int& num_elems)
{
  auto& c = Context::singleton();

  Elements& e = c.create<Elements> ();
  const SimulationParams& params = c.get<SimulationParams>();

  const bool consthv = (params.hypervis_scaling==0.0);
  e.init (num_elems, consthv, /* alloc_gradphis = */ true);

  // Init also the tracers structure
  Tracers& t = c.create<Tracers> ();
  t.init(num_elems,params.qsize);

  // In the context, we register also Elements[Geometry|State|DerivedState|Forcing],
  // making sure they store the same views as in the subobjects of Elements.
  // This allows objects that need only a piece of Elements, to grab it from the Context,
  // while still knowing that what they grab contains the same views as the object stored in the
  // Elements inside the Context
  // WARNING: after this point, you should NOT do things like
  //             e.m_geometry.m_phis = ...
  //          since they would NOT be reflected into the ElementsGeometry stored in the Context.
  //          In other words, you cannot reset the views. If you *really* need to do it,
  //          you must reset the view in both c.get<Elements>().m_geometry AND
  //          c.get<ElementsGeometry>()
  c.create_ref<ElementsGeometry>(e.m_geometry);
  c.create_ref<ElementsState>(e.m_state);
  c.create_ref<ElementsDerivedState>(e.m_derived);
#ifndef HOMMEXX_GB_CONFIG
  // For GB submission, skip forcing
  c.create_ref<ElementsForcing>(e.m_forcing);
#endif
}

void init_functors_c ()
{
  // We init all the functors in the Context, so that every call to
  // Context::singleton().get<[FunctorName]>()
  // will return a functor already initialized.
  // This avoids the risk of having a class doing
  //   FunctorName f = Context::singleton().get<FunctorName>();
  //   f.init(some_args);
  // and then, somewhere else, we find
  //   FunctorName f = Context::singleton().get<FunctorName>(some_args);
  // The problem is that the first call created an uninitialized functor,
  // *copied it*, and initialized the copy. The second call to get,
  // sees that there is *already* an object of type FunctorName in the
  // Context, and therefore does *not* create a FunctorName object.
  // A solution would be to use a reference in the first call, or to always
  // use the get method with the initialization arguments. However, if the functor
  // is already init-ed, what do we do? Ignore: there's the risk that the user thinks
  // that the functor was created with the inputs he/she provided Throw: we create
  // the headache of first finding out if a functor already exists, and, if not,
  // create it and init it.
  // It is way easier to create all functors here once and for all,
  // so that the user does not have to initialize them when calling them.

  auto& elems   = Context::singleton().get<Elements>();
  auto& tracers = Context::singleton().get<Tracers>();
  auto& ref_FE  = Context::singleton().get<ReferenceElement>();
  auto& hvcoord = Context::singleton().get<HybridVCoord>();
  auto& params  = Context::singleton().get<SimulationParams>();

  // Check that the above structures have been inited
  Errors::runtime_check(elems.inited(),    "Error! You must initialize the Elements structure before initializing the functors.\n", -1);
  Errors::runtime_check(tracers.inited(),  "Error! You must initialize the Tracers structure before initializing the functors.\n", -1);
  Errors::runtime_check(ref_FE.inited(),   "Error! You must initialize the ReferenceElement structure before initializing the functors.\n", -1);
  Errors::runtime_check(hvcoord.m_inited,  "Error! You must initialize the HybridVCoord structure before initializing the functors.\n", -1);
  Errors::runtime_check(params.params_set, "Error! You must initialize the SimulationParams structure before initializing the functors.\n", -1);

  // First, sphere operators, then the others
  auto& sph_op = Context::singleton().create<SphereOperators>(elems.m_geometry,ref_FE);
  auto& caar = Context::singleton().create<CaarFunctor>(elems,tracers,ref_FE,hvcoord,sph_op,params);
  auto& esf  = Context::singleton().create<EulerStepFunctor>();
  auto& hvf  = Context::singleton().create<HyperviscosityFunctor>();
  auto& fbm  = Context::singleton().create<FunctorsBuffersManager>();
  auto& ff   = Context::singleton().create<ForcingFunctor>();
  auto& diag = Context::singleton().create<Diagnostics> (elems.num_elems());
  auto& vrm  = Context::singleton().create<VerticalRemapManager>();

  const bool need_dirk = (params.time_step_type==TimeStepType::IMEX_KG243 ||   
                          params.time_step_type==TimeStepType::IMEX_KG254 ||
                          params.time_step_type==TimeStepType::IMEX_KG255 ||
                          params.time_step_type==TimeStepType::IMEX_KG355);

  if (need_dirk) {
    // Create dirk functor only if needed
    Context::singleton().create<DirkFunctor>(elems.num_elems());
  }

  // Make the functor request their buffer to the buffers manager
  // Note: diagnostics also needs buffers
  fbm.request_size(caar.requested_buffer_size());
  fbm.request_size(esf.requested_buffer_size());
  fbm.request_size(hvf.requested_buffer_size());
  fbm.request_size(diag.requested_buffer_size());
  fbm.request_size(ff.requested_buffer_size());
  fbm.request_size(vrm.requested_buffer_size());
  if (need_dirk) {
    const auto& dirk = Context::singleton().get<DirkFunctor>();
    fbm.request_size(dirk.requested_buffer_size());
  }

  // Allocate the buffers in the FunctorsBuffersManager, then tell the functors to grab their buffers
  fbm.allocate();

  caar.init_buffers(fbm);
  esf.init_buffers(fbm);
  hvf.init_buffers(fbm);
  diag.init_buffers(fbm);
  ff.init_buffers(fbm);
  vrm.init_buffers(fbm);
  if (need_dirk) {
    auto& dirk = Context::singleton().get<DirkFunctor>();
    dirk.init_buffers(fbm);
  }
}

void init_elements_2d_c (const int& ie,
                         CF90Ptr& D, CF90Ptr& Dinv, CF90Ptr& fcor,
                         CF90Ptr& spheremp, CF90Ptr& rspheremp,
                         CF90Ptr& metdet, CF90Ptr& metinv,
                         CF90Ptr& phis, CF90Ptr& gradphis,
                         CF90Ptr &tensorvisc, CF90Ptr &vec_sph2cart)
{
  Elements& e = Context::singleton().get<Elements> ();
  const SimulationParams& params = Context::singleton().get<SimulationParams>();

  const bool consthv = (params.hypervis_scaling==0.0);
  e.m_geometry.set_elem_data(ie,D,Dinv,fcor,spheremp,rspheremp,metdet,metinv,phis,tensorvisc,vec_sph2cart,consthv);

  // Note: yes, we *could* compute gradphis from grad, but at the time of this call,
  //       we do not yet have the SphereOperators functor initialized. For simplicity,
  //       we just require gradphis as input, and copy it manually.
  //       We also create the view here (rather than in ElementsGeometry init function),
  //       so that the view is not allocated in ElementsGeometry for the preqx model.
  HostViewUnmanaged<const Real [2][NP][NP]> h_gradphis(gradphis);
  sync_to_device(h_gradphis,Homme::subview(e.m_geometry.m_gradphis,ie));
}

void init_elements_states_c (CF90Ptr& elem_state_v_ptr,       CF90Ptr& elem_state_w_i_ptr, CF90Ptr& elem_state_vtheta_dp_ptr,
                             CF90Ptr& elem_state_phinh_i_ptr, CF90Ptr& elem_state_dp3d_ptr,
                             CF90Ptr& elem_state_ps_v_ptr,    CF90Ptr& elem_state_Qdp_ptr)
{
  ElementsState& state = Context::singleton().get<ElementsState> ();
  state.pull_from_f90_pointers(elem_state_v_ptr,elem_state_w_i_ptr,elem_state_vtheta_dp_ptr,
                               elem_state_phinh_i_ptr,elem_state_dp3d_ptr,elem_state_ps_v_ptr);
  Tracers &tracers = Context::singleton().get<Tracers>();
  tracers.pull_qdp(elem_state_Qdp_ptr);
}

void init_reference_states_c (CF90Ptr& elem_theta_ref_ptr, 
                              CF90Ptr& elem_dp_ref_ptr,
                              CF90Ptr& elem_phi_ref_ptr)
{
  auto& state = Context::singleton().get<ElementsState> ();
  auto& ref_states = state.m_ref_states;

  const int num_elems = state.m_ref_states.num_elems();
  assert(num_elems>0);

  HostViewUnmanaged<const Real*[NUM_PHYSICAL_LEV][NP][NP]>  theta_ref(elem_theta_ref_ptr,num_elems);
  HostViewUnmanaged<const Real*[NUM_PHYSICAL_LEV][NP][NP]>  dp_ref(elem_dp_ref_ptr,num_elems);
  HostViewUnmanaged<const Real*[NUM_INTERFACE_LEV][NP][NP]> phi_ref(elem_phi_ref_ptr,num_elems);

  sync_to_device(theta_ref, ref_states.theta_ref);
  sync_to_device(dp_ref,    ref_states.dp_ref);
  sync_to_device(phi_ref,   ref_states.phi_i_ref);
}

void init_diagnostics_c (F90Ptr& elem_state_q_ptr, F90Ptr& elem_accum_qvar_ptr,  F90Ptr& elem_accum_qmass_ptr,
                         F90Ptr& elem_accum_q1mass_ptr, F90Ptr& elem_accum_iener_ptr,
                         F90Ptr& elem_accum_kener_ptr, F90Ptr& elem_accum_pener_ptr)
{
  ElementsGeometry& geometry = Context::singleton().get<ElementsGeometry> ();
  ElementsState&    state    = Context::singleton().get<ElementsState> ();
  Diagnostics&      diags    = Context::singleton().get<Diagnostics> ();

  auto& params  = Context::singleton().get<SimulationParams>();
  auto& hvcoord = Context::singleton().get<HybridVCoord>();
  
  diags.init(state, geometry, hvcoord, params.theta_hydrostatic_mode,
             elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr, elem_accum_q1mass_ptr,
             elem_accum_iener_ptr, elem_accum_kener_ptr, elem_accum_pener_ptr);
}

void init_boundary_exchanges_c ()
{
  auto& params       = Context::singleton().get<SimulationParams>();
  auto  connectivity = Context::singleton().get_ptr<Connectivity>();
  assert (connectivity->is_finalized());

  // Create buffers manager map
  auto& bmm = Context::singleton().create<MpiBuffersManagerMap>();
  bmm[MPI_EXCHANGE]->set_connectivity(connectivity);
  bmm[MPI_EXCHANGE_MIN_MAX]->set_connectivity(connectivity);

  // Euler BEs
  auto& esf = Context::singleton().get<EulerStepFunctor>();
  esf.reset(params);
  esf.init_boundary_exchanges();

  // RK stages BE's
  auto& cf = Context::singleton().get<CaarFunctor>();
  cf.init_boundary_exchanges(bmm[MPI_EXCHANGE]);

  // HyperviscosityFunctor's BE's
  auto& hvf = Context::singleton().get<HyperviscosityFunctor>();
  hvf.init_boundary_exchanges();
}

} // extern "C"

} // namespace Homme
