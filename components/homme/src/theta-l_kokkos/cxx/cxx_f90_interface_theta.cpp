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
#include "ComposeTransport.hpp"
#include "ForcingFunctor.hpp"
#include "FunctorsBuffersManager.hpp"
#include "HommexxEnums.hpp"
#include "HybridVCoord.hpp"
#include "HyperviscosityFunctor.hpp"
#include "LimiterFunctor.hpp"
#include "ReferenceElement.hpp"
#include "SimulationParams.hpp"
#include "SphereOperators.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "GllFvRemap.hpp"
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
                               const int& hypervis_order, const int& hypervis_subcycle, const int& hypervis_subcycle_tom, 
                               const double& hypervis_scaling, const double& dcmip16_mu,
                               const int& ftype, const int& theta_adv_form, const bool& prescribed_wind, const bool& moisture, const bool& disable_diagnostics,
                               const bool& use_cpstar, const int& transport_alg, const bool& theta_hydrostatic_mode, const char** test_case,
                               const int& dt_remap_factor, const int& dt_tracer_factor,
                               const double& rearth, const int& nsplit, const bool& pgrad_correction,
                               const double& dp3d_thresh, const double& vtheta_thresh)
{
  // Check that the simulation options are supported. This helps us in the future, since we
  // are currently 'assuming' some option have/not have certain values. As we support for more
  // options in the C++ build, we will remove some checks
  Errors::check_option("init_simulation_params_c","vert_remap_q_alg",remap_alg,{1,3,10});
  Errors::check_option("init_simulation_params_c","prescribed_wind",prescribed_wind,{false});
  Errors::check_option("init_simulation_params_c","hypervis_order",hypervis_order,{2});
  Errors::check_option("init_simulation_params_c","transport_alg",transport_alg,{0,12});
  Errors::check_option("init_simulation_params_c","time_step_type",time_step_type,{1,4,5,6,7,9,10});
  Errors::check_option("init_simulation_params_c","qsize",qsize,0,Errors::ComparisonOp::GE);
  Errors::check_option("init_simulation_params_c","qsize",qsize,QSIZE_D,Errors::ComparisonOp::LE);
  Errors::check_option("init_simulation_params_c","limiter_option",limiter_option,{8,9});
  Errors::check_option("init_simulation_params_c","ftype",ftype, {-1, 0, 2});
  Errors::check_option("init_simulation_params_c","nu_p",nu_p,0.0,Errors::ComparisonOp::GT);
  Errors::check_option("init_simulation_params_c","nu",nu,0.0,Errors::ComparisonOp::GT);
  Errors::check_option("init_simulation_params_c","dp3d_thresh",dp3d_thresh,0.0,Errors::ComparisonOp::GT);
  Errors::check_option("init_simulation_params_c","vtheta_thresh",vtheta_thresh,0.0,Errors::ComparisonOp::GT);
  Errors::check_option("init_simulation_params_c","nu_div",nu_div,0.0,Errors::ComparisonOp::GT);
  Errors::check_option("init_simulation_params_c","theta_advection_form",theta_adv_form,{0,1});
#ifndef SCREAM
  Errors::check_option("init_simulation_params_c","nsplit",nsplit,1,Errors::ComparisonOp::GE);
#else
  if (nsplit<1 && Context::singleton().get<Comm>().root()) {
    printf ("Note: nsplit=%d, while nsplit must be >=1. We know SCREAM does not know nsplit until runtime, so this is fine.\n"
            "      Make sure nsplit is set to a valid value before calling prim_advance_subcycle!\n",nsplit);
  }
#endif

  // Get the simulation params struct
  SimulationParams& params = Context::singleton().create<SimulationParams>();

  if (remap_alg==1) {
    params.remap_alg = RemapAlg::PPM_MIRRORED;
  } else if (remap_alg == 10) {
    params.remap_alg = RemapAlg::PPM_LIMITED_EXTRAP;
  }

  if (theta_adv_form==0) {
    params.theta_adv_form = AdvectionForm::Conservative;
  } else {
    params.theta_adv_form = AdvectionForm::NonConservative;
  }

  params.limiter_option                = limiter_option;
  params.rsplit                        = rsplit;
  params.qsplit                        = qsplit;
  params.dt_remap_factor               = dt_remap_factor;
  params.dt_tracer_factor              = dt_tracer_factor;
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
  params.hypervis_subcycle_tom         = hypervis_subcycle_tom;
  params.hypervis_scaling              = hypervis_scaling;
  params.disable_diagnostics           = disable_diagnostics;
  params.moisture                      = (moisture ? MoistDry::MOIST : MoistDry::DRY);
  params.use_cpstar                    = use_cpstar;
  params.transport_alg                 = transport_alg;
  params.theta_hydrostatic_mode        = theta_hydrostatic_mode;
  params.dcmip16_mu                    = dcmip16_mu;
  params.nsplit                        = nsplit;
  params.rearth                        = rearth;
  params.pgrad_correction              = pgrad_correction;
  params.dp3d_thresh                   = dp3d_thresh;
  params.vtheta_thresh                 = vtheta_thresh;

  if (time_step_type==5) {
    //5 stage, 3rd order, explicit
    params.time_step_type = TimeStepType::ttype5;
  } else if (time_step_type==7) {
    //5 stage, based on 2nd order explicit KGU table
    //1st order (BE) implicit part
    params.time_step_type = TimeStepType::ttype7_imex;
  } else if (time_step_type==9) {
    //5 stage, based on 3rd order explicit KGU table
    //2nd order implicit table
    params.time_step_type = TimeStepType::ttype9_imex;
  } else if (time_step_type==10) {
    //5 stage, based on the 2nd order explicit KGU table
    //2nd order implicit table
    params.time_step_type = TimeStepType::ttype10_imex;
  } else {
    Errors::runtime_abort("Invalid time_step_time" 
                          + std::to_string(time_step_type), Errors::err_not_implemented);
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
    params.ftype = ForcingAlg::FORCING_0;
  } else if (ftype == 2) {
    params.ftype = ForcingAlg::FORCING_2;
  }

  // TODO Parse a fortran string and set this properly. For now, our code does
  // not depend on this except to throw an error in apply_test_forcing.
  std::string test_name(*test_case);
  //TEMP
  params.test_case = TestCase::JW_BAROCLINIC;

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
}

//currently, we do not need FVTheta and FPHI, because they are computed from FT and FQ
//in applycamforcing_tracers inside xx
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

  Tracers &tracers = Context::singleton().get<Tracers>();
  if (tracers.fq.data() == nullptr) {
    tracers.fq = decltype(tracers.fq)("fq", num_elems, tracers.num_tracers());
  }
  HostViewUnmanaged<Real * [QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]> fq_f90(
      elem_derived_FQ, num_elems);
  sync_to_device(fq_f90, tracers.fq);
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
  e.init (num_elems, consthv, /* alloc_gradphis = */ true,
          params.rearth,
          /* alloc_sphere_coords = */ params.transport_alg > 0);

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
  c.create_ref<ElementsForcing>(e.m_forcing);
}

void init_functors_c (const bool& allocate_buffer)
{
  auto& c = Context::singleton();

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

  auto& elems    = c.get<Elements>();
  auto& tracers  = c.get<Tracers>();
  auto& ref_FE   = c.get<ReferenceElement>();
  auto& hvcoord  = c.get<HybridVCoord>();
  auto& params   = c.get<SimulationParams>();
  auto& geometry = c.get<ElementsGeometry>();
  auto& state    = c.get<ElementsState>();
  auto& derived  = c.get<ElementsDerivedState>();

  // Check that the above structures have been inited
  Errors::runtime_check(elems.inited(),    "Error! You must initialize the Elements structure before initializing the functors.\n", -1);
  Errors::runtime_check(tracers.inited(),  "Error! You must initialize the Tracers structure before initializing the functors.\n", -1);
  Errors::runtime_check(ref_FE.inited(),   "Error! You must initialize the ReferenceElement structure before initializing the functors.\n", -1);
  Errors::runtime_check(hvcoord.m_inited,  "Error! You must initialize the HybridVCoord structure before initializing the functors.\n", -1);
  Errors::runtime_check(params.params_set, "Error! You must initialize the SimulationParams structure before initializing the functors.\n", -1);

  // First, sphere operators, then the others
  auto& sph_op = c.create<SphereOperators>(elems.m_geometry,ref_FE);
  auto& limiter = c.create_if_not_there<LimiterFunctor>(elems,hvcoord,params);

  // Some functors might have been previously created, so
  // use the create_if_not_there() function.
  auto& caar = c.create_if_not_there<CaarFunctor>(elems,tracers,ref_FE,hvcoord,sph_op,params);
  if (params.transport_alg == 0) c.create_if_not_there<EulerStepFunctor>();
#ifdef HOMME_ENABLE_COMPOSE
  else                           c.create_if_not_there<ComposeTransport>();
#endif
  auto& hvf     = c.create_if_not_there<HyperviscosityFunctor>();
  auto& ff      = c.create_if_not_there<ForcingFunctor>();
  auto& diag    = c.create_if_not_there<Diagnostics> (elems.num_elems(),params.theta_hydrostatic_mode);
  auto& vrm     = c.create_if_not_there<VerticalRemapManager>(elems.num_elems());

  auto& fbm     = c.create_if_not_there<FunctorsBuffersManager>();

//OG why are if-statement here -- above calls define which constructor is called

  // If any Functor was constructed only partially, setup() must be called.
  // This does not apply to Diagnostics or DirkFunctor since they only
  // contain one constructor.
  if (caar.setup_needed()) {
    caar.setup(elems, tracers, ref_FE, hvcoord, sph_op);
  }
  if (params.transport_alg == 0) {
    auto& esf = c.get<EulerStepFunctor>();
    if (esf.setup_needed()) esf.setup();
  } else {
#ifdef HOMME_ENABLE_COMPOSE	  
    auto& ct = c.get<ComposeTransport>();
    if (ct.setup_needed()) ct.setup();
#endif    
  }
  if (hvf.setup_needed()) {
    hvf.setup(geometry, state, derived);
  }
  if (ff.setup_needed()) {
    ff.setup();
  }
  if (vrm.setup_needed()) {
    vrm.setup();
  }

  const bool need_dirk = (params.time_step_type==TimeStepType::ttype7_imex ||   
                          params.time_step_type==TimeStepType::ttype9_imex ||
                          params.time_step_type==TimeStepType::ttype10_imex  );

  if (need_dirk) {
    // Create dirk functor only if needed
    c.create_if_not_there<DirkFunctor>(elems.num_elems());
  }

  // If memory in the buffer manager was previously allocated, skip allocation here
  if (allocate_buffer) {
    // Make the functor request their buffer to the buffers manager
    // Note: diagnostics also needs buffers
    fbm.request_size(caar.requested_buffer_size());
    if (params.transport_alg == 0)
      fbm.request_size(c.get<EulerStepFunctor>().requested_buffer_size());
#ifdef HOMME_ENABLE_COMPOSE
    else	    
      fbm.request_size(c.get<ComposeTransport>().requested_buffer_size());
#endif
    fbm.request_size(hvf.requested_buffer_size());
    fbm.request_size(diag.requested_buffer_size());
    fbm.request_size(ff.requested_buffer_size());
    fbm.request_size(vrm.requested_buffer_size());
    fbm.request_size(limiter.requested_buffer_size());
    if (need_dirk) {
      const auto& dirk = Context::singleton().get<DirkFunctor>();
      fbm.request_size(dirk.requested_buffer_size());
    }

    // Allocate the buffers in the FunctorsBuffersManager, then tell the functors to grab their buffers
    fbm.allocate();
  }

  caar.init_buffers(fbm);
  if (params.transport_alg == 0)
    Context::singleton().get<EulerStepFunctor>().init_buffers(fbm);
#ifdef HOMME_ENABLE_COMPOSE
  else
    Context::singleton().get<ComposeTransport>().init_buffers(fbm);
#endif
  hvf.init_buffers(fbm);
  diag.init_buffers(fbm);
  ff.init_buffers(fbm);
  vrm.init_buffers(fbm);
  limiter.init_buffers(fbm);
  if (need_dirk) {
    auto& dirk = Context::singleton().get<DirkFunctor>();
    dirk.init_buffers(fbm);
  }
  // The SCREAM-side Hommexx interface will initialize GllFvRemap if it's
  // needed. But it expects the Homme-side Hommexx interface to init buffers, so
  // do that here.
  if (c.has<GllFvRemap>()) {
    auto& gfr = c.get<GllFvRemap>();
    gfr.setup();
    gfr.init_buffers(fbm);
  }
}

void init_elements_2d_c (const int& ie,
                         CF90Ptr& D, CF90Ptr& Dinv, CF90Ptr& fcor,
                         CF90Ptr& spheremp, CF90Ptr& rspheremp,
                         CF90Ptr& metdet, CF90Ptr& metinv,
                         CF90Ptr &tensorvisc, CF90Ptr &vec_sph2cart,
                         double* sphere_cart_vec, double* sphere_latlon_vec)
{
  auto& c = Context::singleton();
  Elements& e = c.get<Elements> ();
  const SimulationParams& params = c.get<SimulationParams>();

  const bool consthv = (params.hypervis_scaling==0.0);
  e.m_geometry.set_elem_data(ie,D,Dinv,fcor,spheremp,rspheremp,metdet,metinv,tensorvisc,
                             vec_sph2cart,consthv,sphere_cart_vec,sphere_latlon_vec);
}

void init_geopotential_c (const int& ie,
                          CF90Ptr& phis, CF90Ptr& gradphis)
{
  auto& geo = Context::singleton().get<ElementsGeometry>();

  geo.set_phis(ie,phis);

  // Note: yes, we *could* compute gradphis from grad, but at the time of this call,
  //       we do not yet have the SphereOperators functor initialized. For simplicity,
  //       we just require gradphis as input, and copy it manually.
  HostViewUnmanaged<const Real [2][NP][NP]> h_gradphis(gradphis);
  sync_to_device(h_gradphis,Homme::subview(geo.m_gradphis,ie));
}

void compute_gradphis_c ()
{
  auto& c = Context::singleton();
  auto& geo = c.get<ElementsGeometry>();
  auto& ref_FE = c.get<ReferenceElement>();
  SphereOperators sph_op(geo,ref_FE);
  auto p = get_default_team_policy<ExecSpace>(geo.num_elems());
  sph_op.allocate_buffers(p);
  Kokkos::parallel_for(p,KOKKOS_LAMBDA(const TeamMember& team) {
    KernelVariables kv(team);
    auto phis     = Homme::subview(geo.m_phis,kv.ie);
    auto gradphis = Homme::subview(geo.m_gradphis,kv.ie);
    sph_op.gradient_sphere_sl(kv,phis,gradphis);
  });
}

void init_elements_states_c (CF90Ptr& elem_state_v_ptr,       CF90Ptr& elem_state_w_i_ptr, CF90Ptr& elem_state_vtheta_dp_ptr,
                             CF90Ptr& elem_state_phinh_i_ptr, CF90Ptr& elem_state_dp3d_ptr,
                             CF90Ptr& elem_state_ps_v_ptr,    CF90Ptr& elem_state_Qdp_ptr)
{
  const auto& c = Context::singleton();
  ElementsState& state = c.get<ElementsState> ();
  state.pull_from_f90_pointers(elem_state_v_ptr,elem_state_w_i_ptr,elem_state_vtheta_dp_ptr,
                               elem_state_phinh_i_ptr,elem_state_dp3d_ptr,elem_state_ps_v_ptr);
  Tracers &tracers = c.get<Tracers>();
  tracers.pull_qdp(elem_state_Qdp_ptr);
  const auto  qdp = tracers.qdp;
  const auto  q = tracers.Q;
  const auto  dp = state.m_dp3d;
  auto& tl = c.get<TimeLevel>();
  const auto n0 = tl.n0;
  const SimulationParams& params = c.get<SimulationParams>();
  tl.update_tracers_levels(params.qsplit);
  const auto n0_qdp = tl.n0_qdp;
  const auto qsize = tracers.num_tracers();
  const auto size = tracers.num_elems()*tracers.num_tracers()*NP*NP*NUM_LEV;
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,size),
                       KOKKOS_LAMBDA(const int idx) {
    const int ie  =  idx / (qsize*NP*NP*NUM_LEV);
    const int iq  = (idx / (NP*NP*NUM_LEV)) % qsize;
    const int igp = (idx / (NP*NUM_LEV)) % NP;
    const int jgp = (idx / NUM_LEV) % NP;
    const int k   =  idx % NUM_LEV;

    q (ie,iq,igp,jgp,k) = qdp (ie,n0_qdp,iq,igp,jgp,k) / dp(ie,n0,igp,jgp,k);
  });
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

  auto& hvcoord = Context::singleton().get<HybridVCoord>();
  
  diags.init(state, geometry, hvcoord,
             elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr, elem_accum_q1mass_ptr,
             elem_accum_iener_ptr, elem_accum_kener_ptr, elem_accum_pener_ptr);
}

void init_boundary_exchanges_c ()
{
  auto& c = Context::singleton();

  auto& params       = c.get<SimulationParams>();
  auto  connectivity = c.get_ptr<Connectivity>();
  assert (connectivity->is_finalized());

  // Create buffers manager map
  auto& bmm = c.create_if_not_there<MpiBuffersManagerMap>();
  if (!bmm[MPI_EXCHANGE]->is_connectivity_set()) {
    bmm[MPI_EXCHANGE]->set_connectivity(connectivity);
  }
  if (!bmm[MPI_EXCHANGE_MIN_MAX]->is_connectivity_set()) {
    bmm[MPI_EXCHANGE_MIN_MAX]->set_connectivity(connectivity);
  }

  if (params.transport_alg == 0) {
    // Euler BEs
    auto& esf = c.get<EulerStepFunctor>();
    esf.reset(params);
    esf.init_boundary_exchanges();
  } else {
#ifdef HOMME_ENABLE_COMPOSE	  
    auto& ct = c.get<ComposeTransport>();
    ct.reset(params);
    ct.init_boundary_exchanges();
#endif    
  }

  // RK stages BE's
  auto& cf = c.get<CaarFunctor>();
  cf.init_boundary_exchanges(bmm[MPI_EXCHANGE]);

  // HyperviscosityFunctor's BE's
  auto& hvf = c.get<HyperviscosityFunctor>();
  hvf.init_boundary_exchanges();

  if (c.has<GllFvRemap>()) {
    auto& gfr = c.get<GllFvRemap>();
    gfr.reset(params);
    gfr.init_boundary_exchanges();
  }
}

} // extern "C"

} // namespace Homme
