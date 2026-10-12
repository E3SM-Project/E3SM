#include "physics/p3/eamxx_p3_process_interface.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/physics/physics_functions.hpp"
#include "share/physics/physics_saturation_impl.hpp"

#ifdef EAMXX_HAS_INFERENCE
#include "create_inference_backend.hpp"
#include "inference_backend.hpp"
#include "tensor.hpp"
#endif

#include <ekat_assert.hpp>

#include <array>
#include <cstdlib>
#include <string>

/*
 * Warm-rain emulator of P3: a model that computes the rates of the four
 * warm-rain collision processes (autoconversion, accretion, droplet and rain
 * self-collection) instead of P3, where it says it applies.
 *
 * The model is a python module, run through the python inference backend of
 * components/emulators, before p3_main, on the state P3's own warm-rain rates
 * are computed from: the state at the start of the step, with the adjustments
 * of p3_main_part1 (see run_warm_rain_emulator). Its rates are inputs of
 * p3_main: part2 merges them into P3's (see warm_rain_emulator_merge), and
 * they go through the same conservation checks.
 *
 * Contract of the model: the module has create_emulator(config), which returns
 * an object with infer(inputs, outputs). Both are dicts of (ncol, nlev) numpy
 * arrays of doubles, by name. They view host memory of P3 (on a GPU, the arrays
 * are copied to the host and back), and the outputs are written in place.
 *   inputs : qc, qr [kg/kg], nc, nr [#/kg] (cell averages, dry), rho [kg/m3]
 *   outputs: qc2qr_autoconv_tend, qc2qr_accret_tend [kg/kg/s], ncautr,
 *            nc2nr_autoconv_tend, nc_accret_tend, nc_selfcollect_tend (<= 0),
 *            nr_selfcollect_tend [#/kg/s], cell averages, with P3's signs;
 *            use_cloud, use_rain: 1 where P3 should use the cloud rates (the
 *            first six) and the rain self-collection, 0 where it should not
 * config is a dict of strings, with model_path, and num_threads: the number of
 * threads the model may use (those of this rank).
 *
 * p3 parameters:
 *   use_warm_rain_emulator                      (false)
 *   warm_rain_emulator_python_module            the module
 *   warm_rain_emulator_python_path              colon-separated directories with the
 *                                               module (the packages of the python EAMxx
 *                                               was configured with are added)
 *   warm_rain_emulator_model_path               given to the module (its model file)
 *   warm_rain_emulator_autoconv_fallback_scale  (1) factor on P3's autoconversion
 *                                               where the emulator does not apply
 *   warm_rain_emulator_cloud_self_collection    (true) use the emulated droplet
 *                                               self-collection (zero in P3)
 *   warm_rain_emulator_replaces_p3              (false) P3 does not compute its own
 *                                               warm-rain rates, and uses the emulator's
 *                                               in every cell, whatever its masks say:
 *                                               for models meant to apply everywhere
 */

namespace scream
{

#ifdef EAMXX_HAS_INFERENCE

namespace {
constexpr int num_inputs  = 5;
constexpr int num_outputs = 9;
const std::array<const char*, num_inputs> input_names = {"qc", "nc", "qr", "nr", "rho"};
const std::array<const char*, num_outputs> output_names = {
  "qc2qr_autoconv_tend", "qc2qr_accret_tend", "ncautr", "nc2nr_autoconv_tend", "nc_accret_tend",
  "nc_selfcollect_tend", "nr_selfcollect_tend", "use_cloud", "use_rain"};
} // anonymous namespace

// The backend, and the (ncol, nlev) arrays of doubles the model reads and
// writes: array n of the inputs is in(n,:,:), and likewise for the outputs. (One
// view each, not a std::array of views: Cuda kernels cannot read those.) The
// model works on their host mirrors (the same memory, when P3 runs on the host).
struct P3Microphysics::WarmRainEmulator {
  using buffer_t = Kokkos::View<double***, Kokkos::LayoutRight, DefaultDevice>;

  std::shared_ptr<emulator::inference::InferenceBackend> backend;
  buffer_t in, out;
  buffer_t::host_mirror_type in_h, out_h;
  emulator::inference::TensorMap inputs, outputs;
};
#endif

// =========================================================================================
void P3Microphysics::init_warm_rain_emulator ()
{
  if (not runtime_options.use_warm_rain_emulator) {
    EKAT_REQUIRE_MSG (not runtime_options.warm_rain_emulator_replaces_p3,
        "[P3Microphysics] Error! warm_rain_emulator_replaces_p3 needs use_warm_rain_emulator.\n");
    return;
  }

#ifndef EAMXX_HAS_INFERENCE
  EKAT_ERROR_MSG ("[P3Microphysics] Error! use_warm_rain_emulator requires EAMXX_ENABLE_INFERENCE=ON.\n");
#else
  using namespace emulator::inference;

  const std::string prefix = "[P3Microphysics] Error! use_warm_rain_emulator: ";
  const auto module = m_params.get<std::string>("warm_rain_emulator_python_module", "");
  EKAT_REQUIRE_MSG (not module.empty(), prefix + "warm_rain_emulator_python_module is not set.\n");

  InferenceConfig config;
  config.model_path      = m_params.get<std::string>("warm_rain_emulator_model_path", "");
  config.input_channels  = num_inputs;
  config.output_channels = num_outputs;
  config.set("python_module", module);
  // The python is the one EAMxx was configured with. Embedded, it needs to be told
  // where its standard library is (unless the environment says so already), and
  // where its packages are (after the directories of the module).
  setenv("PYTHONHOME", EAMXX_INFERENCE_PYTHON_HOME, 0);
  config.set("python_path", m_params.get<std::string>("warm_rain_emulator_python_path", "") +
                            ":" EAMXX_INFERENCE_PYTHON_PACKAGES);
  // Each rank must only use its own threads, not all the cores of the node
  config.set("num_threads", std::to_string(Kokkos::DefaultHostExecutionSpace().concurrency()));

  auto emu = std::make_shared<WarmRainEmulator>();
  try {
    emu->backend = create_backend(BackendType::PYTHON, config);
  } catch (const std::exception& e) {
    EKAT_ERROR_MSG (prefix + "could not load the python module '" + module + "':\n" + e.what() + "\n");
  }

  // The arrays the model reads and writes
  const std::vector<std::int64_t> dims = {m_num_cols, m_num_levs};
  const size_t array_size = static_cast<size_t>(m_num_cols) * m_num_levs;
  emu->in    = WarmRainEmulator::buffer_t("p3_warm_rain_emulator_in",  num_inputs,  m_num_cols, m_num_levs);
  emu->out   = WarmRainEmulator::buffer_t("p3_warm_rain_emulator_out", num_outputs, m_num_cols, m_num_levs);
  emu->in_h  = Kokkos::create_mirror_view(emu->in);
  emu->out_h = Kokkos::create_mirror_view(emu->out);
  for (int n=0; n<num_inputs; ++n) {
    const double* data = emu->in_h.data() + n*array_size;
    emu->inputs.wrap(input_names[n], data, dims);
  }
  for (int n=0; n<num_outputs; ++n) {
    emu->outputs.wrap(output_names[n], emu->out_h.data() + n*array_size, dims);
  }
  m_warm_rain_emulator = emu;

  // What part2 reads (padding stays zero, and is not used)
  const int nk_pack = ekat::npack<Pack>(m_num_levs);
  auto& rates = m_emulated_warm_rain;
  rates.qc2qr_autoconv_tend = view_2d("p3_emu_qc2qr_autoconv_tend", m_num_cols, nk_pack);
  rates.qc2qr_accret_tend   = view_2d("p3_emu_qc2qr_accret_tend",   m_num_cols, nk_pack);
  rates.ncautr              = view_2d("p3_emu_ncautr",              m_num_cols, nk_pack);
  rates.nc2nr_autoconv_tend = view_2d("p3_emu_nc2nr_autoconv_tend", m_num_cols, nk_pack);
  rates.nc_accret_tend      = view_2d("p3_emu_nc_accret_tend",      m_num_cols, nk_pack);
  rates.nc_selfcollect_tend = view_2d("p3_emu_nc_selfcollect_tend", m_num_cols, nk_pack);
  rates.nr_selfcollect_tend = view_2d("p3_emu_nr_selfcollect_tend", m_num_cols, nk_pack);
  rates.use_cloud           = view_2d("p3_emu_use_cloud",           m_num_cols, nk_pack);
  rates.use_rain            = view_2d("p3_emu_use_rain",            m_num_cols, nk_pack);

  m_atm_logger->info("[P3Microphysics] Warm-rain emulator: python module " + module +
                     (config.model_path.empty() ? "" : ", model " + config.model_path));
  m_atm_logger->info("[P3Microphysics]   autoconversion fallback scale: " +
                     std::to_string(runtime_options.warm_rain_emulator_autoconv_fallback_scale) +
                     ", emulated droplet self-collection: " +
                     (runtime_options.warm_rain_emulator_cloud_self_collection ? "on" : "off") +
                     (runtime_options.warm_rain_emulator_replaces_p3 ?
                      ", replaces P3's warm-rain rates everywhere (masks and fallback are not used)" : ""));
#endif
}

#ifdef EAMXX_HAS_INFERENCE
// =========================================================================================
void P3Microphysics::run_warm_rain_emulator ()
{
  using policy_t = Kokkos::MDRangePolicy<KT::ExeSpace, Kokkos::Rank<2>>;
  constexpr int N = Pack::n;

  start_timer("EAMxx::P3::warm_rain_emulator");

  auto& emu = *m_warm_rain_emulator;
  const int ncol = m_num_cols;
  const int nlev = m_num_levs;

  // 1. The inputs, as (ncol, nlev) arrays of doubles: the state P3's process rates
  //    are computed from, i.e., with what p3_main_part1 does to it first (dry
  //    density, clipping of trace cloud and rain, activation or prescription of
  //    nc, melting of trace ice into rain). KEEP IN SYNC with p3_main_part1.
  {
    using PhysF = physics::Functions<Real, DefaultDevice>;
    constexpr Real qsmall     = P3F::C::QSMALL;
    constexpr Real nccnst     = P3F::C::NCCNST;
    constexpr Real g          = P3F::C::gravit.value;
    constexpr Real T_zerodegc = P3F::C::T_zerodegc.value;
    const auto qc  = prog_state.qc;
    const auto nc  = prog_state.nc;
    const auto qr  = prog_state.qr;
    const auto nr  = prog_state.nr;
    const auto qi  = prog_state.qi;
    const auto qv  = prog_state.qv;
    const auto th  = prog_state.th;
    const auto pres           = diag_inputs.pres;
    const auto inv_exner      = diag_inputs.inv_exner;
    const auto dpres          = diag_inputs.dpres;
    const auto dz             = diag_inputs.dz;
    const auto cld_frac_l     = diag_inputs.cld_frac_l;
    const auto nccn           = diag_inputs.nccn;
    const auto nc_nuceat_tend = diag_inputs.nc_nuceat_tend;
    const bool predict_nc     = infrastructure.predictNc;
    const bool prescribed_ccn = infrastructure.prescribedCCN;
    const Real dt             = infrastructure.dt;
    const Real ccn_factor     = runtime_options.spa_ccn_to_nc_factor;
    const Real ccn_exponent   = runtime_options.spa_ccn_to_nc_exponent;
    const auto in  = emu.in;
    const int nk_pack = ekat::npack<Pack>(nlev);
    Kokkos::parallel_for("p3_warm_rain_emulator_inputs", policy_t({0, 0}, {ncol, nk_pack}),
                         KOKKOS_LAMBDA (const int i, const int p) {
      const auto range_mask = ekat::range<IntPack>(p*N) < nlev;
      const Pack rho     = dpres(i,p) / dz(i,p) / g;
      const Pack inv_rho = 1 / rho;
      Pack qc_p = qc(i,p), nc_p = nc(i,p), qr_p = qr(i,p), nr_p = nr(i,p);

      // Trace cloud evaporates; otherwise droplets are activated or prescribed
      const auto dry_c = qc_p < qsmall;
      const auto wet_c = !dry_c && range_mask;
      qc_p.set(dry_c, 0);
      nc_p.set(dry_c, 0);
      if (prescribed_ccn) {
        const Pack inv_cld_frac_l = 1 / cld_frac_l(i,p);
        Pack nccn_scaled = nccn(i,p) / inv_cld_frac_l;
        nccn_scaled = pow(nccn_scaled, ccn_exponent);
        nc_p.set(wet_c, max(nc_p, ccn_factor * nccn_scaled));
      } else if (predict_nc) {
        nc_p.set(wet_c, max(nc_p + nc_nuceat_tend(i,p) * dt, 0.0));
      } else {
        nc_p.set(wet_c, nccnst * inv_rho);
      }

      // Trace rain evaporates
      const auto dry_r = qr_p < qsmall;
      qr_p.set(dry_r, 0);
      nr_p.set(dry_r, 0);

      // Trace ice above freezing melts into rain, unless it sublimates first
      const Pack exner = 1 / inv_exner(i,p);
      const Pack T_atm = th(i,p) * exner;
      const Pack qv_sat_i = PhysF::qv_sat_dry(T_atm, pres(i,p), true, range_mask, PhysF::MurphyKoop,
                                              "p3 warm-rain emulator inputs");
      const Pack qv_supersat_i = qv(i,p) / qv_sat_i - 1;
      const Pack qi_p = qi(i,p);
      const auto sublimates = qi_p < qsmall || (qi_p < 1.e-8 && qv_supersat_i < -0.1);
      const auto melts = !sublimates && qi_p >= qsmall && qi_p < 1.e-8 && T_atm >= T_zerodegc;
      qr_p.set(melts, qr_p + qi_p);

      for (int s=0; s<N; ++s) {
        const int k = p*N + s;
        if (k < nlev) {
          in(0,i,k) = qc_p[s];
          in(1,i,k) = nc_p[s];
          in(2,i,k) = qr_p[s];
          in(3,i,k) = nr_p[s];
          in(4,i,k) = rho[s];
        }
      }
    });
    Kokkos::fence();
  }

  // 2. The model, on the host
  Kokkos::deep_copy(emu.in_h, emu.in);
  try {
    emu.backend->infer(emu.inputs, emu.outputs);
  } catch (const std::exception& e) {
    EKAT_ERROR_MSG ("[P3Microphysics] Error! The warm-rain emulator failed:\n" + std::string(e.what()) + "\n");
  }
  Kokkos::deep_copy(emu.out, emu.out_h);

  // 3. Its rates, as packs for part2. A rate that is not finite must not get into P3.
  {
    const auto out   = emu.out;
    const auto rates = m_emulated_warm_rain;
    int num_bad = 0;
    Kokkos::parallel_reduce("p3_warm_rain_emulator_outputs", policy_t({0, 0}, {ncol, nlev}),
                            KOKKOS_LAMBDA (const int i, const int k, int& bad) {
      const int p = k / N, s = k % N;
      for (int n=0; n<num_outputs; ++n) {
        if (!Kokkos::isfinite(out(n,i,k))) ++bad;
      }
      rates.qc2qr_autoconv_tend(i,p)[s] = out(0,i,k);
      rates.qc2qr_accret_tend(i,p)[s]   = out(1,i,k);
      rates.ncautr(i,p)[s]              = out(2,i,k);
      rates.nc2nr_autoconv_tend(i,p)[s] = out(3,i,k);
      rates.nc_accret_tend(i,p)[s]      = out(4,i,k);
      rates.nc_selfcollect_tend(i,p)[s] = out(5,i,k);
      rates.nr_selfcollect_tend(i,p)[s] = out(6,i,k);
      rates.use_cloud(i,p)[s]           = out(7,i,k);
      rates.use_rain(i,p)[s]            = out(8,i,k);
    }, num_bad);
    Kokkos::fence();
    if (num_bad>0) {
      // Say where the first one is, and what the model was given there
      int first = ncol*nlev;
      Kokkos::parallel_reduce("p3_warm_rain_emulator_first_bad", policy_t({0, 0}, {ncol, nlev}),
                              KOKKOS_LAMBDA (const int i, const int k, int& f) {
        for (int n=0; n<num_outputs; ++n) {
          if (!Kokkos::isfinite(out(n,i,k))) f = Kokkos::min(f, i*nlev + k);
        }
      }, Kokkos::Min<int>(first));
      const int ic = first / nlev, kc = first % nlev;
      std::string where = "column " + std::to_string(ic) + ", level " + std::to_string(kc) + ", with inputs";
      for (int n=0; n<num_inputs; ++n) {
        where += std::string(" ") + input_names[n] + "=" + std::to_string(emu.in_h(n,ic,kc));
      }
      EKAT_ERROR_MSG ("[P3Microphysics] Error! The warm-rain emulator returned " + std::to_string(num_bad) +
                      " values that are not finite (on rank " + std::to_string(get_comm().rank()) + ").\n"
                      "  The first is at local " + where + ".\n");
    }
  }

  stop_timer("EAMxx::P3::warm_rain_emulator");
}
#else
// =========================================================================================
void P3Microphysics::run_warm_rain_emulator ()
{
  EKAT_ERROR_MSG ("[P3Microphysics] Error! use_warm_rain_emulator requires EAMXX_ENABLE_INFERENCE=ON.\n");
}
#endif

} // namespace scream
