#include "catch2/catch.hpp"
#include "physics/mam/eamxx_mam_constituent_fluxes_interface.hpp"
#include "physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
#include "share/data_managers/field_manager.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace {
using namespace scream;

constexpr int ncols = 218;
constexpr int nlevs = 72;
constexpr int nconst = mam4::aero_model::pcnst;
constexpr Real dt = 1800;
constexpr Real dp = Real(99000) / nlevs;
constexpr Real amufac = 1.65979e-23;
const std::string data_dir =
    std::string(SCREAM_DATA_DIR) + "/mam4xx/emissions/ne2np4/";

enum class SourceClass { Both, OnlineOnly, PrescribedOnly, Inactive };

struct TracerExpectation {
  const char* name;
  int mode;
  int species;
  bool is_number;
  SourceClass source_class;
};

// This is an explicit scientific expectation table, not a classification
// inferred from the online-emission implementation under test. Species is the
// within-mode index from mam4::mode_aero_species; it is unused for number.
constexpr std::array<TracerExpectation, 25> aerosol_expectations = {{
    {"so4_a1", 0, 0, false, SourceClass::PrescribedOnly},
    {"pom_a1", 0, 1, false, SourceClass::Inactive},
    {"soa_a1", 0, 2, false, SourceClass::Inactive},
    {"bc_a1", 0, 3, false, SourceClass::Inactive},
    {"dst_a1", 0, 4, false, SourceClass::OnlineOnly},
    {"ncl_a1", 0, 5, false, SourceClass::OnlineOnly},
    {"mom_a1", 0, 6, false, SourceClass::OnlineOnly},
    {"num_a1", 0, -1, true, SourceClass::Both},
    {"so4_a2", 1, 0, false, SourceClass::PrescribedOnly},
    {"soa_a2", 1, 1, false, SourceClass::Inactive},
    {"ncl_a2", 1, 2, false, SourceClass::OnlineOnly},
    {"mom_a2", 1, 3, false, SourceClass::OnlineOnly},
    {"num_a2", 1, -1, true, SourceClass::Both},
    {"dst_a3", 2, 0, false, SourceClass::OnlineOnly},
    {"ncl_a3", 2, 1, false, SourceClass::OnlineOnly},
    {"so4_a3", 2, 2, false, SourceClass::Inactive},
    {"bc_a3", 2, 3, false, SourceClass::Inactive},
    {"pom_a3", 2, 4, false, SourceClass::Inactive},
    {"soa_a3", 2, 5, false, SourceClass::Inactive},
    {"mom_a3", 2, 6, false, SourceClass::Inactive},
    {"num_a3", 2, -1, true, SourceClass::OnlineOnly},
    {"pom_a4", 3, 0, false, SourceClass::PrescribedOnly},
    {"bc_a4", 3, 1, false, SourceClass::PrescribedOnly},
    {"mom_a4", 3, 2, false, SourceClass::Inactive},
    {"num_a4", 3, -1, true, SourceClass::PrescribedOnly},
}};

static_assert(aerosol_expectations.size() ==
              mam_coupling::num_aero_tracers() +
                  mam_coupling::num_aero_modes());

int tracer_slot(const TracerExpectation& tracer) {
  return tracer.is_number
             ? mam4::AeroConfig::numptr_amode(tracer.mode)
             : mam4::AeroConfig::lmassptr_amode(tracer.species, tracer.mode);
}

std::string interstitial_field_name(const TracerExpectation& tracer) {
  return tracer.is_number
             ? mam_coupling::int_aero_nmr_field_name(tracer.mode)
             : mam_coupling::int_aero_mmr_field_name(tracer.mode,
                                                      tracer.species);
}

std::string cloudborne_field_name(const TracerExpectation& tracer) {
  return tracer.is_number
             ? mam_coupling::cld_aero_nmr_field_name(tracer.mode)
             : mam_coupling::cld_aero_mmr_field_name(tracer.mode,
                                                      tracer.species);
}

struct Source {
  std::string name;
  std::vector<std::string> sectors;
  std::string path;
  Real scale_factor = 1;
};

std::vector<Source> sources() {
  std::vector<Source> result = {
      {"dms", {"DMS"}, "", 1},
      {"so2", {"AGR", "RCO", "SHP", "SLV", "TRA", "WST"}, "", 1},
      {"bc_a4", {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"}, "", 1},
      {"num_a1", {"num_a1_SO4_AGR", "num_a1_SO4_SHP", "num_a1_SO4_SLV",
                   "num_a1_SO4_WST"}, "", 1},
      {"num_a2", {"num_a2_SO4_RCO", "num_a2_SO4_TRA"}, "", 1},
      {"num_a4", {}, "", 1},
      {"pom_a4", {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"}, "", 1},
      {"so4_a1", {"AGR", "SHP", "SLV", "WST"}, "", 1},
      {"so4_a2", {"RCO", "TRA"}, "", 1}};
  for (const auto& species : {"BC", "POM"})
    for (const auto& sector : {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"})
      result[5].sectors.push_back(std::string("num_a1_") + species + "_" + sector);
  for (auto& src : result) {
    src.path = data_dir + "surface/";
    if (src.name == "dms")
      src.path +=
          "DMSflux.2010.ne2np4_conserv.POPmonthlyClimFromACES4BGC_"
          "c20260730.nc";
    else
      src.path += "cmip6_mam4_" + src.name +
                  "_surf_ne2np4_2010_clim_c20260730.nc";
  }
  return result;
}

void nc_check(const int status) {
  if (status != NC_NOERR) throw std::runtime_error(nc_strerror(status));
}

struct NcFile {
  int id;
  explicit NcFile(const std::string& path, const int mode = NC_NOWRITE) {
    nc_check(nc_open(path.c_str(), mode, &id));
  }
  ~NcFile() { nc_close(id); }
  NcFile(const NcFile&) = delete;
  NcFile& operator=(const NcFile&) = delete;
  int var(const std::string& name) const {
    int v;
    nc_check(nc_inq_varid(id, name.c_str(), &v));
    return v;
  }
};

// Fixture validation only; this does not add production input guards.
void audit(const Source& src) {
  NcFile file(src.path);
  for (const auto& sector : src.sectors) {
    const int var = file.var(sector);
    int rank;
    nc_check(nc_inq_varndims(file.id, var, &rank));
    if (rank != 2) throw std::runtime_error("Expected time,ncol sector");
    int dims[2];
    nc_check(nc_inq_vardimid(file.id, var, dims));
    for (int d = 0; d < 2; ++d) {
      size_t len;
      char name[NC_MAX_NAME + 1];
      nc_check(nc_inq_dim(file.id, dims[d], name, &len));
      if (len != size_t(d == 0 ? 12 : ncols) ||
          std::string(name) != (d == 0 ? "time" : "ncol"))
        throw std::runtime_error("Unexpected sector dimensions");
    }
    size_t len;
    nc_check(nc_inq_attlen(file.id, var, "units", &len));
    if (len == 0 || len > 256)
      throw std::runtime_error("Invalid fixture units");
    std::string units(len, '\0');
    nc_check(nc_get_att_text(file.id, var, "units", &units[0]));
    if (src.name.rfind("num_", 0) == 0 &&
        units != "(particles/cm2/s) * 6.022e26")
      throw std::runtime_error("Unexpected encoded number units");
    std::vector<double> values(12 * ncols);
    nc_check(nc_get_var_double(file.id, var, values.data()));
    for (const double value : values)
      if (!std::isfinite(value) || value < 0)
        throw std::runtime_error("Invalid fixture source value");
  }
}

// Unique short scratch outside the source tree. Remove only fixture-owned files.
struct Scratch {
  std::string dir;
  std::vector<std::string> files;
  Scratch() {
    char pattern[] = "/tmp/mam-emissions-XXXXXX";
    const auto path = mkdtemp(pattern);
    if (!path) throw std::runtime_error("Cannot create emission test scratch");
    dir = path;
  }
  ~Scratch() {
    for (const auto& path : files) std::remove(path.c_str());
    rmdir(dir.c_str());
  }
  Source copy(const Source& src) {
    Source dst = src;
    dst.path = dir + "/" + std::to_string(files.size()) + ".nc";
    files.push_back(dst.path);
    std::ifstream input(src.path, std::ios::binary);
    std::ofstream output(dst.path, std::ios::binary);
    if (!input || !output)
      throw std::runtime_error("Cannot open emission fixture copy");
    output << input.rdbuf();
    output.close();
    if (!output || input.bad())
      throw std::runtime_error("Cannot copy emission fixture");
    return dst;
  }
};

void set_sources(const Source& src, const double value) {
  NcFile file(src.path, NC_WRITE);
  std::vector<double> values(12 * ncols, value);
  for (const auto& sector : src.sectors)
    nc_check(nc_put_var_double(file.id, file.var(sector), values.data()));
}

double stored_sector_sum_at_origin(const Source& src) {
  NcFile file(src.path);
  const size_t first[2] = {0, 0};
  double sum = 0;
  for (const auto& sector : src.sectors) {
    double value;
    nc_check(nc_get_var1_double(file.id, file.var(sector), first, &value));
    sum += value;
  }
  return sum;
}

ekat::ParameterList parameters(const std::vector<Source>& inputs) {
  ekat::ParameterList p;
  p.set<std::string>("log_level", "warn");
  p.set<std::string>("srf_remap_file", "");
  p.set<std::string>("soil_erodibility_file",
                     data_dir + "dst_ne2np4_c20241028.nc");
  p.set<std::string>(
      "marine_organics_file",
      data_dir +
          "monthly_macromolecules_0.1deg_bilinear_year01_merge_ne2np4_"
          "c20260807.nc");
  p.set("dust_emis_scheme", 2);  // Controlled erodibility is exactly one.
  p.set("srf_emis_scale_factor_for_dust", 1.5);
  p.set("srf_emis_scale_factor_for_seasalt", 0.6);
  for (const auto& src : inputs)
    p.set<std::string>("srf_emis_specifier_for_" + src.name, src.path);
  return p;
}

std::shared_ptr<GridsManager> make_grid(const ekat::Comm& comm) {
  ekat::ParameterList p;
  p.set("grids_names", std::vector<std::string>{"point_grid"});
  auto& grid = p.sublist("point_grid");
  grid.set<std::string>("type", "point_grid");
  grid.set("aliases", std::vector<std::string>{"physics"});
  grid.set("number_of_global_columns", ncols);
  grid.set("number_of_vertical_levels", nlevs);
  auto gm = create_mesh_free_grids_manager(comm, p);
  gm->build_grids();
  return gm;
}

Real comparison_tolerance(const Real actual, const Real expected) {
  const Real scale =
      std::max({std::abs(actual), std::abs(expected), Real(1)});
  return 10 * std::numeric_limits<Real>::epsilon() * scale;
}

void check_close(const Real actual, const Real expected) {
  REQUIRE(std::isfinite(actual));
  CHECK(std::abs(actual - expected) <=
        comparison_tolerance(actual, expected));
}

using Output = std::array<std::vector<Real>, 2>;

Real tracer_baseline(const size_t index,
                     const TracerExpectation& tracer) {
  return tracer.is_number ? Real(1000 + index) : Real(index + 1) * Real(1e-9);
}

void reset_aerosol_fields(FieldManager& fm) {
  for (size_t t = 0; t < aerosol_expectations.size(); ++t) {
    const auto& tracer = aerosol_expectations[t];
    const Real baseline = tracer_baseline(t, tracer);
    fm.get_field(interstitial_field_name(tracer)).deep_copy(baseline);
    fm.get_field(cloudborne_field_name(tracer)).deep_copy(baseline);
  }
}

void check_constituent_consumer(FieldManager& fm,
                                const std::vector<Real>& fluxes) {
  constexpr Real g = physics::Constants<Real>::gravit.value;
  for (size_t t = 0; t < aerosol_expectations.size(); ++t) {
    const auto& tracer = aerosol_expectations[t];
    CAPTURE(tracer.name);
    const Real baseline = tracer_baseline(t, tracer);
    auto& interstitial = fm.get_field(interstitial_field_name(tracer));
    auto& cloudborne = fm.get_field(cloudborne_field_name(tracer));
    interstitial.sync_to_host();
    cloudborne.sync_to_host();
    const auto int_view = interstitial.get_view<const Real**, Host>();
    const auto cld_view = cloudborne.get_view<const Real**, Host>();
    const int slot = tracer_slot(tracer);
    for (int col = 0; col < ncols; ++col) {
      const Real flux = fluxes[col * nconst + slot];
      const Real expected = baseline + flux * dt * g / dp;
      if (flux == 0)
        CHECK(int_view(col, nlevs - 1) == baseline);
      else
        check_close(int_view(col, nlevs - 1), expected);
      for (int lev = 0; lev < nlevs - 1; ++lev)
        CHECK(int_view(col, lev) == baseline);
      for (int lev = 0; lev < nlevs; ++lev)
        CHECK(cld_view(col, lev) == baseline);
    }
  }
}

Output run(const ekat::Comm& comm, const std::vector<Source>& inputs,
           const bool online, const bool check_consumer = false) {
  auto gm = make_grid(comm);
  auto grid = gm->get_grid("physics");
  auto fm = std::make_shared<FieldManager>(grid);
  MAMSrfOnlineEmiss emission(comm, parameters(inputs));
  MAMConstituentFluxes consumer(comm, ekat::ParameterList());
  std::vector<AtmosphereProcess*> procs = {&emission};
  if (check_consumer) procs.push_back(&consumer);

  std::set<std::string> names;
  for (auto* proc : procs) {
    proc->set_grids(gm);
    for (const auto& req : proc->get_field_requests()) {
      fm->register_field(req);
      names.insert(req.fid.name());
    }
  }
  fm->registration_ends();

  const util::TimeStamp t0({2021, 10, 12}, {12, 30, 0});
  for (const auto& name : names) {
    auto& field = fm->get_field(name);
    field.deep_copy(Real(0));
    field.get_header().get_tracking().update_time_stamp(t0);
  }
  fm->get_field("T_mid").deep_copy(Real(280));
  fm->get_field("pseudo_density").deep_copy(dp);
  fm->get_field("pbl_height").deep_copy(Real(1000));
  fm->get_field("phis").deep_copy(Real(1));
  fm->get_field("sst").deep_copy(Real(300));

  for (const auto& name : {"p_int", "p_mid", "ocnfrac", "dstflx",
                           "horiz_winds"}) {
    auto& field = fm->get_field(name);
    field.sync_to_host();
    if (std::string(name) == "horiz_winds") {
      auto view = field.get_view<Real***, Host>();
      for (int col = 0; col < ncols; ++col)
        for (int lev = 0; lev < nlevs; ++lev) {
          view(col, 0, lev) = 8 + col % 5;
          view(col, 1, lev) = 3;
        }
    } else if (std::string(name) == "ocnfrac") {
      auto view = field.get_view<Real*, Host>();
      for (int col = 0; col < ncols; ++col)
        view(col) = online ? Real(col % 3) / 2 : 0;
    } else {
      auto view = field.get_view<Real**, Host>();
      for (int col = 0; col < ncols; ++col)
        for (int index = 0; index < int(view.extent(1)); ++index) {
          if (std::string(name) == "dstflx")
            view(col, index) =
                online ? -Real(1e-10) * (index + 1) *
                             (1 - Real(col % 3) / 2)
                       : 0;
          else
            view(col, index) =
                1000 + dp *
                           (index +
                            (std::string(name) == "p_mid" ? Real(0.5)
                                                           : Real(0)));
        }
    }
    field.sync_to_dev();
  }

  std::array<ATMBufferManager, 2> buffers;
  for (size_t j = 0; j < procs.size(); ++j) {
    auto* proc = procs[j];
    for (const auto& req : proc->get_field_requests()) {
      auto field = fm->get_field(req.fid.name());
      if (req.usage & Required) proc->set_required_field(field.get_const());
      if (req.usage & Computed) proc->set_computed_field(field);
    }
    buffers[j].request_bytes(proc->requested_buffer_size_in_bytes());
    buffers[j].allocate();
    proc->init_buffers(buffers[j]);
    proc->initialize(t0, RunType::Initial);
  }

  Output output;
  for (int step = 0; step < 2; ++step) {
    auto& flux = fm->get_field("constituent_fluxes");
    flux.sync_to_host();
    auto flux_view = flux.get_view<Real**, Host>();
    // A stale nonzero sentinel must be discarded on every invocation.
    for (int col = 0; col < ncols; ++col)
      for (int slot = mam4::utils::gasses_start_ind(); slot < nconst; ++slot)
        flux_view(col, slot) = Real(100 + step);
    flux.sync_to_dev();

    emission.run(dt);
    flux.sync_to_host();
    output[step].resize(ncols * nconst);
    for (int col = 0; col < ncols; ++col)
      for (int slot = 0; slot < nconst; ++slot)
        output[step][col * nconst + slot] = flux_view(col, slot);

    if (check_consumer) {
      reset_aerosol_fields(*fm);
      consumer.run(dt);
      check_constituent_consumer(*fm, output[step]);
    }
  }

  for (auto* proc : procs) proc->finalize();
  return output;
}

void check_expectation_table() {
  std::set<int> slots;
  std::array<int, 4> class_counts = {{0, 0, 0, 0}};
  for (const auto& tracer : aerosol_expectations) {
    CAPTURE(tracer.name);
    const auto field_name = interstitial_field_name(tracer);
    if (!tracer.is_number &&
        mam4::mode_aero_species(tracer.mode, tracer.species) ==
            mam4::AeroId::NaCl) {
      // The surface-flux symbols use ncl_a*, while the established EAMxx
      // prognostic-field mapping uses nacl_a*. Keep the requested flux labels
      // explicit and verify this repository alias at the consumer boundary.
      const std::string mode = std::to_string(tracer.mode + 1);
      REQUIRE(std::string(tracer.name) == "ncl_a" + mode);
      REQUIRE(field_name == "nacl_a" + mode);
    } else {
      REQUIRE(field_name == tracer.name);
    }
    REQUIRE(slots.insert(tracer_slot(tracer)).second);
    ++class_counts[static_cast<int>(tracer.source_class)];
  }
  CHECK(class_counts[static_cast<int>(SourceClass::Both)] == 2);
  CHECK(class_counts[static_cast<int>(SourceClass::OnlineOnly)] == 8);
  CHECK(class_counts[static_cast<int>(SourceClass::PrescribedOnly)] == 5);
  CHECK(class_counts[static_cast<int>(SourceClass::Inactive)] == 10);
}

void check_additivity(const Output& combined, const Output& online,
                      const Output& prescribed) {
  for (int step = 0; step < 2; ++step) {
    for (const auto& tracer : aerosol_expectations) {
      CAPTURE(step, tracer.name);
      const int slot = tracer_slot(tracer);
      Real max_online = 0;
      Real max_prescribed = 0;
      int overlap_columns = 0;
      for (int col = 0; col < ncols; ++col) {
        const int index = col * nconst + slot;
        const Real o = online[step][index];
        const Real p = prescribed[step][index];
        const Real c = combined[step][index];
        REQUIRE(std::isfinite(o));
        REQUIRE(std::isfinite(p));
        REQUIRE(std::isfinite(c));
        REQUIRE(o >= 0);
        REQUIRE(p >= 0);
        max_online = std::max(max_online, o);
        max_prescribed = std::max(max_prescribed, p);
        if (o > 0 && p > 0) ++overlap_columns;

        check_close(c, o + p);
        if (tracer.source_class == SourceClass::OnlineOnly ||
            tracer.source_class == SourceClass::Inactive)
          CHECK(p == 0);
        if (tracer.source_class == SourceClass::PrescribedOnly ||
            tracer.source_class == SourceClass::Inactive)
          CHECK(o == 0);
        if (tracer.source_class == SourceClass::Inactive) CHECK(c == 0);
      }

      if (tracer.source_class == SourceClass::Both ||
          tracer.source_class == SourceClass::OnlineOnly)
        CHECK(max_online > 0);
      if (tracer.source_class == SourceClass::Both ||
          tracer.source_class == SourceClass::PrescribedOnly)
        CHECK(max_prescribed > 0);
      if (tracer.source_class == SourceClass::Both)
        CHECK(overlap_columns > 0);
    }
  }
}

void check_conversion_oracles(const ekat::Comm& comm,
                              const std::vector<Source>& zero,
                              Scratch& scratch) {
  constexpr std::array<const char*, 7> prescribed_targets = {{
      "so4_a1", "so4_a2", "bc_a4", "pom_a4", "num_a1", "num_a2",
      "num_a4"}};
  constexpr int offset =
      mam4::aero_model::pcnst - mam4::gas_chemistry::gas_pcnst;

  for (size_t target_index = 0; target_index < prescribed_targets.size();
       ++target_index) {
    const std::string target_name = prescribed_targets[target_index];
    CAPTURE(target_name);
    auto controlled = zero;
    auto source_it = std::find_if(
        controlled.begin(), controlled.end(),
        [&target_name](const Source& src) { return src.name == target_name; });
    REQUIRE(source_it != controlled.end());
    *source_it = scratch.copy(*source_it);
    set_sources(*source_it, 1e20 + target_index * 1e19);
    audit(*source_it);

    const auto actual = run(comm, controlled, false);
    const auto tracer_it = std::find_if(
        aerosol_expectations.begin(), aerosol_expectations.end(),
        [&target_name](const TracerExpectation& tracer) {
          return tracer.name == target_name;
        });
    REQUIRE(tracer_it != aerosol_expectations.end());
    const int target_slot = tracer_slot(*tracer_it);
    const Real stored_sum = stored_sector_sum_at_origin(*source_it);
    const Real expected =
        stored_sum * source_it->scale_factor * amufac *
        mam4::gas_chemistry::adv_mass[target_slot - offset];
    REQUIRE(expected > 0);

    for (int step = 0; step < 2; ++step)
      for (int col = 0; col < ncols; ++col)
        for (const auto& tracer : aerosol_expectations) {
          const Real value =
              actual[step][col * nconst + tracer_slot(tracer)];
          if (tracer.name == target_name)
            check_close(value, expected);
          else
            CHECK(value == 0);
        }
  }
}

struct ScorpioGuard {
  explicit ScorpioGuard(const ekat::Comm& comm) {
    scorpio::init_subsystem(comm);
  }
  ~ScorpioGuard() { scorpio::finalize_subsystem(); }
};

}  // namespace

TEST_CASE("mam_surface_emissions_all_aerosol_sources",
          "[mam][emissions]") {
  using namespace scream;
  ekat::Comm comm(MPI_COMM_WORLD);
  // The complete 218-column fixture is intentionally local to one rank.
  REQUIRE(comm.size() == 1);
  REQUIRE(mam4::nlev == nlevs);
  check_expectation_table();
  ScorpioGuard guard(comm);

  Scratch scratch;
  const auto original = sources();
  auto zero = original;
  for (size_t index = 0; index < original.size(); ++index) {
    audit(original[index]);
    zero[index] = scratch.copy(original[index]);
    set_sources(zero[index], 0);
    audit(zero[index]);
  }

  const auto online = run(comm, zero, true);
  const auto prescribed = run(comm, original, false);
  const auto combined = run(comm, original, true, true);
  check_additivity(combined, online, prescribed);

  const auto no_source = run(comm, zero, false);
  for (int step = 0; step < 2; ++step)
    for (int col = 0; col < ncols; ++col)
      for (int slot = mam4::utils::gasses_start_ind(); slot < nconst; ++slot)
        CHECK(no_source[step][col * nconst + slot] == 0);

  check_conversion_oracles(comm, zero, scratch);
}
