#include <catch2/catch.hpp>

#include "share/data_managers/field_manager.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field.hpp"
#include "share/field/field_reader.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_comm.hpp>
#include <ekat_fpe.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_units.hpp>

#include <cfenv>
#include <fstream>
#include <string>
#include <vector>

/*
 * Fields sampled at cloud top/bottom with aerocom(..) reach the file, and
 * requests at the same position and phase share one weights diag.
 */

namespace scream {

namespace {

constexpr int ncols = 4;
constexpr int nlevs = 4;

// One cloudy layer per column, at level 1 or 2, with cloud fraction 0.5 and
// three quarters of the condensate liquid. Every other level has no cloud,
// so a division by cldfrac_tot there would be a division by zero.
int cloudy_level (const int gid) { return 1 + gid%2; }
Real T_at (const int k) { return 250 + 10*k; }
constexpr Real cld = 0.5;
constexpr Real phi = 0.75;
constexpr Real pden = 2000;
constexpr Real qv_val = 0.01;
constexpr Real nc_val = 1e8;
Real p_at (const int k) { return 20000 + 20000*k; }

// The names of the diags in the repo that all output streams share, by diag class
struct DiagRepo : public AtmosphereOutput {
  static std::map<std::string,std::string> classes () {
    std::map<std::string,std::string> c;
    for (const auto& it : m_diag_repo) {
      c[it.first] = it.second->name();
    }
    return c;
  }
};

} // anonymous namespace

TEST_CASE ("io_dexpr_aerocom")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  register_diagnostics();

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ncols);
  gm->build_grids();
  auto grid = gm->get_grid("point_grid");
  const auto gname = grid->name();
  const auto layout3d = grid->get_3d_scalar_layout(LEV);
  const auto layout2d = grid->get_2d_scalar_layout();
  const int nlocal = grid->get_num_local_dofs();
  const auto gids = grid->get_dofs_gids().get_view<const int*,Host>();

  util::TimeStamp t0({2023,1,1},{0,0,0});

  auto fm = std::make_shared<FieldManager>(gm);
  const std::vector<std::pair<std::string,ekat::units::Units>> inputs = {
    {"T_mid",K}, {"p_mid",Pa}, {"pseudo_density",Pa}, {"qv",kg/kg},
    {"qc",kg/kg}, {"qi",kg/kg}, {"nc",1/kg}, {"cldfrac_tot",none}
  };
  for (const auto& [n,u] : inputs) {
    Field f (FieldIdentifier(n,layout3d,u,gname));
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    auto v = f.get_view<Real**,Host>();
    for (int i=0; i<nlocal; ++i) {
      for (int k=0; k<nlevs; ++k) {
        v(i,k) = n=="T_mid" ? T_at(k)
               : n=="p_mid" ? p_at(k)
               : n=="pseudo_density" ? pden
               : n=="qv" ? qv_val
               : n=="qc" ? phi*1e-4
               : n=="qi" ? (1-phi)*1e-4
               : n=="nc" ? nc_val
               : (k==cloudy_level(gids(i)) ? cld : 0);
      }
    }
    f.sync_to_dev();
    fm->add_field(f);
  }

  ekat::ParameterList params;
  params.set<std::string>("filename_prefix","io_dexpr_aerocom");
  params.set<std::string>("averaging_type","INSTANT");
  params.set<std::string>("floating_point_precision","real");
  params.sublist("fields").sublist(gname).set<std::vector<std::string>>("field_names",{
    "T_top      := T_mid.aerocom('top')",
    "T_bot      := T_mid.aerocom('bot')",
    "cld_top    := mask.aerocom('top')",
    "cld_liq    := mask.aerocom('top', phase='liq')",
    "nc_top     := nc.aerocom('top', phase='liq')",
    // In-cloud droplet number per volume: the cloud fraction is zero at all
    // but one level, and masked there, so nothing is divided by zero
    "cdnc_top   := (nc*pseudo_density/dz/gravit/cldfrac_tot.where(cldfrac_tot>0)).aerocom('top', phase='liq')",
  });
  auto& ctrl_pl = params.sublist("output_control");
  ctrl_pl.set<std::string>("frequency_units","nsteps");
  ctrl_pl.set<int>("frequency",1);
  ctrl_pl.set<bool>("save_grid_data",false);
  ctrl_pl.set<bool>("skip_t0_output",true);

  // NOTE: scoped, so the stream (and the diags it owns) are gone before Kokkos
  //       is finalized
  {
    OutputManager om;
    om.initialize(comm,params,t0,false);
    om.setup(fm,gm->get_grid_names());

    // Six requests at three (position,phase) pairs: three weights diags
    std::map<std::string,int> nclass;
    for (const auto& [name,cls] : DiagRepo::classes()) {
      ++nclass[cls];
    }
    REQUIRE (nclass["AeroComCld"]==6);
    REQUIRE (nclass["AeroComWeights"]==3);
    const auto repo = DiagRepo::classes();
    for (const std::string w : {"cldfrac_tot.aerocom_weights('top')",
                                "cldfrac_tot.aerocom_weights('top',phase='liq')",
                                "cldfrac_tot.aerocom_weights('bot')"}) {
      REQUIRE (repo.count(w)==1);
    }

    auto t = t0 + 1;
    om.init_timestep(t0,1);
    for (const auto& [n,u] : inputs) {
      fm->get_field(n).get_header().get_tracking().update_time_stamp(t);
    }
    // Trap a division by zero (where enabled: not in device kernels), so that
    // the cdnc request is shown to divide by no zero
    const auto fpes = ekat::get_enabled_fpes();
    ekat::enable_fpes(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    om.run(t);
    ekat::disable_all_fpes();
    ekat::enable_fpes(fpes);
    om.finalize();
  }

  const auto filename = "io_dexpr_aerocom.INSTANT.nsteps_x1.np" +
                        std::to_string(comm.size()) + "." + (t0+1).to_string() + ".nc";
  std::ifstream file_check(filename);
  REQUIRE (file_check.good());
  file_check.close();

  {
    std::vector<Field> out;
    for (const std::string n : {"T_top","T_bot","cld_top","cld_liq","nc_top","cdnc_top"}) {
      out.emplace_back(FieldIdentifier(n,layout2d,none,gname));
      out.back().allocate_view();
    }
    FieldReader reader;
    reader.set_file_specs(filename);
    reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
    reader.set_fields(out);
    reader.read(0);
    reader.clean_up();

    auto v = [&](const int idx, const int i) {
      out[idx].sync_to_host();
      return out[idx].get_view<const Real*,Host>()(i);
    };
    using PF = PhysicsFunctions<HostDevice>;
    using C  = physics::Constants<Real>;
    for (int i=0; i<nlocal; ++i) {
      // A single cloudy layer: seen from either side, it has weight 0.5
      const int k = cloudy_level(gids(i));
      const Real dz = PF::calculate_dz(pden,p_at(k),T_at(k),qv_val);
      REQUIRE (v(0,i)==Approx(cld*T_at(k)));
      REQUIRE (v(1,i)==Approx(cld*T_at(k)));
      REQUIRE (v(2,i)==Approx(cld));
      REQUIRE (v(3,i)==Approx(cld*phi));
      REQUIRE (v(4,i)==Approx(cld*phi*nc_val));
      REQUIRE (v(5,i)==Approx(cld*phi*nc_val*pden/dz/C::gravit.value/cld));
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
